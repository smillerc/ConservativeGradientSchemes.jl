using ConservativeGradientSchemes
using KernelAbstractions: @kernel, synchronize
using Metal: MetalBackend, MtlArray
using Test

@kernel function limited_face_kernel!(left, right, limiter, q, first_derivative, second_derivative)
    reconstruct_edge!(
        left, right, MEG{2}(), limiter, q, CartesianIndex(2), 1,
        nothing, nothing, first_derivative, second_derivative
    )
end

@testset "MEG Superbee and Van Leer on Metal" begin
    backend = MetalBackend()
    q = MtlArray(Float32[0, 1, 4, 5])
    first_derivative = MtlArray(fill(4.0f0, 4))
    second_derivative = MtlArray(zeros(Float32, 4))
    left, right = similar(q), similar(q)
    for (limiter, expected) in ((SuperbeeLimiter(), (2.0f0, 3.0f0)), (VanLeerLimiter(), (1.75f0, 3.25f0)))
        limited_face_kernel!(backend)(left, right, limiter, q, first_derivative, second_derivative; ndrange = 1)
        synchronize(backend)
        @test (Array(left)[2], Array(right)[2]) === expected
    end
end

@testset "MEG fluid derivative orders on Metal" begin
    backend = MetalBackend()
    x = collect(-1.0f0:0.125f0:1.0f0)
    centers = (MtlArray(x),)
    samples = MtlArray(x .+ x .^ 3)
    cells = CartesianIndices(size(samples))
    first_domain = CartesianIndices((8:11,))
    second_domain = CartesianIndices((9:10,))
    first_derivative, second_derivative = similar(samples), similar(samples)

    for order in (2, 4, 6, 8)
        scheme = MEG{order}()
        prepare_first_derivative!(first_derivative, scheme, samples, 1, centers, first_domain, cells, backend)
        prepare_second_derivative!(second_derivative, scheme, samples, first_derivative, 1, centers, second_domain, cells, backend)
        synchronize(backend)
        expected = 1 .+ 3 .* x[8:11] .^ 2
        order == 2 && (expected .+= 0.125f0^2)
        @test Array(first_derivative)[8:11] ≈ expected
        @test Array(second_derivative)[9:10] ≈ 6 .* x[9:10]
    end
end

@testset "MEG$order Float32 metric sweeps on Metal" for order in (2, 4, 6, 8)
    scheme = MEG{order}()
    backend = MetalBackend()
    centers = (0.0f0:0.25f0:2.75f0, 0.0f0:0.5f0:5.5f0, 0.0f0:0.125f0:1.375f0)
    faces = map(c -> range(first(c) - step(c) / 2; step = step(c), length = length(c) + 1), centers)
    samples = MtlArray([x + 2y - 3z for x in centers[1], y in centers[2], z in centers[3]])
    device_centers = map(c -> MtlArray(collect(c)), centers)
    device_faces = map(f -> MtlArray(collect(f)), faces)
    cells = CartesianIndices(size(samples))
    first_derivative = similar(samples)
    second_derivative = similar(samples)
    derivative = similar(samples)

    for axis in 1:3
        dims = ntuple(d -> size(samples, d) + (d == axis), 3)
        edges = similar(samples, dims)
        face_centers = ntuple(d -> d == axis ? faces[d] : centers[d], 3)
        expected_edges = [x + 2y - 3z for x in face_centers[1], y in face_centers[2], z in face_centers[3]]

        prepare_metric_derivatives!(
            first_derivative, second_derivative, scheme, samples, axis, device_centers,
            cells, cells, cells, backend
        )
        interpolate_metric_edges!(
            edges, scheme, samples, first_derivative, second_derivative,
            axis, device_centers, device_faces, CartesianIndices(dims), cells, backend
        )
        metric_cell_derivative!(
            derivative, scheme, samples, first_derivative, second_derivative,
            axis, device_centers, device_faces, cells, cells, backend
        )
        synchronize(backend)

        expected_derivative = (1.0f0, 2.0f0, -3.0f0)[axis]
        @test maximum(abs, Array(edges) .- expected_edges) < 2.0f-4
        @test maximum(abs, Array(derivative) .- expected_derivative) < 2.0f-3
    end
end

@testset "Central Float32 metric sweeps on Metal" begin
    backend = MetalBackend()
    centers = (collect(0.0f0:0.25f0:2.75f0), collect(1.0f0:3.0f0), collect(1.0f0:2.0f0))
    x, y, z = centers
    q = MtlArray([a + 2b - 3c for a in x, b in y, c in z])
    device_centers = map(MtlArray, centers)
    cells = CartesianIndices(size(q))
    derivative = similar(q)
    faces = similar(q, (length(x) + 1, length(y), length(z)))
    expected = [a + 2b - 3c for a in (-0.125f0:0.25f0:2.875f0), b in y, c in z]

    for order in (4, 6, 8)
        scheme = Central{order}()
        interpolate_metric_edges!(faces, scheme, q, 1, CartesianIndices(size(faces)), cells, backend)
        metric_cell_derivative!(derivative, scheme, q, 1, device_centers, cells, cells, backend)
        synchronize(backend)
        @test maximum(abs, Array(faces) .- expected) < 4.0e-5
        @test maximum(abs, Array(derivative) .- 1) < 2.0e-4
    end
end
