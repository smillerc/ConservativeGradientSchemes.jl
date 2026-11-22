using ConservativeGradientSchemes
using KernelAbstractions: synchronize
using Metal: MetalBackend, MtlArray
using Test

@testset "MEG6 Float32 metric sweeps on Metal" begin
    scheme = MEG{6}()
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
