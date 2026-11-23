using ConservativeGradientSchemes
using KernelAbstractions: CPU, synchronize
using Test

struct SampleFaceScheme <: AbstractScheme end

function ConservativeGradientSchemes.reconstruct_edge(::SampleFaceScheme, q::AbstractVector, I, axis, centers, faces)
    return q[I] + 2, q[I + CartesianIndex(1)] - 2
end

@testset "Shared limiter interface for any face reconstruction" begin
    q = [0.0, 0.0, 1.0, 4.0, 5.0, 5.0]
    I = CartesianIndex(3)
    scheme = SampleFaceScheme()
    left, right = similar(q), similar(q)
    for (limiter, expected) in (
            (MinmodLimiter(), (1.5, 3.5)),
            (SuperbeeLimiter(), (2.0, 3.0)),
            (VanLeerLimiter(), (1.75, 3.25)),
            (MP5Limiter(), (3.0, 2.0)),
        )
        @test limit_face_pair(limiter, 3.0, 2.0, q, I, 1) == expected
        @test reconstruct_edge(scheme, limiter, q, I, 1, nothing, nothing) == expected
        @test reconstruct_edge!(left, right, scheme, limiter, q, I, 1, nothing, nothing) === nothing
        @test (left[I], right[I]) == expected
    end
end

@testset "MEG fluid derivative orders" begin
    backend = CPU()
    cells = CartesianIndices((17,))
    I = CartesianIndex(9)
    first_domain = CartesianIndices((8:11,))
    second_domain = CartesianIndices((9:10,))

    for T in (Float32, Float64), order in (2, 4, 6, 8)
        scheme = MEG{order}()
        x = T.(-1:0.125:1)
        centers = (x,)
        faces = (T.(-1.0625:0.125:1.0625),)
        # At x = 0 the exact derivative is 1. The next odd degree detects
        # accidental use of a wider first derivative stencil.
        q = x .+ x .^ order
        @test only(curvilinear_cell_gradient(scheme, q, I, centers)) ≈ one(T)
        q = x .+ (8 .* x) .^ (order + 1) ./ T(8)
        @test !isapprox(only(curvilinear_cell_gradient(scheme, q, I, centers)), one(T))

        q = x .+ x .^ 3
        first_derivative, second_derivative = similar(q), similar(q)
        prepare_first_derivative!(first_derivative, scheme, q, 1, centers, first_domain, cells, backend)
        synchronize(backend)
        expected = 1 .+ 3 .* x[8:11] .^ 2
        order == 2 && (expected .+= T(0.125)^2)
        @test first_derivative[8:11] ≈ expected
        prepare_second_derivative!(second_derivative, scheme, q, first_derivative, 1, centers, second_domain, cells, backend)
        synchronize(backend)
        @test second_derivative[9:10] ≈ 6 .* x[9:10]

        for limiter in ((), (MinmodLimiter(),), (SuperbeeLimiter(),), (VanLeerLimiter(),), (MP5Limiter(),))
            direct = reconstruct_edge(scheme, limiter..., q, I, 1, centers, faces)
            prepared = reconstruct_edge(scheme, limiter..., q, I, 1, centers, faces, first_derivative, second_derivative)
            @test all(isapprox.(direct, prepared))
        end
    end
end

@testset "MEG Superbee and Van Leer slope bounds" begin
    scheme = MEG{2}()
    I = CartesianIndex(2)
    for T in (Float32, Float64), direction in (-1, 1)
        q = T[0, 1, 4, 5] .* direction
        first_derivative = fill(T(4 * direction), 4)
        second_derivative = zeros(T, 4)
        left, right = similar(q), similar(q)
        for (limiter, half_slope) in ((SuperbeeLimiter(), T(1)), (VanLeerLimiter(), T(0.75)))
            expected = (direction * (1 + half_slope), direction * (4 - half_slope))
            @test reconstruct_edge(scheme, limiter, q, I, 1, nothing, nothing, first_derivative, second_derivative) === expected
            reconstruct_edge!(left, right, scheme, limiter, q, I, 1, nothing, nothing, first_derivative, second_derivative)
            @test (left[I], right[I]) === expected
            @test reconstruct_edge(scheme, limiter, T[0, 1, 2, 3, 4, 5], CartesianIndex(3), 1) === (T(2.5), T(2.5))
            for flat_or_extremum in (T[1, 1, 4, 4], T[2, 1, 4, 3])
                @test reconstruct_edge(scheme, limiter, flat_or_extremum, I, 1, nothing, nothing, abs.(first_derivative), second_derivative) === (T(1), T(4))
            end
        end
    end
end

@testset "MEG$order metric boundaries without halos" for order in (2, 4, 6, 8)
    backend = CPU()
    scheme = MEG{order}()
    centers = (collect(0.0:0.25:2.75),)
    faces = (collect(-0.125:0.25:2.875),)
    samples = collect(only(centers))
    cells = CartesianIndices(samples)
    first_derivative = similar(samples)
    second_derivative = similar(samples)
    metric_edges = similar(samples, length(samples) + 1)
    derivative = similar(samples)

    prepare_metric_derivatives!(first_derivative, second_derivative, scheme, samples, 1, centers, cells, cells, cells, backend)
    synchronize(backend)
    interpolate_metric_edges!(
        metric_edges, scheme, samples, first_derivative, second_derivative,
        1, centers, faces, CartesianIndices(metric_edges), cells, backend
    )
    synchronize(backend)
    metric_cell_derivative!(
        derivative, scheme, samples, first_derivative, second_derivative,
        1, centers, faces, cells, cells, backend
    )
    synchronize(backend)

    @test first_derivative ≈ ones(length(samples))
    @test maximum(abs, second_derivative) < 1.0e-10
    @test metric_edges ≈ collect(only(faces))
    @test derivative ≈ ones(length(samples))
    @test sum(derivative) * 0.25 ≈ last(metric_edges) - first(metric_edges)

    # A tile edge inside the valid geometry must use the interior stencil.
    tile = CartesianIndices((5:8,))
    empty_domain = CartesianIndices((1:0,))
    prepare_metric_derivatives!(derivative, second_derivative, scheme, samples, 1, centers, tile, empty_domain, cells, backend)
    synchronize(backend)
    @test derivative[5:8] == first_derivative[5:8]
end

@testset "MEG$order metric first derivative polynomial exactness" for order in (2, 4, 6, 8)
    scheme = MEG{order}()
    backend = CPU()
    width = order == 6 ? 6 : order + 1
    degree = width - 1
    radius = order ÷ 2
    for n in (width, 17)
        x = collect(range(-1.0, 1.0; length = n))
        q = x .^ degree
        cells = CartesianIndices(q)
        first_derivative, second_derivative = similar(q), similar(q)
        prepare_metric_derivatives!(first_derivative, second_derivative, scheme, q, 1, (x,), cells, cells, cells, backend)
        synchronize(backend)
        expected = degree .* x .^ (degree - 1)
        h = x[2] - x[1]
        @test first_derivative ≈ expected
        @test all(isfinite, second_derivative)
        @test collect(metric_boundary_first(scheme, Val(:low), Tuple(q[1:width]))) ≈ h .* expected[1:radius]
        @test collect(metric_boundary_first(scheme, Val(:high), Tuple(q[(end - width + 1):end]))) ≈ h .* expected[(end - radius + 1):end]
        if n == 17
            i = 12
            stencil = Tuple(q[[((i - radius):(i - 1)); ((i + 1):(i + radius))]])
            @test metric_centered_first(scheme, stencil) ≈ h * expected[i] atol = 1.0e-12
        end
    end
end

@testset "MEG6 face conventions" begin
    q = collect(1.0:12.0)
    i = CartesianIndex(6)
    first_derivative = ones(length(q))
    second_derivative = fill(2.0, length(q))
    centers, faces = (1:12,), (0.5:12.5,)
    scheme = MEG{6}()

    fluid_pair = reconstruct_edge(scheme, q, i, 1, centers, faces, first_derivative, second_derivative)
    @test fluid_pair == (6.5 + 1 / 6, 6.5 + 1 / 6)
    @test reconstruct_boundary_edge(scheme, q, CartesianIndex(1), 1, centers, faces, Val(:low), first_derivative, second_derivative) == 1 - (1 / 2 + 1 / 6)
    @test reconstruct_edge(scheme, MinmodLimiter(), q, i, 1) == (6.5, 6.5)
    @test reconstruct_edge(scheme, MP5Limiter(), q, i, 1) == (6.5, 6.5)
    @test reconstruct_edge(scheme, q, i, 1, first_derivative) == (6.5, 6.5)
    @test reconstruct_edge(scheme, q, i, 1) == (6.5, 6.5)
end

@testset "MEG6 MP5 corrects an overshoot" begin
    scheme = MEG{6}()
    q = ones(12)
    first_derivative = zeros(12)
    second_derivative = fill(12.0, 12)
    i = CartesianIndex(6)
    centers, faces = (1:12,), (0.5:12.5,)

    @test reconstruct_edge(scheme, q, i, 1, centers, faces, first_derivative, second_derivative) == (2.0, 2.0)
    @test reconstruct_edge(scheme, MP5Limiter(), q, i, 1, centers, faces, first_derivative, second_derivative) == (1.0, 1.0)
end

@testset "MEG6 logical and physical gradients" begin
    scheme = MEG{6}()
    centers = (0.0:0.5:5.0, 0.0:0.25:2.5)
    q = [2x - 3y for x in centers[1], y in centers[2]]
    i = CartesianIndex(6, 6)
    inverse_jacobian = [2.0 1.0; 0.0 4.0]
    dest = (similar(q), similar(q))

    @test curvilinear_cell_gradient(scheme, q, i, centers) == (2.0, -3.0)
    @test physical_cell_gradient(scheme, q, i, centers, inverse_jacobian) == (4.0, -10.0)
    @test curvilinear_cell_gradient!(dest, scheme, q, i, centers) === nothing
    @test (dest[1][i], dest[2][i]) == (2.0, -3.0)
    @test physical_cell_gradient!(dest, scheme, q, i, centers, inverse_jacobian) === nothing
    @test (dest[1][i], dest[2][i]) == (4.0, -10.0)
end

include("central.jl")
