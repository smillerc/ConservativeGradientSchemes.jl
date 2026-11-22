using ConservativeGradientSchemes
using KernelAbstractions: CPU, synchronize
using Test

@testset "MEG6 metric boundaries without halos" begin
    backend = CPU()
    scheme = MEG{6}()
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

    # A tile edge inside the valid geometry must use the interior stencil.
    tile = CartesianIndices((5:8,))
    empty_domain = CartesianIndices((1:0,))
    prepare_metric_derivatives!(derivative, second_derivative, scheme, samples, 1, centers, tile, empty_domain, cells, backend)
    synchronize(backend)
    @test derivative[5:8] == first_derivative[5:8]
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
