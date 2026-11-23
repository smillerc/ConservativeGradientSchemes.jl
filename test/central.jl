@testset "Central gradients and conservative faces" begin
    backend = CPU()
    for order in (4, 6, 8)
        @testset "Order $order" begin
            scheme = Central{order}()
            x = collect(range(-1.0, 1.0; length = 19))
            centers = (x,)
            cells = CartesianIndices(x)
            edge = Vector{Float64}(undef, length(x) + 1)
            derivative = similar(x)

            # Monomial reproduction includes the complete one-sided region.
            for degree in 0:order
                q = x .^ degree
                expected = degree == 0 ? zeros(length(x)) : degree .* x .^ (degree - 1)
                metric_cell_derivative!(derivative, scheme, q, 1, centers, cells, cells, backend)
                synchronize(backend)
                @test maximum(abs, derivative .- expected) < 3.0e-12
                I = CartesianIndex(7)
                @test first(curvilinear_cell_gradient(scheme, q, I, centers)) ≈ expected[I] atol = 3.0e-12
            end

            q = exp.(x)
            interpolate_metric_edges!(edge, scheme, q, 1, CartesianIndices(edge), cells, backend)
            prepare_first_derivative!(derivative, scheme, q, 1, centers, cells, cells, backend)
            synchronize(backend)
            @test derivative ≈ diff(edge) ./ (x[2] - x[1]) atol = 3.0e-12
            @test sum(derivative) * (x[2] - x[1]) ≈ last(edge) - first(edge) atol = 3.0e-12
            tile = CartesianIndices((2:7,))
            tiled = fill(NaN, length(x))
            metric_cell_derivative!(tiled, scheme, q, 1, centers, tile, cells, backend)
            synchronize(backend)
            @test tiled[2:7] == derivative[2:7]
            @test isnan(tiled[1]) && all(isnan, tiled[8:end])

            errors = map((16, 32)) do n
                nodes = collect(range(0.0, 1.0; length = n))
                values = exp.(3 .* nodes)
                d = similar(values)
                indices = CartesianIndices(values)
                metric_cell_derivative!(d, scheme, values, 1, (nodes,), indices, indices, backend)
                synchronize(backend)
                maximum(abs, d .- 3 .* values)
            end
            @test log2(errors[1] / errors[2]) > order - 0.8

            logical = (similar(x),)
            physical = (similar(x),)
            i = CartesianIndex(7)
            Jinv = reshape([2.0], 1, 1)
            curvilinear_cell_gradient!(logical, scheme, x, i, centers)
            physical_cell_gradient!(physical, scheme, x, i, centers, Jinv)
            @test logical[1][i] ≈ 1.0
            @test physical[1][i] ≈ 2.0
        end
    end
end

@testset "Central 3D conservative metric identity" begin
    backend = CPU()
    centers = map(n -> collect(range(0.0, 1.0; length = n)), (11, 12, 13))
    ξ, η, ζ = centers
    warp(a, b, c) = sinpi(a) * sinpi(b) * sinpi(c)
    coordinates = (
        [a + 0.03 * warp(a, b, c) for a in ξ, b in η, c in ζ],
        [b - 0.02 * warp(a, b, c) for a in ξ, b in η, c in ζ],
        [c + 0.01 * warp(a, b, c) for a in ξ, b in η, c in ζ],
    )
    cells = CartesianIndices(first(coordinates))

    for order in (4, 6, 8)
        scheme = Central{order}()
        function D(q, axis)
            dest = similar(q)
            metric_cell_derivative!(dest, scheme, q, axis, centers, cells, cells, backend)
            synchronize(backend)
            return dest
        end

        # The conservative metric columns are discrete curls. Their divergence
        # must vanish when the same directional operators are used throughout.
        residuals = map(((2, 3), (3, 1), (1, 2))) do (j, k)
            a, b = coordinates[j], coordinates[k]
            p = ntuple(axis -> a .* D(b, axis), 3)
            m1 = D(p[3], 2) - D(p[2], 3)
            m2 = D(p[1], 3) - D(p[3], 1)
            m3 = D(p[2], 1) - D(p[1], 2)
            D(m1, 1) + D(m2, 2) + D(m3, 3)
        end
        @test maximum(r -> maximum(abs, r), residuals) < 1.0e-10
        # Contract the metric divergence with a constant physical flux.
        @test maximum(abs, 0.3 * residuals[1] - 0.2 * residuals[2] + 0.1 * residuals[3]) < 1.0e-10
    end
end
