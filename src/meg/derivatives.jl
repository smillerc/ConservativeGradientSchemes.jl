@inline shift(I::CartesianIndex{N}, axis, offset) where {N} =
    CartesianIndex(ntuple(d -> I[d] + ifelse(d == axis, offset, 0), Val(N)))

@inline spacing(q::T, ::Nothing, axis) where {T} = one(T)
@inline spacing(q::T, centers, axis) where {T} =
    T(centers[axis][firstindex(centers[axis]) + 1] - first(centers[axis]))

@inline difference(::Fluid, a, b) = a - b

@inline function difference(::Metric, a::T, b::T) where {T}
    return (a - b) * !isapprox(a, b; rtol = sqrt(eps(T)), atol = 10 * eps(T))
end

@inline compensated_sum(values) = compensated_sum(values, zero(eltype(values)))
@inline compensated_sum(values::NTuple{N, T}) where {N, T} = compensated_sum(values, zero(T))

@inline function compensated_sum(values, initial)
    total = initial
    correction = initial
    for value in values
        adjusted = value - correction
        next = total + adjusted
        correction = (next - total) - adjusted
        total = next
    end
    return total + correction
end

@inline function centered_difference(tag, q, I, axis, offset)
    return difference(tag, q[shift(I, axis, offset)], q[shift(I, axis, -offset)])
end

@inline centered_first(::MEG{2}, q, I, axis) = centered_difference(Fluid(), q, I, axis, 1) / 2

# The higher order fluid stencils match Eddy's reconstruction/meg.jl.
@inline function centered_first(::MEG{4}, q, I, axis)
    d1 = centered_difference(Fluid(), q, I, axis, 1)
    d2 = centered_difference(Fluid(), q, I, axis, 2)
    return (8d1 - d2) / 12
end

@inline function centered_first(::MEG{6}, q, I, axis)
    d1 = centered_difference(Fluid(), q, I, axis, 1)
    d2 = centered_difference(Fluid(), q, I, axis, 2)
    d3 = centered_difference(Fluid(), q, I, axis, 3)
    return (45d1 - 9d2 + d3) / 60
end

@inline function centered_first(::MEG{8}, q, I, axis)
    d1 = centered_difference(Fluid(), q, I, axis, 1)
    d2 = centered_difference(Fluid(), q, I, axis, 2)
    d3 = centered_difference(Fluid(), q, I, axis, 3)
    d4 = centered_difference(Fluid(), q, I, axis, 4)
    return (672d1 - 168d2 + 32d3 - 3d4) / 840
end

@inline function centered_metric_first(scheme::MEG, q, I, axis)
    weights, divisor = metric_centered_weights(scheme)
    T = eltype(q)
    terms = ntuple(Val(length(weights))) do k
        (T(weights[k]) / T(divisor)) * centered_difference(Metric(), q, I, axis, k)
    end
    value = compensated_sum(terms)
    return value * (abs(value) >= eps(T))
end

@inline function centered_metric_first(::MEG{6}, q, I, axis)
    d1 = centered_difference(Metric(), q, I, axis, 1)
    d2 = centered_difference(Metric(), q, I, axis, 2)
    d3 = centered_difference(Metric(), q, I, axis, 3)
    return centered_first_value(Metric(), d1, d2, d3)
end

@inline function centered_first_value(::Metric, a::T, b::T, c::T) where {T}
    value = compensated_sum((T(3) / 4 * a, -T(3) / 20 * b, c / 60))
    return value * (abs(value) >= eps(T))
end

# The six-point metric closures are reflected at the high boundary.
@inline function inward_first(::MEG{6}, q, I, axis, depth, direction)
    J = shift(I, axis, -direction * depth)
    a, b, c, d, e, f = ntuple(k -> q[shift(J, axis, direction * (k - 1))], Val(6))
    T = typeof(a)
    if depth == 0
        value = compensated_sum((-T(137) / 60 * a, 5b, -5c, T(10) / 3 * d, -T(5) / 4 * e, f / 5))
    elseif depth == 1
        value = compensated_sum((-a / 5, -T(13) / 12 * b, 2c, -d, e / 3, -f / 20))
    else
        value = compensated_sum((a / 20, -b / 2, -c / 3, d, -e / 4, f / 30))
    end
    return direction * value * (abs(value) >= 50 * eps(T))
end

@inline function inward_first(scheme::MEG{Order}, q, I, axis, depth, direction) where {Order}
    J = shift(I, axis, -direction * depth)
    samples = ntuple(k -> q[shift(J, axis, direction * (k - 1))], Val(Order + 1))
    T = eltype(samples)
    return direction * metric_inward_value(scheme, samples, depth, Val(T), 50 * eps(T))
end

@inline function directional_first(scheme::MEG{Order}, q, I, axis, lo, hi) where {Order}
    lower_depth = I[axis] - lo
    upper_depth = hi - I[axis]
    if lower_depth < Order ÷ 2
        return inward_first(scheme, q, I, axis, lower_depth, 1)
    elseif upper_depth < Order ÷ 2
        return inward_first(scheme, q, I, axis, upper_depth, -1)
    end
    return centered_metric_first(scheme, q, I, axis)
end

@inline curvature(::Fluid, a, b, c, dm, dp) = 2 * (c - 2b + a) - (dp - dm) / 2

@inline function curvature(scheme::Metric, a, b, c, dm, dp)
    return compensated_sum((2 * (c + a), -4b, -difference(scheme, dp, dm) / 2))
end

@inline function directional_second(scheme::MEG{Order}, q, first, I, axis, h, lo, hi) where {Order}
    if I[axis] - lo < Order ÷ 2
        return inward_first(scheme, first, I, axis, I[axis] - lo, 1) / h
    elseif hi - I[axis] < Order ÷ 2
        return inward_first(scheme, first, I, axis, hi - I[axis], -1) / h
    end
    return interior_curvature(Metric(), q, first, I, axis, h) / h^2
end

@inline function interior_curvature(scheme, q, first, I, axis, h)
    m, p = shift(I, axis, -1), shift(I, axis, 1)
    return curvature(scheme, q[m], q[I], q[p], h * first[m], h * first[p])
end
