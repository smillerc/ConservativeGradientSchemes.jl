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

@inline function centered_first(scheme, q, I, axis)
    d1 = difference(scheme, q[shift(I, axis, 1)], q[shift(I, axis, -1)])
    d2 = difference(scheme, q[shift(I, axis, 2)], q[shift(I, axis, -2)])
    d3 = difference(scheme, q[shift(I, axis, 3)], q[shift(I, axis, -3)])
    return centered_first_value(scheme, d1, d2, d3)
end

@inline centered_first_value(::Fluid, a, b, c) = (45a - 9b + c) / 60

@inline function centered_first_value(::Metric, a::T, b::T, c::T) where {T}
    value = compensated_sum((T(3) / 4 * a, -T(3) / 20 * b, c / 60))
    return value * (abs(value) >= eps(T))
end

# The six-point metric closures are reflected at the high boundary.
@inline function inward_first(q, I, axis, depth, direction)
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

@inline function directional_first(scheme, q, I, axis, lo, hi)
    lower_depth = I[axis] - lo
    upper_depth = hi - I[axis]
    if lower_depth < 3
        return inward_first(q, I, axis, lower_depth, 1)
    elseif upper_depth < 3
        return inward_first(q, I, axis, upper_depth, -1)
    end
    return centered_first(scheme, q, I, axis)
end

@inline curvature(::Fluid, a, b, c, dm, dp) = 2 * (c - 2b + a) - (dp - dm) / 2

@inline function curvature(scheme::Metric, a, b, c, dm, dp)
    return compensated_sum((2 * (c + a), -4b, -difference(scheme, dp, dm) / 2))
end

@inline function directional_second(scheme, q, first, I, axis, h, lo, hi)
    if I[axis] - lo < 3
        return inward_first(first, I, axis, I[axis] - lo, 1) / h
    elseif hi - I[axis] < 3
        return inward_first(first, I, axis, hi - I[axis], -1) / h
    end
    return interior_curvature(scheme, q, first, I, axis, h) / h^2
end

@inline function interior_curvature(scheme, q, first, I, axis, h)
    m, p = shift(I, axis, -1), shift(I, axis, 1)
    return curvature(scheme, q[m], q[I], q[p], h * first[m], h * first[p])
end
