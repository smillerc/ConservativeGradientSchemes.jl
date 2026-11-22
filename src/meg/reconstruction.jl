@inline function snap(q::T, face) where {T}
    threshold = max(10 * eps(T), sqrt(eps(T)) * max(abs(q), abs(face)))
    return ifelse(abs(face - q) <= threshold, q, face)
end

@inline face_value(::Fluid, q, d1, d2, sign) = snap(q, (q + sign * d1 / 2) + d2 / 12)
@inline face_value(::Metric, q, d1, d2, sign) = q + sign * (d1 / 2 + d2 / 12)

@inline fluid_first(q, I, axis, h, ::Nothing) = centered_first(Fluid(), q, I, axis)
@inline fluid_first(q, I, axis, h, first) = h * first[I]

@inline function fluid_face_terms(q, I, axis, h, first, ::Nothing)
    m, j, p = shift(I, axis, -1), shift(I, axis, 1), shift(I, axis, 2)
    dm = fluid_first(q, m, axis, h, first)
    di = fluid_first(q, I, axis, h, first)
    dj = fluid_first(q, j, axis, h, first)
    dp = fluid_first(q, p, axis, h, first)
    return di, curvature(Fluid(), q[m], q[I], q[j], dm, dj),
        dj, curvature(Fluid(), q[I], q[j], q[p], di, dp)
end

@inline function fluid_face_terms(q, I, axis, h, first, second)
    J = shift(I, axis, 1)
    h2 = h * h
    return h * first[I], h2 * second[I], h * first[J], h2 * second[J]
end

@inline minmod(a, b) = ifelse(a * b > zero(a), sign(a) * min(abs(a), abs(b)), zero(a))

@inline function reconstruct_edge(::MEG{6}, q, I, axis, centers, faces, first = nothing, second = nothing)
    J = shift(I, axis, 1)
    h = spacing(q[I], centers, axis)
    l1, l2, r1, r2 = fluid_face_terms(q, I, axis, h, first, second)
    return face_value(Fluid(), q[I], l1, l2, 1), face_value(Fluid(), q[J], r1, r2, -1)
end

@inline function reconstruct_edge(scheme::MEG{6}, limiter::MinmodLimiter, q, I, axis, centers, faces, first = nothing, second = nothing)
    J = shift(I, axis, 1)
    h = spacing(q[I], centers, axis)
    l1, l2, r1, r2 = fluid_face_terms(q, I, axis, h, first, second)
    delta = q[J] - q[I]
    ls = minmod(q[I] - q[shift(I, axis, -1)], delta)
    rs = minmod(delta, q[shift(J, axis, 1)] - q[J])
    left = q[I] + l1 / 2 + l2 / 12
    right = q[J] - r1 / 2 + r2 / 12
    return snap(q[I], q[I] + minmod(left - q[I], ls / 2)),
        snap(q[J], q[J] - minmod(q[J] - right, rs / 2))
end

# Compute MP5 bounds with ordinary second differences.
@inline function _minmod(a::T, b::T) where {T}
    half = T(1) / T(2)
    return half * (sign(a) + sign(b)) * min(abs(a), abs(b))
end

@inline function _minmod(a::T, b::T, c::T, d::T) where {T}
    eighth = T(1) / T(8)
    return eighth *
        (sign(a) + sign(b)) *
        abs((sign(a) + sign(c)) * (sign(a) + sign(d))) *
        min(abs(a), abs(b), abs(c), abs(d))
end

@inline function _mp5_left_face_value(
        ϕLᵢ₊½::T, ϕᵢ₋₂::T, ϕᵢ₋₁::T, ϕᵢ::T, ϕᵢ₊₁::T, ∂²ϕᵢ₋₁::T, ∂²ϕᵢ::T, ∂²ϕᵢ₊₁::T
    ) where {T}
    ϵ = T <: AbstractFloat ? eps(one(T))^4 : zero(T)
    two = T(2)
    three = T(3)
    four = T(4)
    α = four
    β = four / three

    ϕ_mp = ϕᵢ + _minmod(ϕᵢ₊₁ - ϕᵢ, α * (ϕᵢ - ϕᵢ₋₁))
    activate_limiter = (ϕLᵢ₊½ - ϕᵢ) * (ϕLᵢ₊½ - ϕ_mp) >= ϵ

    dᵢ₊½ = _minmod(four * ∂²ϕᵢ - ∂²ϕᵢ₊₁, four * ∂²ϕᵢ₊₁ - ∂²ϕᵢ, ∂²ϕᵢ, ∂²ϕᵢ₊₁)
    dᵢ₋½ = _minmod(four * ∂²ϕᵢ - ∂²ϕᵢ₋₁, four * ∂²ϕᵢ₋₁ - ∂²ϕᵢ, ∂²ϕᵢ, ∂²ϕᵢ₋₁)

    ϕ_md = (ϕᵢ + ϕᵢ₊₁ - dᵢ₊½) / two
    ϕ_ul = ϕᵢ + α * (ϕᵢ - ϕᵢ₋₁)
    ϕ_lc = (three * ϕᵢ - ϕᵢ₋₁) / two + β * dᵢ₋½

    ϕ_min = max(min(ϕᵢ, ϕᵢ₊₁, ϕ_md), min(ϕᵢ, ϕ_ul, ϕ_lc))
    ϕ_max = min(max(ϕᵢ, ϕᵢ₊₁, ϕ_md), max(ϕᵢ, ϕ_ul, ϕ_lc))
    candidate = ϕLᵢ₊½ + _minmod(ϕ_min - ϕLᵢ₊½, ϕ_max - ϕLᵢ₊½)
    return ifelse(activate_limiter, candidate, ϕLᵢ₊½)
end

@inline function reconstruct_edge(::MEG{6}, ::MP5Limiter, q, I, axis, centers, faces, first = nothing, second = nothing)
    m2, m1, i, p1, p2, p3 = ntuple(k -> q[shift(I, axis, k - 3)], Val(6))
    h = spacing(i, centers, axis)
    l1, l2, r1, r2 = fluid_face_terms(q, I, axis, h, first, second)
    left = i + l1 / 2 + l2 / 12
    right = p1 - r1 / 2 + r2 / 12
    dm1 = m2 - 2m1 + i
    di = m1 - 2i + p1
    dp1 = i - 2p1 + p2
    dp2 = p1 - 2p2 + p3
    return _mp5_left_face_value(left, m2, m1, i, p1, dm1, di, dp1),
        _mp5_left_face_value(right, p3, p2, p1, i, dp2, dp1, di)
end

@inline reconstruct_edge(scheme::MEG{6}, q, I, axis, first = nothing) =
    reconstruct_edge(scheme, q, I, axis, nothing, nothing, first)

@inline reconstruct_edge(scheme::MEG{6}, limiter::Union{MinmodLimiter, MP5Limiter}, q, I, axis, first = nothing) =
    reconstruct_edge(scheme, limiter, q, I, axis, nothing, nothing, first)

@inline function reconstruct_edge!(left, right, scheme::MEG{6}, q, I, axis, centers, faces, first = nothing, second = nothing)
    left[I], right[I] = reconstruct_edge(scheme, q, I, axis, centers, faces, first, second)
    return nothing
end

@inline function reconstruct_edge!(left, right, scheme::MEG{6}, limiter::Union{MinmodLimiter, MP5Limiter}, q, I, axis, centers, faces, first = nothing, second = nothing)
    left[I], right[I] = reconstruct_edge(scheme, limiter, q, I, axis, centers, faces, first, second)
    return nothing
end

@inline function reconstruct_boundary_edge(scheme::MEG{6}, q, I, axis, centers, faces, side::Union{Val{:low}, Val{:high}}, first, second)
    q_i = q[I]
    sign = side isa Val{:low} ? -1 : 1
    j = I[axis]
    face = j + (side isa Val{:high})
    h = typeof(q_i)(2 * sign * (faces[axis][face] - centers[axis][j]))
    return face_value(Metric(), q_i, h * first[I], h^2 * second[I], sign)
end

@inline function interior_face(q, I, axis, centers, first, second)
    J = shift(I, axis, 1)
    h = spacing(q[I], centers, axis)
    h2 = h * h
    left = face_value(Metric(), q[I], h * first[I], h2 * second[I], 1)
    right = face_value(Metric(), q[J], h * first[J], h2 * second[J], -1)
    return (left + right) / 2
end

@inline function metric_difference(high, low)
    delta = high - low
    return delta * (abs(delta) >= 50 * eps(typeof(delta)))
end

@inline function curvilinear_cell_gradient(scheme::MEG{6}, q, I::CartesianIndex{N}, centers) where {N}
    return ntuple(axis -> centered_first(Fluid(), q, I, axis) / spacing(q[I], centers, axis), Val(N))
end

@inline function physical_cell_gradient(scheme::MEG{6}, q, I::CartesianIndex{N}, centers, inverse_jacobian) where {N}
    logical = curvilinear_cell_gradient(scheme, q, I, centers)
    return ntuple(Val(N)) do b
        sum(ntuple(a -> inverse_jacobian[a, b] * logical[a], Val(N)))
    end
end

@inline function curvilinear_cell_gradient!(dest, scheme::MEG{6}, q, I, centers)
    values = curvilinear_cell_gradient(scheme, q, I, centers)
    for d in eachindex(values)
        dest[d][I] = values[d]
    end
    return nothing
end

@inline function physical_cell_gradient!(dest, scheme::MEG{6}, q, I, centers, inverse_jacobian)
    values = physical_cell_gradient(scheme, q, I, centers, inverse_jacobian)
    for d in eachindex(values)
        dest[d][I] = values[d]
    end
    return nothing
end
