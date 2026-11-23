@inline function snap(q::T, face) where {T}
    threshold = max(10 * eps(T), sqrt(eps(T)) * max(abs(q), abs(face)))
    return ifelse(abs(face - q) <= threshold, q, face)
end

@inline face_value(::Fluid, q, d1, d2, sign) = snap(q, (q + sign * d1 / 2) + d2 / 12)
@inline face_value(::Metric, q, d1, d2, sign) = q + sign * (d1 / 2 + d2 / 12)

@inline fluid_first(scheme, q, I, axis, h, ::Nothing) = centered_first(scheme, q, I, axis)
@inline fluid_first(scheme, q, I, axis, h, first) = h * first[I]

@inline function fluid_face_terms(scheme, q, I, axis, h, first, ::Nothing)
    m, j, p = shift(I, axis, -1), shift(I, axis, 1), shift(I, axis, 2)
    dm = fluid_first(scheme, q, m, axis, h, first)
    di = fluid_first(scheme, q, I, axis, h, first)
    dj = fluid_first(scheme, q, j, axis, h, first)
    dp = fluid_first(scheme, q, p, axis, h, first)
    return di, curvature(Fluid(), q[m], q[I], q[j], dm, dj),
        dj, curvature(Fluid(), q[I], q[j], q[p], di, dp)
end

@inline function fluid_face_terms(scheme, q, I, axis, h, first, second)
    J = shift(I, axis, 1)
    h2 = h * h
    return h * first[I], h2 * second[I], h * first[J], h2 * second[J]
end

@inline function reconstruct_edge(scheme::MEG, q, I, axis, centers, faces, first = nothing, second = nothing)
    J = shift(I, axis, 1)
    h = spacing(q[I], centers, axis)
    l1, l2, r1, r2 = fluid_face_terms(scheme, q, I, axis, h, first, second)
    return face_value(Fluid(), q[I], l1, l2, 1), face_value(Fluid(), q[J], r1, r2, -1)
end

# Limit the raw candidates before applying the existing MEG snapping rule.
@inline function reconstruct_edge(scheme::MEG, limiter::Union{SlopeLimiter, MP5Limiter}, q, I, axis, centers, faces, first = nothing, second = nothing)
    J = shift(I, axis, 1)
    h = spacing(q[I], centers, axis)
    l1, l2, r1, r2 = fluid_face_terms(scheme, q, I, axis, h, first, second)
    left, right = limit_face_pair(limiter, q[I] + l1 / 2 + l2 / 12, q[J] - r1 / 2 + r2 / 12, q, I, axis)
    return limiter isa MP5Limiter ? (left, right) : (snap(q[I], left), snap(q[J], right))
end

@inline reconstruct_edge(scheme::MEG, q, I, axis, first = nothing) =
    reconstruct_edge(scheme, q, I, axis, nothing, nothing, first)

@inline reconstruct_edge(scheme::MEG, limiter::Union{SlopeLimiter, MP5Limiter}, q, I, axis, first = nothing) =
    reconstruct_edge(scheme, limiter, q, I, axis, nothing, nothing, first)

@inline function reconstruct_edge!(left, right, scheme::MEG, q, I, axis, centers, faces, first = nothing, second = nothing)
    left[I], right[I] = reconstruct_edge(scheme, q, I, axis, centers, faces, first, second)
    return nothing
end

@inline function reconstruct_edge!(left, right, scheme::MEG, limiter::Union{SlopeLimiter, MP5Limiter}, q, I, axis, centers, faces, first = nothing, second = nothing)
    left[I], right[I] = reconstruct_edge(scheme, limiter, q, I, axis, centers, faces, first, second)
    return nothing
end

@inline function reconstruct_boundary_edge(scheme::MEG, q, I, axis, centers, faces, side::Union{Val{:low}, Val{:high}}, first, second)
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

@inline function curvilinear_cell_gradient(scheme::MEG, q, I::CartesianIndex{N}, centers) where {N}
    return ntuple(axis -> centered_first(scheme, q, I, axis) / spacing(q[I], centers, axis), Val(N))
end

@inline function physical_cell_gradient(scheme::MEG, q, I::CartesianIndex{N}, centers, inverse_jacobian) where {N}
    logical = curvilinear_cell_gradient(scheme, q, I, centers)
    return ntuple(Val(N)) do b
        sum(ntuple(a -> inverse_jacobian[a, b] * logical[a], Val(N)))
    end
end

@inline function curvilinear_cell_gradient!(dest, scheme::MEG, q, I, centers)
    values = curvilinear_cell_gradient(scheme, q, I, centers)
    for d in eachindex(values)
        dest[d][I] = values[d]
    end
    return nothing
end

@inline function physical_cell_gradient!(dest, scheme::MEG, q, I, centers, inverse_jacobian)
    values = physical_cell_gradient(scheme, q, I, centers, inverse_jacobian)
    for d in eachindex(values)
        dest[d][I] = values[d]
    end
    return nothing
end
