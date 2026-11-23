@inline face_cell(I::CartesianIndex{N}, axis, offset) where {N} =
    CartesianIndex(ntuple(d -> I[d] + ifelse(d == axis, offset, 0), Val(N)))

"""
    limit_face_pair(limiter, left, right, samples, I, axis)

Bound the two scalar candidates at the midpoint face between cell `I` and
its next neighbor along `axis`. Supply samples on a uniform logical axis.
Slope limiters read offsets -1 through +2 from `I`. MP5 reads offsets -2
through +3. Return the bounded pair without modifying the samples.
"""
@inline function limit_face_pair(limiter::SlopeLimiter, left, right, q, I, axis)
    J = face_cell(I, axis, 1)
    delta = q[J] - q[I]
    dl = limit_displacement(limiter, left - q[I], q[I] - q[face_cell(I, axis, -1)], delta)
    dr = limit_displacement(limiter, q[J] - right, delta, q[face_cell(I, axis, 2)] - q[J])
    return q[I] + dl, q[J] - dr
end

@inline function limit_face_pair(::MP5Limiter, left, right, q, I, axis)
    m2, m1, i, p1, p2, p3 = ntuple(k -> q[face_cell(I, axis, k - 3)], Val(6))
    dm1 = m2 - 2m1 + i
    di = m1 - 2i + p1
    dp1 = i - 2p1 + p2
    dp2 = p1 - 2p2 + p3
    return mp5_left_face_value(left, m2, m1, i, p1, dm1, di, dp1),
        mp5_left_face_value(right, p3, p2, p1, i, dp2, dp1, di)
end

@inline function reconstruct_edge(scheme::AbstractScheme, limiter::Union{SlopeLimiter, MP5Limiter}, q, I, axis, centers, faces, derivatives...)
    left, right = reconstruct_edge(scheme, q, I, axis, centers, faces, derivatives...)
    return limit_face_pair(limiter, left, right, q, I, axis)
end

@inline function reconstruct_edge!(left, right, scheme::AbstractScheme, limiter::Union{SlopeLimiter, MP5Limiter}, q, I, axis, centers, faces, derivatives...)
    left[I], right[I] = reconstruct_edge(scheme, limiter, q, I, axis, centers, faces, derivatives...)
    return nothing
end
