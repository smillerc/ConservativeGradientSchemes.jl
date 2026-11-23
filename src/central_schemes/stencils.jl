@inline shift(I::CartesianIndex{N}, axis, offset) where {N} =
    CartesianIndex(ntuple(d -> I[d] + ifelse(d == axis, offset, 0), Val(N)))

@inline logical_spacing(::Type{T}, ::Nothing, axis) where {T} = one(T)
@inline function logical_spacing(::Type{T}, centers, axis) where {T}
    @inbounds coordinates = centers[axis]
    j = firstindex(coordinates)
    @inbounds return T(coordinates[j + 1] - coordinates[j])
end

@inline divisor(::Central{4}) = 12
@inline divisor(::Central{6}) = 60
@inline divisor(::Central{8}) = 840

@inline derivative_weights(::Central{4}) = (8, -1)
@inline derivative_weights(::Central{6}) = (45, -9, 1)
@inline derivative_weights(::Central{8}) = (672, -168, 32, -3)

# Differencing two consecutive interior faces recovers derivative_weights.
@inline face_weights(::Central{4}) = (-1, 7, 7, -1)
@inline face_weights(::Central{6}) = (1, -8, 37, 37, -8, 1)
@inline face_weights(::Central{8}) = (-3, 29, -139, 533, 533, -139, 29, -3)

# Rows are indexed by face distance from the low boundary, starting at zero.
# Solve the polynomial moment equations through degree Order using the
# interior face rule as the target functional. Reflect these rows at the high
# boundary. Matching the degree-Order moment also cancels its constant face
# error when taking differences across the transition to the interior.
@inline function boundary_weights(::Central{4}, row)
    rows = ((27, -31, 25, -11, 2), (2, 17, -11, 5, -1))
    @inbounds return rows[row]
end

@inline function boundary_weights(::Central{6}, row)
    rows = (
        (156, -267, 372, -343, 197, -64, 9),
        (9, 93, -78, 57, -28, 8, -1),
        (-1, 16, 72, -43, 22, -7, 1),
    )
    @inbounds return rows[row]
end

@inline function boundary_weights(::Central{8}, row)
    rows = (
        (2375, -5173, 9899, -13509, 12783, -8217, 3431, -841, 92),
        (92, 1547, -1861, 2171, -1917, 1191, -489, 119, -13),
        (-13, 209, 1079, -769, 533, -279, 99, -21, 2),
        (2, -31, 281, 911, -517, 281, -111, 27, -3),
    )
    @inbounds return rows[row]
end

@inline function centered_derivative(scheme, samples, I, axis, spacing)
    T = eltype(samples)
    result = zero(T)
    weights = derivative_weights(scheme)
    @inbounds for k in eachindex(weights)
        delta = samples[shift(I, axis, k)] - samples[shift(I, axis, -k)]
        result += T(weights[k]) * delta
    end
    return (result / T(divisor(scheme))) / spacing
end

@inline function weighted_face(scheme, samples, start, axis, direction, weights)
    T = eltype(samples)
    @inbounds reference = samples[start]
    result = zero(T)
    @inbounds for j in 2:length(weights)
        result += T(weights[j]) * (samples[shift(start, axis, direction * (j - 1))] - reference)
    end
    return reference + result / T(divisor(scheme))
end

@inline function metric_face(scheme::Central{Order}, samples, face, axis, low, high) where {Order}
    radius = Order ÷ 2
    j = face[axis]
    low_depth = j - low
    high_depth = high + 1 - j
    if low_depth < radius
        weights = boundary_weights(scheme, low_depth + 1)
        return weighted_face(scheme, samples, shift(face, axis, -low_depth), axis, 1, weights)
    elseif high_depth < radius
        weights = boundary_weights(scheme, high_depth + 1)
        return weighted_face(scheme, samples, shift(face, axis, high_depth - 1), axis, -1, weights)
    end
    return weighted_face(scheme, samples, shift(face, axis, -radius), axis, 1, face_weights(scheme))
end
