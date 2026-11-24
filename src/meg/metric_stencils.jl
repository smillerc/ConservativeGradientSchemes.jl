@inline metric_boundary_width(::MEG{Order}) where {Order} = Order + 1
@inline metric_boundary_width(::MEG{6}) = 6

# Centered coefficients multiply the paired differences in order of increasing offset.
@inline metric_centered_weights(::MEG{2}) = ((1,), 2)
@inline metric_centered_weights(::MEG{4}) = ((8, -1), 12)
@inline metric_centered_weights(::MEG{6}) = ((45, -9, 1), 60)
@inline metric_centered_weights(::MEG{8}) = ((672, -168, 32, -3), 840)

# Each row differentiates the Lagrange interpolant at its cell index.
# We store integer numerators and a common denominator for each row.
@inline function metric_inward_weights(::MEG{2}, depth)
    rows = (
        ((-3, 4, -1), 2),
    )
    return rows[depth + 1]
end

@inline function metric_inward_weights(::MEG{4}, depth)
    rows = (
        ((-25, 48, -36, 16, -3), 12),
        ((-3, -10, 18, -6, 1), 12),
    )
    return rows[depth + 1]
end

@inline function metric_inward_weights(::MEG{8}, depth)
    rows = (
        ((-2283, 6720, -11760, 15680, -14700, 9408, -3920, 960, -105), 840),
        ((-105, -1338, 2940, -2940, 2450, -1470, 588, -140, 15), 840),
        ((15, -240, -798, 1680, -1050, 560, -210, 48, -5), 840),
        ((-5, 60, -420, -378, 1050, -420, 140, -30, 3), 840),
    )
    return rows[depth + 1]
end

@inline function metric_inward_value(scheme, samples::NTuple{N}, depth, ::Val{C}, ϵ) where {N, C}
    weights, divisor = metric_inward_weights(scheme, depth)
    center = samples[depth + 1]
    terms = ntuple(k -> (C(weights[k]) / C(divisor)) * (samples[k] - center), Val(N))
    value = compensated_sum(terms)
    return value * (abs(value) >= ϵ)
end

"""
    metric_centered_first(scheme::MEG{Order}, samples::NTuple{Order}, ϵ, coefficient_type)

Evaluate the centered metric first derivative in unit logical spacing. Supply
samples in increasing index order and omit the center sample. The tolerance
defaults to machine epsilon of the sample type. Pass `Val(Float32)` or
`Val(Float64)` to select the coefficient type. The default is `Val(Float64)`.
"""
@inline function metric_centered_first(scheme::MEG{Order}, samples::NTuple{Order, T}, ϵ = eps(T), ::Val{C} = Val(Float64)) where {Order, T, C}
    weights, divisor = metric_centered_weights(scheme)
    radius = Order ÷ 2
    terms = ntuple(Val(radius)) do k
        (C(weights[k]) / C(divisor)) * difference(Metric(), samples[radius + k], samples[radius + 1 - k])
    end
    value = compensated_sum(terms)
    return value * (abs(value) >= ϵ)
end

"""
    metric_boundary_first(scheme::MEG, side, samples, coefficient_type; ϵ)

Return the first derivatives at the low or high end of the stencil in increasing
index order. Pass `Val(:low)` or `Val(:high)` and a tuple in increasing index
order. Supply `Order + 1` samples, or six samples for MEG6.
"""
@inline function metric_boundary_first(scheme::MEG{Order}, ::Val{:low}, samples::NTuple{N, T}, coeff::Val = Val(Float64); ϵ = 50 * eps(T)) where {Order, N, T}
    width = metric_boundary_width(scheme)
    N == width || throw(ArgumentError("MEG$Order requires $width boundary samples"))
    return ntuple(k -> metric_inward_value(scheme, samples, k - 1, coeff, ϵ), Val(Order ÷ 2))
end

@inline function metric_boundary_first(scheme::MEG, ::Val{:high}, samples::Tuple, coeff::Val = Val(Float64); kwargs...)
    inward = metric_boundary_first(scheme, Val(:low), reverse(samples), coeff; kwargs...)
    return map(-, reverse(inward))
end

@inline function metric_boundary_first(scheme::MEG{6}, side::Val{:low}, samples::NTuple{6, T}, coeff::Val = Val(Float64); ϵ::T = 50 * eps(T)) where {T}
    return metric_boundary_first(scheme, side, samples..., coeff; ϵ)
end

@inline function metric_boundary_first(scheme::MEG{6}, side::Val{:high}, samples::NTuple{6, T}, coeff::Val = Val(Float64); ϵ::T = 50 * eps(T)) where {T}
    return metric_boundary_first(scheme, side, samples..., coeff; ϵ)
end

@inline metric_upper_coefficient(::Val{Float64}, q) = (137 // 60) * q
@inline metric_upper_coefficient(::Val{C}, q) where {C} = C(137 / 60) * q

"""Evaluate the centered MEG6 metric first stencil in unit logical spacing."""
@inline function metric_centered_first(
        ::MEG{6},
        ϕᵢ₋₃::T, ϕᵢ₋₂::T, ϕᵢ₋₁::T, ϕᵢ₊₁::T, ϕᵢ₊₂::T, ϕᵢ₊₃::T, ϵ = eps(T), ::Val{C} = Val(Float64)
    ) where {T, C}
    ∂ϕ = compensated_sum(
        (
            C(3 / 4) * difference(Metric(), ϕᵢ₊₁, ϕᵢ₋₁),   #
            C(-3 / 20) * difference(Metric(), ϕᵢ₊₂, ϕᵢ₋₂), #
            C(1 / 60) * difference(Metric(), ϕᵢ₊₃, ϕᵢ₋₃),  #
        )
    )

    ∂ϕ = ∂ϕ * (abs(∂ϕ) >= ϵ)
    return ∂ϕ
end

"""Return the three first derivatives at the low end of a six-cell metric stencil."""
@inline function metric_boundary_first(
        ::MEG{6}, ::Val{:low},
        ϕ₁::T, ϕ₂::T, ϕ₃::T, ϕ₄::T, ϕ₅::T, ϕ₆::T, ::Val{C} = Val(Float64); ϵ::T = 50eps(T)
    ) where {T, C}
    # FiniteDifferenceMethod(0:5, 1)
    ∂ϕ₁ = compensated_sum(
        (
            C(-137 / 60) * ϕ₁,  #
            5ϕ₂,               #
            -5ϕ₃,              #
            C(10 / 3) * ϕ₄,     #
            -C(5 / 4) * ϕ₅,     #
            C(1 / 5) * ϕ₆,      #
        )
    )

    # FiniteDifferenceMethod(-1:4, 1)
    ∂ϕ₂ = compensated_sum(
        (
            C(-1 / 5) * ϕ₁,    #
            -C(13 / 12) * ϕ₂,  #
            2ϕ₃,              #
            -ϕ₄,              #
            C(1 / 3) * ϕ₅,     #
            -C(1 / 20) * ϕ₆,   #
        )
    )

    # FiniteDifferenceMethod(-2:3, 1)
    ∂ϕ₃ = compensated_sum(
        (
            C(1 / 20) * ϕ₁,  #
            -C(1 / 2) * ϕ₂,  #
            -C(1 / 3) * ϕ₃,  #
            ϕ₄,             #
            -C(1 / 4) * ϕ₅,  #
            C(1 / 30) * ϕ₆,  #
        )
    )

    return (
        ∂ϕ₁ * (abs(∂ϕ₁) >= ϵ), #
        ∂ϕ₂ * (abs(∂ϕ₂) >= ϵ), #
        ∂ϕ₃ * (abs(∂ϕ₃) >= ϵ), #
    )
end

"""Return the three first derivatives at the high end in increasing index order."""
@inline function metric_boundary_first(
        ::MEG{6}, ::Val{:high},
        ϕ₁::T, ϕ₂::T, ϕ₃::T, ϕ₄::T, ϕ₅::T, ϕ₆::T, ::Val{C} = Val(Float64); ϵ::T = 50eps(T)
    ) where {T, C}
    # FiniteDifferenceMethod(-5:0, 1)
    ∂ϕ₆ = compensated_sum(
        (
            C(-1 / 5) * ϕ₁,   #
            C(5 / 4) * ϕ₂,    #
            -C(10 / 3) * ϕ₃,  #
            5ϕ₄,             #
            -5ϕ₅,            #
            metric_upper_coefficient(Val(C), ϕ₆),  #
        )
    )

    # FiniteDifferenceMethod(-4:1, 1)
    ∂ϕ₅ = compensated_sum(
        (
            C(1 / 20) * ϕ₁,   #
            -C(1 / 3) * ϕ₂,   #
            ϕ₃,              #
            -2ϕ₄,            #
            C(13 / 12) * ϕ₅,  #
            C(1 / 5) * ϕ₆,    #
        )
    )

    # FiniteDifferenceMethod(-3:2, 1)
    ∂ϕ₄ = compensated_sum(
        (
            C(-1 / 30) * ϕ₁,  #
            C(1 / 4) * ϕ₂,    #
            -ϕ₃,             #
            C(1 / 3) * ϕ₄,    #
            C(1 / 2) * ϕ₅,    #
            -C(1 / 20) * ϕ₆,  #
        )
    )

    return (
        ∂ϕ₄ * (abs(∂ϕ₄) >= ϵ), #
        ∂ϕ₅ * (abs(∂ϕ₅) >= ϵ), #
        ∂ϕ₆ * (abs(∂ϕ₆) >= ϵ), #
    )
end

"""Evaluate the centered MEG metric second stencil in unit logical spacing."""
@inline metric_centered_second(::MEG, a, b, c, first_lo, first_hi) =
    curvature(Metric(), a, b, c, first_lo, first_hi)

"""
    metric_edge_state(scheme::MEG, side, q, first, second)

Return one cell's contribution to a conservative finite-difference metric face.
Supply a center point sample and its derivatives in unit logical spacing.
Use `Val(:low)` for `q - first/2 + second/12` and `Val(:high)` for
`q + first/2 + second/12`. The curvature term has the same sign on both faces.
The `1/12` coefficient corrects the face difference for use as a cell-center
derivative. Evaluate the mapping at the face when you need its point coordinates.
"""
@inline metric_edge_state(::MEG, ::Val{:low}, q, first, second, ::Val{C} = Val(Float64)) where {C} =
    q + (-C(1 / 2) * first + C(1 / 12) * second)

@inline metric_edge_state(::MEG, ::Val{:high}, q, first, second, ::Val{C} = Val(Float64)) where {C} =
    q + (C(1 / 2) * first + C(1 / 12) * second)
