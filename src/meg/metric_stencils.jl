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

"""Evaluate the centered MEG6 metric second stencil in unit logical spacing."""
@inline metric_centered_second(::MEG{6}, a, b, c, first_lo, first_hi) =
    curvature(Metric(), a, b, c, first_lo, first_hi)

"""
    metric_edge_state(scheme::MEG{6}, side, q, first, second)

Reconstruct the specified face of one cell from unit-logical-spacing metric
derivatives. Pass `Val(:low)` or `Val(:high)`. No adjacent cell is read.
"""
@inline metric_edge_state(::MEG{6}, ::Val{:low}, q, first, second, ::Val{C} = Val(Float64)) where {C} =
    q - (C(1 / 2) * first + C(1 / 12) * second)

@inline metric_edge_state(::MEG{6}, ::Val{:high}, q, first, second, ::Val{C} = Val(Float64)) where {C} =
    q + (C(1 / 2) * first + C(1 / 12) * second)
