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

"""
Bound a candidate on the high face of a cell using MP5. Supply the candidate,
four consecutive cell samples from offsets -2 through +1, and ordinary second
differences from offsets -1 through +1. Reverse the stencil for the low face.
"""
@inline function mp5_left_face_value(
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
