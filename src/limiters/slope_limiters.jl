@inline minmod(a, b) = ifelse(a * b > zero(a), sign(a) * min(abs(a), abs(b)), zero(a))

@inline limited_slope(::MinmodLimiter, a, b) = minmod(a, b)

@inline function limited_slope(::SuperbeeLimiter, a, b)
    sign(a) == sign(b) || return zero(a)
    lo, hi = minmax(abs(a), abs(b))
    return copysign(min(2lo, hi), a)
end

@inline function limited_slope(::VanLeerLimiter, a, b)
    sign(a) == sign(b) || return zero(a)
    lo, hi = minmax(abs(a), abs(b))
    iszero(lo) && return zero(a)
    # Scale the harmonic mean to avoid multiplying the two slopes.
    return copysign(lo / ((1 + lo / hi) / 2), a)
end

"""
Limit a candidate displacement toward the high face of a cell. Use the two
adjacent forward differences to bound the magnitude by half the limited
slope. For the low face, pass the cell value minus the candidate face value.
"""
@inline function limit_displacement(limiter::SlopeLimiter, displacement, backward, forward)
    return minmod(displacement, limited_slope(limiter, backward, forward) / 2)
end
