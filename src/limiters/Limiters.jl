"""
Scalar slope limiters and monotonicity preserving face bounds.
"""
module Limiters

using ..Interfaces: AbstractLimiter, AbstractScheme
import ..Interfaces: reconstruct_edge, reconstruct_edge!

"""
Limit a reconstructed face displacement with the minmod slope of its parent cell.
"""
struct MinmodLimiter <: AbstractLimiter end

"""Limit a reconstructed face displacement with the Superbee slope of its parent cell."""
struct SuperbeeLimiter <: AbstractLimiter end

"""Limit a reconstructed face displacement with the Van Leer slope of its parent cell."""
struct VanLeerLimiter <: AbstractLimiter end

"""
Apply the MP5 monotonicity-preserving bounds to a reconstructed face value.

Suresh and Huynh describe the MP5 bounds in *Journal of Computational Physics*.
[DOI 10.1006/jcph.1997.5745](https://doi.org/10.1006/jcph.1997.5745).
"""
struct MP5Limiter <: AbstractLimiter end

"""Limiters that bound a face displacement with adjacent cell differences."""
const SlopeLimiter = Union{MinmodLimiter, SuperbeeLimiter, VanLeerLimiter}

include("slope_limiters.jl")
include("mp5.jl")
include("face_pairs.jl")

export MinmodLimiter, SuperbeeLimiter, VanLeerLimiter, MP5Limiter
export SlopeLimiter, limit_face_pair

end
