"""
Sixth-order MEG derivative preparation and scalar face reconstruction.
"""
module MEGSchemes

using KernelAbstractions: @index, @kernel
using ..Interfaces: AbstractScheme, AbstractLimiter
import ..Interfaces: reconstruct_edge, reconstruct_edge!, reconstruct_boundary_edge,
    prepare_first_derivative!, prepare_second_derivative!, prepare_metric_derivatives!, interpolate_metric_edges!,
    metric_cell_derivative!, curvilinear_cell_gradient, curvilinear_cell_gradient!,
    physical_cell_gradient, physical_cell_gradient!

"""
    MEG{6}()

Sixth-order MEG with operation-specific fluid and metric formulas. Supports
floating-point samples on uniform logical axes with midpoint faces. Point
operations require an interior stencil. Metric preparation requires at least
six valid cells along the sweep axis.
"""
struct MEG{Order} <: AbstractScheme
    MEG{6}() = new{6}()
end

"""
Limit the fluid MEG face displacement with the minmod slope of its parent cell.
"""
struct MinmodLimiter <: AbstractLimiter end

"""Apply the MP5 monotonicity-preserving bounds to the fluid MEG face pair."""
struct MP5Limiter <: AbstractLimiter end

const Fluid = Val{:fluid}
const Metric = Val{:metric}

include("derivatives.jl")
include("metric_stencils.jl")
include("reconstruction.jl")
include("sweeps.jl")

export MEG, MinmodLimiter, MP5Limiter
export metric_centered_first, metric_centered_second, metric_boundary_first, metric_edge_state

end
