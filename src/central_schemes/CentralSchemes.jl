"""
Centered scalar gradients and conservative metric differentiation.
"""
module CentralSchemes

using KernelAbstractions: @index, @kernel
using ..Interfaces: AbstractScheme
import ..Interfaces: curvilinear_cell_gradient, curvilinear_cell_gradient!,
    physical_cell_gradient, physical_cell_gradient!, prepare_first_derivative!,
    interpolate_metric_edges!, metric_cell_derivative!

"""
    Central{Order}()

Centered finite differences of order 4, 6, or 8 on uniformly spaced logical
axes. Point gradients require an interior stencil. Metric sweeps work on
rectangular index domains and use inward face stencils at both ends.
Each valid line must contain at least `Order + 1` cells.
"""
struct Central{Order} <: AbstractScheme
    Central{4}() = new{4}()
    Central{6}() = new{6}()
    Central{8}() = new{8}()
end

include("stencils.jl")
include("gradients.jl")
include("metric_sweeps.jl")

export Central

end
