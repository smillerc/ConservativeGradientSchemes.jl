"""
MEG derivative preparation and scalar face reconstruction.
"""
module MEGSchemes

using KernelAbstractions: @index, @kernel
using ..Interfaces: AbstractScheme
using ..Limiters: MinmodLimiter, MP5Limiter, SlopeLimiter, limit_face_pair
import ..Interfaces: reconstruct_edge, reconstruct_edge!, reconstruct_boundary_edge,
    prepare_first_derivative!, prepare_second_derivative!, prepare_metric_derivatives!, interpolate_metric_edges!,
    metric_cell_derivative!, curvilinear_cell_gradient, curvilinear_cell_gradient!,
    physical_cell_gradient, physical_cell_gradient!

"""
    MEG{Order}()

MEG with a centered first derivative of order 2, 4, 6, or 8 for fluid
reconstruction, gradients, and metric preparation. Supports floating-point samples on uniform
logical axes with midpoint faces. Point gradients need `Order ÷ 2` samples
on each side. For metric preparation, each valid line needs at least
`Order + 1` cells, with six cells sufficient for the MEG6 closure.

Chamarthi describes MEG6 reconstruction and derivative reuse in *Computers &
Fluids*. [DOI 10.1016/j.compfluid.2022.105706](https://doi.org/10.1016/j.compfluid.2022.105706).

Chandravamsi, Chamarthi, Hoffmann, and Frankel describe MEG6 and metric
interpolation on curvilinear meshes in *Computers & Fluids*.
[DOI 10.1016/j.compfluid.2023.105859](https://doi.org/10.1016/j.compfluid.2023.105859).
"""
struct MEG{Order} <: AbstractScheme
    MEG{2}() = new{2}()
    MEG{4}() = new{4}()
    MEG{6}() = new{6}()
    MEG{8}() = new{8}()
end

const Fluid = Val{:fluid}
const Metric = Val{:metric}

include("derivatives.jl")
include("metric_stencils.jl")
include("reconstruction.jl")
include("sweeps.jl")

export MEG, MinmodLimiter, MP5Limiter
export metric_centered_first, metric_centered_second, metric_boundary_first, metric_edge_state

end
