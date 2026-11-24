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

MEG for floating-point samples on uniform logical axes with midpoint faces.
`Order` is the order of the centered first-derivative stencil. It needs
`Order ÷ 2` samples on each side. Metric preparation needs at least
`Order + 1` cells on each valid line, or six cells for MEG6.

- **Advection.** With center point samples, use `reconstruct_edge` to build
  face values for conservative finite differences. For constant-speed
  advection, reconstruct the upwind state. For variable-speed conservative
  advection, reconstruct the point flux or its split components. The `1/12`
  curvature coefficient corrects the face difference for a center derivative.
  The same coefficient reconstructs a quadratic from uniform one-dimensional
  cell averages. On a mapped finite-volume grid, account for the averaging
  measure and face quadrature.
- **Grid metrics.** Supply center point samples to the metric sweeps for a
  conservative finite-difference derivative with matching MEG flux differencing.
  For physical point coordinates at a face, evaluate the mapping there. In
  Taylor interpolation from point samples with exact derivatives, use `1/8`
  for the quadratic term at the midpoint.
- **Cell-centered gradients.** Use `curvilinear_cell_gradient` for direct
  logical derivatives of center point samples. Use `physical_cell_gradient`
  to transform those derivatives with the inverse Jacobian at the same cell.
  Both operations use the centered first-derivative stencil directly.

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
