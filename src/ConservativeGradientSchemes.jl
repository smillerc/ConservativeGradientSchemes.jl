"""
Interfaces for scalar edge reconstruction, limiting, and cell-center gradients.
"""
module ConservativeGradientSchemes

include("interfaces/Interfaces.jl")
include("limiters/Limiters.jl")
include("meg/MEGSchemes.jl")
include("central_schemes/CentralSchemes.jl")
include("weno/WENO.jl")
include("ppm/PPM.jl")
include("bvd/BVD.jl")
include("mlp/MLP.jl")

using .CentralSchemes: Central
using .Limiters: MinmodLimiter, SuperbeeLimiter, VanLeerLimiter, MP5Limiter, limit_face_pair
using .MEGSchemes:
    MEG, metric_centered_first, metric_centered_second,
    metric_boundary_first, metric_edge_state
using .Interfaces:
    AbstractScheme, AbstractLimiter, reconstruct_edge, reconstruct_edge!,
    reconstruct_boundary_edge, prepare_first_derivative!, prepare_second_derivative!, prepare_metric_derivatives!,
    interpolate_metric_edges!, metric_cell_derivative!, curvilinear_cell_gradient, curvilinear_cell_gradient!, physical_cell_gradient, physical_cell_gradient!

export Central
export MEG
export MinmodLimiter, SuperbeeLimiter, VanLeerLimiter, MP5Limiter, limit_face_pair
export metric_centered_first, metric_centered_second, metric_boundary_first, metric_edge_state
export AbstractScheme, AbstractLimiter, reconstruct_edge, reconstruct_edge!
export reconstruct_boundary_edge, prepare_first_derivative!, prepare_second_derivative!, prepare_metric_derivatives!
export interpolate_metric_edges!, metric_cell_derivative!
export curvilinear_cell_gradient, curvilinear_cell_gradient!, physical_cell_gradient, physical_cell_gradient!

end
