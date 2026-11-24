"""
Public contracts for scalar reconstruction, limiting, and gradients.

Point operations must support calls from device kernels with device-compatible
arguments. They must not launch kernels, allocate heap storage, or access host
data on the device path. Sweeps launch through the supplied backend. Their
implementations must use caller-owned backend-compatible buffers and must not
use host scalar indexing or stage field data through the host. Callers must
order dependent launches before consuming the output.
"""
module Interfaces

"""
Supertype for numerical schemes. Implement the operations each scheme supports.
Document the supported sample representation
(point values or cell averages), stencil footprint, and coordinate assumptions
for each concrete scheme.
"""
abstract type AbstractScheme end

"""
Supertype for optional reconstruction limiters. Document the required stencil
footprint, including transverse samples, and the bounds each limiter enforces.
"""
abstract type AbstractLimiter end

"""
    reconstruct_edge(scheme::AbstractScheme, samples, I, axis, centers, faces)
    reconstruct_edge(scheme::AbstractScheme, limiter::AbstractLimiter,
                     samples, I, axis, centers, faces)

    reconstruct_edge(scheme::AbstractScheme, samples, I, axis, centers, faces, first)
    reconstruct_edge(scheme::AbstractScheme, samples, I, axis, centers, faces, first, second)

Return the scalar pair `(left, right)` at the interface between cell `I` and
its neighbor in the positive `axis` direction. The left state comes from the
lower-index side. Both states have the same units as the samples.

Supply cell-indexed scalar samples and tuples of one-dimensional logical
coordinate arrays. On each axis, `centers[d][j]` locates cell `j` and
`faces[d][j + 1]` locates its high face. Supply valid samples and coordinates
over the entire stencil, including halos where needed. The call reads its inputs
without mutating them.

The first form applies the reconstruction scheme as defined, including any
intrinsic weighting or monotonicity rules. With a built-in limiter, the default
overload computes that pair and calls `limit_face_pair`. Use it for midpoint
faces on a uniform logical axis. A scheme can specialize the overload when
it needs to limit intermediate face candidates.

With MEG and center point samples, use the returned states to form a
conservative finite-difference flux. The `1/12` curvature coefficient corrects
the flux difference for a cell-center derivative. For midpoint point values,
use a point interpolant.

`first` and `second` are read-only cell-indexed arrays for the same samples and
axis, with the units of the first and second logical derivatives, respectively.
With only `first`, a scheme can form its second derivatives locally. With both,
it uses the supplied derivatives. The caller must prepare every entry read by
the selected stencil before reconstruction. MEG also supports unit logical
spacing via `reconstruct_edge(scheme, samples, I, axis[, first])` and its
explicit limiter form.
"""
function reconstruct_edge end

"""
    reconstruct_edge!(left, right, scheme::AbstractScheme, samples, I, axis, centers, faces)
    reconstruct_edge!(left, right, scheme::AbstractScheme, limiter::AbstractLimiter,
                      samples, I, axis, centers, faces)

Write the two interface values from `reconstruct_edge` to `left[I]` and
`right[I]`. Index the caller-owned output arrays by the cell on the lower-index
side of the interface. Each call writes one interface along `axis` and returns
`nothing`. Append `first` or `first, second` for the prepared overloads, in
either limiter form.

Use the sample, coordinate, and stencil contract of `reconstruct_edge`.
The output arrays must have writable entries at `I` and must be compatible with
the selected backend. They must not alias each other or the read-only inputs.
"""
function reconstruct_edge! end

"""
    reconstruct_boundary_edge(scheme::AbstractScheme, samples, I, axis, centers, faces,
                              side, first, second)

Return one scalar face value from the in-domain boundary cell `I`. For
`Val(:low)`, use the low face of `I` at index `I[axis]` in `faces[axis]`.
For `Val(:high)`, use the high face at index `I[axis] + 1`.
Read only the sample and prepared directional derivatives at `I` and the
coordinates for that cell and face. No exterior geometry is required.
"""
function reconstruct_boundary_edge end

"""
    prepare_first_derivative!(first, scheme::AbstractScheme, samples, axis, centers,
                              domain, valid_cells, backend)
    prepare_second_derivative!(second, scheme::AbstractScheme, samples, first, axis, centers,
                               domain, valid_cells, backend)

Write one directional derivative at each index in `domain` and return
`nothing`. Write first derivatives in field units per logical coordinate and
second derivatives in field units per logical coordinate squared. The second
sweep borrows the prepared first derivative at each location its stencil reads.

`valid_cells` identifies the available stencil region, independently of the
output domain. Supply valid halo samples where the centered stencil requires
them. All arrays must use compatible indexing and storage on `backend`.
The caller owns the output and must keep it disjoint from the inputs.

For `Central`, use the same derivative and inward closures as
`metric_cell_derivative!`. Use `centers = nothing` for unit logical spacing.
"""
function prepare_first_derivative! end

"""
    prepare_second_derivative!(second, scheme::AbstractScheme, samples, first, axis, centers,
                               domain, valid_cells, backend)

Write the second directional derivative using the prepared first derivative.
See [`prepare_first_derivative!`](@ref) for units, halo coverage, and storage.
"""
function prepare_second_derivative! end

"""
    prepare_metric_derivatives!(first, second, scheme::AbstractScheme, samples, axis,
                                centers, first_domain, domain, valid_cells, backend)

Prepare directional derivatives for metric reconstruction. Write first
derivatives on `first_domain` and second derivatives on `domain`. Both output
domains are subsets of the rectangular `valid_cells` region. Each line must
have enough valid entries for the scheme's inward stencils. A tile edge inside
the valid region does not trigger a closure.

The caller supplies disjoint, backend-compatible output buffers. Each entry
read by second derivative preparation must be ready after the first sweep.
An empty output domain skips its sweep. Return `nothing`.
"""
function prepare_metric_derivatives! end

"""
    interpolate_metric_edges!(dest, scheme::Central, samples, axis,
                              face_domain, valid_cells, backend)
    interpolate_metric_edges!(dest, scheme::AbstractScheme, samples, first, second,
                              axis, centers, faces, face_domain, valid_cells, backend)

Write single-valued metric faces in `face_domain` and return `nothing`.
For `Central`, evaluate each face directly from the samples. Its face weights
produce the centered derivative after differencing in the interior. Use the
inward closure where the interior stencil would cross `valid_cells`.

For MEG, supply center point samples. Use the arithmetic mean of the two
reconstructions on interior faces. On a low or high boundary of `valid_cells`,
call `reconstruct_boundary_edge` with the adjacent in-domain cell. Read no
exterior geometry or derivative values. Use these face values for conservative
finite differences. Evaluate the mapping at the face for physical point coordinates.

Index `dest` by the face coordinate index along `axis` and by cell indices
in transverse directions. The low and high faces of cell `j` along `axis`
are at face indices `j` and `j + 1`. The caller supplies writable storage
for all requested faces, including boundary faces. The storage must be
compatible with `backend` and disjoint from the inputs.
"""
function interpolate_metric_edges! end

"""
    metric_cell_derivative!(dest, scheme::Central, samples, axis, centers,
                            domain, valid_cells, backend)
    metric_cell_derivative!(dest, scheme::AbstractScheme, samples, first, second,
                                 axis, centers, faces, domain, valid_cells, backend)

Write one logical derivative component at each cell in `domain` and return
`nothing`. For `Central`, difference the same face values produced by
`interpolate_metric_edges!` and divide by the constant logical spacing along
`axis`. Supply one coordinate array per logical axis, or `nothing` for unit
logical spacing. The valid domain must contain every sample read by the sweep.

For the MEG metric operator, use the arithmetic mean of the left
and right reconstructions as the single value of each interior face. At the
low or high end of `valid_cells`, use the corresponding in-domain boundary
face value from `reconstruct_boundary_edge`. For each cell, subtract the low
face value from the high face value and divide by its logical face spacing.

The caller supplies read-only samples on the stencil region. For MEG, also
supply prepared derivatives with one-sided preparation where geometry halos
are unavailable.
The sweep needs no geometry beyond those boundary faces. It can evaluate face
values as needed without allocating face arrays. `dest` is caller-owned,
indexed by cell, compatible with `backend`, and disjoint from the inputs.
"""
function metric_cell_derivative! end

"""
    curvilinear_cell_gradient(scheme::AbstractScheme, samples, I, centers)

Return the cell-center derivatives of the scalar field with respect to the
logical coordinates, in axis order. Each component has units of the scalar
field per unit of its logical coordinate.

Supply cell-indexed scalar samples and a tuple of one-dimensional logical
coordinate arrays with `centers[d][j]` locating cell `j` on axis `d`.
Supply valid samples and coordinates over the entire stencil, including halos
where needed. `Central` also accepts `centers = nothing` for unit logical
spacing. The call reads its inputs without mutating them.
"""
function curvilinear_cell_gradient end

"""
    curvilinear_cell_gradient!(dest, scheme::AbstractScheme, samples, I, centers)

Write the logical gradient at cell `I` to caller-owned component arrays in
`dest`. Write component `d` to `dest[d][I]` in logical axis order. Return
`nothing`.

Use the sample, coordinate, and stencil contract of `curvilinear_cell_gradient`.
Supply one writable array per component with an entry at `I` and storage
compatible with the selected backend. The arrays must not alias each other or
the read-only inputs.
"""
function curvilinear_cell_gradient! end

"""
    physical_cell_gradient(scheme::AbstractScheme, samples, I, centers, inverse_jacobian)

Return the scalar gradient at cell `I` in Cartesian physical coordinates and in
field units per physical length. Use the same samples, coordinates, and
stencil contract as [`curvilinear_cell_gradient`](@ref).

Supply the inverse Jacobian at `I` with logical rows and physical columns,
such that entry `(a, b)` is the derivative of logical coordinate `a` with
respect to physical coordinate `b`. The physical gradient is the transpose
of this matrix times the logical gradient, with both gradients treated as
column vectors. The call reads its inputs without mutating them.
"""
function physical_cell_gradient end

"""
    physical_cell_gradient!(dest, scheme::AbstractScheme, samples, I, centers, inverse_jacobian)

Write the physical gradient at cell `I` to caller-owned component arrays in
`dest`. Write component `d` to `dest[d][I]` in Cartesian physical axis order.
Return `nothing`.

Use the sample, coordinate, stencil, and inverse Jacobian contract of
`physical_cell_gradient` and the output storage contract of
`curvilinear_cell_gradient!`.
"""
function physical_cell_gradient! end

export AbstractScheme, AbstractLimiter, reconstruct_edge, reconstruct_edge!
export reconstruct_boundary_edge, prepare_first_derivative!, prepare_second_derivative!
export prepare_metric_derivatives!, interpolate_metric_edges!, metric_cell_derivative!
export curvilinear_cell_gradient, curvilinear_cell_gradient!, physical_cell_gradient, physical_cell_gradient!

end
