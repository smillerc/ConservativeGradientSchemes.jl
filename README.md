# ConservativeGradientSchemes.jl

`ConservativeGradientSchemes` is designed to provide finite-difference gradients and scalar face reconstruction on structured grids. These grids must use a constant step along each logical axis.

- **Cell-center gradients.** Pass logical cell coordinates as `centers` and
  sample the field at the mapped physical cell centers. For a physical gradient,
  also supply the inverse Jacobian at the cell.
- **Geometry derivatives.** Use physical coordinates as samples and differentiate
  them along a logical axis. For example, differentiate the physical $x$
  coordinates along $\xi$ to estimate $\partial x / \partial \xi$.
- **Face reconstruction.** Reconstruct the left and right states at the logical
  midpoint between adjacent cell centers. On a stretched grid, the mapped face
  can lie away from the arithmetic midpoint of the physical cell centers.

The stencil uses one logical step per axis. On a stretched grid, pass the
logical coordinates in `centers` and keep the physical coordinates in the grid
mapping or geometry samples.

## Schemes

| Scheme | Cell-center gradients | Grid Metrics | Face reconstruction |
| --- | --- | --- | --- |
| `MEG{2,4,6,8}` | Yes | Yes | Yes |
| `Central{4,6,8}` | Yes | Yes | No |

## What the schemes do

Face reconstruction estimates the left and right values of a scalar field at a cell face. Use these values when you calculate a numerical flux.

Metric sweeps interpolate geometry samples to shared face values and compute
cell-center derivatives from their differences. Both schemes use one-sided
stencils at the edge of the valid geometry region when a centered stencil would
extend past it.

## Halo cells

Halo cells provide neighboring samples for stencils near a tile boundary. They are not part of the valid geometry region. When geometry halos are available, include them in the sample and coordinate arrays. At a physical or geometry boundary without halos, metric sweeps use one-sided derivatives. For MEG, prepare the metric derivatives with one-sided stencils at that boundary.

## Limiters

Use these limiters with MEG face reconstruction:

- `MinmodLimiter()`
- `SuperbeeLimiter()`
- `VanLeerLimiter()`
- `MP5Limiter()`

Pass a limiter after the scheme in `reconstruct_edge`. You can also call `limit_face_pair` when you already have the two face values.

## Example

```julia
using ConservativeGradientSchemes

scheme = MEG{6}()
samples = collect(1.0:12.0)
I = CartesianIndex(6)
left, right = reconstruct_edge(scheme, samples, I, 1)
limited = reconstruct_edge(scheme, MinmodLimiter(), samples, I, 1)

centers = (collect(1.0:12.0),)
gradient = curvilinear_cell_gradient(scheme, samples, I, centers)
```

`left` and `right` are the reconstructed values at the face between cells 6 and 7. Use `Central{4}()`, `Central{6}()`, or `Central{8}()` for centered differences.

## Jacobian and grid metrics

On a curvilinear grid, the mapping takes logical coordinates
$(\xi, \eta, \zeta)$ to Cartesian coordinates $(x, y, z)$. The forward
Jacobian matrix $\mathbf{J}$ has physical rows and logical columns. For
example, $x_\xi$ denotes $\partial x / \partial \xi$.

$$
\mathbf{J} =
\begin{bmatrix}
x_\xi & x_\eta & x_\zeta \\
y_\xi & y_\eta & y_\zeta \\
z_\xi & z_\eta & z_\zeta
\end{bmatrix}
$$

The Jacobian determinant $J = \det(\mathbf{J})$ gives the local volume scaling.

The inverse Jacobian matrix has logical rows and physical columns.

$$
\mathbf{J}^{-1} =
\begin{bmatrix}
\xi_x & \xi_y & \xi_z \\
\eta_x & \eta_y & \eta_z \\
\zeta_x & \zeta_y & \zeta_z
\end{bmatrix}
$$

Its determinant is $J^{-1} = \det(\mathbf{J}^{-1})$.

Pass the inverse matrix $\mathbf{J}^{-1}$ at cell `I` to
`physical_cell_gradient(scheme, samples, I, centers, inverse_jacobian)`.
The function computes the logical gradient of the scalar field $q$ and
applies the chain rule.

$$
\nabla_{\mathbf{x}} q = (\mathbf{J}^{-1})^{\mathsf{T}} \nabla_{\boldsymbol{\xi}} q
$$

For grid metric calculations, use `metric_cell_derivative!` to take
logical derivatives of geometry samples. For example, sweep the physical
$x$ coordinates along the logical $\xi$ axis to estimate
$\partial x / \partial \xi$. The sweep subtracts the low face value from
the high face value and divides by the logical spacing between the faces.
Adjacent cells share one value at their common face.

The metric sweep writes one derivative component at a time. The caller or grid
package assembles the Jacobian matrix, computes its determinant and inverse,
and constructs the grid metrics. Supply the inverse Jacobian when you call
`physical_cell_gradient`.
