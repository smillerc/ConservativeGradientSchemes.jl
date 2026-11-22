# ConservativeGradientSchemes.jl

Scalar face reconstruction and cell-center gradients on structured grids,
including curvilinear grids with nonuniform physical spacing.

## Install

Requires Julia 1.13 and KernelAbstractions 0.10.4 or later in the 0.10 series.
Once registered, install with `import Pkg; Pkg.add("ConservativeGradientSchemes")`.

## Example

```julia
using ConservativeGradientSchemes

samples = collect(1.0:12.0)
left, right = reconstruct_edge(MEG{6}(), samples, CartesianIndex(6), 1)
```

This reconstructs the two states at the face between cells 6 and 7.
Pass a limiter after the scheme to limit the reconstruction.

## Tests

From the package root:

```sh
julia --project -e 'using Pkg; Pkg.test(; julia_args=["--check-bounds=auto"])'
```

For the optional 3D Metal test, run this in Julia from the same directory:

```julia
using Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=pwd())
Pkg.add(["KernelAbstractions", "Metal"])
include("test/metal.jl")
```
