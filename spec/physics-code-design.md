# Physics-Code Design

<!-- BEGIN COMMON PHYSICS-CODE DESIGN revision=3 -->

## 1. Purpose, authority, and terminology

This standard governs the architecture and engineering of Julia physics packages.
It is written for domain scientists, developers, and coding agents. Its priorities
are correct physics, measured performance, and implementation that a scientist can
review directly. It does not select equations, discretizations, solver defaults,
or supported execution modes.

**PCD-GOV-001 — Authority.** This document MUST govern code organization and
engineering conventions. Package physics RFCs MUST retain authority over equations
and numerical methods; package acceptance specifications MUST retain their oracles,
tolerances, test tiers, and capability qualifications. Conflicts MUST be identified
and resolved explicitly, never by silently changing a numerical contract. An
architecture target or a passing structural check MUST NOT grant physics/backend
support or establish that a migration is complete.

**PCD-GOV-002 — Requirement language.** MUST and MUST NOT identify requirements;
SHOULD and SHOULD NOT identify recommendations whose departure needs a recorded
reason; MAY identifies a permitted choice. Requirement IDs are stable references.
The package appendix MUST distinguish observed implementation, proposed targets,
and known nonconformance, with an inspection date and supporting evidence.

A **theme** is a coherent scientific or algorithmic responsibility. A **module
entry** declares its Julia module and maps its vocabulary, dependencies, and files.
An **implementation file** defines methods or private types inside that namespace.
A **separable module** loads through one entry, declares its dependencies, owns its
files, and can be tested through named interfaces without accidental parent names
or unrelated include order. Separability does not require a separate Julia package.

**PCD-GOV-003 — Matching copies.** The text between the common-body markers MUST be
identical in all maintained copies and carry the same design revision. Changes to
common requirements MUST update all maintained copies together. Package names,
local links, ownership registries, and migration records MUST remain outside these markers.
Local appendices MUST NOT silently weaken the common rules. Existing violations
are migration debt, not exceptions established by current code.

## 2. Scientific integrity and readable mathematics

**PCD-SCI-001 — Numerical contracts.** Changes MUST preserve the applicable
conservation laws, flux signs, units, variable centering, basis conventions,
boundary treatment, nullspaces, symmetry, admissibility, residual definitions,
and convergence properties unless a separate numerical design deliberately changes
them. Record the applicable discrete formula or invariant before a nontrivial
numerical change. Do not invent physics to complete an implementation.

**PCD-SCI-002 — Visible formulas.** Implementations MUST keep the governing operation,
index/halo assumptions, and geometry or basis transformations reviewable. Keep
geometry/metric access, local arithmetic, and writeback distinguishable. Prefer
direct equations and loops over deep helper chains, metaprogramming, or generic
abstractions that obscure the method. Reuse the owning geometry/model interface
instead of maintaining a competing formula or cache of the same physical data.

**PCD-SCI-003 — Justified changes.** Repairs, clipping, regularization, floors,
fallbacks, altered tolerances, and changed boundary or singularity treatment MUST
have a mathematical reason and an independent validation oracle. They MUST NOT be
introduced merely to remove a NaN or make a test pass. A numerical bug fix MUST
exercise the former failure and demonstrate the correct mathematical behavior.

**PCD-SCI-004 — Honest simplification.** Organizational changes MUST preserve
formulas, defaults, tolerances, and public behavior. Do not mix mechanical movement
with numerical changes. Dead code, speculative frameworks, unexplained constants,
swallowed errors, and duplicated formulas without a justified execution need MUST
NOT substitute for a clear implementation.

## 3. Modules, directories, and file ownership

**PCD-ORG-001 — Module boundaries.** Major themes MUST be actual Julia modules with
a distinct responsibility, a supported interface, and an acyclic dependency
relationship. A feature SHOULD remain in its existing owner when it shares that
owner's responsibility and interfaces. File count, a family of private types, or
CPU/GPU variants alone MUST NOT require another module. A generic `Utils` module
MUST NOT substitute for deciding ownership.

**PCD-ORG-002 — Layout.** Each themed module MUST own one directory and one entry
whose filename matches its module name, including capitalization. Module names
MUST use descriptive CamelCase; directories and implementation filenames MUST use
lowercase snake case. Family names SHOULD be plural and MUST avoid collisions
with natural type names or existing bindings. In the target layout the only root
source file is the package facade; all implementation files belong under themes.
Optional package extensions follow Section 7, not this root-source layout.

**PCD-ORG-003 — Feature folders and names.** A multi-file feature SHOULD use a named
subfolder inside its owner. A folder does not create a namespace. Judge filenames
in their directory context: `line_relaxation/workspace.jl` and feature-local
`types.jl` are acceptable when their responsibility is clear. Bare theme-level
`utils.jl`, `helpers.jl`, `common.jl`, or `misc.jl` MUST NOT collect unrelated
behavior; use a specific responsibility such as `residual_norms.jl`. A prefix or
folder alone does not justify a dumping ground. README material MUST live in
`docs/` or the top-level README, not under `src/`.

**PCD-ORG-004 — Entry contents.** Entries MUST provide, in this order where
practical: a concise module docstring; the module declaration; explicit imports;
owned exported abstract/concrete types and interface declarations; lightweight
validating constructors; implementation includes in dependency order; explicit
exports; and the closing `end`. This vocabulary includes supported internal
interfaces, not just package API. Related types and lightweight constructors
SHOULD remain together. Substantial construction, full algorithms, backend kernels,
and unrelated convenience behavior MUST live in implementation files. Small policy
validation MAY stay with its vocabulary. The entry MUST remain a readable map.

**PCD-ORG-005 — Single inclusion.** Each implementation file MUST implement one
coherent algorithm or narrowly related method family and be included exactly once,
directly by its owning entry, including files in feature subfolders. Includes MUST
be literal, unconditional top-level statements and MUST stay within the owner's
directory tree. Implementations MUST NOT contain `module`, `using`, `import`,
`export`, or `include` declarations, load their owner, or assume unrelated includes
have established names. Private helpers and types MAY live beside the
implementation that uses them. Document real internal ordering dependencies in
the entry; do not duplicate definitions to avoid declaring a dependency.

**PCD-ORG-006 — Facade.** The package entry MUST contain only package documentation,
one include per themed entry in dependency order, selected explicit imports, and
the package export list. It MUST NOT implement numerical behavior, rebind private
helpers, or expose names accidentally through include order. Replace automatic
reexports during migration with a deliberate public name list. Preserve the
identity of reexported functions and types rather than introducing substitutes.

Illustrative layout (names are not a required package API):

```text
src/
├── PhysicsPackage.jl
└── relaxations/
    ├── Relaxations.jl
    ├── point_relaxation.jl
    └── line_relaxation/
        ├── workspace.jl
        ├── sweeps.jl
        ├── cpu_cached.jl
        └── gpu_cached.jl
```

`Relaxations.jl` declares `module Relaxations` and directly includes all five
implementation files. Splitting a feature into these files does not require a
`LineRelaxation` module.

## 4. Interfaces, imports, and dependencies

**PCD-API-001 — Three interface levels.** The following distinctions MUST be explicit:

| Level | Declaration and consumers | Obligation |
|---|---|---|
| Package API | Facade exports and any other deliberately documented public surface | Preserve documented behavior and compatibility |
| Supported internal interface | Explicit exports of its owning themed module | Document purpose, inputs, outputs, invariants, mutation, and intended consumers; coordinate changes across those consumers |
| Private implementation | Names not part of either declared supported surface | Only the owner consumes them in production |

Module exports do not automatically become package API. An underscore prefix alone
does not determine privacy. Sibling consumers MUST use exported interfaces; merely
exporting an implementation helper does not establish a sound boundary. Each
cross-module interaction MUST have a scientific or algorithmic purpose and a
documented contract. Tests MAY inspect private state through its owner for a
focused invariant, without promoting that state to a supported interface.

**PCD-API-002 — Dependency direction.** The documented module graph MUST be acyclic.
Lower-level algorithms MUST NOT import solver orchestration for a type, helper,
or storage access, nor reach through the package facade. Move a genuinely shared
contract to the lowest coherent owner, or use a higher-level adapter under
PCD-STATE-005. Review ownership before introducing another shared module.
Dependency maps MUST define arrow direction and distinguish current from target
edges, including geometry, topology, I/O, and diagnostics services.

**PCD-API-003 — Import forms.** Sibling interfaces MUST use explicit symbol-scoped
relative imports. External dependencies MAY use a module import for qualified
access. Unqualified names MUST use symbol-scoped imports. Broad `using` is
prohibited. A module adding methods to another module's function MUST explicitly
import that function. All these declarations belong in the owning entry.

```julia
import LinearAlgebra                  # qualified LinearAlgebra.norm(...)
using LinearAlgebra: dot              # unqualified dot(...)
using ..Operators: apply_operator!    # consume a sibling's exported interface
import ..Operators: operator_inputs   # extend a declared interface
```

**PCD-API-004 — Method ownership.** Numerical methods belong to the theme owning
the behavior, including mapped-grid, multispecies, and backend variants. Multiple
dispatch does not justify duplicating a kernel or reversing dependencies. Methods
adapting higher-level owned types to lower-level interfaces follow PCD-STATE-005;
optional-dependency methods follow PCD-EXT-001. Type piracy is prohibited.

**PCD-API-005 — Public compatibility.** Reorganization MUST preserve public names,
binding identity, method behavior, constructors, defaults, and supported extension
points. Deliberate API changes require documentation, consumer/compatibility
analysis, examples when useful, and changelog consideration in a separate change.
Internal interface changes MUST update all affected siblings and extensions.

## 5. State, workspaces, caches, and adapters

**PCD-STATE-001 — Ownership and lifetime.** Algorithm modules MUST own their workspace
types, storage requirements, and preparation/refresh/invalidation operations.
Solver orchestration MUST compose and retain the reusable storage for the solve.
This distinguishes code ownership from the runtime object's lifetime. Policy
objects SHOULD describe algorithms without carrying mutable execution scratch.
Keep shared physical/geometry data with its conceptual owner and feature-specific
factors or scratch with the algorithm that uses them.

**PCD-STATE-002 — Concrete boundaries.** Prefer explicit arguments or small,
concretely typed borrowed views containing the data an operation requires. A
borrowed view MUST NOT take ownership of buffers, hide copies, or require a
lower-level import of the solver. It MUST state lifetime and mutation assumptions.
Avoid passing an entire solver or growing a universal options/context object as
a substitute for designing the interface.

**PCD-STATE-003 — Storage contract.** Each reusable workspace/cache interface MUST
identify who allocates, retains, reads, and mutates its storage; scalar/backend and
shape requirements; centering and interior/halo coverage; permitted aliasing;
lifetime; and concurrency/reentrancy limitations. Repeated operations MUST reuse
valid storage. A resizing or backend change MUST use an explicit reconstruction
or reallocation path rather than a hidden hot-path allocation or transfer.

**PCD-STATE-004 — Cache validity.** Cache owners MUST define what makes data ready,
what dependencies invalidate it, and which operation refreshes it before use.
Consider coefficient, geometry, boundary, topology, hierarchy, and backend changes
when applicable. Distinguish allocated storage from valid derived contents. Track
readiness or generations only where needed; do not mandate one universal cache
framework. Refresh MAY reuse storage when its shape/backend contract remains
valid. Consumers MUST NOT silently use stale data or mutate another owner's
private cache fields to repair readiness.

**PCD-STATE-005 — Adapters.** A solver-specific method MAY extend an exported
interface declared by the lower-level algorithm. The adapter belongs with the
solver-owned type and MUST only prepare the required data/services and delegate
to the algorithm. The numerical implementation stays with the algorithm owner.
Adapters MUST preserve type stability, storage reuse, and dependency direction.
Do not add a wrapper when a direct call already expresses the same boundary.

**PCD-STATE-006 — Contexts.** Algorithm-specific execution contexts MAY be used when
they express a demonstrated family of execution needs more clearly than arguments
or borrowed views. The algorithm owner MUST declare the narrow exported protocol;
the solver owns adapters for solver-specific contexts. A universal context or
hook framework is not required. Contexts MUST NOT conceal unrelated solver policy,
allocation, communication, or dynamic dispatch in local numerical operations.

## 6. Kernel design and measured performance

**PCD-KERN-001 — Setup and execution.** Keep validation, flexible configuration,
backend selection, scheduling, and user policy outside local kernels. Use function
barriers to enter concrete computation. Hot paths MUST avoid abstractly typed
fields, `Any` containers, type-changing locals, and avoidable dynamic dispatch;
prefer parametric structs and concrete storage. Global configuration bindings
MUST be `const`; kernels MUST NOT depend on mutable global state.

**PCD-KERN-002 — Memory and arithmetic.** Hot loops, sweeps, transfers, halo paths,
and repeated cycle/step execution MUST avoid unnecessary allocations, copies,
temporary arrays, allocating slices, captured heap state, logging, and string
construction. Use explicit loops or nonallocating views as appropriate. Favor
Julia's contiguous column-major access, regular indexing, disjoint writes, and
coarse enough task/tiling granularity to amortize scheduling. Hoist repeated terms
when useful without obscuring the formula. Hidden allocating broadcasts and
comprehensions MUST NOT enter kernel bodies.

**PCD-KERN-003 — Safe optimizations.** `@inbounds`, SIMD assumptions, and fast-math
transformations require established safety invariants and numerical validation.
Keep common mathematical operations shared. Specializations or duplicated kernels
require a concrete backend/throughput need and the same mathematical oracle.
Explain less-obvious optimizations, their benefit, and their assumptions locally.

**PCD-KERN-004 — Measurement.** Changes that can affect hot-path performance MUST
provide relevant type-stability, allocation, runtime, and compilation evidence.
Use representative probes such as `@code_warntype`, `@allocated`, profiling, and
benchmarks; explain any omitted measurement. Compare the same problem, scalar
type, Julia/package environment, hardware, backend, thread/rank layout, tolerance,
and warm-up state. Separate compilation/first-call effects from warmed execution;
report noise and reproducible measurements. This standard introduces no universal
timing or allocation threshold. Do not claim speed from shorter source code.

**PCD-KERN-005 — Simplicity and throughput.** Prefer the simplest implementation
meeting the measured requirements. A cosmetic refactor MUST NOT regress measured
hot-path performance; a less-readable optimization needs a meaningful measured
benefit. Establish correctness before benchmarking. Retain clear equations when
choosing between abstractions with comparable performance.

## 7. Parallel execution and optional extensions

**PCD-EXEC-001 — Ownership in parallel execution.** Distinguish global, rank-local,
block-local, interior, and halo indices. Preserve halo width and field ownership
at every level or stage. Make the owner of communication, completion, and boundary
application explicit; keep communication outside local kernels where possible.
Reductions MUST use the intended local/global domain without double-counting halos.

**PCD-EXEC-002 — Scheduling and failure.** Shared-buffer lifetimes and writes MUST
remain safe across threads, tasks, asynchronous device work, and halo exchange.
Do not insert iterative barriers, collectives, or inner-loop synchronization
without necessity and performance evidence. Preserve valid overlap opportunities.
Distributed failure paths MUST respect collective participation so one rank's
error does not strand its peers. Test affected error paths as well as success.

**PCD-EXEC-003 — Backend boundaries.** Preserve backend-generic local kernels where
supported. Isolate necessary specializations with the algorithm owner. Avoid
device scalar indexing, accidental host/device transfers, and silent host
fallbacks. Serial code MUST work without initialized MPI or an available GPU.
Unsupported modes MUST follow an explicit setup/error contract, not masquerade as
successful execution on another backend.

**PCD-EXEC-004 — Execution cases.** Applicable validation MUST consider scalar type,
dimension, geometry/basis, physical boundaries, periodicity, multiblock ownership,
uneven/thin partitions, coordinate singularities, and multigrid levels. Test empty
domains if supported, or their documented rejection otherwise. Verify decomposition
independence across valid rank layouts and backend agreement using justified
tolerances; do not infer untested support from shared source code.

**PCD-EXT-001 — Extension ownership.** Optional-dependency methods MUST live in the
package extension when that dependency is required. Extensions MUST consume core
code through declared interfaces, never include `src/` files or rely on include
order. For a core-owned operation the core owner MUST declare/export its interface;
extensions explicitly import it when adding methods. For an external operation
(such as plotting), import that operation and specialize on the supported core
types. Core modules MUST NOT depend on optional extension implementation details.

**PCD-EXT-002 — Extension verification.** An affected module migration MUST check
extension imports, activation, dispatch, and behavior. Check core loading without
optional packages, loading with affected dependencies in applicable orders, and
device execution when supported hardware is available. Missing hardware is a
reported validation gap, not a passing backend test. Extension entries MUST follow
explicit import and ownership rules; implementation files included by an extension
belong solely to that extension.

## 8. Documentation and coding style

**PCD-DOC-001 — Formatting and names.** Format Julia source with
[Runic](https://github.com/fredrikekre/Runic.jl), declared as a dependency in the
package's `Project.toml`. Use its fixed formatting rules, including
four-space indentation. Do not add a separate indentation or line-width
configuration. Runic controls mechanical layout; PCD-DOC-005 through PCD-DOC-007
control logical groups, comments, and spacing between definitions. The final
source MUST pass the Runic check after that style review.

Do not reformat unrelated code during a focused change. A deliberate package-wide
formatting migration MAY format all maintained Julia files together. Keep that
migration separate from numerical changes. Generated output and dependency caches
are not maintained Julia source.

Struct fields MUST communicate physical or algorithmic meaning. Use conventional indices
such as `i`, `j`, `k`, `d` and meaningful intermediates such as `J`, `Δx`, and
`∂u∂x`. Standard mathematical Unicode MAY follow the file's established style
when it improves clarity. Keep public names descriptive and stable; shorter local
names and unambiguous local aliases MAY aid compact kernels.

**PCD-DOC-005 — Logical groups and vertical spacing.** Functions MUST use one
blank line to separate adjacent logical groups. Each group SHOULD perform one
clear step, such as input validation, geometry access, state loading, local
arithmetic, or result storage. Keep statements that form one equation,
calculation, or tightly related set together. A function with several distinct
steps MUST NOT present those steps as one uninterrupted block.

Blank lines MUST identify meaningful boundaries, not a fixed line count. Do not
put a blank line after every statement. A short function that performs one
operation MAY remain compact. Do not add blank lines immediately after a function
or control-flow opening, or immediately before its matching `end`.

Separate the final return from preceding work when it forms a distinct result or
packaging step. Keep an immediate return beside the condition or calculation it
completes.

Apply the same grouping rules inside loops and branches. In numerical kernels,
keep state or geometry access, numerical calculation, and writeback visually
distinct when these are separate steps. Keep a short loop or branch compact when
its body performs one operation. Do not split tightly coupled component
calculations merely to create space. Keep `elseif`, `else`, `catch`, and `finally`
attached to the preceding block without an intervening blank line.

This example separates setup, cell-local scaling, state updates, and return. It
illustrates layout only; it does not define a numerical update contract.

```julia
function update_cells!(U, rhs, volumes, dt, domain)
    nvariables = size(U, 1)

    for I in domain
        volume = volumes[I]
        scale = dt / volume

        for variable in 1:nvariables
            U[variable, I] += scale * rhs[variable, I]
        end
    end

    return U
end
```

Apply these rules to new code and code changed for the task. Do not reformat
unrelated functions. Do not extract helpers, reorder calculations, change loop
structure, or alter numerical expressions solely to create visual groups. Use
blank lines and concise comments when they are sufficient.

**PCD-DOC-006 — Comments and section labels.** Put a comment immediately above
the statements it explains. When the comment starts a new logical group, put the
blank line before the comment. Add a short section comment when the purpose of a
substantial group is not clear from its code. Explain the equation, invariant, or
reason as required by PCD-DOC-004. Do not label routine statements with comments
such as "Set variables" or "Run loop." Functions MUST NOT use decorative banners
or numbered step comments. Use ordinary comments and blank lines.

**PCD-DOC-007 — Spacing between definitions.** Separate top-level definitions
with one blank line. Keep a docstring and any associated macros directly attached
to their definition. Related one-line methods MAY form a compact group without
blank lines between them. Separate different method families with one blank line.

**PCD-DOC-002 — Function documentation.** Supported package and internal interfaces
MUST document their contract. New/updated functions SHOULD have docstrings where
their behavior needs explanation; function docstrings MUST follow this shape,
omitting empty or inapplicable sections. Include returns when informative, and
explain mutation, ownership, units, centering, edge cases, and failure behavior
where relevant.

```julia
"""
    function_signature(arg1, arg2; kwarg=value)

Brief description of the operation.

Additional contract, invariants, mutation, and failure behavior as needed.

# Arguments
  - `arg1`: Meaning and relevant storage/units assumptions.
  - `arg2`: Meaning and relevant storage/units assumptions.

# Keywords
  - `kwarg`: Meaning and default.

# Returns
Return value, including whether it aliases an input, when useful.
"""
```

**PCD-DOC-003 — Type and module documentation.** Supported types MUST explain their
meaning, invariants, and ownership; type docstrings MUST use a concise description
and a `# Fields` section for relevant fields. Private types SHOULD document
non-obvious assumptions. Module entry docstrings MUST summarize the theme and its
boundary. A type template is:

```julia
"""
Brief description of the type's purpose.

Additional invariants and ownership/lifetime information as needed.

# Fields
  - `field1`: Physical or algorithmic meaning.
  - `field2`: Meaning, including mutation or readiness when relevant.
"""
```

**PCD-DOC-004 — Explain reasons.** Comments SHOULD explain the equation, invariant,
index convention, or optimization reason instead of paraphrasing syntax. Cite the
owning numerical specification when the derivation is not evident. Documentation
MUST distinguish intended behavior from observed implementation and tested
capability. Update relevant docs/examples when their contract changes.

## 9. Tests and acceptance evidence

**PCD-TEST-001 — Oracles.** Numerical tests MUST check independent mathematical
behavior: analytic/manufactured solutions, conservation, polynomial reproduction,
symmetry/energy properties, or justified convergence bands. State intent, oracle,
failure/mutation caught, and size/runtime budget. Assertions that only prove code
runs or mirror its implementation do not establish numerical correctness. Smoke,
loading, and structure checks MAY test their own explicitly limited contracts.

**PCD-TEST-002 — Comparisons.** Floating-point comparisons MUST be scaled to the
problem and scalar type with justified tolerances; absolute comparisons need a
known magnitude bound. Exact comparisons are for exact invariants. Use controlled
random seeds and robust norm/band comparisons for unavoidable execution-order
differences. Do not weaken tolerances, replace an oracle with a snapshot, or remove
coverage to accommodate a refactor. Package-specific thresholds remain in the
acceptance specifications.

**PCD-TEST-003 — Test organization.** Unit-test directories SHOULD mirror themed
source ownership with focused runners. Tests MUST load production code through
the package or module boundary, never by directly including production files.
Integration/regression suites MAY group by physical problem, geometry, or
verification objective. Preserve package regression-module lifecycle rules.
Retain meaningful tests of interfaces, invariants, dispatch, and failure modes
through migrations; add characterization where coverage is missing.

**PCD-TEST-004 — Structural enforcement.** Introduce checks with each migrated
theme for entry names/docstrings, owned files, literal direct includes, prohibited
namespace declarations, explicit imports, supported sibling interfaces, the
declared acyclic graph, and facade exports/binding identity. Check extensions when
affected. Derive dependency edges from declarations, including qualified accesses;
do not rely on an outdated hand-maintained list alone. An exact registry MAY
constrain architecture, but changes MUST reconcile the spec, registry, and docs.
Existing unmigrated code needs an explicit debt record, not a weakened final rule.

**PCD-TEST-005 — Validation selection.** Use the narrowest meaningful checks first,
then the affected package acceptance tiers. Mathematical invariants, integration,
CPU/device comparisons, multiple MPI layouts, and allocation/runtime/compilation
checks answer different questions; none substitutes for the others. Run expensive
suites when required by the change or package gates. Report precisely what was
and was not run. A spec-only change needs document parity, references, requirement
coverage, and scope verification; it does not require numerical benchmarks.

## 10. Migration and review

**PCD-MIG-001 — Inventory first.** Before moving production code, record the include
tree, symbol/type ownership, public exports, private root aliases, duplicate
includes, internal consumers, extension methods, and include-order dependencies.
Propose the acyclic target graph and identify interfaces/storage boundaries needed
before extraction. Add missing characterization and capture representative
performance baselines before moving affected hot paths. Inventory itself does
not move production files.

**PCD-MIG-002 — Bounded phases.** Migrate one coherent theme or dependency layer at
a time. Establish foundational contracts before their consumers; never require
partially initialized facade bindings. Each phase MUST leave the package
functional with appropriate validation. A package appendix records sequencing
and debt; it need not prescribe the entire source migration in advance.

**PCD-MIG-003 — Compatibility and cleanup.** Temporary explicit facade imports MAY
preserve supported API. Temporary cyclic dependencies, duplicated implementations,
runtime/conditional includes, and migration-only forwarding shims between themes
are prohibited. Purposeful adapters under PCD-STATE-005 are permitted. Remove
obsolete aggregators, duplicate definitions, and migration-only shims as their
owners migrate; preserve intentional public compatibility separately.

**PCD-MIG-004 — Scope and evidence.** Preserve unrelated user edits and avoid
unrelated formatting or renaming. If a boundary change requires different physics,
stop that migration and handle the numerical change separately. Each phase MUST
report changed ownership/interfaces, correctness and performance evidence, required
unexecuted checks, and remaining nonconformance. Update the relevant specs,
architecture docs, examples, and structural checks in an implementation migration.
A spec-only revision instead records the follow-up discrepancies explicitly.

Review acceptance requires answers to all applicable questions:

| Question | Governing requirements |
|---|---|
| Are the formulas, physical contracts, and public behavior preserved? | PCD-SCI-001 through PCD-SCI-004; PCD-API-005 |
| Is each theme/file clearly owned and its entry a readable map? | PCD-ORG-001 through PCD-ORG-006 |
| Are supported interfaces purposeful, explicit, and acyclic? | PCD-API-001 through PCD-API-004 |
| Are storage lifetime, borrowing, mutation, and refresh responsibilities explicit? | PCD-STATE-001 through PCD-STATE-006 |
| Are execution, extensions, and measured performance covered? | PCD-KERN-001 through PCD-KERN-005; PCD-EXEC-001 through PCD-EXEC-004; PCD-EXT-001 and PCD-EXT-002 |
| Are documentation, code grouping and spacing, mathematical tests, structural checks, and gaps accounted for? | PCD-DOC-001 through PCD-DOC-007; PCD-TEST-001 through PCD-TEST-005 |
| Is this a bounded phase with an accurate handoff and synchronized spec copies? | PCD-MIG-001 through PCD-MIG-004; PCD-GOV-003 |

<!-- END COMMON PHYSICS-CODE DESIGN -->


## Appendix A. Package application

### Current organization

Observed on 2026-11-23.

| Module | Owns | Depends on |
| --- | --- | --- |
| [Package entry](../src/ConservativeGradientSchemes.jl) | Module wiring and public exports | Themed modules |
| [Interfaces](../src/interfaces/Interfaces.jl) | Abstract scheme and limiter types and operation contracts | No scheme modules |
| [Limiters](../src/limiters/Limiters.jl) | Scalar slope limiters and MP5 face bounds | Interfaces |
| [MEGSchemes](../src/meg/MEGSchemes.jl) | MEG stencils, reconstruction, and sweeps | Interfaces, Limiters, KernelAbstractions |
| [CentralSchemes](../src/central_schemes/CentralSchemes.jl) | Centered stencils, gradients, and shared-face metric sweeps | Interfaces, KernelAbstractions |
| WENO, PPM, BVD, MLP | Names reserved for future schemes | No scheme modules |

### Verification

- Run the [CPU suite](../test/runtests.jl) through the package test entry point.
- Run the [Metal suite](../test/metal.jl) for affected device sweeps on
  compatible hardware.
- Check changed Julia files with the Runic 1.13.0 version pinned in [CI](../.github/workflows/CI.yml).
  Follow Section 8 for logical grouping and docstrings. CI uses a separate
  tooling environment, which differs from the dependency placement in PCD-DOC-001.
