# AGENTS.md

Follow the [Physics-Code Design standard](spec/physics-code-design.md).
Use short sentences, active voice, and consistent technical terms.

- Trace the public call path before editing. Check sample representation,
  stencil reach, logical spacing, derivative units, precision, and face indexing.
- For sweeps, document the caller-owned buffers, valid stencil domain, output
  domain, aliasing, and launch order. Select closures from the valid domain.
- Test numerical changes through the public API against an independent
  reference. Check affected backends and measure allocations and runtime
  when a change can affect a hot path.
- Follow the [package verification instructions](spec/physics-code-design.md#verification)
  for tests and formatting. Report checks you could not run.
