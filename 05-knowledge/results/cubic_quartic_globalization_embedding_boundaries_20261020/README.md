# Cubic/quartic globalization and embedding boundaries

20 October 2026. **Independently reviewed PASS.**

- A saturated relative plane algebra with monic original coordinate,
  immersed primitive noninjective slices, and generic semigroup <3,4>
  need not have any global monic cubic.
- Its exact conductor drops from degree six generically to degree two at
  a nodal special fiber.
- A primitive-cubic lattice identifies the precise degree-drop divisor.
  Excluding low-degree special-fiber elements suffices to globalize the
  entire monic cubic/quartic presentation.
- The immersed rational triple-point algebra C[y³,y⁴+y] has two explicitly
  inequivalent plane embeddings. Collision periods give a self-contained
  obstruction to ambient equivalence.

Read `CUBIC_QUARTIC_GLOBALIZATION_AND_EMBEDDING_BOUNDARIES.md`.
Primary-source scope is in `PRIMARY_SOURCES.md`; exact computational
controls are in `verify_cubic_quartic_boundaries.py`. The independent proof
and primary-source audit is in `review/AUDIT.md`: 74 producer checks and
56 independent exact controls, with matching normal and optimized runs.

Reproduce from this directory with:

    python3 verify_cubic_quartic_boundaries.py
    python3 -O verify_cubic_quartic_boundaries.py

These are failure boundaries and a conditional positive result, not a
Keller counterexample, an exclusion of the generic conductor-six sector,
or a resolution of the planar Jacobian conjecture.
