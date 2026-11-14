# Three positive Chan–Pak routes: terminal division, asymmetric reservoirs, and symmetry breaking

**Chan–Pak remains open.** These are physical constructions and
exact target inverses, not a larger interval of numerical coverage.

- **[Seven without a seven-path component](AXIS_FRONTIER.md).**
  The nonsingular length-7 corridor has
  \(C_7(1,4)^t=7(47,56)^t\). At any terminal width \(4\cdot8^j\)
  it gives an exact 7-divided ordinary count on a positive primitive
  source row. If the *right* side does not divide, a separate
  *left* lattice condition does. The prescribed-target decoder
  checks both cases and constructs the actual word or rules out
  that lane. Explicit examples include \(7\cdot271\)
  and \(7\cdot3631\). [Cofactor congruence](PROJECTIVE_COFACTORS.md)
  remains complete for these fixed rows without making the
  absolute image dense.
  [A logarithmic inverse chain](OUTER_DEPTH.md) takes every large
  odd cofactor to a relaxed source with coordinates O(sqrt(p));
  the final state still needs a physical prefix.
- **[The ±2 reload is now classified](UNIT_PROBES.md).**
  For the old physical same-reset pair, across every possible common
  corridor suffix there are just two unbounded ±2 rails and four
  isolated additional repairs. Moreover no positive common loop
  can return any nonzero difference vector even up to sign. This
  redirects the positive search toward changing both sources,
  rather than a longer shared reload.
- **[One apex, two independently adjustable reservoirs](NONUNIFORM_APEX.md).**
  For arbitrary nonuniform block contacts there is an exact
  positive all-color interval transfer and bilateral seam law.
  At the four-block pattern \(1,0,1,0\) it becomes
  \(2H-3hq=(2\ell+h)(2L-q)\) in the two profiles
  \((h,\ell),(q,L)\). Factoring the shifted target inverts both
  independent dyadic widths of a fixed core pair. Literal prime and
  \(7p\) outputs check both uses; picking the correct core pair
  remains the open part.
- **[A genuine order-7 route](SEVEN_ORBIT.md).**
  For primitive \(7p\) at \(p\ge11\), any nontrivial symmetry
  is exactly order 7 and all moving orbits must be nonuniformly
  connected. A homogeneous Paley 7-module forces 63, so by
  itself it rules out those targets. But two coupled Paley
  orbits give the actual \(7\cdot2163467\), with an exact
  larger padding rail too. Ordinary prime-count witnesses
  at \(p\ge5\) are rigid by the same sieve.
- **[Sharp two-contact threshold](MINIMAL_PALEY_BREAK.md).**
  One nonuniform apex on a Paley module still forces a
  factor 3, even in arbitrary otherwise uniform surroundings.
  Two contacts give a literal prime \(2621\) and a literal
  \(7\cdot353\). The latter yields the rigid cofactor rail
  \(5171\cdot2^r-4465\), with a complete padding inverse.
  **[The 43-coordinate contact menu](PALEY_CONTACT_MENU.md)**
  replaces all \(2^{14}\) literal contact probes on each fixed
  outside reservoir by 43 canonical probes and makes
  its exact inverse explicit.

The actionable common goal is now more precise. The terminal-seven
route needs an attained state at an exact affine point, not just
a congruence. The shared-apex route needs to choose a core pair
that makes the shifted target factor in its two actual width
channels. The Paley route provides a finite modifiable contact
palette, but needs a varying outside bank to move its baseline.
These are independent and physically tested obligations; none is
silently promoted to full coverage.

Run the Python replays alongside the notes. The axis, nonuniform-apex,
orbit, and Paley-contact inverse tests are separate and reproducible.
