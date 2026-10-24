# Independent review of the separator compiler

Reviewer: complement/trijection lane,2026-10-24. Status: **VERIFIED** proof
and implementation review; finite executions reproduced in normal and
optimized Python. Reviewed `../theory/SEPARATOR_SECTIONS.md`, `compiler.py`,
`check.py`, and `check_reset.py`. No producer edits were needed.

The running-intersection condition makes the owned-slot sets partition the
carrier. This is the exact reason charges and costs may be summed once even
though full separator records occur in several bags. The count recurrence
is a bijective decomposition of assignments; the two-arrival recurrence is
the same decomposition on a pair of assignments. In particular the
implementation merges literal witness dictionaries, checking common records,
rather than merely matching their charge.

The subgroup/coset section rule is sound as a sufficient rule conditional on
one retained separator record and one actual local bag tuple. It neither
assumes all feasible images are cosets nor forgets a common separator when
summing independent charge directions. Integer margin coordinates are
correctly fixed by offsets because the certified subgroup is finite.

The exact pair-fibre minimum is deliberately not advertised as a route
metric. The A={0}, B={0,1}, C={1} hostile demonstrates why independent bridge
minima cannot be concatenated. Source-fixed tables or independent
within-fibre connectivity are necessary. Actual occurrence costs remain
additive per labelled slot; counting distinct shared owner variables would
require another explicit ownership convention.

Reproduction:

* `theory/check.py`:307,767 checks in both normal and `-O`;2048 complete
  path-three count instances,512 two-arrival path instances,4096 branched
  integer-times-C3 pair instances, coset controls, and native C3-kernel
  boundary examples.
* `theory/check_reset.py`:30,886 checks in both modes. It independently
  compares source-fixed min-plus witnesses with literal triples, then
  optimizes the fixed native route with layer sizes1,3,8,3,3,1 (216 routes).
  The exact peak103,total415 improves the inherited114,433. The separate
  bottleneck pass followed by additive minimization under the attained cap
  avoids the invalid one-pass lexicographic partial-path optimization.

The native application stays finite and conditional on complete actual
factor tables. Neither compiler nor audit manufactures missing record
domains, capacity, a global Hamiltonian observer, or new native microedges.
Bounded separator width does not imply bounded event support: the actual
kernel-lift chain retains its3m occurrence cost. This matches the present
lane's distinction between safe individual pair relations and their unsafe
three-helper conjunction.
