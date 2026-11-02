# Proof and computation audit

Baseline fb7c88537a4a4c57e1e85e5394aa085f08005da5. Audit is same-session,
with independent implementations, not an independent researcher or Lean audit.

## Mathematical checks

- Both signs of det A are retained, with full-cover multipliers0 and2.
- Selected-cycle multiplicities are compiled from every partition with parts
  at least2, rather than an assumed list of negative terms.
- D6 adds back M46 because its negative coefficient is already consumed by D4.
- The four-pair/triangle payment divides by FOUR selected submatchings,
  giving M22223<=M2223/16. The triangle-pair payment has a five-vertex
  complement with15 two-edge matchings.
- The triangle reserve is53/128. Only7/40 is spent on W, leaving153/640
  times tau. M23 remains available for P21. No allocation is reused.
- All endpoint graph polynomials are separately affine on[0,1/4] per edge.
- The augmentation pruning deletes only coordinates with nonnegative linear
  coefficient; every quadratic coefficient is at least56. The automatic
  bound drops only nonnegative quadratic contributions.
- Exact controls include transitivity, complete fairness, the equality
  witness, perturbations and P3 plus4K2. The latter has final slack367/1280,
  so graph-payment equality is not called tournament equality.

## Universes and independent paths

Generator: `/opt/homebrew/Cellar/nauty/2.9.3/bin/geng`, SHA256
`a755568c69fc6f3d5c3601381c5776ce01c85171d80bf898559b89428cd9305f`.
Completeness of its unfiltered enumeration is an explicit dependency;
matching a known count alone is not offered as a completeness proof.

`triangle8.json`: all12,346 simple graphs up to isomorphism, minimum212,
with equality only at2K4. An independent explicit edge-subset kernel
checks the same whole universe.

`extension11.json`: all12,005,168 ten-vertex simple bases, no degree,
connectivity, edge-count or tournament-realizability filters. Every labelled
eleven-vertex graph is isomorphic to one of their1,024 extensions. The
proof covers all extensions even where the program safely avoids listing
them. There are143,360 direct comparisons with a second matching-count
kernel. The independent Python kernel checks all19,456 extensions of the
19 exceptional bases. Normal and optimized Python outputs agree.

`ledger.json` gives exact rational cycle allocations. The two
`audit_matrices*.json` files check18 rational matrices and the actual
theorem consequence. These are controls, not the universal proof.

## Reproduction

Run in this directory, with nauty's `geng` on PATH:

```bash
clang++ -O3 -std=c++17 graph_certificate.cpp -o /tmp/eleven-wedge-certificate
set -o pipefail
geng -q 8 | /tmp/eleven-wedge-certificate triangle8 12346 > /tmp/triangle8.replay.json
geng -q 10 | /tmp/eleven-wedge-certificate extension11 12005168 > /tmp/extension11.replay.json
cmp triangle8.json /tmp/triangle8.replay.json
cmp extension11.json /tmp/extension11.replay.json
python3 ledger.py > /tmp/ledger.replay.json
python3 audit.py graphs > /tmp/audit.graphs.json
python3 audit.py matrices > /tmp/audit.matrices.json
cmp ledger.json /tmp/ledger.replay.json
cmp audit_graphs.json /tmp/audit.graphs.json
cmp audit_matrices.json /tmp/audit.matrices.json
```

Replace `python3` by `python3 -O` for assertion-independent replays. Explicit
exceptions enforce checks. The matrix consumer imports the predecessor's
exact determinant/path/cycle routines; it is independent of the partition
compiler, not a separately authored replacement for those inherited routines.

## Geometry and twenty-object extension

`geometry_controls.py` independently reconstructs the root5 native inverse
tree, honeycomb and triangular balls through26. It checks8,001 pairs at
each base2,3,4,5 through binary depth6 using exact rational cosh distances;
the infinite metric statement rests on the analytic proof. It enumerates
all6,431 K4-free eight-vertex graphs, verifies the two-diamond unit realization
over Q(sqrt(3)), and checks three exact rational consequences of the improved
weighted inequality. Actual divisor lists test the user's arithmetic inputs.

`twenty_controls.py` enumerates the complete square25 path bank from its
forced leaf18, verifies both transported operations on all20 objects,
and reconstructs the polynomial from actual union sizes. Invariant color
refinement proves graph rigidity with12,23,25 classes. Six-slot rotation
is an explicit failed transport. The finite-field counts in TWENTY_SHADOW
are elementary fixed-field counts, not computational conjectures.

`extension26.py` checks all4,864 active G25 orientations and all four new
port orientations. It verifies old-state transition closure, the36/38/38
refinement, and the same final equivalence starting from current legality
rather than the old quotient. Literal arc membership independently verifies
the exhibited one-edge source pair's counts1 and2. The old dead fibre is
not enumerated and is excluded from the minimality claim.

`next12_probe.py` checks37,038 padded-core cases and77 clique partitions.
Those are declared restricted universes; the all-order12 claim stays OPEN.
All four side-lane checkers pass in normal and optimized Python with identical output:

```bash
python3 geometry_controls.py > /tmp/geometry.replay.json
python3 twenty_controls.py > /tmp/twenty.replay.json
python3 next12_probe.py > /tmp/next12.replay.json
python3 extension26.py > /tmp/extension26.replay.json
cmp geometry_controls.json /tmp/geometry.replay.json
cmp twenty_controls.json /tmp/twenty.replay.json
cmp next12_probe.json /tmp/next12.replay.json
cmp extension26.json /tmp/extension26.replay.json
```

The historical THM-461 repair has elementary exact witnesses: four points
at spacing2 have no unit distances; two equilateral triangles sharing an
edge realize u(4)=5, and K4 cannot be a planar unit-distance graph. The
second real-root correction follows from243<15625/64. These do not depend
on a current unverified extremal table. The source correction and its
mechanism are logged in MISTAKES.md.
