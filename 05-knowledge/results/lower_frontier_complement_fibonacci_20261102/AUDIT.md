# Proof and computation audit

Status: **VERIFIED** same-agent audit, 2026-11-02. No independent-agent audit
is claimed. The exact graph check plus analytic reduction proves the scoped
order-ten theorem; numerical tests alone would not do so.

## Exact coverage and independent checks

- All 42 unordered cycle collections of total size at most ten were compiled
  from determinant signs and the inherited odd-cycle Hamilton identity.
  Subtracting every proof group leaves nonnegative coefficients in BOTH
  signed determinant branches. Selected equal-sized cycles retain their
  multiplicities; no sign of det(A) is guessed.
- The seven-vertex graph orbit census covers all 2,097,152 labeled graphs.
  Exhaustive augmentation and exact canonicalization give 12,346 eight-
  and 274,668 nine-vertex representatives. Every ten-vertex graph appears
  up to relabeling in the 140,630,016 tested augmentations. The graph counts
  are obtained by the program, not imported as a completeness assumption.
- The quadratic extension formula is checked on 6,200 augmentations against
  direct vertex-subset matching DP and disjoint-wedge counting. Canonical
  relabeling is checked on 200 seeded controls. All arithmetic is integer.
- The analytic completed square uses exact fractions. Separate direct
  Hamilton DP, integer Bareiss determinants, and cycle DP verify the FULL
  inequality on 50 seeded quarter-integral tournaments and six ordinal
  fair-matching controls. They are additional controls, not exhaustive
  weighted enumeration.
- The square25 program enumerates all 4,864 active guard orientations,
  every legal complement and every target fibre. It verifies all 32 active
  quotient classes, their weights, the inherited path-count polynomial,
  and the stable partition after one refinement. The general matching-bank
  theorem separately explains minimality. The dead class has no enabled
  command. Full tournament fibres have 268 unused bits, preserved in weights.
- At square15, same restricted state and ordinary H=1 versus H=9 is an
  explicit scope hostile. The 81/27 prior tensors remain distinct and are
  cited, not recomputed or identified by this quotient.
- Both normal and optimized Python use active exceptions, never assertions.
  No claim relies on a random positive test, approximate rank, or hash-only
  graph equivalence. Fibonacci identities are checked exactly through n=20;
  their all-size proof is the recurrence with its initial values.

## Reproduction

From the repository root, compile outside the tracked tree if preferred:

```text
g++ -O3 -std=c++17 05-knowledge/results/lower_frontier_complement_fibonacci_20261102/graph10.cpp -o graph10
./graph10 > 05-knowledge/results/lower_frontier_complement_fibonacci_20261102/graph10.json
python -B -X utf8 05-knowledge/results/lower_frontier_complement_fibonacci_20261102/audit.py
python -B -O -X utf8 05-knowledge/results/lower_frontier_complement_fibonacci_20261102/audit.py
```

On Windows use `graph10.exe`; the session executable was kept in the user
Temp directory. `--scout` runs seeded controls only and is explicitly marked
as such; it is not the exhaustive certificate.

## Publication limits

The low-order and all-large-order stability proofs have separate sources.
This packet advances the lower cutoff to ten, leaving 11..41 open. Its
five-vertex identity comes from Session138, its graph augmentation method
extends that packet's order-nine construction, and its path bank comes
from the joint-incidence packet. The quadratic graph update and matching-bank
complement closure are the additional mechanisms used here.

The packing source remains an inherited accepted premise with a specified
contact family. There is no new primary-literature claim about optimal square
packing. Q8 occurs in a different frame transport action. The equality
classification in stability is not a complete metagraph classification;
the 33-state carrier has an explicitly restricted path alphabet.
