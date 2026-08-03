# Weighted-backbone dominance sections and a selector-sign representative

**Status:** PROVED SYNTHESIS NOTE.  The finite statement is promoted as
[THM-3287](../01-canon/theorems/THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut.md)
after an independently implemented hostile audit.  Both primary and audit
companions agree with their stored transcripts in normal and optimized mode.

## Inheritance and definition

The closest proved mechanisms are
[THM-3254 -- first-shell two-row clutch](../01-canon/theorems/THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go.md),
[THM-3269 -- scale-invariant weighted polarization](../01-canon/theorems/THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization.md),
[THM-3277 -- critical-phase geodesic backbone](../01-canon/theorems/THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas.md),
and the newer
[THM-3278 -- selector-origin bipartition](../01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md).

For an ordered core edge `u -> v`, write the blend as
`lambda f_u+f_v`.  Its exact clutch strength is attained by pairs

```text
(small-ratio trap state, large-ratio trap state).
```

Define the oriented maximizing-witness relation by reversing this pair:

```text
R_(u->v) = {(large-ratio trap, small-ratio trap)}.
```

Thus the relation follows dominance from row `u` at large `lambda` to row
`v` at small `lambda`. Reversing the row edge reverses the relation exactly.
A lift of a row path is a sequence of reset-link states whose adjacent pairs
belong to these relations. It is a compatibility section of static
maximizers, not a response trajectory.

## Exact signal

Independent reconstruction from THM-3238's coefficient formulas gives:

- all `36` nonempty oriented vertex-simple paths in THM-3277's seven-edge
  backbone lift;
- all `14` oriented phase minimizers lift, including both reversal ties at
  targets zero and six;
- only `68/132` oriented paths in the selected tree and `2442/21226` in the
  full core lift;
- the backbone has exactly four global sections on row order
  `(2,3,7,10,11,18,19,21,22)`.

The four sections are

```text
(Q+4,Q-1,Q-7,Q-1,Q+4,X,Q-1,Q-1,Q+4),
X in {Q+1,Q+2,Q+4,Q+5}.
```

The mechanism is therefore not merely path-by-path luck: one global section
restricts to every backbone path, and the only free fibre is row 18.

The clean hostile extension is the canonical-root path `2-17-16`.  The first
edge requires the middle witness `Q-1`, while the second requires `Q-8`; the
fibres are disjoint.  There are exactly four shortest oriented nonlifts in
the selected tree:

```text
2-17-16, 16-17-2, 11-3-16, 16-3-11.
```

This identifies the first failure as a length-two fibre-gluing obstruction,
not a late product-weight effect.

## The common state and the newer selector cut

The physical link state

```text
B=Q+{4}=(1,3,3,4,4,5,6,7,8)
```

occurs in a maximizing pair on every backbone edge. More strongly, its exact
response signs on all twelve core rows are

```text
positive: {2,11,16,18,22},
negative: {3,7,10,13,17,19,21}.
```

These are exactly THM-3278's small/full selector-origin classes. Hence the
availability cut has one literal analytic representative: evaluation at the
single physical state `Q+{4}`.

There is an explicit primitive positive integral weighting of all twelve
rows for which evaluation at `B` is `+2M` on each small-side row and `-M` on
each full-side row. Every one of the 22 core edges therefore has source
margin `M>0`, and the full twelve-row blend has margin `3M>0`. The companion
prints the weights and the exact value

```text
M=612782438251008238327605167443444079001600.
```

Restricting the same construction to the nine backbone rows gives equal
edge margin

```text
355648542223452256719445831365899059200
```

and full nine-row margin three times that value.

## Type and loss ledger

```text
source: oriented simple row paths in the weighted response core
target: compatible sections of maximizing reset-link trap states
map: each edge -> its large-dominance-to-small-dominance witness relation
preserved: exact clutch maximum and equality of the shared middle witness
lost: actual ratio endpoint, response magnitude, and chronological state
sidecar: reset Q plus the ordered row-dominance gauge
hostile: 2-17-16 with disjoint Q-1/Q-8 middle fibres
```

Across both row orientations there are eighteen distinct dominance arrows.
Every one changes two pole-count coordinates and has count-vector `L1`
distance two; a physical one-pole transition has distance one.  Every
Q-directed reset-link edge ends at `Q`.  Consequently the result is a
Cech-style dominance handoff only. It is not response composition, a positive
current, a Gaussian moment functional, an owner phase, or an LRC decrement.
This boundary also obeys MISTAKE-354: an abstract compatibility/realization
statement cannot be read as original-response-compatible sequential dynamics.

## Reproduction

Run

```text
python3 04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py
python3 -O 04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py
```

and compare with
`05-knowledge/results/gmc_backbone_maximizing_witness_section_scout_20260803.out`.
The LF-normalized hashes are

```text
script aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6
output 89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3
```

The exact companion has no assertion nodes, floating literals, randomness,
or imported maximizing-witness certificate.

The independent audit builds the reset link by admissible add/delete
operations, solves the source-blend inequalities directly in both
orientations, and imports or executes neither the primary companion nor its
relation bank.  Run

```text
python3 04-computation/gmc_backbone_maximizing_witness_section_independent_audit_20260803.py
python3 -O 04-computation/gmc_backbone_maximizing_witness_section_independent_audit_20260803.py
```

and compare with
`05-knowledge/results/gmc_backbone_maximizing_witness_section_independent_audit_20260803.out`.
Its LF-normalized source/output hashes are

```text
bab55f108a2dd1695465b4195b6565b15c8b8646028ebe756972ea0230621a98
123ddb9ed988907fd758e6470c1474eaa36c5d52e394ea6efb42d9d114cdf332.
```

## Exact sequence closure

[THM-3288](../01-canon/theorems/THM-3288-maximizing-witness-lifted-walk-rational-series.md)
turns the same static relations into finite symmetric decorated adjacency
matrices.  Their walk counts have rational generating functions with minimal
prefix orders `10`, `14`, and `15` on the backbone, selected tree, and full
core.  The full core has a genuine one-step initial atom: its even
degree-fourteen relation has scalar residual `-1392` at index fourteen and
annihilates only after one adjacency step.  This supplies efficient exact
sequence evaluation while retaining the same nonchronological scope.
