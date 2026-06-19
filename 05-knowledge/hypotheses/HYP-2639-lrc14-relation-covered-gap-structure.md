---
id: HYP-2639
title: LRC(14) relation-covered GAP structure - high small-relation coverage must be typed by summand shell and sign
status: OPEN
source: codex-2026-06-19-S26
depends_on:
  - HYP-2637
  - HYP-2638
  - HYP-2635
  - HYP-2636
  - HYP-2634
  - HYP-2067
  - HYP-2118
related:
  - THM-538
  - HYP-2607
  - HYP-2604
  - HYP-2633
  - OPEN-Q-108
---

# HYP-2639 - LRC(14) Relation-Covered GAP Structure

## Claim

The HYP-2635 live lead, as sharpened by HYP-2637's weighted relation-fiber
split,

```text
no dissociated stranger -> every element in a small relation
-> high additive energy -> Freiman GAP
```

needs one more typed visibility/sign layer.  Relation coverage alone does not imply the
elementary Freiman `3k-4` small-doubling window.  The useful object is instead a
small-relation hypergraph whose edges are labelled by:

```text
summand shell          C=a+b
multiplicand test      C | w
relation sign type     balanced vs observer-coupled
visibility             C in S vs C not in S
```

The theorem target should be:

```text
dissociated stranger peel
OR
relation-covered row has enough low visible/hidden summand-shell payload
to channelize the signed reciprocal tail before absolute values.
```

This refines rather than rejects the additive-energy/GAP route.  HYP-2637 shows
that bounded weighted summand fibers and low relation nullity are the right
Freiman-side object, and HYP-2638 reserves the Freiman `3k-4` small-excess
finite pocket.  HYP-2639 adds the LRC-side guardrail: a global energy
count is too coarse unless the retained summand node, observer visibility, and
sign/parity type survive the quotient.

## Evidence

The scout `04-computation/lrc14_relation_covered_gap_scout_codex_s26.py`
compares AP, shifted AP, `V*`, the two KPS third-pocket examples, and sampled
wide relation-covered rows.

| row | `|S+S|` | `3k-4?` | energy | visible folds | all nonzero covered | `M*n` |
|---|---:|---:|---:|---:|---:|---:|
| `AP13` | 25 | yes | 1469 | 36 | yes | 1.000 |
| `shiftAP13` | 25 | yes | 1469 | 0 | yes | 4.789 |
| `Vstar` | 36 | no | 1169 | 31 | yes | 1.000 |
| `KPS_third_pocket_A` | 31 | no | 168 | 10 | yes | 2.571 |
| `KPS_third_pocket_B` | 31 | no | 152 | 8 | yes | 2.152 |

The decisive guardrail is the AP/shifted-AP pair: they have the same additive
energy and the same sumset size, but one is floor-tight and the other is very
safe.  The difference is not energy; it is that AP's pair-sum collisions land as
observer-coupled visible folds, while shifted AP pushes the same balanced
collision profile into hidden summand shells.

The KPS third-pocket examples and sampled wide rows have every nonzero element
in a small motif, yet `|S+S|` is far above `3k-4` for `k=8`.  KPS S12's
sumset-excess scan is consistent with this: excess `0` is the AP/Freiman line,
whereas two-dimensional GAP pockets and the third pocket are much looser in
`L_y`.  Thus a direct one-dimensional Freiman route is too strong for this
pocket.  What survives is relation coverage with low-shell labels.

## Integration With HYP-2637, HYP-2638, And KPS S12

HYP-2637 supplies the upstream inverse-combinatorics split:

```text
uncovered element -> peel
full bounded relation-fiber coverage + small nullity -> Freiman/GAP pocket
```

HYP-2638 then reserves the low-dimensional Freiman small-excess pocket:
when `|E+E| <= 3k-4`, the row is finite after AP normalization.

The present refinement says that the LRC proof cannot stop at the scalar GAP
address.  It must also retain whether a collision is an observer-visible fold
or a hidden balanced shell.  The AP/shifted-AP calibration is the minimal
example: identical pair energy and sumset size, opposite LRC hardness.

KPS S12 adds the complementary sumset-excess signal: excess `0` is AP, small
positive excess still stays far below the cap, and the named third-pocket rows
sit at larger excess with large slack.  So the live proof obligation is not to
force every relation-covered row into a one-dimensional GAP.  It is to show
that every non-AP relation-covered pocket either has enough Freiman dimension
slack or enough signed shell cancellation once summand visibility and
multiplicand parity are retained.

The later KPS S12 Freiman-dimension recursion sharpens this into a cleaner
partition: `d=1` is the AP pocket and apparent unique top; `d>=2` GAP pockets
show a large `L_y` dimension penalty; the remaining analytic work is proving
that penalty rigorously enough and then handling the signed relation-rich tail.
HYP-2639 is the sign/visibility interface for that tail: dimension tells how
many GAP directions exist, while summand-shell labels say which relations are
actually visible to the LRC observer.

## Addition, Multiplication, Even/Odd, Positive/Negative

The S560 summand/multiplicand graph dictionary is the right carrier:

- Addition supplies candidate pinch denominators `C in S+S`.
- Multiplication tests clearance: at `t=m/C` with `(m,C)=1`, runner `w` is not
  cleared exactly when `C | w`.
- Balanced relations such as `a+b=c+d` have coefficient sum `0`; they are
  translation-invariant even-scalar data.
- Folds such as `a+b=c` have coefficient sum `1`; they are
  observer-coupled, translation-sensitive, odd-marked data.
- Midpoint relations `a+c=2b` are balanced but carry the 2-adic midpoint
  channel.

So "positive/negative" in a relation is not decoration.  It decides whether the
relation is a balanced energy shadow or a signed observer-coupled fold that can
move the lonely threshold.

## Proof Route

1. Replace the failed dissociation dichotomy by a three-way split:
   dissociated stranger, Freiman-small GAP, and relation-covered non-GAP.
2. For the relation-covered non-GAP pocket, build the labelled relation
   hypergraph over pair-sum nodes `C`.
3. Delete exact low-height defect-zero motifs as in HYP-2634.
4. Apply HYP-2636's block-frequency transfer by summand channel before any
   absolute value.
5. Use multiplicand/divisibility data only as the clearance sieve on each
   candidate denominator.

## Tournament Analysis

The pairwise observable ranks proof lanes by certificate power, decoy
resistance, sign preservation, graph-bridge content, and maturity.

Hamiltonian path:

```text
observer_coupled_visible_folds
> low_hidden_summand_shells
> multiplicand_clearance_sieve
> relation_coverage_hypergraph
> freiman_small_doubling_GAP
> balanced_pair_energy
> raw_runner_vertices
```

Fingerprint: score histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, no directed
3-cycles, one Hamiltonian path.

Assumption challenge: the vertices are not runners.  I considered runners,
pair-sum nodes, relation hyperedges, divisibility obstructions, sign patterns,
visible folds, hidden summand shells, Fourier modes, and proof obligations.
This quotient preserves whether a row has low observer-coupled relation payload;
it destroys exact time geometry until the multiplicand sieve and block-frequency
transfer are reattached.  The challenged assumption is that high additive energy
alone is a valid LRC hardness certificate.
