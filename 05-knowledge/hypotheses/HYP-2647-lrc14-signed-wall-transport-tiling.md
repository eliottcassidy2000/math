---
id: HYP-2647
title: LRC(14) signed wall transport tiling - weighted positives and negatives are a quotient transport ledger
status: OPEN; conceptual synthesis
source: codex-2026-06-19-S30
depends_on:
  - HYP-2646
  - HYP-2643
  - HYP-2642
  - HYP-2640
  - HYP-2639
  - HYP-2638
related:
  - HYP-2249
  - HYP-2039
  - HYP-2628
  - HYP-2627
  - THM-013
  - THM-116
  - THM-353
  - PROP-001
  - OPEN-Q-108
---

# HYP-2647 - LRC(14) Signed Wall Transport Tiling

## Claim

The weighted positive and weighted negative terms in the k=9 AP-to-defect
certificate should be regarded as a signed transport matrix on a common wall
tiling, not as two unrelated scalar totals.

For two rows `E` and `E'`, form the common refinement of all sector walls for
the two lonely-runner arrangements.  On each atom `I`, keep a bucket such as

```text
alpha_E(I) = (missed-sector set, empty count N, visible-fold target, sign type).
```

Then the perturbation `E -> E'` defines a finite transport ledger

```text
T_{a,b}(E,E') = meas{x in I : alpha_E(I)=a and alpha_E'(I)=b}.
```

The HYP-2642 wall-transfer totals are the signed valuation of this matrix:

```text
positive_weight = sum_{a,b} max(0, g(b)-g(a)) T_{a,b}
negative_weight = sum_{a,b} max(0, g(a)-g(b)) T_{a,b}.
```

The useful object is the addressed matrix `T`, not the two scalar sums.  The
scalars become meaningful only after the quotient has retained enough address:
which wall moved, which fold target moved, whether the sign is
observer-coupled or balanced, and which residue/coimage shell carries the
transition.

In the binding comparison

```text
A = (0,1,2,3,4,5,6,7,8)
D = (0,1,2,3,4,5,6,7,9),
```

HYP-2642 computes

```text
weighted positive transfers = 9749/158760
weighted negative transfers = 2659/39690
negative - positive         = 887/158760 = L_y(A)-L_y(D).
```

HYP-2643 identifies a fold-level coordinate behind the same loss: raw fold
count stays `12`, but three top folds move from target `8` to target `9`,
giving reciprocal transport `3/8 -> 3/9`, i.e. `1/24`.

KPS HYP-2646 adds a compatible signed quotient fact on the reciprocal-tail
side: the support-6 kernel factorizes as a finite mod-7 coefficient times a
pure reciprocal product.  That is the analytic analogue of the warning here:
keep the signed/coset address before scalarizing.

The synthesis is that the near-AP row is an arc/wall flip in a hidden signed
tiling.  It does not destroy the tile count.  It reattaches mass to a neighboring
tile with worse valuation.

## Evidence

Script:
`04-computation/lrc14_signed_wall_transport_tiling_codex_s30.py`

Stored output:
`05-knowledge/results/lrc14_signed_wall_transport_tiling_codex_s30.out`

The scout keeps the moving endpoint sector address for the perturbation
`8 -> 9`.  On the common wall refinement, each atom records:

```text
(old sector of 8x, new sector of 9x, AP missed set, defect missed set).
```

The endpoint sector transport has exact row and column balance:

```text
total mass = 1
old-sector row masses = 1/7 in every sector 0..6
new-sector column masses = 1/7 in every sector 0..6
```

Thus the signed surplus is not a raw density imbalance in the moving runner's
sector map.  It is a valuation imbalance after the balanced transport is passed
through missed-sector state.

The script reproduces the HYP-2642 scalar shadow:

```text
weighted positive = 9749/158760
weighted negative = 2659/39690
signed D-AP       = -887/158760
AP-D surplus      = 887/158760
```

Before applying weights, the common-wall atoms split as

```text
positive mass = 274/2205
negative mass = 2269/17640
neutral mass  = 4393/5880.
```

So about three quarters of the wall tiling is neutral under `g9`; the proof
pressure lives in a small addressed off-diagonal submatrix.  This supports the
HYP-2647 claim that future bounds should pair addressed transitions, not only
compare total positive and negative mass.

## Geometry

The LRC object is not just the existence of an empty sector.  By pigeonhole,
some large gap exists in many of the reduced models.  The difficult assertion
is that enough of that gap can be transported to the observer after all
projections and quotienting.

This is the HYP-2039 "hole transport" view with a sharper local model:

```text
hole worldline in the common wall tiling
  -> wall flip changes its address
  -> signed valuation measures whether the address became more or less useful
  -> fold/coimage/residue data record the information lost by coarse projection.
```

The AP is a critical confinement state.  The near-AP perturbation is a small
boundary twist: it keeps the same number of visible fold events but shifts their
target address outward.  In topological language, the scalar `L_y` sees only a
weighted boundary integral; the hidden transport matrix remembers the cellular
chain that produced it.

This matches three older repo patterns.

1. Arc flips in tournaments are local but not scalar-local.  PROP-001 and
   THM-013 say the lost/gained odd cycles through a flipped arc must be kept
   with complement corrections; THM-116 says the delta lives naturally in a
   contraction.  Likewise, an LRC wall flip should be measured after retaining
   its contraction address: missed-sector bucket, fold target, and sign shell.
2. Tiling transport checksums are bookkeeping identities.  THM-353 and the
   tiling-bucket-balance variable say off-diagonal transport is conserved by
   row.  The content is where the off-diagonal mass goes and how it is valued.
3. Euler-copy and crossing-number notes warn against quotienting too early.
   HYP-2628 says exact-period `phi` packets should be kept before squarefree
   compression.  HYP-2627 says the raw `7*6*6*5=1260` carrier matters before
   dividing to the crossing number.  Here the analogous carrier is the full
   wall/fold transport matrix before projection to total weighted sign.

## Flipping a Negative Arc to Positive

The phrase "flip an arc from negative to positive" should be read as a change
of gauge on a signed adjacency object.

One possible object is the directed graph whose vertices are wall atoms and
whose edges are atom transitions under the perturbation `E -> E'`.  Each edge
has:

- orientation: the perturbation direction;
- sign: whether `g(b)-g(a)` is positive or negative;
- weight: atom length;
- address: fold target, residue shell, and observer-coupled/balanced type.

Flipping a negative edge to positive without changing the address is not a
legal proof operation; it changes the signed valuation.  The legal operation is
closer to a tournament arc flip: compare the two orientations, then account for
the cycles or contractions whose addresses changed.  The sign scalar is only
the shadow of a local reattachment.

In the AP9-to-nearAP9 row, the flip is especially small.  The fold ledger says
the tile count is preserved and three units are moved from target `8` to `9`.
The wall ledger says the resulting atom transitions have more negative than
positive `L_y` weight by exactly `887/158760`.

## Proposed Proof Obligation

The next sharp target is a signed wall-transport lemma.

For the k=9 high-`L_y` non-AP envelope, prove that every row either:

1. has a transport matrix whose signed valuation is at most the endpoint
   near-AP valuation, or
2. exits the near-AP transport class and is already bounded by HYP-2638
   small-excess/Freiman dimension or HYP-2639 signed-shell cancellation.

A more concrete first version:

```text
Among one-gap endpoint rows F_s=(0,1,2,3,4,5,6,7,7+s),
the addressed wall-transport matrix from AP9 has maximal L_y at s=2.
```

This should be stronger than a discrepancy scalar because it can pair positive
and negative entries by address.  The row-balance identity should come from the
common wall tiling; the strict inequality should come from the target shift
`8 -> 9` and from residue/coimage phases that prevent later gaps from
recovering the lost valuation.

## Tournament Analysis

The tournament vertices here are not runners.  They are possible proof
quotients:

```text
signed_wall_transport_matrix
fold_target_transport
arc_flip_contraction_address
tiling_bucket_checksum
signed_pair_sum_gauge
exact_period_phi_packet
raw_weighted_positive_negative_totals
raw_runner_vertices
```

Pairwise observable: quotient `Q1` beats quotient `Q2` if it preserves a
strictly finer part of the LRC predicate relevant to the AP-to-defect delta:
missed-sector identity, sign of valuation, fold target, residue/coimage shell,
or contraction address.

Switch/gauge: choose the quotient that still admits a finite transport checksum
and can project to `L_y(A)-L_y(D)`.

Tie Hamiltonian path:

```text
signed_wall_transport_matrix
> fold_target_transport
> arc_flip_contraction_address
> tiling_bucket_checksum
> signed_pair_sum_gauge
> exact_period_phi_packet
> raw_weighted_positive_negative_totals
> raw_runner_vertices
```

Fingerprint: the proposed priority tournament is transitive, with score
histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}`, no directed 3-cycles, and one
Hamiltonian path.  This is not an empirical tournament result; it is the
assumption challenge recorded as a proof-design tournament.

## Assumption Challenge

Alternate vertex sets considered: runners, runner arcs, gaps, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, fold targets, matroid circuits, exact-period packets, and proof
obligations.

The chosen quotient preserves the LRC predicate only after valuation: it knows
which atoms gained or lost `L_y` weight and which hidden address caused that
gain/loss.  It destroys literal time intervals, full runner identities, and
some high-tail arithmetic.  That destruction is acceptable only if the retained
transport matrix can be checked by a row-balance identity and then projected
back to `p0` or `L_y`.

The challenged assumption is that weighted positive and weighted negative
transfers are primary.  They are probably projections of a richer object.  The
proof should first prove conservation and address-pairing in the wall tiling,
then take the signed weighted shadow.
