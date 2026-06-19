---
id: HYP-2638
title: LRC(14) Freiman small-excess certificate - finite 3k-4 pocket for the additive-energy route
status: PARTIAL certificate; Freiman 3k-4 pocket PASS for k=8,9,10
source: codex-2026-06-19-S27
depends_on:
  - HYP-2637
  - HYP-2635
  - HYP-2607
  - HYP-2604
related:
  - THM-531
  - THM-535
  - OPEN-Q-108
---

# HYP-2638 - LRC(14) Freiman Small-Excess Certificate

## Claim Being Tested

The additive-energy lead from HYP-2635/HYP-2637 has a sharp first finite
pocket: ordinary sumset excess

```text
exc(E) = |E+E| - (2|E|-1)
```

should be enough to close the low-dimensional side up to the Freiman `3k-4`
threshold.  If `|E+E| <= 3k-4`, equivalently `exc(E) <= k-3`, Freiman's
`3k-4` theorem places `E` inside an arithmetic progression of length at most
`|E+E|-k+1 = k+exc(E)`.  Since the LRC sector functional is invariant under
translation and dilation of the offset set, this turns the infinite
small-excess pocket into a finite normalized enumeration.

The intended certificate is:

1. Enumerate normalized `k`-subsets of `[0, k+e-1]` containing `0`, grouped by
   exact excess `e <= k-3`.
2. Compute exact `L_y(E)` and the sector distribution for each row.
3. Verify that `e=0` is precisely the AP row and every `e>=1` row has a
   strict margin below the AP/cap frontier for `k=8,9,10` (the dangerous
   LRC14 small-cluster sizes).
4. Use this as the finite theorem-shaped pocket inside the broader
   relation-fiber/GAP split: low excess is finite; high excess must either
   peel an uncovered relation vertex or enter a higher-rank GAP/tail estimate.

## Computation

Script:

- `04-computation/lrc14_freiman_small_excess_certificate_codex_s27.py`
- output: `05-knowledge/results/lrc14_freiman_small_excess_certificate_codex_s27.out`

The script enumerates primitive affine normal forms after translation and
dilation.  This matters because a dilated AP such as `{0,2,4,...}` has excess
zero but large ambient span; its primitive normal form is the consecutive AP.

For each `k=8,9,10`, it enumerates all primitive `k`-subsets of
`[0,2k-4]` containing `0`, filters exact `exc(E)<=k-3`, computes exact
sector distributions, and checks the Freiman hull inequality

```text
max(E)+1 <= k+exc(E).
```

No hull failures occur.

## Exact Table

```text
k=8:
  primitive forms in [0,12]: 792
  small-excess forms: 225
  AP L_y = 2633/7350 = 0.358231293
  cap_8 = 2243/5880 = 0.381462585
  best positive-excess L_y = 297/980 = 0.303061224
      at exc=1, E=(0,2,3,4,5,6,7,8)
  AP gap = 811/14700 = 0.055170068
  cap gap = 461/5880 = 0.078401361

k=9:
  primitive forms in [0,14]: 3003
  small-excess forms: 620
  AP L_y = 26083/52920 = 0.492876039
  cap_9 = 1979/4004 = 0.494255744
  best positive-excess L_y = 38681/79380 = 0.487288990
      at exc=1, E=(0,1,2,3,4,5,6,7,9)
  AP gap = 887/158760 = 0.005587050
  cap gap = 39541/5675670 = 0.006966755

k=10:
  primitive forms in [0,16]: 11440
  small-excess forms: 1644
  AP L_y = 45253/79380 = 0.570080625
  cap_10 = 55/91 = 0.604395604
  best positive-excess L_y = 3307/5880 = 0.562414966
      at exc=1, E=(0,1,2,3,4,5,6,7,8,10)
  AP gap = 1217/158760 = 0.007665659
  cap gap = 3209/76440 = 0.041980638
```

The same table checks `p0=meas(S7)`.  Its positive-excess maxima also stay
strictly below the AP row:

```text
k=8:  AP p0 - best positive p0 = 3/56
k=9:  AP p0 - best positive p0 = 71/5880
k=10: AP p0 - best positive p0 = 187/17640
```

## Reading

This closes the first finite low-dimensional pocket, assuming the standard
Freiman `3k-4` theorem.  It is not a monotonicity theorem in excess:
for example, in the table the maxima can rise again at later positive-excess
layers.  The correct object is a finite layer certificate, not an assertion
that `L_y` decreases with exact excess.

The tightest non-AP row is `k=9`, excess `1`, with only
`0.006966755` margin to the cap.  This shows why the near-AP finite certificate
must remain exact.  By contrast, the HYP-2635 third-pocket examples have much
larger excess and sit far below the cap; they belong to the higher-rank
relation-fiber/GAP tail, not this `3k-4` pocket.

## Integration With KPS S12

After the initial HYP-2638 claim, KPS S12 pushed a Freiman-dimension pocket
reframe under a colliding `HYP-2637` detail filename.  That note is compatible
with this result: its broad dimension program says AP dimension `1` is the top
pocket and higher-dimensional GAPs are penalized.  HYP-2638 supplies the exact
finite subcertificate for the first rigorous threshold, `exc<=k-3`, where
Freiman `3k-4` converts the infinite problem to normalized finite tables.

## What Remains

1. Formalize the Freiman `3k-4` import in the project proof ledger.
2. Treat small-excess rows as done by this finite table for `k=8,9,10`.
3. For `exc>k-3`, use the HYP-2637 relation-fiber split: either peel an
   uncovered relation vertex, or prove a higher-rank GAP / signed-lattice tail
   bound.  This is where the HYP-2635 wide no-dissociated-stranger examples
   still live.

## Assumption Challenge

The first challenged assumption is that the useful tournament vertices are raw
runners or ordinary pair sums.  For this pocket the vertices should be proof
obligations:

```text
Freiman_3k4_pocket
> exact_excess_layer
> AP_invariance
> relation_fiber_cover
> higher_rank_GAP_tail
> raw_pair_sum_energy
> raw_runner_vertices
```

This quotient preserves the LRC predicate only after recording the exact sector
functional on each finite normal form; it destroys the individual multiplicand
sign and reciprocal-tail data, which must remain in HYP-2636/HYP-2634.
