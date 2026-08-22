---
id: THM-3258
title: "Depth-two affine Farkas clutch and complete reset-distance gauge no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the complete support-(1,3), bank-I2 physical state space, none of
  THM-3244's 31 lawful two-row covering pairs can be flattened to a positive
  gauge whose two coefficients depend only on reset distance.  THM-3254
  closes 23 pairs on the distance-one link.  For the remaining eight, exact
  target-envelope decomposition of all 55 distance-two states gives 109
  rational ratio cells, compressed to fourteen affine Farkas circuits with
  at most four inequalities each.  Thus every graded two-row gauge already
  fails by distance two.  This does not address three-row or state-dependent
  gauges and has no GMC consequence.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-03
audit: >
  The assertion-independent exact companion pins promoted THM-3238 and
  THM-3254; loads only THM-3238's immutable definition prefix; and directly
  reconstructs all 22 lawful response values on the 67 states used by the
  proof.  It verifies the 11/55 shell census, the eight exact first-link
  gaps, all 109 target-envelope cells, fourteen declared circuit supports,
  their exact positive null weights, strictly negative affine constants,
  target-dominance domains, and complete interval coverage.  No discovered
  response value or circuit weight is trusted.  Normal, optimized and stored
  transcript replay and the LF hashes are required.  An independent hostile
  audit replayed both interpreter modes against stored output, verified the
  declared hashes, and separately rederived the arbitrary-target maximum,
  all 109 cells including the unbounded ones, envelope seams, circuit signs,
  and exact Farkas contradiction.  No defect remained.
depends_on:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
related:
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
script: 04-computation/gmc_depth_two_graded_pair_clutch_thm3258.py
output: 05-knowledge/results/gmc_depth_two_graded_pair_clutch_thm3258.out
script_sha256: 713d66641043c1846c5cfeeee2537c594d92c08f72393af43a95e894ce8af37f
output_sha256: 1d64f904e0110ea66ac06e8dc4445cdfce9ff611857c082295b814a8f06424db
hash_basis: LF-normalized bytes
---

# THM-3258 -- depth-two affine Farkas clutch and complete reset-distance gauge no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3244 constructs 31 two-row local atlases on the complete physical state
bank: at every nonreset state, at least one of the two rows has a strict
one-pole ascent toward the unique reset.  THM-3254 proves that no fixed
positive blend of either pair works globally.  It also proves the stronger
reset-distance-graded no-go for 23 pairs whose clutch is already visible on
the distance-one link.

This theorem closes the eight residual pairs.  The obstruction is already on
the next shell and is exact affine convex geometry, not a large-state or
asymptotic phenomenon.

## 1. Statement

Let

```text
Q=(1,3,3,4,5,6,7,8)                                  (1)
```

be THM-3238's unique reset, let `f_1,...,f_22` be its lawful response rows,
and write `d(sigma)=dist(sigma,Q)` for physical edit distance.  Fix any of
THM-3244's 31 covering pairs `(i,j)`.  There do not exist positive numbers

```text
a_d>0, b_d>0                                           (2)
```

depending only on `d` such that

```text
F(sigma)=a_d f_i(sigma)+b_d f_j(sigma)                 (3)
```

has, at every nonreset state `sigma`, some one-pole neighbor `tau` satisfying

```text
d(tau)=d(sigma)-1,       F(tau)>F(sigma).              (4)
```

For every pair, a contradiction occurs using only states with `d<=2`.

The 23 first-link pairs are already handled by THM-3254.  It remains to treat

```text
(2,7), (3,9), (7,14), (11,17),
(11,21), (12,13), (12,19), (14,19).                   (5)
```

## 2. The first-shell ratio and the second-shell envelope

Because every distance-one state has unique target `Q` and every response
vanishes at `Q`, scale by `b_1` and put

```text
lambda=a_1/b_1>0.                                     (6)
```

Strict ascent on all eleven first-shell states forces `lambda` into the open
interval `J_(i,j)` in the following table.  The last two columns give the
number of exact second-shell target-envelope cells and the number of affine
circuits needed after compression.

| pair | exact open interval `J_(i,j)` | raw cells | circuits |
|---|---|---:|---:|
| `(2,7)` | `(221630501005/672431853909, 1)` | 14 | 3 |
| `(3,9)` | `(1/8, 5576212/26946933)` | 15 | 2 |
| `(7,14)` | `(1/55, 21657540264/221630501005)` | 14 | 2 |
| `(11,17)` | `(554766960/270430481, 28/9)` | 15 | 2 |
| `(11,21)` | `(453787039/270430481, 3)` | 15 | 2 |
| `(12,13)` | `(1420574831/3464404479, infinity)` | 18 | 1 |
| `(12,19)` | `(38075529020/3464404479, infinity)` | 14 | 1 |
| `(14,19)` | `(351232551587/161618731130, 133113695/47074217)` | 4 | 1 |

Now put

```text
x=a_2/b_1>0,       y=b_2/b_1>0.                       (7)
```

For a distance-two state `sigma` and a directed neighbor
`tau in N_Q(sigma)`, ascent is the strict affine inequality

```text
lambda f_i(tau)+f_j(tau)
  -x f_i(sigma)-y f_j(sigma) > 0.                     (8)
```

There are exactly 55 distance-two states.  Three have one `Q`-directed
target and 52 have two.  Since an arbitrary edge may be chosen, the necessary
and sufficient condition at `sigma` is

```text
G_sigma(lambda,x,y)
 :=max_(tau in N_Q(sigma))
      {lambda f_i(tau)+f_j(tau)}
      -x f_i(sigma)-y f_j(sigma) > 0.                 (9)
```

Each maximum in `(9)` is the upper envelope of at most two rational affine
lines.  Their exact crossing points split the eight intervals in the table
into

```text
14+15+14+15+15+18+14+4=109                            (10)
```

open cells.  On each cell the maximizing target in every `(9)` is fixed.

## 3. Four-row affine Farkas circuits

On one such target domain, every required condition has the form

```text
v_h dot (lambda,x,y)+c_h > 0.                         (11)
```

The rows include selected instances of `(8)`, positivity of `x` or `y`, and,
when needed, one ratio-domain bound.  An affine Farkas circuit is a collection
of positive integers `w_h` with

```text
sum_h w_h v_h=0,
sum_h w_h c_h<0.                                      (12)
```

If all inequalities `(11)` held, their weighted sum would be strictly
positive; `(12)` says it is a negative rational number.  This is impossible.

The exact companion reconstructs fourteen such circuits.  Twelve use four
rows and two use three.  Their counts over the pairs in `(5)` are

```text
3,2,2,2,2,1,1,1.                                     (13)
```

The target-dominance interval of every selected edge is checked exactly.
The resulting circuit domains form an overlapping chain across every
`J_(i,j)`.  Thus the 109 raw cells are all covered by fourteen certificates.

At an internal envelope crossing, two targets have equal value.  Either
target is therefore a lawful maximizer.  Likewise, when a circuit ratio-bound
row becomes zero, its positively weighted strict shell inequalities remain.
Hence shared endpoints introduce no gap.  The outer endpoints are excluded
already by the strict first-shell inequalities.

## 4. A displayed three-row circuit

The pair `(12,13)` is especially transparent.  Its entire unbounded gap is
killed by the three rows

```text
x>0,
G_(1,1,1,3,3,4,5,6,7,8)>0
  using target (1,1,3,3,4,5,6,7,8),
G_(1,3,3,4,4,5,5,6,7,8)>0
  using target (1,3,3,4,4,5,6,7,8).                  (14)
```

The respective primitive positive weights are

```text
4452809708057982010202880,
4277107928796831,
371631008.                                            (15)
```

Their three `(lambda,x,y)` coefficient vectors sum to zero, while their
weighted constants sum to

```text
-9167765286556241877685416.                           (16)
```

Thus `(14)` cannot hold.  The other thirteen certificates are the same
mechanism with different exact target domains; their complete support,
weights, domains and constants have digest

```text
6e0cbf932f99cea022dcf9fa81bbaa2adec0f6251e27b5a72fcbd32aca6a9263. (17)
```

Equations `(10)--(17)` close the eight pairs in `(5)`.  Together with
THM-3254, this proves the statement for all 31 pairs.

## 5. Holotopy interpretation

THM-3254 identifies an oriented positive-cone clutch between two local
response charts.  The present theorem shows that allowing one gauge per
distance shell does not trivialize that clutch.  For the residual pairs, its
next obstruction is carried by rank-at-most-four affine circuits in the
three parameters `(lambda,x,y)`.

This is best viewed as a finite oriented-matroid shadow of a two-stage
clutch, not as topological monodromy.  The first shell restricts the transition
ratio; the second shell supplies a target-envelope cocircuit that prevents
the next local section from extending it.  A successful repair must retain
more than radial depth--for example a state label, an edge selector, a larger
row support, or a genuinely non-scalar transition object.

## 6. Scope

The theorem is exact only for THM-3238's complete support-(1,3), bank-I2
physical state space, its 22 lawful rows, THM-3244's 31 covering pairs, and
one-pole edges strictly decreasing distance to `Q`.

It does not rule out:

1. a positive combination of three or more lawful rows;
2. coefficients depending on the full state, chosen edge, pole label, or a
   richer flag than reset distance;
3. a nonlinear, matrix-valued, or signed transition clutch; or
4. another support/bank, arbitrary-radial NC2, the Gaussian Moment
   Conjecture, or any `LRC(14)` decrement.

THM-3244's successful state-dependent two-row atlas remains intact.  The
sharp conclusion is that no one-dimensional depth grading can flatten any of
its 31 two-chart atlases.

## 7. Exact verification

Run

```text
python 04-computation/gmc_depth_two_graded_pair_clutch_thm3258.py
python -O 04-computation/gmc_depth_two_graded_pair_clutch_thm3258.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_depth_two_graded_pair_clutch_thm3258.out.
```

The companion uses exact integer and rational arithmetic only.  It has no
floating point, random search, discovery cache, optimization-sensitive
assertion, or embedded response-value certificate.  The fourteen support
choices are treated only as proposed witnesses: every target, weight,
constant and domain is rebuilt and checked independently from THM-3238's
original formulas.

QED.
