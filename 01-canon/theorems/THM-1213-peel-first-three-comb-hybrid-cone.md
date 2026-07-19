---
id: THM-1213
title: Exact fractional discrepancy after one peel closes the d/a at least 25/4 cone in r=5
status: PROVED analytically.  A saturated first-killer transfer followed by exact phase-aware three-comb mass/incidence accounting gives a scale-free gate.  A complete sawtooth-endpoint audit closes every proportional shape with 4d>=25a, uniformly in core, scale, and fifth killer.  The former coarse 49/6 cone remains a short corollary.  THM-1214 subsequently closes the whole ambient clustered r=5 stratum by a different owner-hypergraph quotient
source: codex-2026-07-18-S77
depends_on: [THM-1097, THM-1137, THM-1148]
related: [THM-1159, THM-1168, THM-1214, THM-1140, THM-1167]
script:
  - 04-computation/lrc14_r5_peel_threecomb_hybrid_referee_codex_S77.py
  - 04-computation/lrc14_r5_exact_discrepancy_25_4_referee_codex_S77audit.py
  - 04-computation/lrc14_r5_sharp_family_address_repair_referee_codex_S77.py
output:
  - 05-knowledge/results/lrc14_r5_peel_threecomb_hybrid_referee_codex_S77.out
  - 05-knowledge/results/lrc14_r5_exact_discrepancy_25_4_referee_codex_S77audit.out
  - 05-knowledge/results/lrc14_r5_sharp_family_address_repair_referee_codex_S77.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCPeelThreeCombHybrid.lean
---

# THM-1213 -- exact discrepancy after peeling the first comb

For a positive integer `k`, put

```text
D_k={t in R/Z: ||kt||<1/14}.
```

Let `P` be an eight-element subset of `{1,...,12}`, let `M=max(P)`, and
write `ell(P)` for the length of a largest closed component of its safe set.
Fix integers

```text
1<=a<b<c<d,       m>=1,       ma>13M,                 (1)
```

and take the first four killers to be

```text
(k1,k2,k3,k4)=m(a,b,c,d).                              (2)
```

The fifth killer is arbitrary subject to `k5>md`; it need not lie on the
same ray.  The goal is a closed component longer than `1/(7md)`, which no
open tooth of `D_(k5)` can cover.

## 1. The first removal is saturated

The exact eight-core atlas in THM-1148 proves

```text
ell(P)(13M+1)>=72/35.                                  (3)
```

Since `ma` and `13M` are integers, (1) gives `ma>=13M+1`, and therefore

```text
ma ell(P)>=72/35>13/7.                                 (4)
```

THM-1137's one-period transfer applied to a longest core-safe interval is in
its saturated branch.  It supplies a closed interval `I`, safe for the core
and `D_(ma)`, of length at least `6/(7ma)`.  Shrink it to

```text
L=|I|=6/(7ma).                                         (5)
```

This converts the first comb into a complete safe gap.  Only the remaining
three combs pay discrepancy and cut fees.

## 2. Keep the exact fractional discrepancy

For `q in {b,c,d}`, the comb `D_(mq)` sees the scaled interval length

```text
x_q=mqL=6q/(7a),       theta_q={x_q}.                  (6)
```

Define

```text
epsilon(theta)=min(theta,1/7)-theta/7,
nu_a(q)=floor(6q/(7a)+8/7).                            (7)
```

The periodic danger arc has length `1/7`.  Splitting a lift of `I` into
`floor(x_q)` full periods and a final interval of length `theta_q` gives the
phase-aware bound

```text
|I intersect D_(mq)|
 <=L/7+epsilon(theta_q)/(mq).                          (8)
```

Unlike the generic estimate in THM-1097, (8) does not replace every
`epsilon(theta_q)` by its worst possible value `6/49`.

A tooth meeting `I` has its centre in the closed radius enlargement of `I`.
After scaling, that centre interval has length `x_q+1/7`, and hence contains
at most

```text
floor(x_q+1/7)+1=floor(x_q+8/7)=nu_a(q)                (9)
```

integer centres.  Formula (9) includes the integral endpoint case.  Since
the danger teeth are open, a merely tangent tooth is harmless; using a closed
centre enlargement only overcounts.

After deleting the three remaining combs, the surviving mass `A` and number
of positive-length components `C` therefore satisfy

```text
A>=4L/7-(1/m) sum_{q=b,c,d} epsilon(theta_q)/q,
C<=1+sum_{q=b,c,d} nu_a(q).                            (10)
```

Consequently a component is longer than `1/(7md)` whenever

```text
24d/(7a)-7d sum_{q=b,c,d} epsilon(theta_q)/q
  >1+sum_{q=b,c,d} nu_a(q).                            (11)
```

Indeed, the lower and upper bounds in (10) show that (11) forces
`7md A>C`.  A longest component then has length
strictly greater than `1/(7md)>1/(7k5)`, so `D_(k5)` cannot cover it.

> **Exact one-peel gate.**  Under (1)--(2), inequality (11) implies a
> `1/14`-lonely time for every fifth killer `k5>md`.

There is no height cutoff, Covering hypothesis, or sampled scaling step.

## 3. The old gate and why peeling `a` is optimal

Replacing `epsilon` by `6/49` and `nu_a(q)` by `6q/(7a)+8/7` recovers
THM-1097's coarse three-comb gate

```text
Q_hyb=6(3d-b-c)/a-6d/b-6d/c-37>0.                    (12)
```

It has the exact positive-dispersion identity

```text
Q_hyb
 =6d/a-49
  +6(b-a)(d-b)/(ab)
  +6(c-a)(d-c)/(ac).                                  (13)
```

Thus `6d>=49a` implies strict positivity, even at equality.  The exact gate
(11) dominates (12), because it retains the true fractional fee and the
integer incidence floor.

The smallest shape coordinate is also the best single comb to peel before
applying this three-comb average.  Direct cleared-denominator subtraction
gives

```text
Q_a-Q_b = positive multiple of
          (b-a)(4d-a-b-c),
Q_a-Q_c = positive multiple of
          (c-a)(4d-a-b-c),                            (14)
```

and `4d-a-b-c>d>0`.  Peeling `d` gives a gate already below `18-37<0`.
The primary referee checks (14) on all 8,214,570 ordered shapes through
height 120.  The new threshold below is therefore not hiding a better choice
among the other three one-peel orders.

## 4. The sawtooth cost and its five champions

No `Q4` slope restriction is needed in this section.  Put

```text
(r,x,y)=(d/a,b/a,c/a).                                 (15)
```

The shape order alone gives the complete domain needed below:

```text
1<x<y<r.                                               (16)
```

For a phase coordinate `z`, define its exact obstruction cost

```text
kappa_r(z)
 =floor(6z/7+8/7)
  +7r epsilon({6z/7})/z.                               (17)
```

Then (11) becomes

```text
24r/7>1+kappa_r(x)+kappa_r(y)+kappa_r(r).              (18)
```

On each cell between a zero, peak, or incidence jump of (17), write
`n=floor(6z/7)`.  On the rising `epsilon` branch, `kappa_r` has the form

```text
n+constant+36r/7-6nr/z,
```

and is strictly increasing.  On the falling branch it has the form

```text
n+constant+r(n+1)/z-6r/7,
```

and is strictly decreasing.  Thus on `(1,r)` only the fixed left boundary,
cell endpoints, and the moving right endpoint can maximize the cost.

For

```text
25/4<=r<49/6,                                          (19)
```

the complete champion list is

```text
6  ->  41/6  ->  moving endpoint r  ->  43/6  ->  8.  (20)
```

The relevant fixed values are

```text
K_6 =6+r/7,
K_41=7+6r/287,
K_43=7+36r/301,
K_8 =8+r/56.                                           (21)
```

For the proof it is convenient to use the following piecewise upper envelope:

```text
[25/4,41/6):       K_6,
[41/6,246/35]:     K_41,
[246/35,43/6):     kappa_r(r),
[43/6,8):          K_43,
[8,49/6):          K_8.                               (22)
```

At the isolated entry values `r=41/6` and `r=8`, the displayed `K_41` and
`K_8` are conservative one-sided jump bounds because their fixed point is
the excluded endpoint `z=r`; this only strengthens the upper bound.  Away
from those two entry points, (22) is the actual supremum envelope.  The switch
`K_41=kappa_r(r)` occurs exactly at `r=246/35`.

For completeness, the cells below `17/3` cannot be discarded merely because
of a `Q4` slope condition.  On the full interval `(1,r)`, the additional
peak/jump candidates and their exact costs are

| type | `z` | `kappa_r(z)` |
|---|---:|---:|
| left jump bound | `1` | `2+r/7` |
| peak | `4/3` | `2+9r/14` |
| jump | `13/6` | `3+6r/91` |
| peak | `5/2` | `3+12r/35` |
| jump | `10/3` | `4+3r/70` |
| peak | `11/3` | `4+18r/77` |
| jump | `9/2` | `5+2r/63` |
| peak | `29/6` | `5+36r/203` |
| jump | `17/3` | `6+3r/119` |

Zeros of the sawtooth are local minima.  Every entry in this table is
strictly below the envelope (22).  Since all costs and all five pieces of
(22) are affine in `r`, it suffices to compare their rational endpoint
values.  In the order displayed above, the infima of
`(envelope)-(earlier candidate)` over (19) are

```text
4, 7/12, 181/52, 49/30, 21/8,
91/66, 61/36, 133/174, 25/34.                         (22a)
```

All are positive.  Equivalently, the worst earlier candidate on the five
successive envelope pieces is the peak at `z=4/3`, and the corresponding
piecewise gap infima are

```text
7/12, 22/35, 22/35, 35/43, 43/48.                    (22b)
```

This completes the analytic enumeration on the full ordered domain (16).
The denominator-120 champion lattice and the denominator-4096 boundary hunt
in the referee cross-check the later champion cells.  The nine earlier rows
are independently checked against every envelope piece in
`LRCPeelThreeCombHybrid.lean`.  These checks are not substitutes for the
analytic enumeration.

As additional non-proof telemetry, an independent exact integer sweep of the
formerly misrouted branch `2d>a+b+c` through `d<=160` checked 10,586,985
ordered cone shapes with no failure.  Its smallest exact gate margin was
`91/116` at `(a,b,c,d)=(12,58,72,75)`.  Neither this finite sweep nor its
height bound is used in the analytic argument.

Using the champion twice for `x,y`, the lower margin in (18) is, on the six
successive cells,

```text
[25/4,41/6):       4r-25,
[41/6,7]:          1218r/287-28,
[7,246/35]:        14-504r/287,
[246/35,43/6]:     86-12r,
[43/6,8):          1218r/301-29,
[8,49/6):          119r/28-32.                        (23)
```

Their endpoint values are respectively

```text
(0,7/3), (1,70/41), (70/41,58/35),
(58/35,0), (0,145/43), (2,65/24).                     (24)
```

All are nonnegative.  The two zero cases are still strict where needed:

- at `r=25/4`, `K_6` has the unique maximizer `z=6`, while `x!=y`, so
  `kappa_r(x)+kappa_r(y)<2K_6`;
- at `r=43/6`, the unique maximizing moving endpoint is `z=r`, excluded by
  `x,y<r`.

Therefore (18) is strict throughout (19).  For `r>=49/6`, the coarse
positive-dispersion identity (13) is already strict.

> **Uniform 25/4 cone.**  If the proportional shape in (2) satisfies
>
> ```text
> 4d>=25a,
> ```
>
> then every legal scale and every larger fifth killer close.  Equations
> (16)--(24), including the complete earlier-candidate audit (22a), prove the
> exact gate (11) directly; no asymptotic `Q4` dispatch or additional scale
> threshold is used.

## 5. Exact telemetry and method sharpness

On THM-1168's exact height-64 primitive residual bank, the progression is

```text
residual shapes                         95,336
old Q_hyb successes                         484
exact gate (11) successes                 4,028
clean 25/4-cone successes                   920
```

The exact gate adds 3,544 rows beyond (12).  Its smallest positive bank
margin is

```text
1/1120 at (a,b,c,d)=(11,32,35,39).                    (25)
```

The clean threshold improves the old `49/6` cone from 351 to 920 height-64
residual shapes.  These finite rows were already closed by THM-1168's exact
shape erosion; they serve as an independent cross-check.  The theorem here
is uniform over unbounded height.

The threshold is asymptotically sharp for this phase-blind mass/incidence
envelope.  For every `n>=3`, the Q4-residual family

```text
(a,b,c,d)=(4n,24n,24n+1,25n-1)                        (26)
```

lies immediately below `25/4`, and its exact-gate margin is

```text
-(71n+5)/[4n(24n+1)]<0.                               (27)
```

At the boundary `d=25n`, the margin becomes

```text
25/[4(24n+1)]>0.                                      (28)
```

Thus improving `25/4` requires phase placement, overlap, or component-address
information beyond the separate mass/incidence bounds in (8)--(10).  In fact,
the nominal sharp family itself becomes the first example where retaining an
address strictly improves the theorem.

### The phase-blind sharp family is repaired by its first-gap address

The saturated transfer does not merely return an interval of the correct
length: it returns a complete `ma`-safe gap.  Up to integer translation in
the `ma`-coordinate it is

```text
I_j=[(14j+1)/(14ma),(14j+13)/(14ma)].                  (29)
```

For the below-threshold family (26), `b=6a`, and multiplication by `mb`
sends (29) exactly to

```text
[6j+3/7,6j+5+4/7].                                    (30)
```

This interval meets exactly the five complete `mb`-danger teeth centred at
`6j+1,...,6j+5`.  Consequently

```text
|I_j intersect D_(mb)|=5/(7mb)
                      =|I_j|/7-1/(49mb).              (31)
```

The phase-blind bound used the discrepancy `+6/(49mb)` here, so (31) saves
one full tooth.  For `c=24n+1` and `d=25n-1`, the exact fractional fees and
incidence costs are

```text
epsilon_c=6/49-3/(98n),     epsilon_d=9/98+3/(98n),
nu_b=nu_c=nu_d=6.                                      (32)
```

Substituting the actual `b` occupancy and (32) into `7md A-C` gives

```text
(15n+1)(40n-31)/[24n(24n+1)]>0                       (33)
```

for every `n>=3`.  Hence the entire family (26), at every scale and for every
fifth killer above `md`, closes.  The genuine next obstruction below `25/4`
must therefore evade not only the scalar sawtooth envelope but also the
first-gap address orbit.  This is the first concrete gain from coupling a
metric comb estimate back to its wall-cell address.

## 6. Recursive audit and the later owner closure

Two further natural compositions were checked:

1. peel `ma,mb` exactly and apply THM-1094 to `mc,md`;
2. peel `ma,mb,mc` and apply the sharp one-comb component threshold to `md`.

They add zero rows inside the height-64 THM-1148 residual; their successes
are already owned by older transfer/ratio gates.  The one-peel exact
three-comb gate is the first new rung of this recursion.

THM-1214 subsequently closes the entire eight-core/five-killer clustered
stratum with a carrier-owner hypergraph, including nonproportional killers.
Accordingly THM-1213 is no longer an open proof-map branch.  It remains an
independent metric theorem: it identifies the exact fractional obstruction,
gives an unbounded shape cone without residue-carrier casework, cross-checks
the owner proof on 4,028 previously hard rays, and records the sharp stopping
point of the phase-blind mass/incidence method.

## 7. Tournament Analysis and assumption challenge

The naked runner/shape-order tournament on `(a,b,c,d)` is transitive, with
scores `(3,2,1,0)`, no directed cycles, four singleton SCCs, one Hamiltonian
path, and six flips under reversal.  It cannot see (11).

The proof-bearing vertices are instead the **phase obligations** `z`, labelled
by the pair

```text
(nu_a(z), epsilon({6z/7})).
```

The gauge is `r=d/a`, and an edge is oriented by comparison of the weighted
cost `kappa_r`.  Edge flips occur at the rational switches in (20), producing
the champion path

```text
6 -> 41/6 -> r -> 43/6 -> 8.
```

This challenges runners, gaps, wall events, residues, shape coordinates, and
proof obligations as possible vertices.  The faithful carrier is the first
safe gap plus the labelled sawtooth cost envelope.  A numerical-order
tournament destroys the incidence floor, fractional fee, core-entry width,
and strict endpoint data.

## 8. Reproducibility and formal boundary

The primary referee verifies the coarse identity, optimal peel order, recursive
negative audit, and height-64 telemetry.  The independent exact referee checks
406,350 local epsilon/incidence rows, 1,104,200 rational champion rows,
1,158,880 complete cone rows through `d<=400`, 9,634,045 adversarial rational
boundary rows through `a<=4096`, and 9,998 sharp-family rows.  Normal,
`python -O`, and frozen outputs are byte-identical.  A third exact referee
checks the address repair (29)--(33) on all 9,998 rows `3<=n<=10000`.

```text
04-computation/lrc14_r5_peel_threecomb_hybrid_referee_codex_S77.py
SHA-256 6084835f9e82b848415cfb2bcc067be4bac5b01c1745305473a3d2aa74e67fc8

05-knowledge/results/lrc14_r5_peel_threecomb_hybrid_referee_codex_S77.out
SHA-256 8ee249356f673d86468f2e003131f3b6ff626b6e92d60732ec95a556e81ccb15

04-computation/lrc14_r5_exact_discrepancy_25_4_referee_codex_S77audit.py
SHA-256 5f13c676fec2494fd5e584506e15a8ceb29824ddcc5d907a5ff28ad908c48cc7

05-knowledge/results/lrc14_r5_exact_discrepancy_25_4_referee_codex_S77audit.out
SHA-256 a84edb543f57c3d92c66bef9321e819d566a927f0c446d70c2ea246a220195dc

04-computation/lrc14_r5_sharp_family_address_repair_referee_codex_S77.py
SHA-256 ec824e5a5e58740f925e15b3b2c56e984fcaad5451a0f83ff73e232a9d355865

05-knowledge/results/lrc14_r5_sharp_family_address_repair_referee_codex_S77.out
SHA-256 dfd21e4bb599b696803c358a84ea574f0e5f161d4a552de3c12d6d83d9de4839
```

`LRCPeelThreeCombHybrid.lean` kernel-checks the saturation inequality, scale
cancellation, coarse dispersion identity, strict `49/6` fallback, exact-score
scale cancellation, champion assembly, all nine earlier candidates against
each of the five affine envelope pieces, all six gate-margin pieces, and the
strict address-repair numerator with no `sorry`, `native_decide`, or new axiom.
The fractional-period decomposition, sawtooth monotonicity and endpoint
exhaustiveness, strict-maximizer endpoint cases, and interval/component
assembly are the explicit next formalization layer.
