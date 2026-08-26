---
id: THM-4235
title: "Finite four-dimensional torus Kakeya quadratic carrier"
status: >
  FINITE-EXACT + VERIFIED.  Over F_61, the THM-4035 nonzero torus
  direction atlas has an exact four-step affine-placement ladder and a
  quadratic carrier with maximum multiplicity eight.  Four quadratic boundary
  charts give a genuine one-line-per-direction Kakeya set of size 1,868,641
  and maximum multiplicity twelve.  This is not a Euclidean Kakeya estimate or
  a finite-field-to-Euclidean transfer theorem.
source: codex-kakeya4d-broadness-multiscale-20260826 + finite_spine_upgrade
depends_on:
  - THM-4035-sixty-clock-separation-and-finite-kakeya-spine
related:
  - THM-4236-four-dimensional-kakeya-monic-cubic-focal-selector
  - HYP-2235-finite-field-kakeya-falconer-carrier
script: 04-computation/kakeya4d_quadratic_carrier_thm4235.py
output: 05-knowledge/results/kakeya4d_quadratic_carrier_thm4235.out
source_reference: 05-knowledge/reference/CORE-PAPERS-KAKEYA-4D-2026-08-24.md
script_sha256: 54d6fcae9f2c2527473b11a4f54a159c93ff56481c227b1363cc5386ce1e360f
output_sha256: 916b8c03e54a7e4039f7bb51d4463e38e8ac00d1059ab45a13f145e9fbb74dfc
---

# THM-4235 -- finite four-dimensional torus Kakeya quadratic carrier

**FINITE-EXACT + VERIFIED.**  This theorem attaches affine basepoints,
point multiplicities, local direction ranks, and the missing projective
boundary charts to THM-4035.  It is a theorem over `F_61`; it does not change
the open status of Euclidean Kakeya in `R^4`.

## 1. The inherited torus and four affine placements

Put `q=61`.  THM-4035 identifies `F_q^*` with one `C_60` clock and gives the
projective direction torus

```text
D* = {[1:u:v:w] : u,v,w in F_q^*},       |D*|=60^3=216,000.
```

For `k=0,1,2,3`, choose one affine line in each direction by

```text
L^(k)_(u,v,w)
  = {(t,f_1(t,u),f_2(t,v),f_3(t,w)) : t in F_q},

f_i(t,s)=s^2+t s  for i<=k,
f_i(t,s)=t s      for i>k.                         (1)
```

All four families have the same directions and hence the same determinant
and projective-plane statistics.  Their exact union sizes are

```text
k=0: 12,960,001
k=1:  6,696,030
k=2:  3,460,500
k=3:  1,814,460.                                  (2)
```

Thus changing only the basepoint law changes the union by the factor

```text
12,960,001 / 1,814,460 > 7.                         (3)
```

The `k=0` hostile is the concurrent star
`{0} union (F_q^*)^4`.  Direction-only broadness cannot distinguish it from
the fully quadratic carrier.

## 2. Fibre mechanism and exact multiplicity

For `s in F_q^*`, the map

```text
s |-> s^2+t s
```

has thirty doubleton fibres when `t=0`; for every `t!=0` it has two singleton
and twenty-nine doubleton fibres.  Once the first point coordinate fixes `t`,
the three spatial fibres are independent.  Consequently the `k=3` carrier has
the exact positive multiplicity histogram

```text
m=1:       480 points
m=2:    20,880 points
m=4:   302,760 points
m=8: 1,490,340 points.                              (4)
```

The weighted check is

```text
480+2(20,880)+4(302,760)+8(1,490,340)
  =13,176,000=216,000*61.                           (5)
```

The script derives `(4)` from the one-coordinate fibres and independently
enumerates all `13,176,000` incidences in the `61^4=13,845,841` point ambient
space.  The multiplicity histograms agree exactly.

## 3. Global and pinned broadness

Write a projective plane as

```text
a_0+a_1u+a_2v+a_3w=0.
```

If `r` of the three variable coefficients are nonzero, the number of nonzero
`r`-tuples with prescribed zero or nonzero sum is respectively

```text
((q-1)^r+(q-1)(-1)^r)/q,
((q-1)^r-(-1)^r)/q.                                 (6)
```

Counting coefficient-support types gives the full dual-space histogram

```text
|H intersect D*|=0:       4 planes
|H intersect D*|=3540: 14,400 planes
|H intersect D*|=3541:216,000 planes
|H intersect D*|=3600:    360 planes.               (7)
```

In particular, every plane contains at most `3,600=|D*|/60` torus directions.

At a point with `j` doubleton coordinate fibres, the incident directions form
an affine `j`-cube and span exactly `j+1` vector dimensions.  Every
multiplicity-eight point is therefore rank four, not plany.  After an
invertible affine shear and diagonal scaling, its eight directions are

```text
(1,+/-1,+/-1,+/-1).
```

Among the `C(8,4)=70` quartets, exactly `58` have nonzero determinant and `12`
have rank three.  Thus `(4)` contains exactly

```text
58*1,490,340=86,439,720                              (8)
```

pinned transverse-quartet incidences.  Strong global plane avoidance and
abundant local rank-four packets still do not determine the affine union.

## 4. Completing every projective direction

The torus omits `230,764-216,000=14,764` directions.  Partition all of
`P^3(F_q)` by the first nonzero homogeneous coordinate.  For `j=0,1,2,3`, put

```text
D_j={[0:...:0:1:s_(j+1):...:s_3] : s_i in F_q},
|D_j|=q^(3-j).                                       (9)
```

These four charts are disjoint and their sizes sum to `q^3+q^2+q+1`.  Attach
the line

```text
L_(j,s)={(0,...,0,t,
          s_(j+1)^2+t s_(j+1),...,
          s_3^2+t s_3):t in F_q}.                   (10)
```

Its direction is the corresponding element of `D_j`, so the union `K_quad`
is a genuine finite-field Kakeya set: it contains one affine line in every
projective direction.

For the all-field quadratic map define

```text
n_t(y)=#{s in F_q:s^2+t s=y}.
```

If `x=(x_0,x_1,x_2,x_3)`, its line multiplicity in `(10)` is exactly

```text
n_(x_0)(x_1)n_(x_0)(x_2)n_(x_0)(x_3)
+[x_0=0] n_(x_1)(x_2)n_(x_1)(x_3)
+[x_0=x_1=0] n_(x_2)(x_3)
+[x_0=x_1=x_2=0].                                  (11)
```

Exhausting all `61^4` points through `(11)` gives

```text
m=1:       150 points
m=2:     9,450 points
m=4:   210,151 points
m=6:       120 points
m=7:        30 points
m=8: 1,641,540 points
m=9:        60 points
m=10:    1,260 points
m=12:    5,880 points.                              (12)
```

Hence

```text
|K_quad|=1,868,641,
max_x m(x)=12,
sum_x m(x)=14,076,604=230,764*61.                  (13)
```

A second path directly enumerates every line incidence from `(10)` and gives
the same histogram.  The crude chart-sum upper bound is `1,877,824`.  For
every odd prime power `q`, the same hierarchy has size at most

```text
q sum_(r=0)^3 ((q+1)/2)^r = q^4/8+O(q^3).          (14)
```

This is consistent with the asymptotically sharp finite-field density scale
of Bukh--Chao.  It is not a new finite-field asymptotic bound.

## 5. The Euclidean shadow and the exact stopping boundary

The quadratic torus chart is the reduction modulo `61` of the intercept
`a(s)=(s_1^2,s_2^2,s_3^2)` and incidence map

```text
F(s,t)=(t,s_1^2+t s_1,s_2^2+t s_2,s_3^2+t s_3),
P_F(s,t):=det(D_s a+tI_3)=product_i(2s_i+t).         (15)
```

Thus the finite multiplicity-eight cube is a threefold quadratic fold, and
the absolute Jacobian is governed by a monic focal cubic with three focal
times (`det DF=-P_F` in `(s,t)` coordinate order).  THM-4236 isolates that
Euclidean mechanism.  Equation `(15)` is a source-target map, not a transfer
theorem: reduction modulo a prime preserves the polynomial and fibre degree
but destroys Euclidean angle, tube volume, shading distribution, two-ends,
and scale.

The full projective completion `(10)` closes the boundary-direction debt and
adds affine lines and exact multiplicity.  It still has no Euclidean tube
thickness, angular separation, multiscale parents, shadings, polynomial-Wolff
nonconcentration, or Hausdorff-dimension consequence.  Euclidean Kakeya in
`R^4` remains **OPEN**.

## 6. Replay

From the repository root:

```text
python3 -B 04-computation/kakeya4d_quadratic_carrier_thm4235.py
python3 -B -O 04-computation/kakeya4d_quadratic_carrier_thm4235.py
```

The output streams are byte-identical.  **QED.**
