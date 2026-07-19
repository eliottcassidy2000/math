---
id: THM-1177
title: The balanced (1,2,3) continuum mechanism and the exact four-point gap carrier
status: CORRECTED PROOF SKETCH plus PROVED carrier equivalence.  The exactly balanced configuration and the 2/21 proportional-family value are exact.  The former claim that positive bad measure forces exact balance is false: bad gap vectors fill a polytope, and nonproportional (2,3,5) contains two exact bad intervals of total measure 4/105.  The general 2/21 ceiling and uniform r=5 therefore remain OPEN
source: kind-pasteur-2026-07-18-S128 (cont.77; owner: prove d ∝ (1,2,3) is the maximiser)
depends_on:
  - THM-1174    # the measured ceiling this seeks to explain
  - THM-1173    # the continuum limit and exact endpoints
related: [THM-1141, MISTAKE-171]
script: 04-computation/maximiser_proof_kps_S128c77.py (+ .out)
---

# THM-1177 — the balanced mechanism and its exact gap carrier

> **Absorbing codex-S74 / MISTAKE-171.** My THM-1141 replacement target, "max gap ≥ (4/3)·mean",
> is false: their exact row P = {1,…,8}, J = [1/14,13/112], killers (108,109,110,111) gives
> L/μ_actual = 638/573 < 4/3. The error was that my measured ratio 3.34 was taken against the
> *uniform-interleaving benchmark* m₀ = 3/(7Σk), not against the *actual component mean* —
> different denominators. Their D·B decomposition (D = L_max/μ_actual, B = μ_actual/m₀) keeps
> the two effects apart, and the row succeeds through baseline gain B despite D < 4/3. Nothing
> below depends on that retracted target; this is the continuum-measure line, not the ratio line.

## Referee correction: bad does not force exact balance

Three disjoint full teeth have total width `1/2`, so their four complementary
gaps also total `1/2`.  The condition that every gap is at most `1/6` does
**not** force the four gaps to equal `1/8`: it leaves the polytope

```text
p_i>=0,             p_i<=1/6,             sum_i p_i=1/2.       (1)
```

Its vertices include the permutations of `(0,1/6,1/6,1/6)`.  Exact balance
is the centroid of the bad polytope, not a consequence of badness.  This is
the precise logical gap in the original maximizer sketch.

## 1. Exact balance sits at ratio 1 : 2 : 3

Laying out `1/8, tooth, 1/8, tooth, 1/8, tooth, 1/8` puts the tooth right
edges at

```text
h=(7/24,7/12,7/8),                                    (2)
```

which is exactly in ratio `1:2:3`.  The four pieces are all `1/8`, so this
configuration is bad.

## 2. What the affine-ratio argument really proves

Since `h_i=(7/6) frac(-d_i u)`, **exact** balance requires the three
fractions to have ratio `1:2:3`.  On any interval with no wrapping,
`frac(-d_i u)=c_i-d_i u` is affine, and

```text
h_3/h_2=(c_3-d_3u)/(c_2-d_2u)
```

is nonconstant unless `c_3/c_2=d_3/d_2`.  Consequently an interval of
exactly balanced configurations forces `d_3/d_2=2` and `d_4/d_2=3`.
Positive bad measure does not require such an interval, because it may move
through the interior of (1).  The affine-ratio argument identifies the
balanced ray but does not prove that this ray maximizes the full bad set.

The failure is concrete.  For the nonproportional drift triple `(2,3,5)`,
every

```text
5/21 <= u <= 9/35                                      (3)
```

is bad.  The cyclic order of the four points is

```text
0 < 5u-1 < 2u < 3u < 1,
```

and its cyclic gaps are

```text
5u-1,        1-3u,        u,        1-3u.              (4)
```

Each is at most `2/7` throughout (3).  Reflection `u -> 1-u` supplies a
second interval of the same width, so this nonproportional row has bad
measure at least

```text
2(9/35-5/21)=4/105>0.                                  (5)
```

## 3. The ratio probe, correctly scoped

On a 2520-grid, the fraction of `u` whose tooth-edge ratio lies within one
percent of `1:2:3` is `0.3329` for `(1,2,3)`, `0.3325` for `(2,4,6)`, and
zero for the sampled nonproportional rows `(1,2,4)`, `(1,3,5)`, `(2,3,4)`,
`(1,2,5)`, and `(3,5,7)`.  This probes exact balance only; it is not a
bad-set census, as (3)--(5) demonstrate.

## 4. On the proportional family the value is exactly 2/21

For `d=(m,2m,3m)`, THM-1173 gives `2m` reflection-paired runs of width
`1/(21m)`, hence total bad measure `2/21`.  The first run is

```text
[5/(21m),2/(7m)],
```

with `F=1/6` at both endpoints and `F=5/36` at its midpoint.  Multiplication
by `m` is measure-preserving on the circle, so the calculation is uniform in
`m`, not just a check of the first three rows.

## 5. Exact four-point gap carrier

Put `x_i=frac(-d_i u)` and order the three values as
`x_(1)<=x_(2)<=x_(3)`.  The corresponding tooth is

```text
[7x_i/6-1/6,7x_i/6] intersect [0,1].                  (6)
```

The four possible complementary gaps therefore have lengths

```text
max(0,(7x_(1)-1)/6),
max(0,(7(x_(2)-x_(1))-1)/6),
max(0,(7(x_(3)-x_(2))-1)/6),
max(0,1-7x_(3)/6).                                    (7)
```

Every quantity in (7) is at most `1/6` if and only if

```text
x_(1)<=2/7,
x_(2)-x_(1)<=2/7,
x_(3)-x_(2)<=2/7,
1-x_(3)<=2/7.                                         (8)
```

Reflection of the circle removes the minus sign.  Hence the continuum row is
bad exactly when

```text
the largest cyclic gap of {0,{d_2u},{d_3u},{d_4u}} is at most 2/7.  (GC4)
```

This is the faithful carrier: four moving points and their cyclic gap word,
not the ratio of three tooth edges.  On a fixed cyclic-order chamber
`0 -> p -> q -> r -> 0`, `(GC4)` says that the four directed Hamilton-cycle
increments

```text
{pu}, {(q-p)u}, {(r-q)u}, {-ru}
```

all lie in `(0,2/7]`; their sum is automatically one.  The chamber endpoints
are beat frequencies of the four labels.  This turns the open ceiling into a
one-dimensional rational-polytope intersection problem and aligns it with
THM-1170's beat clock.

## Honest status

The balanced-ray calculation, proportional `2/21` value, explicit
nonproportional intervals (3), and carrier equivalence `(GC4)` are exact.  The
general inequality

```text
measure{u:(GC4)} <= 2/21                               (9)
```

is still open.  The original affine-ratio sketch does not prove (9), and a
wrapping case split cannot repair the missing implication from badness to
exact balance.  A proof must bound how an arbitrary integer Kronecker line
crosses the full gap polytope (8), summed over its cyclic-order chambers.
Uniform `r=5` remains open.

## Named next

- Partition by the six cyclic orders of the four labels and bound the total
  length of their Hamilton-cycle window intersections.
- Test a rearrangement or majorization theorem for four-point Kronecker gap
  words; `(0,1,2,3)` is the equality target and `(0,2,3,5)` is the first
  exact nonproportional stress row.
