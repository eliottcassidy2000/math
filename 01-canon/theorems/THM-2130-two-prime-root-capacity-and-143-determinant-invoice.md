---
id: THM-2130
title: "Two-prime root capacity and the 143-determinant invoice"
status: >
  PROVED. Let p be 13 or 11 and let a rank-two torus have a guard nonzero
  modulo p and at most twelve terminal characters. Partition the terminals
  into p-content blockers, nonzero reductions on the guard's projective line,
  and transverse reductions, with counts (b,r,t). If b>=1, all blockers can
  be made safe on the divided base. A two-coordinate p-root count then forbids
  a cover whenever r<=4,t<=6 for p=13, or r<=3,t<=5 for p=11. For rank eight,
  the mod-13 survivors have r>=5. Modulo 11 the only alternatives to r>=4 are
  (b,r,t)=(1,0,7),(2,0,6),(1,1,6). Hence, when the guard is nonzero modulo both
  primes and the mod-11 sparse exceptions fail, some terminal satisfies
  143|det(g,c_i). This is a necessary two-prime invoice, not LRC(14).
source: codex-2026-07-22-LRC-two-prime-root-capacity
depends_on:
  - THM-2114
  - THM-2120
  - THM-2122
  - THM-2123
  - THM-2125
related:
  - THM-2104
  - THM-2119
  - THM-2124
external: Settled Lonely Runner Conjecture for at most twelve integer speeds (through thirteen total runners).
---

# THM-2130 -- two-prime root capacity

## 1. The general `p`-root setup

Let `p in {11,13}`, let `Gamma` be a saturated rank-two character lattice,
and put

```text
K=Hom(Gamma,R/Z).                                      (1)
```

Take a nonzero guard `g` and `n<=12` nonzero terminal characters `c_i`.
Assume `g mod p` is nonzero and partition the terminal labels as

```text
B={i:c_i in p Gamma},                                  b=|B|;
R={i:c_i mod p in F_p^* (g mod p)},                    r=|R|;
T={1,...,n}\(B union R),                               t=|T|.       (2)
```

Suppose `b>=1`. We prove the capacity implication

```text
p=13, r<=4, t<=6  -> no almost-everywhere cover;
p=11, r<=3, t<=5  -> no almost-everywhere cover        (3)
```

of the guard-safe region

```text
C={X in K:||g.X||>1/7}                                 (4)
```

by the terminal danger bands `{||c_i.X||<1/14}`.

## 2. Make all true blockers safe

For `i in B`, write `c_i=p u_i`. Discard duplicate sign classes and choose a
generic integral cocharacter `d` on which the remaining `u_i` have distinct
nonzero absolute values. There are at most `b<=12` such integer speeds.
Settled LRC through twelve speeds gives `s in R/Z` with

```text
||u_i.(s d)||>=1/(number of sign classes+1)
             >=1/13>1/14.                             (5)
```

Thus some nonempty open neighborhood `U` of `Y_0=s d` satisfies

```text
||u_i.Y||>1/14                  for i in B, Y in U.    (6)
```

This uses LRC only on divided blocker directions. It is not a shifted LRC
claim about the original terminals.

## 3. Coordinates on the root fiber

Multiplication by `p` is a finite covering `[p]:K->K`. Choose a local inverse
sheet `sigma:U->K`. Since `g mod p` is nonzero, take a basis `z,v` of `K[p]`
such that

```text
g.z=1/p,                         g.v=0.                (7)
```

The `p^2` local roots of `Y` are

```text
X_jk(Y)=sigma(Y)+jz+kv,             (j,k) in F_p^2.   (8)
```

For `i in B`, equations (6) and (8) give

```text
c_i.X_jk(Y)=u_i.Y,                                      (9)
```

so every true blocker is safe on every root sheet.

## 4. Column capacity

Fix `Y in U`. As `j` varies, the guard values form a translated `p`-grid,
independent of `k`. A closed radius-`1/7` arc has length `2/7`; on either the
eleven- or thirteen-grid it contains at most four points. Hence the guard
forbids at most four columns.

For `i in R`, the value `c_i.X_jk` is also independent of `k`. Its open danger
arc has length `1/7`, so it forbids at most two columns. Therefore the guard
and all aligned nonblockers forbid at most

```text
4+2r.                                                  (10)
```

This is below `p` when `r<=4` for `p=13`, and when `r<=3` for `p=11`. Choose a
column `j` on which the guard and every label in `R` are safe.

## 5. Capacity inside a surviving column

For `i in T`, transversality gives `c_i.v!=0`. As `k` varies in the chosen
column, its values form a complete translated `p`-grid, and its danger band
hits at most two points. All transverse terminals therefore forbid at most

```text
2t.                                                     (11)
```

This is below `p` for `t<=6` at `p=13`, and for `t<=5` at `p=11`. Some `k`
survives every terminal band as well as the guard.

The selected pair `(j,k)` may depend on `Y`. Define one measurable subset of
`U` for each of the finitely many root indices on which that index is safe.
These subsets cover `U`, so one has positive measure. Its local-sheet image is
a positive-measure strict escape set after the finitely many equality
hypersurfaces are removed. This contradicts an almost-everywhere cover and
proves (3). QED.

## 6. Consequences at thirteen, including higher guarded ranks

At `p=13`, equation (3) contains the one-orbit extension of THM-2122:

```text
r=0 and n-b<=6  -> no cover.                           (12)
```

Thus, under residual transversality, it excludes

```text
(n,b)=(9,b>=3),(10,b>=4),(11,b>=5),(12,b>=6).          (13)
```

It also extends THM-2125 to every `n<=12`: any branch with at least one true
blocker, at most four guard-parallel nonblockers, and at most six transverse
labels is empty.

For rank eight, THM-2114 supplies `b>=1` when the guard is nonzero modulo
thirteen. THM-2120/2122/2123 exclude `r=0`, while (3) excludes `1<=r<=4`.
Hence `r>=5`, and the complete count list is

```text
(b,r,t)=(1,5,2),(1,6,1),(1,7,0),
        (2,5,1),(2,6,0),(3,5,0).                      (14)
```

No count in (14) is asserted empty.

## 7. The prime-eleven split

For rank eight, THM-2114 likewise supplies a terminal `11`-blocker whenever
the guard is nonzero modulo eleven. If `0<=r<=3`, equation (3) forces `t>=6`.
Together with `b>=1` and `b+r+t=8`, the only possibilities are

```text
(b,r,t)=(1,0,7),(2,0,6),(1,1,6).                      (15)
```

Every other non-guard-blocker branch has

```text
r_11>=4.                                               (16)
```

This is a projective residue statement, not rational proportionality.

## 8. The `143` determinant invoice

Assume now that the rank-eight guard is nonzero modulo both eleven and
thirteen, and that none of the three mod-eleven sparse patterns (15) occurs.
Equations (14) and (16) give terminal-label sets

```text
|R_13|>=5,                         |R_11|>=4.           (17)
```

They are subsets of the same eight labels, so `R_13 intersect R_11` is
nonempty. For a label `i` in the intersection,

```text
det(g,c_i)=0 mod 13,                det(g,c_i)=0 mod 11.
```

Consequently

```text
143 divides det(g,c_i).                                  (18)
```

Divisibility in (18) is independent of the chosen integral basis of `Gamma`.
Combining the guard-content alternatives gives the exact rank-eight invoice:

```text
13|cont(g), or 11|cont(g), or one of (15),
or there exists i with 143|det(g,c_i).                 (19)
```

## 9. Hostile control and assumption challenge

The determinant invoice is not a closure. For example, take

```text
g=(1,0),
c_i=(1,143i),                         1<=i<=5,
c_6=13(1,1),             c_7=11(1,2),             c_8=(2,3).       (20)
```

Modulo thirteen the count is `(1,5,2)`; modulo eleven it is also `(1,5,2)`.
The first five determinants pay (18), but their exact rational directions are
distinct. Under the specialization `d=(1,1)` all speeds are positive and
distinct, and the row also pays the separate `11/13` content invoices. This is
a compatibility witness for the necessary filters, not a torus cover.

The challenged assumption was that the forced primes should be used one at a
time. Their lawful intersection is label-level determinant divisibility.
Scalar counts alone would lose which terminal pays both primes.

Candidate tournament vertices were terminals, prime labels, root columns,
root points, and capacity obligations. No orientation preserves both stages of
the proof: the faithful object is the labelled bipartite incidence carrier

```text
prime -> aligned terminal label -> forbidden column/point,                (21)
```

together with the blocker-safe open set and finite-sheet null-set sidecar.
There is no intrinsic binary relation or honest tie Hamiltonian path, so
tournament fingerprints are not proof invariants here. QED.
