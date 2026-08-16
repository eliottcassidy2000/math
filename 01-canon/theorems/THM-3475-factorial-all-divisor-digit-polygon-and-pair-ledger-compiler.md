---
id: THM-3475
title: "All-divisor digit polygons and the factorial pair-ledger compiler"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED.  If a
  prime p divides N=d-1,
  the complete common Newton-slope ledger of the two resonant factorial
  rows A_(N-1)^(N+1) and A_N^(N+1) is computable from base-p digit weights,
  without constructing their large integer coefficients and without a
  Kurepa hypothesis.  For odd p the first row has an exact two-plateau
  recursion; for p=2 only exact odd anchors can enter its positive-slope
  hull.  This is a necessary factor-degree compiler, not a factor
  construction or a proof of SFC or FC.  On the 38 seven-exit residuals
  2501<=d<=2600, all-divisor ledgers close exactly 33 and leave five explicit
  degree packets; in particular the contiguous exact quadratic boundary
  reaches r=2513.
audit: >
  The odd and binary recurrences, exact coefficient anchors, convex-hull
  reduction, and separation from the unresolved nonpositive prefix were
  independently proof-audited.  A pinned coefficient engine verifies 209
  odd-prime cells, 100 binary cells, and all 309 pair ledgers.  A second,
  formula-only implementation independently reconstructs the seven-exit
  universe, all divisor hulls, the 33 closures, and the same five survivor
  sets.  Normal and optimized replays are byte-identical.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
related:
  - THM-3138-left-factorial-kurepa-bridge-and-leading-coefficient-nonvanishing
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3467-factorial-seven-exit-newton-barcode-extension-after-r2498
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
scripts:
  - 04-computation/factorial_all_divisor_digit_pair_compiler_thm3475.py
  - 04-computation/factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py
outputs:
  - 05-knowledge/results/factorial_all_divisor_digit_pair_compiler_thm3475.out
  - 05-knowledge/results/factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.out
script_sha256:
  - 834d0913eb5cd5b15684c7fb88af60e42d2a6ef36feb821e261c0498f55027ab
  - 9330ca1b991b9d5875779b9975fc88701ab36855a6527e1865e821e6cd3ea665
output_sha256:
  - 1fee475f1a09e2f191d295817bb20af1c055185fc16dd621091d55173bc87ad5
  - 6aeb576763a63412fad02da33bc777569e69bce87b936afcb0624b144acaec9f
semantic_sha256:
  - 886f2bcd66711a44668b003717dd4f39643fc0e9cb0b694708b8670eeaf21499
  - e019eb61019620cfaffa8b5bb5769e8d171f08fb5ffeb2153f209c0128d42115
hash_basis: raw bytes
---

# THM-3475 -- all-divisor digit polygons and the factorial pair-ledger compiler

**PROVED + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Let

```text
L(x^r)=r!,
A_n^(d)(v)=L((d-x+v x^2)^n),
N=d-1,
F=A_(N-1)^(N+1),                 G=A_N^(N+1).              (1)
```

Fix a prime `p|N`.  Define the two digital weights

```text
u_N(j)=nu_p C(N-1,j)+nu_p((2j)!),       0<=j<N,
w_N(j)=nu_p C(N,j)+nu_p((2j)!),         0<=j<=N.            (2)
```

THM-3161 proves that the complete actual lower Newton polygon of `G` is

```text
NP_p(G)=lower-hull{(j,w_N(j)):0<=j<=N}.                    (3)
```

The present theorem supplies the missing companion row:

- if `p` is odd and `h=(p-1)/2`, then the part of `NP_p(F)` at and to the
  right of `(h,0)` is exactly the lower hull of the raw points
  `(j,u_N(j))` in that range;
- if `p=2`, then the part of `NP_2(F)` at and to the right of `(1,1)` is
  exactly the lower hull of its odd-index raw points.

Every omitted prefix edge of `F` has nonpositive slope, whereas every edge
of `G` has positive slope.  Therefore the **complete** common finite-slope
ledger of `(F,G)` is obtained by intersecting (3) with the displayed
digit-only suffix polygon of `F`.  Its edge capacities and reduced slope
denominators give the complete local necessary factor-degree barcode by
THM-3152.  No coefficient of `F` or `G` need be constructed.

Because the constant coefficient of `G` is a unit, the common coordinate-root
capacity is zero.  A nonempty barcode is only a possible degree address; it
does not construct a rational factor.

## 2. Coefficient lower weight

Write

```text
[v^j]F=C(N-1,j)(2j)! Z_j,                                 (4)

Z_j=sum_(ell=0)^(N-1-j)
 C(N-1-j,ell)(N+1)^(N-1-j-ell)(-1)^ell(2j+1)_ell.          (5)
```

Thus

```text
nu_p([v^j]F)>=u_N(j),                                      (6)
```

with equality whenever `Z_j` is a `p`-adic unit.  The issue is not merely
to bound coefficients, but to prove equality at every point that can be a
vertex of the raw lower hull.

## 3. Odd-prime two-plateau recursion

Suppose `p` is odd.  Write

```text
N=pM,                    h=(p-1)/2,
j=pq+a,                  0<=a<p,                           (7)
```

and use the same definition (2) for `u_M(q)`, still with `nu_p`.  Put

```text
K_q=2q+u_M(q).                                             (8)
```

Legendre's formula, together with the base-`p` binomial recursion, gives

```text
nu_p C(pM-1,pq+a)=nu_p C(M-1,q),                           (9)
```

and hence

```text
u_N(pq+a)=
  K_q                              if 0<=a<=h,
  K_q+1+nu_p(2q+1)                if h<a<p.                (10)
```

For `0<=q<=M-2`, the next lower plateau satisfies

```text
K_(q+1)-K_q
 =2+nu_p(M-1-q)+nu_p(2q+1).                               (11)
```

Consequently the two successive gaps are

```text
1+nu_p(2q+1),                 1+nu_p(M-1-q),               (12)
```

and are strictly positive.  The raw weight word is therefore a strictly
rising sequence of two constant plateaux in every `p`-block.

## 4. Exact right anchors and the odd suffix polygon

At the right endpoint `a=h` of the lower plateau,

```text
2j+1=p(2q+1).                                              (13)
```

Every positive-`ell` term of (5) contains this factor, while its `ell=0`
term is `1` modulo `p`.  Hence `Z_j==1 (mod p)`.

At the right endpoint `a=p-1` of the upper plateau,

```text
N-1-j=p(M-q-1),                 p | (2j+2).                (14)
```

For `1<=ell<p`, the binomial coefficient in (5) is divisible by `p`; for
every `ell>=2`, the rising product contains `2j+2`.  Again

```text
Z_j==1 (mod p).                                            (15)
```

This includes the terminal case, where there are no positive-`ell` terms.
Thus the right endpoint of every plateau is an actual coefficient point at
the raw height.

In a strictly rising step word, a lower-hull vertex can occur only at a
right endpoint: an earlier point of a plateau lies strictly above the chord
from any earlier lower plateau to that endpoint.  Equations (6), (12), and
the exact endpoint equalities therefore imply

```text
NP_p(F)|_(j>=h)=lower-hull{(j,u_N(j)):h<=j<N}.              (16)
```

The point `(h,0)` is exact.  Every coefficient valuation to its left is
nonnegative, while every later plateau has positive height.  Any chord from
a prefix point to a later suffix point therefore lies strictly above
`(h,0)`, so this anchor cannot be bypassed.  Every possible prefix edge
ending there has nonpositive slope.  This is precisely where a
left-factorial/Kurepa issue may remain, and precisely why it is harmless for
the common ledger.

## 5. Binary suffix polygon

Now let `p=2`, write `N=2M`, and put

```text
K_q=2q+u_M(q).                                             (17)
```

For `0<=q<M`, the first two formulas below hold; for `0<=q<=M-2`, the
increment formula also holds:

```text
u_N(2q)=K_q,                   u_N(2q+1)=K_q+1,
K_(q+1)-K_q=3+nu_2(M-1-q)>0.                               (18)
```

If `j` is odd, `N-1-j` is even.  The `ell=1` term in (5) is killed by its
binomial factor and every `ell>=2` term by an even factor in the rising
product.  Thus

```text
Z_j==1 (mod 2).                                            (19)
```

Every odd point is exact.

If `j` is even, both `N-1-j` and `2j+1` are odd.  Modulo two the `ell=0`
and `ell=1` terms cancel, while every `ell>=2` term vanishes.  Hence

```text
Z_j==0 (mod 2),                                            (20)
```

unless the coefficient itself is zero, and in either case the actual even
point is absent or has height at least `K_q+1`, equal to the exact odd point
immediately to its right.  Since the odd anchors rise strictly by (18), no
even index after zero enters the lower hull.  Therefore

```text
NP_2(F)|_(j>=1)=lower-hull{(2q+1,K_q+1):2q+1<N}.            (21)
```

At `j=0`, (20) gives height at least one, while all later odd anchors have
height strictly greater than one.  Thus `(1,1)` cannot be bypassed and the
possible prefix edge into it again has nonpositive slope.

## 6. Why the prefix cannot meet G

The constant coefficient of `G` is a `p`-adic unit.  Moreover

```text
w_N(j)>0                         for every 0<j<=N.          (22)
```

Indeed, if `2j>=p`, the factorial term is positive; if `2j<p`, then
`0<j<p` and `p|N` makes `C(N,j)` divisible by `p`.  Thus every lower edge
of `G` has positive slope.  Sections 4--5 show that the only uncompiled
part of `F` has nonpositive slope.  It cannot contribute a common slope.

Intersecting the exact digit hull (3) with (16) or (21) therefore gives all
and only the common slopes, with their exact common capacities and reduced
denominators.  This proves the compiler.  QED.

## 7. Exact verification and the 2501--2600 divisor census

The primary companion hash-pins THM-3152's exact coefficient engine and
THM-3161's digital skeleton.  It checks:

```text
odd:    every odd p|N, 3<=N<=160,             209 pairs;
binary: every even 2<=N<=200,                 100 cells;
pair:   digit versus actual blocks/barcodes,  309 ledgers.             (23)
```

Every odd plateau recursion, both endpoint-unit congruences, every binary
parity congruence, and every actual positive `F`-hull agree with the theorem.
The respective semantic digests are

```text
odd     1fe62de92a499ecaa9eefd4aaca90ae77330f4857fd33103f15f9b848c116142,
binary  c1a56e649a584069678c9e275a8a28b69d80e359f27cbc64aef36f078daaeed8,
pair    02915c06dcb6569bbd7a19193e38b964a4e35320ec3368ff6177958492d5c202. (24)
```

The binary audit includes the exact zero coefficient at `(N,j)=(4,2)`.
The hostile `(N,p)=(6,3)` has `G` blocks of slopes `2/3` and `1`, but the
pair shares only the `2/3` block; this refutes any unclipped-principal-block
shortcut.  A raised-prefix hostile changes only the nonpositive `F` prefix
and leaves its positive suffix unchanged.  Planted factors `v` and `v+1`
retain degree one.

For the inherited seven-exit universe `2501<=d<=2600`, the primary compiler
intersects the exact pair barcode over every distinct prime dividing
`N=d-1`.  It closes `33` of the `38` residual rows.  The exact survivors are

```text
d=2516: {503,1006,1509,2012};
d=2564: {466,699,1165,1631,1864,2097,2330};
d=2571: {2056};
d=2576: {103,206,...,2472}  (all 24 positive multiples of 103);
d=2586: {47,141,188,235,282,329,
         2209,2256,2303,2350,2397,2444,2491,2538}.         (25)
```

Its full-record digest is

```text
27c1c59959a58054a3420d0fed7944da7a3ce480ad91a00e9e753862af83efc7. (26)
```

A second script imports no theorem engine: it recomputes valuations by
Legendre/Kummer formulas, uses an integer-cross-product hull and its own
denominator dynamic program, reconstructs the seven-exit counts
`(89,78,69,60,52,44,38)`, and returns the same `33/5` split and the literal
sets (25).  Its independently shaped semantic digest is

```text
e019eb61019620cfaffa8b5bb5769e8d171f08fb5ffeb2153f209c0128d42115. (27)
```

Together with the inherited seven exits and THM-3467, the first unresolved
row is now `d=2516`, or `r=d-2=2514`.  Hence every exact-support quadratic
three-moment window beginning at

```text
1<=r<=2513                                                  (28)
```

has a nonzero moment.  The five rows in (25) are not failures of the theorem;
they are exact local survivor packets requiring a nondivisor place or another
sidecar.

Reproduce with

```text
python3 04-computation/factorial_all_divisor_digit_pair_compiler_thm3475.py
python3 04-computation/factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py
```

and repeat both commands with `python3 -O`.  Each output is byte-identical to
its declared transcript, and neither script contains a Python `assert` gate.

## 8. Scope and failure boundaries

- The divisor condition `p|N` is load-bearing in both digit recurrences.
- The pair compiler does not say that every allowed degree is realized by a
  factor.  It is a necessary local obstruction only.
- A surviving pair degree may require another divisor place, a nondivisor
  place, or a Euclidean row to eliminate it.
- A hypothetical failure of the relevant left-factorial nonvanishing changes
  only the nonpositive prefix and does not change the common ledger.
- The old unclipped principal-capacity rule is false already at `(N,p)=(6,3)`;
  the full digit hull, not one preferred edge, is essential.
- The theorem concerns the exact-support quadratic resonance pair.  It does
  not imply arbitrary-support `SFC(3)` or the Factorial Conjecture.

The source-to-target map is now explicit: THM-3161's exact digital polygon
for `G` is paired with the present exact suffix polygon for `F`; their common
native Newton edges produce the local degree barcode.  The destroyed data are
factor existence and global compatibility.  The retained sidecars are the
prime label, full edge capacities, reduced denominators, and the unresolved
prefix sign, rather than the large coefficient arrays.

**QED.**
