---
id: THM-3467
title: "Prime-power carry faces in the first factorial Euclidean flag and exact closure at r=2499"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For p>=5, d=a p^k+1 with 2<=a<p, the two resonant factorial
  moment rows and their first full Euclidean row have an explicit common
  Newton face of reduced denominator p^k and exact common capacity
  min(a-1,(p-1)/2)p^k.  This face is the complete common p-adic ledger
  whenever 2<=a<=(p-1)/2 or a=p-1.  At d=2501 the complete 2-, 3-, and 5-adic
  degree barcodes are incompatible, closing the first residual beyond
  THM-3201 and extending the exact quadratic three-moment boundary to
  r=2499.  This is not an all-height theorem or FC(3).
audit: >
  The carry-face proof derives coefficientwise supporting inequalities from
  Lucas, Kummer, Legendre, and base-p digit sums, proves exactness at every
  anchor, and proves that one term uniquely controls the first Euclidean row.
  In the small-multiplier range, G's proved face already spans its full
  degree, immediately forcing singleton completeness.  A second pair of
  supporting lines gives the complete two-edge polygon of G when a=p-1,
  while a unique digit-sum maximum excludes G's second slope from F.
  A 48-cell exact hostile audit found the predicted face to be the complete
  common ledger in every tested cell, but completeness is not claimed
  uniformly.  At d=2501 two separately pinned Fraction-hull and
  determinant-hull engines reconstruct the rows, agree on every raw local
  degree ledger and common-slope ledger, retain planted factors v and v+1,
  and produce the same progressive trace and digest.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3201-seven-exit-factorial-newton-euclidean-closure-through-r2403
related:
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
script: 04-computation/factorial_seven_exit_extension_2501_thm3467.py
output: 05-knowledge/results/factorial_seven_exit_extension_2501_thm3467.out
script_sha256: 8ba583ae08a91bfd70dd9911fca927ba75e7376e7b49887c3915ac3c905e02e3
output_sha256: 43ef5f55c4b88bc27fca312f8512879f906acfd85d8efd99af2deee1821143c8
hash_basis: LF-normalized bytes
---

# THM-3467 -- prime-power carry faces and closure at r=2499

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and scope

Let

```text
L(x^r)=r!,                       A_n(v)=L((d-x+v x^2)^n).    (1)
```

Fix a prime `p>=5`, integers `k>=1` and `1<=a<p`, and put

```text
H=p^k,       N=aH,       d=N+1,
F=A_(N-1),   G=A_N.                                        (2)
```

The first full Euclidean row of THM-3152 is

```text
E=(2N-3)(G-2N(2N-1)vF)+2N(N+1)F.                          (3)
```

Set

```text
m=(p-1)/2,                 B=(H-1)/2,
U=2(H-1)/(p-1),            sigma=U/H,
V=(H-1)/(p-1)-k,
T=min(a-1,m),              T_G=min(a,m).                   (4)
```

Writing `f_j`, `g_j`, and `e_j` for the coefficients of `v^j` in
`F`, `G`, and `E`, respectively, and writing `nu=nu_p`, one has globally

```text
nu(g_j) >= sigma j,
 equality exactly at j=tH,             0<=t<=T_G,
 nu(g_(tH))=tU;                                             (5)

nu(f_j) >= V+sigma(j-B),
 equality exactly at j=B+tH,           0<=t<=T,
 nu(f_(B+tH))=V+tU;                                        (6)

nu(e_j) >= V+k+sigma(j-B-1),
 equality exactly at j=B+1+tH,         0<=t<=T,
 nu(e_(B+1+tH))=(H-1)/(p-1)+tU.                            (7)
```

Here a zero coefficient has valuation infinity.  Consequently, when
`a>=2`, the three lower Newton polygons have a common face of slope
`sigma`, and its exact common horizontal capacity is

```text
T H=min(a-1,(p-1)/2)p^k.                                  (8)
```

Since

```text
U=2(1+p+...+p^(k-1)) == 2 (mod p),                         (9)
```

`sigma` is reduced and its denominator is exactly `H=p^k`.  When `a=1`,
`T=0`: the `G` face remains, but `F` and `E` have no positive-length common
`sigma` face.

The theorem proves this explicit face for every displayed parameter and proves
that it is the only common face when `2<=a<=m` or `a=p-1`.  It does **not**
assert uniform completeness for the remaining middle multipliers
`m<a<p-1`; that stronger statement requires a separate ledger audit.

The finite-exact corollary is

```text
d=2501,                         r=d-2=2499.                 (10)
```

The complete local barcodes of `F,G,E` at `p=2,3,5` have empty positive
intersection.  Hence `F` and `G` are coprime over `Q`, so no exact quadratic
factorial window beginning at `r=2499` has three consecutive zero moments.
Together with THM-3201, every such window beginning at

```text
1<=r<=2499                                                   (11)
```

contains a nonzero moment.  This remains a bounded support-`{0,1,2}` result;
it is not arbitrary one-variable support, `SFC(3)`, or the three-variable
Factorial Conjecture.

## 2. Coefficient formula

For `0<=j<=n`, direct expansion gives

```text
[v^j]A_n
 = C(n,j) sum_(ell=0)^(n-j)
     C(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!
 = C(n,j)(2j)! Z_(n,j),                                   (12)

Z_(n,j)=sum_(ell=0)^(n-j)
 C(n-j,ell)d^(n-j-ell)(-1)^ell(2j+1)_ell in Z.             (13)
```

The last factor is an integer, so it may raise a valuation but never lower
one.  Lucas' theorem gives

```text
nu C(N-1,j)=0                         (0<=j<=N-1),          (14)
```

because the base-`p` digits of `N-1` are
`(a-1,p-1,...,p-1)`.  Kummer's theorem gives, for `0<j<=N`,

```text
nu C(N,j)=k-nu(j).                                         (15)
```

Thus

```text
nu(f_j)>=nu((2j)!),
nu(g_j)>=k-nu(j)+nu((2j)!).                                (16)
```

## 3. The G supporting line

For `j>0`, put `h=nu(j)` and write

```text
2j=cH+s,                         0<=s<H.                    (17)
```

If `s>0`, then `h=nu(s)<k` and

```text
s_p(2j)=s_p(c)+s_p(s)
       <=c+(p-1)(k-h)
       < 2j/H+(p-1)(k-h).                                  (18)
```

If `s=0`, then `c=2t`, `h=k`, and

```text
s_p(2j)=s_p(c)<=c,                                         (19)
```

with equality exactly when `c<p`, equivalently `0<=t<=m`.
Legendre's formula

```text
nu((2j)!)=(2j-s_p(2j))/(p-1)                              (20)
```

combined with (15), (18), and (19) proves the inequality in
(5), with possible equality exactly at `j=tH`, `0<=t<=min(a,m)`.
The case `t=0` is included by the exactness argument below.

## 4. The F supporting line

The proposed affine bound can be rewritten as

```text
V+sigma(j-B)
 =(2j+1)(H-1)/((p-1)H)-k.                                 (21)
```

By (16) and Legendre, it is enough to prove

```text
s_p(2j)+1 <= k(p-1)+(2j+1)/H.                              (22)
```

Write `2j+1=cH+s`, `0<=s<H`.  If `s=0`, then

```text
2j=(c-1)H+(H-1),
s_p(2j)=s_p(c-1)+k(p-1)<=c-1+k(p-1).                      (23)
```

Equality holds exactly when `c<=p`.  Since `c` is odd, write
`c=2t+1`; the equality indices are `j=B+tH` with `t<=m`, and the
degree bound imposes `t<=a-1`.

If `s>0`, then

```text
2j=cH+(s-1),
s_p(s-1)<=k(p-1)-1,
s_p(2j)+1<=c+k(p-1)<c+k(p-1)+s/H.                          (24)
```

This proves (6), apart from exactness at the displayed anchors.

## 5. Exactness at the F and G anchors

For `G` at `j=tH`, the factor `n-j=(a-t)H`.  In (13), every term
with `1<=ell<p` has a binomial coefficient divisible by `p`, while every
term with `ell>=p` has `(2tH+1)_ell` divisible by `p`.  Hence

```text
Z_(N,tH)==1 (mod p),
C(N,tH)==C(a,t) !=0 (mod p).                               (25)
```

For `F` at `j=B+tH`, one has

```text
2j+1=(2t+1)H.                                              (26)
```

Thus every `ell>=1` summand of (13) is divisible by `p`, while (14)
makes the binomial factor a unit.  Therefore `Z_(N-1,j)==1 (mod p)`.
Legendre now gives exactly

```text
nu(g_(tH))=tU,
nu(f_(B+tH))=V+tU.                                        (27)
```

This completes (5) and (6).

## 6. Survival through the Euclidean projection

Coefficientwise, (3) is, with `f_(-1)=0`,

```text
e_j=alpha g_j+N beta f_(j-1)+N gamma f_j,                  (28)

alpha=2N-3,
beta=-2(2N-3)(2N-1),
gamma=2(N+1).                                               (29)
```

Modulo `p`, these three scalars are `-3,-6,2`, respectively, so they are
units when `p>=5`.  Put

```text
Lambda_E(j)=V+k+sigma(j-B-1).                              (30)
```

Equations (5)--(6) imply

```text
nu(N beta f_(j-1)) >= Lambda_E(j),
nu(N gamma f_j)     >= Lambda_E(j)+sigma,
nu(alpha g_j)       >= Lambda_E(j)+delta,                  (31)

delta=(H-1)/((p-1)H)>0.                                   (32)
```

At `j=B+1+tH`, `0<=t<=T`, the first term of (28) has exact integral
valuation `Lambda_E(j)`, while the other two terms have strictly larger
valuation.  It therefore dominates uniquely and cannot cancel.  Off those
indices all three terms are strictly above the line.  This proves (7).

The last `E` anchor is within `deg E=N-2`.  If `T=a-1`, its gap below that
degree is `(H-5)/2`; if `T=m`, the gap is

```text
((2a-p)H-5)/2 >= (H-5)/2 >=0.                              (33)
```

The horizontal lengths of the three faces are `TH`, `T_G H`, and `TH`.
Their exact common capacity is therefore (8), completing the structural
proof.  QED.

## 7. Complete ledgers in two multiplier ranges

### 7.1 Small multipliers

Suppose first

```text
2<=a<=m=(p-1)/2.
```

Then `T_G=a`, so the equality face (5) of `G` runs from `(0,0)` to
`(aH,aU)`.  Since `deg G=aH`, this one segment spans the entire coefficient
range and is the whole finite lower Newton polygon of `G`.  Hence `G` has no
finite slope other than `sigma`.  Both `F` and `E` contain their proved
`sigma` faces of common width `(a-1)H`, and `g_0` is a unit.  Therefore the
complete triple ledger and local degree address are

```text
((sigma, capacity (a-1)H, denominator H)),
D_p(F,G,E)={0,H,2H,...,(a-1)H}.
```

No claim about the other slopes of `F` or `E` is needed: intersecting with
the single-slope polygon of `G` already proves completeness.  QED.

### 7.2 The endpoint multiplier a=p-1

There is one uniform completeness cell.  Suppose now

```text
a=p-1,                 m=(p-1)/2,                 N=(p-1)H,
tau=2/(p-1).                                                (34)
```

The `sigma` face of `G` from (5) runs from `(0,0)` to
`(mH,mU)=(mH,H-1)`.  A second global supporting inequality follows from
(15)--(20).  Write `2j=cH+s`, `0<=s<H`.  On the full range
`0<=j<=N`,

```text
s_p(c)<=p-1,              s_p(s)<=(p-1)(k-nu(j)),
```

where `j=0` is treated by its unit anchor.  Therefore

```text
nu(g_j)>=tau j-1.                                           (35)
```

Equality holds at `j=mH` and `j=N`: the first is the last `sigma` anchor,
while

```text
g_N=(2N)!,                  nu(g_N)=2H-1=tau N-1.           (36)
```

Thus `G` has a `tau` face from `(mH,H-1)` to `(N,2H-1)`.  Since
`sigma<tau`, the two supporting lines meet at `mH`; the `sigma` line is the
larger one to the left and the `tau` line is the larger one to the right.
Every coefficient point lies above both, and the displayed anchors lie on
their piecewise maximum.  This two-segment convex chain covers the entire
coefficient range, so these are all finite lower Newton slopes of `G`.

The second slope is absent from `F`.  From (14), (16), and Legendre,

```text
nu(f_j)-tau j >= -s_p(2j)/(p-1).                            (37)
```

On

```text
0<=2j<=2(p-1)H-2,
```

the digit sum has the unique maximum `(k+1)(p-1)` at

```text
2j=pH-1,                 j=B+mH.                            (38)
```

Indeed, below `pH` equality requires all `k+1` digits to be `p-1`.  Above
`pH`, write `2j=pH+r` with `0<=r<=(p-2)H-2`; then
`1+s_p(r)<=(k+1)(p-1)-1`.  At the unique index (38), (27) is exact and gives

```text
nu(f_j)-tau j=-(k+1).                                      (39)
```

Hence the affine functional `y-tau x` has a unique minimum on the
coefficient points of `F`.  No lower Newton edge of `F` can have slope
`tau`.  Since all slopes of `G` are `sigma,tau`, while `F` and `E` both
have the proved `sigma` face, their complete common finite-slope ledger is

```text
((sigma, capacity mH, denominator H)).                     (40)
```

The common coordinate-root capacity is zero because `g_0` is a unit.
Therefore

```text
D_p(F,G,E)={0,H,2H,...,mH}.                                (41)
```

This endpoint-family completeness is structural; no bounded hull scan is
used in its proof.  QED.

## 8. Exact d=2501 corollary

For `p=5`, `k=4`, `H=625`, `a=4`, the structural theorem gives

```text
B=312,          U=312,          V=152,          T=2,
sigma=312/625.                                                (42)
```

The predicted face endpoints are reconstructed independently as

```text
F: (312,152)  -> (1562,776),
G: (0,0)      -> (1250,624),
E: (313,156)  -> (1563,780).                               (43)
```

The complete local common-slope ledgers are

```text
p=2: (511/256, cap 256,  den 256),
     (4095/2048,cap 2048,den 2048);

p=3: (2/3,cap 3,den 3), (8/9,cap 9,den 9),
     (26/27,cap 27,den 27), (242/243,cap 243,den 243),
     (2186/2187,cap 2187,den 2187), (1,cap 27,den 1);

p=5: (312/625,cap 1250,den 625).                           (44)
```

In particular, the raw positive degree barcodes at `p=2` and `p=5` are

```text
D_2^+={256,2048,2304},              D_5^+={625,1250}.      (45)
```

The progressive intersection is

```text
after p=2: {256,2048,2304},
after p=3: {256},
after p=5: empty.                                          (46)
```

By THM-3152's degree-barcode lemma, `F,G,E` have no nonconstant common
rational factor.  Since `gcd(F,G)=gcd(F,E)`, the resonant pair is coprime.
THM-3124 converts this back to the three-moment statement, proving
(10)--(11).  QED.

The reason for (46) is bounded-capacity incompatibility.  Coprime slope
denominators alone do not force disjoint degree sets: sums and products of
their blocks may still meet.  Completeness of each local ledger is the
load-bearing sidecar.

The exact seven-exit factorization invoice is

```text
2501=41*61,              2500=2^2*5^4,
2499=3*7^2*17,           2498=2*1249,
2497=11*227,             2496=2^6*3*13,
2495=5*499.                                                  (47)
```

Thus `d=2501` is genuinely the first residual after THM-3201's seven
uniform exits, rather than a row already covered by one of them.

## 9. Audit, controls, and failure boundaries

The companion hash-pins the primary THM-3201 Fraction-hull engine and its
independent determinant-hull referee.  Both reconstruct `F,G,E`, whose exact
degrees are `(2499,2500,2498)`, and require literal equality of every raw
degree ledger, common-slope ledger, and progressive trace.  The semantic
trace digest is

```text
365533925519a4d8d44db78394f0785e87be5f4cc03e0a98d759f93609fb09ee.
                                                                    (48)
```

Planted common factors `v+1` and `v` retain degree `1` in both observers.
The script contains no Python `assert`, so optimized replay checks the same
conditions.

An exact 48-cell audit used

```text
(p,k)=(5,1..3),(7,1..2),(11,1..2),(13,1),
2<=a<p.                                                       (49)
```

In every cell the complete common ledger was the singleton predicted by
(4)--(8).  This is strong finite evidence for the sharper completeness
conjecture, not a proof of it.

The hypotheses in the structural statement are sharp for the proof:

- at `p=3`, `alpha` and `beta` in (29) cease to be units; for `d=55`, the
  same slope `26/27` occurs but the `E` face begins at `B`, not `B+1`;
- `p=2` has no nontrivial one-digit range `2<=a<p`;
- dropping `a<p` introduces additional quotient digits and carry faces:
  `(p,k,a,d)=(5,1,6,31)` has common slope `12/25`, capacity `25`, and
  denominator `25`, not `5`;
- when `p` divides `a`, `k` is not the exact reset depth; at `a=5,d=26`
  the first-flag common ledger is empty.

The endpoint condition `a=p-1` in Section 7 is also load-bearing for its
two-edge description.  For `(p,k,a)=(5,1,3)`, the tail of `G` has additional
slopes `1/2` and `2/3`; only the proved `sigma` face survives from the general
statement.

Finally, the theorem produces a necessary local degree address, not a
factor, and the finite closure at `d=2501` does not supply an all-height
prime bank.  The first next seven-exit residual is `d=2502`, `r=2500`.
