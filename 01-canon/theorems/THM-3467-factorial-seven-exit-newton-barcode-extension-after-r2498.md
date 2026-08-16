---
id: THM-3467
title: "Prime-power carry ledgers and adaptive exact closure through r=2500"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  For p>=5, d=a p^k+1 with 2<=a<p, the two resonant factorial
  moment rows and their first full Euclidean row have an explicit common
  Newton face of reduced denominator p^k and exact common capacity
  min(a-1,(p-1)/2)p^k.  This face is the complete common p-adic ledger for
  every 2<=a<p.  At d=2501 the complete 2-, 3-, and 5-adic
  degree barcodes are incompatible.  At d=2502 the 61-adic singleton
  reset and THM-3161's 41-adic digital polygon force a hypothetical common
  factor degree to be divisible by 2501, beyond the Euclidean-row degree.
  These close the first two residuals beyond THM-3201 and extend the exact
  quadratic three-moment boundary through r=2500.  This does not close every
  height or FC(3).
audit: >
  The carry-face proof derives coefficientwise supporting inequalities from
  Lucas, Kummer, Legendre, and base-p digit sums, proves exactness at every
  anchor, and proves that one term uniquely controls the first Euclidean row.
  In the small-multiplier range, G's proved face spans its full degree.  In
  the remaining range, two more digit-sum supporting lines give the complete
  G and F tails, whose slopes are strictly separated.
  A 48-cell exact hostile audit found the predicted face to be the complete
  common ledger in every tested cell.  At d=2501 and again at d=2502,
  separately pinned Fraction-hull and determinant-hull engines reconstruct
  the rows, agree on every raw local degree ledger and common-slope ledger,
  retain planted factors v and v+1, and produce the declared traces and
  digests.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3201-seven-exit-factorial-newton-euclidean-closure-through-r2403
related:
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
scripts:
  - 04-computation/factorial_seven_exit_extension_2501_thm3467.py
  - 04-computation/factorial_adaptive_divisor_extension_2502_thm3467.py
outputs:
  - 05-knowledge/results/factorial_seven_exit_extension_2501_thm3467.out
  - 05-knowledge/results/factorial_adaptive_divisor_extension_2502_thm3467.out
script_sha256:
  - 8ba583ae08a91bfd70dd9911fca927ba75e7376e7b49887c3915ac3c905e02e3
  - 8964f7366329c469366cd5bd18ec7a176a9aac5c3bf074c633e9dd27d4067a23
output_sha256:
  - 43ef5f55c4b88bc27fca312f8512879f906acfd85d8efd99af2deee1821143c8
  - f87c6b4dd1765bcc7a2aabecf2e3ca885ded851c8524794d1175eb39f0ae81eb
hash_basis: LF-normalized bytes
---

# THM-3467 -- prime-power carry ledgers and closure through r=2500

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

For every `2<=a<p`, Section 7 proves that this is the only common finite
slope.  Thus the complete local address is

```text
D_p(F,G,E)={0,H,2H,...,T H}.                               (9a)
```

The case `a=1` remains the zero-capacity degeneration.

The first finite-exact corollary is

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

contains a nonzero moment.  Section 9 proves a second, adaptive-divisor
corollary at `d=2502`, extending (11) through `r=2500`.  This remains a
bounded support-`{0,1,2}` result; it is not arbitrary one-variable support,
`SFC(3)`, or the three-variable Factorial Conjecture.

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

### 7.2 The remaining multipliers m<a<p

Assume now `m<a<p` and put

```text
S=(H-1)/(p-1),             q=2a-p,
j_0=mH,                    j_*=(pH-1)/2,
tau=2/(p-1),
rho_G=2(qS+1)/(qH+1),      rho_F=2(qS+1)/(qH-1).           (34)
```

Here `q` is odd and `1<=q<=p-2`.  The proved `sigma` face of `G` ends at

```text
(j_0,H-1).                                                  (35)
```

The global bound

```text
nu(g_j)>=tau j-1                                            (36)
```

follows from Kummer and Legendre exactly as in Section 3: if
`2j=cH+r`, then `s_p(c)<=p-1`, and the lower `k-nu(j)` digits contribute at
most `(p-1)(k-nu(j))`.  Equality is exact at `j_0` and `j_*`, with heights
`H-1` and `pS`; at `j_*`, the first rising factor is
`2j_*+1=pH`, so `Z_(N,j_*)` is a unit.  This gives the `tau` edge

```text
(j_0,H-1) -> (j_*,pS).                                     (37)
```

It remains to control the tail.  For `j>j_*`, write

```text
2j=pH+v,                    1<=v<=qH,
c=floor(v/H),               h=nu_p(v).                      (38)
```

The coefficient lower skeleton satisfies

```text
w_G(j)-pS=k-h+1+(v-s_p(v))/(p-1)
          >=1+(v-c)/(p-1).                                 (39)
```

Set

```text
D=p-q,              M=qH+1,              A=qS+1,
(p-1)A=qH+D-1.                                             (40)
```

The chord from `(j_*,pS)` to the leading point has rise
`A(v+1)/M`.  After multiplying the difference between (39) and that chord by
the positive number `(p-1)M`, it is bounded below by

```text
Phi=M(p-1+v-c)-(qH+D-1)(v+1).                              (41)
```

If `c<=q-1`, the coefficient of `v` in `Phi` is `2-D<=0`, so its minimum on
that digit block is at the right endpoint.  With `u=q-1-c`, the endpoint
value is

```text
Hq+1+u(H(p-2)+1)>0.                                        (42)
```

For `c=q`, necessarily `v=qH`, the leading endpoint.  There

```text
g_N=(2N)!,             nu(g_N)=(p+q)S+1,                   (43)
```

so equality holds.  Hence the complete `G` polygon has slopes

```text
sigma, tau, rho_G,                                           (44)
```

with the last two merging when `q=p-2`.

For `F`, the `sigma` face ends at `(j_*,pS-k)`.  Its leading point has height
`pS-k+qS+1`.  At an interior tail index `2j=pH+v`,
`1<=v<=qH-2`, the digit bound

```text
s_p(v)<=k(p-1)+q-2                                         (45)
```

holds: if `c<=q-2`, bound the leading digit by `q-2`; if `c=q-1`, use
`v<=(qH-2)` to save one unit in the lower digits.  Since
`s_p(pH-1)=(k+1)(p-1)`, Legendre gives

```text
nu((2j)!)-nu((pH-1)!) >= (v+p-q+1)/(p-1).                 (46)
```

The chord increase of slope `rho_F` is `(qS+1)(v+1)/(qH-1)`.  The
cross-multiplied difference between (46) and this chord is exactly

```text
(p-q)(qH-2-v),                                             (47)
```

which is positive before the endpoint.  The leading coefficient of `F` is
exactly `(2N-2)!`, so `F` has one post-`sigma` edge, of slope `rho_F`.

Finally,

```text
sigma<tau<=rho_G<rho_F.                                    (48)
```

By convexity, every pre-`sigma` slope of `F` is below `sigma`; its only
post-`sigma` slope is `rho_F`.  The complete `G` slope list is (44), so the
two rows share only `sigma`.  The row `E` retains the `sigma` face proved in
Section 6.  Therefore, for every `2<=a<p`, not merely the two boundary
ranges,

```text
common ledger=((sigma, capacity T H, denominator H)),
D_p(F,G,E)={0,H,2H,...,T H},       T=min(a-1,m).            (49)
```

The coordinate-root capacity is zero because `g_0` is a unit.  This proves
uniform singleton completeness.  QED.

## 8. Exact d=2501 corollary

For `p=5`, `k=4`, `H=625`, `a=4`, the structural theorem gives

```text
B=312,          U=312,          V=152,          T=2,
sigma=312/625.                                                (50)
```

The predicted face endpoints are reconstructed independently as

```text
F: (312,152)  -> (1562,776),
G: (0,0)      -> (1250,624),
E: (313,156)  -> (1563,780).                               (51)
```

The complete local common-slope ledgers are

```text
p=2: (511/256, cap 256,  den 256),
     (4095/2048,cap 2048,den 2048);

p=3: (2/3,cap 3,den 3), (8/9,cap 9,den 9),
     (26/27,cap 27,den 27), (242/243,cap 243,den 243),
     (2186/2187,cap 2187,den 2187), (1,cap 27,den 1);

p=5: (312/625,cap 1250,den 625).                           (52)
```

In particular, the raw positive degree barcodes at `p=2` and `p=5` are

```text
D_2^+={256,2048,2304},              D_5^+={625,1250}.      (53)
```

The progressive intersection is

```text
after p=2: {256,2048,2304},
after p=3: {256},
after p=5: empty.                                          (54)
```

By THM-3152's degree-barcode lemma, `F,G,E` have no nonconstant common
rational factor.  Since `gcd(F,G)=gcd(F,E)`, the resonant pair is coprime.
THM-3124 converts this back to the three-moment statement, proving
(10)--(11).  QED.

The reason for (54) is bounded-capacity incompatibility.  Coprime slope
denominators alone do not force disjoint degree sets: sums and products of
their blocks may still meet.  Completeness of each local ledger is the
load-bearing sidecar.

The exact seven-exit factorization invoice is

```text
2501=41*61,              2500=2^2*5^4,
2499=3*7^2*17,           2498=2*1249,
2497=11*227,             2496=2^6*3*13,
2495=5*499.                                                  (55)
```

Thus `d=2501` was genuinely the first residual after THM-3201's seven
uniform exits, rather than a row already covered by one of them.

## 9. Adaptive divisor closure at d=2502

Put

```text
d=2502,        N=d-1=2501=41*61,
F=A_2500^(2502),       G=A_2501^(2502),       deg E=2499. (55a)
```

At `p=61`, use `H=61` and `a=41` in the uniform theorem above.  Since
`T=min(40,30)=30`, the complete common ledger is

```text
(slope 2/61, capacity 1830, denominator 61),
D_61(F,G,E)={0,61,122,...,1830}.                           (55b)
```

Thus every positive common-factor degree is divisible by `61`.

At `p=41`, THM-3161's exact digital polygon for the single row `G` has the
two maximal blocks

```text
(slope 2/41,    capacity 820,  denominator 41),
(slope 84/1681, capacity 1681, denominator 1681).          (55c)
```

Every factor degree of `G` is consequently divisible by `41`: the first
block contributes a multiple of `41`, while the second contributes either
zero or `1681=41^2`.  A common factor of `F` and `G` must therefore have
degree divisible by

```text
lcm(41,61)=2501.                                          (55d)
```

But it also divides the Euclidean combination `E`, whose degree is `2499`.
No positive degree is possible, so `gcd_Q(F,G)=1`.  By THM-3124, every exact
quadratic window beginning at `r=2500` has a nonzero member.  Together with
THM-3201 and Section 8, every exact-support quadratic window beginning at

```text
1<=r<=2500                                                 (55e)
```

contains a nonzero moment.  This proof is structural: the `41`-divisibility
comes from THM-3161 and the `61`-divisibility from the all-`a` reset theorem,
not from the finite replay.

The hostile explains why both divisor places matter.  Starting with the
`61` reset and intersecting only the old small observers `p=2,3` leaves the
spurious address `{61}`.  Adding the divisor place `p=41` removes it.  The
same row is also closed finite-exactly by the fixed trace `p=2,3,5`; that is
an independent control, not the mechanism of the proof.

## 10. Audit, controls, and failure boundaries

The `d=2501` companion hash-pins the primary THM-3201 Fraction-hull engine and
its independent determinant-hull referee.  Both reconstruct `F,G,E`, whose exact
degrees are `(2499,2500,2498)`, and require literal equality of every raw
degree ledger, common-slope ledger, and progressive trace.  The semantic
trace digest is

```text
365533925519a4d8d44db78394f0785e87be5f4cc03e0a98d759f93609fb09ee.
                                                                    (56)
```

Planted common factors `v+1` and `v` retain degree `1` in both observers.
The script contains no Python `assert`, so optimized replay checks the same
conditions.

An exact 48-cell audit used

```text
(p,k)=(5,1..3),(7,1..2),(11,1..2),(13,1),
2<=a<p.                                                       (57)
```

In every cell the complete common ledger is the singleton predicted by
(4)--(8), independently supporting the uniform proof in Section 7.

The `d=2502` companion independently reconstructs the three rows with the
same two engine families and requires equality of all raw profiles and local
degree sets at `p=2,3,5,41,61`.  It freezes the full `41`- and `61`-adic
profiles, the hostile survivor `{61}`, the empty adaptive trace, both planted
degree-one controls, and semantic digest

```text
c767afd684dbf63910db607a742c530ed9f03a8d61c97b6e2e4dd5cb22bddf98. (57a)
```

Normal and optimized replays are byte-identical to the stored transcript.

The hypotheses in the structural statement are sharp for the proof:

- at `p=3`, `alpha` and `beta` in (29) cease to be units; for `d=55`, the
  same slope `26/27` occurs but the `E` face begins at `B`, not `B+1`;
- `p=2` has no nontrivial one-digit range `2<=a<p`;
- dropping `a<p` introduces additional quotient digits and carry faces:
  `(p,k,a,d)=(5,1,6,31)` has common slope `12/25`, capacity `25`, and
  denominator `25`, not `5`;
- when `p` divides `a`, `k` is not the exact reset depth; at `a=5,d=26`
  the first-flag common ledger is empty.

The tail analysis in Section 7 is load-bearing beyond the first face.  For
`(p,k,a)=(5,1,3)`, `G` has additional slopes `1/2` and `2/3`; the proof must
show that neither occurs in `F`, rather than treating the `sigma` face alone
as a complete ledger.

Finally, the theorem produces a necessary local degree address, not a
factor, and its two boundary closures do not supply a contiguous all-height
prime bank.  The first next seven-exit residual is `d=2510`, `r=2508`.
