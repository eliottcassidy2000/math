---
id: THM-3176
title: "Six-step prime resonance third Euclidean--Newton separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every odd
  prime p, no exact quadratic has three zero factorial moments beginning at
  r=p+4.  Three Euclidean rows separate the first genuine two-row slope
  collision, a complete wall atlas handles every projected degeneration, and
  first-Witt Hensel wedges close all four offset-six zero-face resultant
  exceptions intrinsically.
audit: >
  Two independent immutable audits accepted the theorem.  Normal and
  optimized replay agree with the 33-line transcript and declared hashes;
  the direct coefficient sum agrees with the corrected full recurrence;
  symbolic quotient and endpoint identities, every wall factorization and
  primality check, the p=109 and p=232961 normalized polygons, and all four
  first-Witt closures were independently rederived.  The precanonical
  constant-term double-count was repaired before the candidate commit and is
  guarded by an exact direct-sum/full-recurrence equality control.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
  - THM-3175-first-witt-hensel-wedge-and-infinitesimal-pluecker-covariance
related:
  - THM-3170-five-step-prime-resonance-euclidean-newton-holotopy
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
script: 04-computation/factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.py
output: 05-knowledge/results/factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.out
script_sha256: a988112c15153b1925fe5c4fb18b7132fdd1eed998247850fc4b691fbc0c51f7
output_sha256: 71b73fdcf2247d23ef602d1879b55206846687a03c1891efa560f3beb975cd66
hash_basis: LF-normalized bytes
---

# THM-3176 -- six-step prime resonance third Euclidean--Newton separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L(t^j)=j!,                         Q(t)=a+bt+ct^2,
abc!=0,
r=p+4,                            d=r+2=p+6,              (1)
```

where `p` is an odd prime.  Then

```text
L(Q^r), L(Q^(r+1)), L(Q^(r+2))                            (2)
```

cannot all vanish.

The first three natural Newton rows have the same positive slope `2/p`.
A third Euclidean remainder is the first row which separates that collision.
Every arithmetic wall is explicit.  The fixed zero-face resultant has four
exceptional primes, and the first-Witt Hensel wedge excludes all four without
neighboring-prime transport.

## 2. Resonant pair

THM-3124 forces `b/a=-1/d` if `(2)` vanishes.  Divide by `a`, put
`v=d(c/a)`, and define

```text
A_n(v)=L((d-t+vt^2)^n),
A=A_(p+4),                         B=A_(p+5).              (3)
```

It suffices to prove `gcd(A,B)=1` over every `Q_p`.  The coefficient formula
is

```text
[v^j]A_n=binom(n,j) sum_(ell=0)^(n-j)
 binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.              (4)
```

THM-3148 gives the fixed residual rows

```text
A == 6F_4,
B == 6F_5                                             (mod p),

F_4=744+384v+864v^2-2880v^3+40320v^4,
F_5=4056+2160v+1440v^2+57600v^3
       -604800v^4+3628800v^5.                            (5)
```

For `p>=197`, the endpoint coordinates in `(5)` are units and the usual
Lucas/factorial split, with `p=2m+1`, gives

```text
v_p(A_j)>=1, 5<=j<=m;
v_p(A_j)>=2, m+1<=j<=p+4.                                (6)
```

Since `A_(p+4)=(2p+8)!` has valuation two,

```text
NP_p(A):(0,0)--(4,0)--(p+4,2),       positive slope 2/p. (7)
```

## 3. Three Euclidean rows

The full first quotient gives

```text
C=(2p+10)(2p+9),

R=(2p+7)B-[(2p+7)Cv-2(p+5)(p+6)]A.                      (8)
```

Put `H=24p+109` and

```text
D_2=(p+5)(2p+3)H^2,
N_1=-2(2p+5)(2p+7)(2p+3)H,
N_0=4(p+6)(2p+5)(28p+123),

S=D_2 A-(N_1v+N_0)R.                                    (9)
```

The top two coordinates cancel in each step, so

```text
deg R<=p+3,                         deg S<=p+2.           (10)
```

The quotient of `R` by `S` is the third collision.  Define

```text
J=256p^4-27648p^3-365600p^2-1528800p-2096649,

K=5120p^5-810240p^4-14447872p^3-92004672p^2
   -256323456p-265142241,

D_3=(2p-1)(2p+7)J^2,
P_1=-2(2p+1)(2p+3)(24p+109)(2p-1)J,
P_0=4(p+6)(2p+1)K,

T=D_3R-(P_1v+P_0)S.                                    (11)
```

Again the top two coordinates cancel exactly, so `deg T<=p+1`.
Every common factor of `A,B` divides `R,S,T`.

Modulo `p`, complete simplification gives

```text
A ==144(31+16v+36v^2-120v^3+1680v^4),
R ==144(3043-17940v-7500v^2-13080v^3),
S ==15120(-375143+3212382v-2795532v^2),
T ==970059888(-51108355843
               +392547973190v-174357330840v^2).          (12)
```

The valuation propagation through `(8)--(11)` is

```text
R: v_p(R_j)>=1 for j>=4,       >=2 for j>=m+2;
S: v_p(S_j)>=1 for j>=3,       >=2 for j>=m+3;
T: v_p(T_j)>=1 for j>=3,       >=2 for j>=m+4.           (13)
```

Away from coordinate walls, the first two remainders have

```text
NP_p(R):(0,0)--(3,0)--(p+3,2),
NP_p(S):(0,0)--(2,0)--(p+2,2),                            (14)
```

so `A,R,S` all retain the same positive slope `2/p`.  This exact collision
is why two Euclidean rows do not suffice at offset six.

## 4. Midpoint and high sidecars for the third row

Put

```text
ell=m+3=(p+5)/2.                                         (15)
```

Nested use of `(4),(8),(9),(11)` and Wilson's theorem gives

```text
T_ell/p ==
 -474570855293572800 (-1)^m 6^(m+5)
 /[m(m-1)(m-2)(m-3)(m-4)]                       (mod p), (16)

474570855293572800=2^6*3^7*5^2*7^2*109^2*232961.         (17)
```

The endpoint and penultimate coordinates are

```text
T_(p+1)= -32 p(p-3)(p-2)(p-1)(p+1)(p+2)(p+3)(p+4)(p+5)
 (2p-7)(2p-5)(2p+3)(2p+7)(24p+109)^2
 H_8(p)(2p-8)!,                                          (18)

T_p= 16 p(p-3)(p-2)(p-1)(p+1)(p+2)(p+3)(p+4)(p+5)
 (2p-7)(2p+3)(2p+7)(24p+109)^2
 H_9(p)(2p-8)!,                                          (19)
```

where

```text
H_8=327680p^8+91422720p^7-1525587968p^6-58319874048p^5
 -605420542720p^4-3106619397120p^3-8698849881696p^2
 -12762150724080p-7697077854849,

H_9=327680p^9+87818240p^8-2581766144p^7-38490294272p^6
 +127672001792p^5+4526613136640p^4+30758409758112p^3
 +98688343652784p^2+157343795925183p+100505677023531.
                                                                    (20)
```

The load-bearing constant factorizations are

```text
H_8(0)=-3^4*89*1067703961,
H_9(0)= 3^4*13*17*52511*106921.                           (21)
```

## 5. Complete third-row wall atlas

For `p>=197`, outside the walls below, `(12)--(21)` force

```text
NP_p(T):(0,0)--(2,0)--(ell,1)--(p+1,2),
positive slopes             2/(p+1), 2/(p-3).            (22)
```

Both are distinct from `2/p`.

The low residual in `(12)` has scalar factor

```text
970059888=2^4*3^6*7*109^2                                (23)
```

and normalized coordinate factors

```text
51108355843                                      (constant),
2*5*22447*1748777                                  (linear),
2^3*3^5*5*7*11*232961                           (quadratic). (24)
```

Thus:

- at `p=22447,1748777`, the middle coordinate rises, but the constant and
  quadratic coordinates stay units; the horizontal edge and `(22)` remain;
- at `p=51108355843`, the constant rises while the linear and quadratic
  coordinates stay units; the polygon either starts at `x=1` or has an
  initial negative-slope edge, and the positive slopes in `(22)` remain;
- at `p=232961`, the quadratic coordinate and midpoint rise together.  A
  direct coefficient-sum computation modulo `p^4` gives, for
  `ell=(p+5)/2`,

  ```text
  T_(ell-1)/p=121373,             T_ell/p^2=136986 (mod p),

  NP_p(T):(0,0)--(1,0)--((p+3)/2,1)--(p+1,2),
  positive slopes             2/(p+1),2/(p-1).            (25)
  ```

  These are again distinct from `2/p`.

The endpoint wall is

```text
p=1067703961,
H_8(p)/p=114714544,               H_9(p)=567766679 (mod p). (26)
```

It has

```text
NP_p(T):(0,0)--(2,0)--(ell,1)--(p,2)--(p+1,3),
positive slopes             2/(p+1),2/(p-5),1.           (27)
```

The penultimate walls `p=52511,106921` only raise a nonvertex while the
endpoint remains at height two, so `(22)` is unchanged.  The quotient-
numerator wall `p=2767` does not divide a final low, midpoint, endpoint, or
frame coordinate and also leaves `(22)` unchanged.

For structural completeness, at the small frame wall `p=109` the whole
third polygon shifts by two:

```text
NP_109(T):(0,2)--(2,2)--(57,3)--(110,4).                  (28)
```

Its positive slopes are still the generic ones.  All other scalar, frame,
endpoint, and penultimate walls in `(17),(21),(23),(24)` are below `197`.
They are in the bounded lane.

## 6. The zero face and first-Witt separation

THM-3148 gives

```text
delta_6=Res(F_4,F_5)
 =44965855750876894470144533873830133760000
 =2^47*3^15*5^4*7*139*3767*12041*807241.                 (29)
```

Outside the four displayed odd exceptional primes above `7`, the original
pair `A,B` has no common unit root.  At those four primes, the monic residual
gcd and the exact first-Witt data are

```text
p=139:    gcd=v+26,       root 113,    omega=55;
p=3767:   gcd=v-1641,     root 1641,   omega=131;
p=12041:  gcd=v+3037,     root 9004,   omega=10951;
p=807241: gcd=v-22489,    root 22489,  omega=439007.       (30)
```

Here, for a lift `v_0`,

```text
omega=(A(v_0)/p)B'(v_0)-(B(v_0)/p)A'(v_0)        (mod p). (31)
```

THM-3175 says `(31)` is lift-independent and
covariant by the determinant under every polynomial row-frame change.  A
nonzero value excludes a common factor above a simple linear residual gcd.
For completeness, replacing `v_0` by `v_0+ph` adds
`h(A'(v_0),B'(v_0))` to the value-over-`p` pair, leaving its determinant
with the derivative pair unchanged.  A common lift would make those two
pairs proportional by Taylor expansion, forcing `(31)` to vanish.
Thus every exceptional zero face in `(29)` is excluded intrinsically.  At
the largest prime the complete first jet is

```text
(A/p,B/p,A',B')=(341729,683313,364516,409175),
(z_A,z_B)=(414823,557512).                                (32)
```

This is stronger than a neighboring-prime transport: it proves that the two
local branches already separate in their first infinitesimal residue layer.

## 7. Completion

Every odd `p<=193` has `r=p+4<=197`, so THM-3124 closes it exactly.  For
`p>=197`, Sections 3--5 separate every nonzero Newton slope of a hypothetical
common factor of `A,B`, including every coordinate and frame wall.  Section 6
excludes every zero-slope factor.  Therefore `A,B` are coprime over `Q_p`,
contradicting `(2)`.  QED.

## 8. Three-frame holotopy

The row changes

```text
(A,B)->(A,R)->(R,S)->(S,T)                                (33)
```

have determinants

```text
2p+7,                    -D_2,                    -D_3.   (34)
```

Consequently

```text
A wedge R=(2p+7)(A wedge B),
R wedge S=-(2p+7)D_2(A wedge B),
S wedge T=(2p+7)D_2D_3(A wedge B).                        (35)
```

The full exterior row transports flatly while its Newton projection passes
through two consecutive slope collisions before separating on the third row.
At a common residual root, the first-Witt wedge obeys the parallel law

```text
omega(MX)=det(Mbar(a))omega(X).                            (36)
```

Equations `(35),(36)` are the two layers of the same holotopy: the generic
Pluecker row over the DVR and its first infinitesimal residue.  A determinant
wall may destroy the projected chart, which is why `p=232961` needed the
direct normalized polygon `(25)` even though common-factor divisibility
remained lawful.

## 9. Exact evidence and scope

Run

```text
python 04-computation/factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.py
python -O 04-computation/factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.out.
```

The companion reconstructs all three quotient rows, fixed residuals,
midpoint and high formulas, generic and exceptional polygons, the offset-six
resultant and modular gcds, and all seven Hensel tuples.  It computes the
`p=232961` midpoint coordinates by a targeted coefficient-sum method
independent of the polynomial recurrence, and replays normally and under
`-O`.

Together with the earlier offset theorems and THM-3161, any still-open exact
quadratic factorial window has `r>=1999` and

```text
d,d-1,d-2,d-3,d-4,d-5,d-6 all composite.                 (37)
```

This is one exact `{0,1,2}` / `SFC(1)` offset.  It does not prove arbitrary
fixed offset, arbitrary support, `SFC(3)`, `NC(2)`, or the Gaussian Moment
Conjecture in two dimensions.  The Euclidean flag preserves common factors
but destroys a distinguished root and phase; Newton polygons preserve
valuations but destroy residual roots; the first-Witt sidecar restores only
the first infinitesimal zero-face mismatch.
