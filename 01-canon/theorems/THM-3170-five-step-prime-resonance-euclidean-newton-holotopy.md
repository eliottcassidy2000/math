---
id: THM-3170
title: "Five-step prime resonance Euclidean--Newton holotopy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  odd prime p, no exact quadratic has three zero factorial moments beginning
  at r=p+3.  Two Euclidean frame changes separate all positive Newton slopes;
  the offset-five resultant and two proved neighboring-prime transports close
  the zero face.  Beyond the small p=5 boundary, the frame determinants
  explain the unique p=17 slope wall.
audit: >
  A separate exact derivation independently reconstructed both Euclidean
  quotients, every valuation tier, the midpoint and high-coefficient formulas,
  all arithmetic-wall polygons, the lawful A/R zero-face splice, both
  prime-neighbor coordinates, and the exterior/resultant transport.  A final
  immutable file-level hostile audit accepted the proof and scope.  The
  canonical script replays normally and under optimization exactly against
  the stored transcript, with both declared LF-normalized hashes verified.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
  - THM-3159-right-neighbor-prime-resonance-zero-face-separation
related:
  - THM-3153-four-step-prime-resonance-second-euclidean-newton-separation
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
script: 04-computation/factorial_five_step_prime_resonance_euclidean_newton_thm3170.py
output: 05-knowledge/results/factorial_five_step_prime_resonance_euclidean_newton_thm3170.out
script_sha256: ef76af79a109b7d4b7da39fa51177af107c9c14797b9518744a2e2393719f2ba
output_sha256: a3b286417a6b8f62824d3619bd4573621dd4a5dbd495ed859df826bba79279d8
hash_basis: LF-normalized bytes
---

# THM-3170 -- five-step prime resonance Euclidean--Newton holotopy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(t^j)=j!,                         Q(t)=a+bt+ct^2           (1)
```

with `abc!=0`, and let `p` be an odd prime.  Put

```text
r=p+3,                             d=r+2=p+5.               (2)
```

Then the three consecutive moments

```text
L(Q^r), L(Q^(r+1)), L(Q^(r+2))                              (3)
```

cannot all vanish.

The proof uses two exact Euclidean rows for positive slopes, THM-3148's
offset-five resultant for the zero face, THM-3146 for the exceptional prime
`20747`, and THM-3159 for the exceptional prime `249721`.  The two remainder
steps also form an exact exterior-frame transport.  Beyond the small `p=5`
boundary, its determinant wall is precisely the exceptional midpoint prime
`17`.

## 1. Resonant pair and first polygon

THM-3124 says that `(3)` would force `b/a=-1/d`.  Divide by `a`, put
`v=d(c/a)`, and define

```text
A_n(v)=L((d-t+vt^2)^n) in Z[v].                            (4)
```

It is enough to prove coprimality of

```text
A=A_(p+3),                         B=A_(p+4).                (5)
```

Their coefficients are

```text
[v^j]A_n
 =binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.               (6)
```

Write `p=2m+1`.  THM-3148 gives the complete low face

```text
A == 10(360v^3+21v+37)                             (mod p). (7)
```

For `p>=199`, its constant, linear, and cubic coordinates are units.  Lucas
reduction of the outer binomial in `(6)`, followed by the factorial
valuation, gives

```text
v_p(A_j)>=1 for 4<=j<=m,
v_p(A_j)>=2 for m+1<=j<=p+3.                               (8)
```

Indeed, in the first range the binomial supplies one `p`; in the second
range either that binomial and the factorial each supply one `p`, or
`j>=p` and the factorial supplies two.  The leading coefficient
`A_(p+3)=(2p+6)!` has valuation exactly two.  Therefore

```text
NP_p(A):(0,0)--(3,0)--(p+3,2),       slopes 0,2/p.          (9)
```

All primes below `199` will be discharged by the bounded theorem, so no
small-prime assertion is hidden in `(7)--(9)`.

## 2. Two exact Euclidean rows

The full linear quotient of `B` by `A` has coefficients

```text
C=(2p+8)(2p+7),
q=-2(p+4)(p+5)/(2p+5).                                    (10)
```

Clear the denominator and define

```text
R=(2p+5)B-[(2p+5)Cv-2(p+4)(p+5)]A.                        (11)
```

Then `deg R<=p+2`.  Direct simplification of the next two coefficients gives

```text
R_(p+2)=-2(p+2)(p+3)(p+4)(24p+85)(2p+2)!,

R_(p+1)= 2(p+1)(p+2)(p+3)(p+4)
          (24p^2-19p-380)(2p)!.                           (12)
```

The quotient of `A` by `R` is `q_1v+q_0`, where

```text
q_1=-2(2p+3)(2p+5)/[(p+4)(24p+85)],

q_0=4(p+5)(2p+3)(28p+95)
    /[(p+4)(2p+1)(24p+85)^2].                              (13)
```

Put

```text
D =(p+4)(2p+1)(24p+85)^2,
N_1=-2(2p+3)(2p+5)(2p+1)(24p+85),
N_0=4(p+5)(2p+3)(28p+95),

S=D A-(N_1v+N_0)R.                                        (14)
```

The top two coefficients cancel exactly, so `deg S<=p+1`.  Every common
factor of `A,B` divides both `R` and `S`.

Modulo `p`, the complete fixed residuals are

```text
A == 10(360v^3+21v+37),
B == 5(40320v^4-5760v^3+720v^2+160v+329),
R == -75(544v^2+1216v-307),
S == 250(2338491v-482198).                         (mod p) (15)
```

In particular,

```text
S_0=-120549500=-2^2*5^3*353*683,
S_1= 584622750= 2*3*5^3*601*1297                 (mod p). (16)
```

Away from those four low-coordinate primes, `S_0,S_1` are units.

## 3. Midpoint and high sidecars

Set

```text
k=m+2=(p+3)/2.                                             (17)
```

Equations `(6),(11),(14)` give the uniform bounds

```text
v_p(R_j)>=1 for j>=3,          v_p(R_j)>=2 for j>=m+2,
v_p(S_j)>=1 for j>=2,          v_p(S_j)>=2 for j>=m+3.     (18)
```

For the first row, use the exact quadratic residual in `(15)` and then apply
`(8)` term by term in `(11)`.  For the second, use its exact linear residual
and the first-row bounds in `(14)`.

At `j=k`, the terms `D A_k` and `N_0R_k` have valuation at least two.  Only
the shifted `R_(m+1)` term survives after division by `p`; inside it only the
shifted `A_m` term survives.  Wilson's theorem and

```text
(m!)^2 == (-1)^(m+1)                              (mod p)
```

give

```text
A_m/p ==
 6(-1)^(m-3)5^(m+4)/[m(m-1)(m-2)(m-3)]             (mod p),

S_k/p ==
 -4284000(-1)^(m-3)5^(m+4)
 /[m(m-1)(m-2)(m-3)]                               (mod p). (19)
```

Here

```text
4284000=2^5*3^2*5^3*7*17.                                  (20)
```

Thus `(19)` is nonzero for every `p>=199`.

The endpoint and penultimate coefficients are

```text
S_(p+1)=4p(p+1)(p+2)(p+3)(p+4)(2p+5)E_4(p)(2p-2)!,

E_4(p)=256p^4-28672p^3-281120p^2-881568p-905545,

S_p=-4p(p-1)(p+1)(p+2)(p+3)(p+4)(2p+5)E_5(p)(2p-4)!,

E_5(p)=256p^5-30720p^4-27936p^3+1633888p^2
       +7109519p+8370560.                                  (21)
```

Since

```text
905545=5*61*2969,
8370560=2^7*5*11*29*41,                                   (22)
```

the endpoint has valuation two for `p>=199`, except at `p=2969`, where it
has valuation three and the penultimate coefficient has valuation two.

## 4. The complete arithmetic-wall atlas

For `p>=199` outside

```text
{353,601,683,1297,2969},                                   (23)
```

equations `(16)--(22)` force

```text
NP_p(S):(0,0)--(1,0)--(k,1)--(p+1,2),
slopes       0,       2/(p+1), 2/(p-1).                    (24)
```

At the two constant-coordinate walls,

```text
p in {353,683}:
NP_p(S):(0,1)--(1,0)--(k,1)--(p+1,2),
slopes       -1,      2/(p+1), 2/(p-1).                    (25)
```

At the two linear-coordinate walls,

```text
p in {601,1297}:
NP_p(S):(0,0)--(k,1)--(p+1,2),
slopes       2/(p+3), 2/(p-1).                             (26)
```

At the raised endpoint,

```text
p=2969:
NP_p(S):(0,0)--(1,0)--(k,1)--(p,2)--(p+1,3),
slopes       0,       2/(p+1), 2/(p-3), 1.                (27)
```

Every positive slope in `(24)--(27)` is distinct from `2/p`, the sole
positive slope of `A`.  The intermediate residual wall `p=307`, where the
constant coordinate of `R mod p` vanishes, has

```text
NP_307(R):(0,1)--(1,0)--(2,0)--(309,2),                    (28)
```

but leaves `(24)` unchanged.  This is an exact hostile showing why the final
flag, rather than every intermediate polygon separately, is the correct
object.

The prime `17` is the sharp midpoint/denominator hostile.  Its exact second
remainder has hull

```text
NP_17(S):(0,0)--(1,0)--(18,2),                              (29)
```

so the positive slope `2/17` overlaps that of `A`.  This is why the bounded
fallback in Section 6 is retained.

## 5. The lawful zero face

The second row can acquire spurious zero-face roots through its quotient
multiplier, so the lawful unit-root test is the original pair `A,B`,
equivalently `A,R` because `2p+5` is a unit for `p>=199`.  THM-3148 gives

```text
delta_5=Res_v(F_(3,5),F_(4,5))
 =2246282972173187481600
 =2^16*3^7*5^2*11^2*20747*249721.                          (30)
```

Thus, outside `p=20747,249721`, `A` and `R` have no common unit root.
Together with the positive-slope separation, `A,B` are coprime.

This use of `A,R` is load-bearing.  At a root of `A`, equation `(11)` gives
`R=(2p+5)B`; by contrast `(14)` gives
`S=-(N_1v+N_0)R`, so `A,S` can share a spurious residual root at a zero of
the displayed multiplier.

## 6. Finite and neighboring-prime completion

Every odd `p<=197` has `r=p+3<=200`, so THM-3124 closes it exactly.  This
includes the genuine low-face exception `p=11` and the midpoint hostile
`p=17`.

The two large resultant exceptions in `(30)` are transported to already
proved neighboring-prime windows:

```text
p=20747:   q=20749 is prime,
            r=p+3=q+1, d=p+5=q+3;
            THM-3146 closes this exact window.

p=249721:  q=249727 is prime,
            r=p+3=q-3, d=p+5=q-1;
            THM-3159 closes this exact window.              (31)
```

All odd primes are now covered, proving `(3)`.  QED.

## 7. Euclidean-frame holotopy

The two changes of frame

```text
(A,B) -> (A,R) -> (R,S)
```

have polynomial matrices

```text
M_1=[[1,0],[-Q_1,2p+5]],
M_2=[[0,1],[D,-(N_1v+N_0)]],                               (32)
```

where `Q_1=(2p+5)Cv-2(p+4)(p+5)`.  Their determinants are `2p+5` and
`-D`.  Hence, away from `{5,17}`, they are invertible changes of the two-row
lattice frame over `Z_p[v]`; the exterior row satisfies

```text
A wedge R=(2p+5)(A wedge B),
R wedge S=-(2p+5)D(A wedge B).                              (33)
```

Equivalently,

```text
Res(A,R)=(2p+5)^deg(A) Res(A,B).                            (34)
```

The Euclidean flag is a flat common-factor transport before projection to
Newton polygons.  Beyond `p=5`, the prime `17` is exactly the nonunit
determinant wall of the second frame, explaining structurally why `(29)` has
a shared slope.  In contrast, `p=307` is only a coordinate-chart wall of the
projected first
remainder; the full frame and final flag remain nondegenerate.

This is the Euclidean--Newton analogue of THM-3160's Pluecker lesson: the
full transported exterior row is flat, while a projected coordinate need
not be.  No GMC selector or positivity statement is transferred between the
two settings.

## 8. Exact companion

Run

```text
python 04-computation/factorial_five_step_prime_resonance_euclidean_newton_thm3170.py
python -O 04-computation/factorial_five_step_prime_resonance_euclidean_newton_thm3170.py
```

and compare byte-for-byte with the declared stored output.  The companion
uses integer, modular, and rational arithmetic only.  It checks the two
quotients, low residuals, midpoint formula, high identities, generic and
all exceptional polygons, including `p=307` and `p=2969` modulo `p^5`;
recomputes `(30)` by Bareiss elimination and all three modular residual gcds;
checks primality and coordinates in `(31)`; and enumerates all `44` odd
primes in the THM-3124 fallback range.  No floating point or fitted
polynomial is used.

## 9. Consequence and scope

Together with THM-3131, THM-3138, THM-3143, THM-3146, THM-3153, and the
bounded THM-3124 range, every still-open bad exact-quadratic resonance has

```text
r>=201,
d=r+2,d-1,d-2,d-3,d-4,d-5 all composite.                   (35)
```

The source object is the resonant pair `(5)`.  The Euclidean flag
`(11),(14)` preserves every common factor and destroys a distinguished root
and phase.  Newton polygons preserve root valuations and destroy residual
roots.  The fixed resultant `(30)` is the missing zero-face sidecar, while
`(31)` is a prime-neighbor transport sidecar.  The wall atlas is finite
because every relevant fixed-offset coordinate reduces to a fixed integer
or to a polynomial in `p` with fixed constant term.

This proves one exact `{0,1,2}` / `SFC(1)` offset.  It is not an arbitrary
fixed-offset induction, arbitrary support, `SFC(3)`, `NC(2)`, or the full
Gaussian Moment Conjecture in two dimensions.
