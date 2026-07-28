---
id: THM-2638
title: "Radial-height graded Wick decoder and Laplace-forgetting boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Before radial integration, the angular constant term of a Gaussian moment
  has a canonical radial-height grading.  Coefficient extraction at height A
  decodes a channel exactly when that channel is alone at height A; all
  channels admit a nonnegative graded decoder exactly when radial height is
  injective on the balanced channel set.  For
  P=a Z^6+b W^2+c W^18, A=6j+12t at moment 4j, so every one of the j+1
  scalar-colliding channels becomes a private radial row.  The scalar moment
  is the pushforward L(s^A)=A!, and its fourth-moment cancellation proves that
  scalar nullity does not imply vanishing of the lifted rows.  More sharply,
  every algebraic source-torus one-parameter grading restricts on balanced
  channels to a multiple of A, while every Gaussian-preserving angular
  grading is zero.
  The quadratic support Z^2,ZW,W^2 has two distinct second-moment channels at
  the same A=2, giving the minimal source-character no-go.  Independent
  coefficient-torus twists always separate a finite channel fibre by a
  mixed-radix DFT, but those are external observables and no nullity
  consequence.  No new NC2 or GMC(2) conclusion follows.
source: wild-holotopy-mining-2026-07-28
depends_on:
  - THM-1645-gmc2-angular-layer-is-dvdk-the-gap-is-purely-radial
  - THM-2631-homogeneous-wick-channel-linear-decoder-and-private-row-no-go
related:
  - THM-2020-gmc2-finite-place-channel-separation
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
script: 04-computation/gmc2_radial_height_graded_decoder_thm2638.py
output: 05-knowledge/results/gmc2_radial_height_graded_decoder_thm2638.out
script_sha256: b63782f6130a4a8b19232a43a38f7e2e7e461b503b09f1a4f07cbfc4266e4d19
output_sha256: f2e2890b9b3e927c263038dcfc094f49d2926580b68a1b274176364ed0ce5606
hash_basis: LF-normalized bytes
---

# THM-2638 -- the radial shell lift restores exactly one forgotten grade

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2631 proves that raw scalar moments cannot linearly decode two Wick
channels in one homogeneous degree.  This theorem identifies a lawful lift in
which those rows can become private, and then determines the precise limit of
the lift.  The missing coordinate is radial height.  It lives before the
Laplace/factorial functional and is genuinely destroyed by that functional.

The result is a decoder theorem for the lifted shell polynomial.  It is not a
decoder from scalar Gaussian moments and therefore does not evade THM-2631.

## 1. The shell polynomial and its scalar pushforward

Let

```text
P=sum_(i=1)^k c_i Z^(a_i) W^(b_i),       c_i!=0,
q_i=a_i-b_i,

B_m={r in N^k: |r|=m and q dot r=0},
A(r)=sum_i a_i r_i=sum_i b_i r_i.                         (1)
```

For `s=|Z|^2` and angular coordinate `u`, put

```text
Lambda_s(u)=P(s^(1/2)u,s^(1/2)u^(-1)).                   (2)
```

Although an individual term of (2) may contain a half-integral power of `s`,
every balanced product has integral exponent `A(r)`.  Its angular constant
term is therefore the honest polynomial

```text
G_m(s)=CT_u Lambda_s(u)^m
      =sum_(r in B_m) [m!/prod_i r_i!] c^r s^A(r).        (3)
```

Let `L` be the factorial/Laplace functional

```text
L(s^A)=A!.                                                (4)
```

Then the Wick expansion is exactly

```text
M_m=E[P^m]=L(G_m)
   =sum_(r in B_m) [m! A(r)!/prod_i r_i!] c^r.            (5)
```

Thus the passage from the shell polynomial to the raw moment is the
pushforward

```text
(m,A) -> m.                                               (6)
```

This is the Gaussian analogue of forgetting a carry coordinate: distinct
radial rows can become one dense scalar row.

## 2. Exact radial-height decoder

For a fixed `m` and height `A`, write

```text
G_(m,A)=[s^A]G_m
       =sum_(r in B_m, A(r)=A) [m!/prod_i r_i!] c^r.      (7)
```

Call the rows `(m,A)` the radial-height graded bank and the columns `(m,r)`
the level-tagged Wick channels.

> **Graded decoder theorem.**  A prescribed channel `r0 in B_m` has a
> nonnegative scalar decoder from the graded rows if and only if
>
> ```text
> {r in B_m:A(r)=A(r0)}={r0}.                             (8)
> ```
>
> In that case
>
> ```text
> c^r0=[prod_i r0_i!/m!] G_(m,A(r0)).                     (9)
> ```
>
> The entire channel fibre `B_m` has a nonnegative left inverse exactly when
> `A:B_m->N` is injective.

Indeed, row `(m,A)` is positive precisely on the height fibre in (7).  The
private-row criterion of THM-2631 therefore gives (8), and (9) is the explicit
inverse.  Notice that no signed cancellation is required: coefficient
extraction is a positive one-row decoder whenever the height is private.

More generally, for any additive grade `g=(g_i)`, the same proof says that
the `g dot r` lift decodes exactly the singleton grade fibres.  The radial
grade is distinguished because it is canonically present in polar
coordinates before `L`.

## 3. Every MISTAKE-211 channel becomes private

For

```text
P=a Z^6+b W^2+c W^18,                                    (10)
```

THM-2631 gives, at `m=4j`,

```text
r(t)=(j+2t,3(j-t),t),                 0<=t<=j.            (11)
```

The radial height is

```text
A(r(t))=6(j+2t)=6j+12t.                                  (12)
```

Hence `A` is injective on every `B_(4j)`, and all `j+1` channels have private
graded rows.  Explicitly,

```text
G_(4j)(s)=sum_(t=0)^j
  (4j)!/[(j+2t)!(3(j-t))!t!]
  a^(j+2t)b^(3(j-t))c^t s^(6j+12t).                      (13)
```

At the first return,

```text
G_4(s)=4ab^3 s^6+4a^3c s^18,                             (14)

L(G_4)=4*6!*ab^3+4*18!*a^3c.                            (15)
```

Taking `a=b=1` and `c=-6!/18!` makes (15) zero while both coefficients in
(14) remain nonzero.  This is the sharp scope control: scalar nullity gives
only `L(G_m)=0`; it does not give `G_m=0` or the vanishing of any private
radial row.  The kernel of `L` is load-bearing.

## 4. The source-character boundary is exactly radial height

The full algebraic source torus scales

```text
(Z,W)->(lambda Z,mu W).                                   (16)
```

An integer one-parameter subgroup
`(lambda,mu)=(t^alpha,t^beta)` assigns monomial `i` the grade

```text
alpha a_i+beta b_i.                                      (17)
```

On a balanced channel `r`, its total grade is

```text
sum_i (alpha a_i+beta b_i)r_i
  =alpha A(r)+beta A(r)
  =(alpha+beta)A(r).                                     (18)
```

Consequently every source-torus character sees at most radial height on the
balanced fibre.  The Gaussian-preserving angular subgroup
`(lambda,mu)=(zeta,zeta^(-1))` has `alpha+beta=0`, so every balanced channel
has character zero.  Angular Fourier refinement produces one trivial row and
cannot improve the scalar rank.

There is a minimal equal-height hostile.  Let

```text
P=aZ^2+bZW+cW^2.                                         (19)
```

At `m=2`, the two channels

```text
(1,0,1),             (0,2,0)                             (20)
```

both have `A=2`, and

```text
G_2(s)=(2ac+b^2)s^2,             M_2=2!(2ac+b^2).        (21)
```

They therefore remain one row under every source-torus character, not only
under angular rotations.  Three support monomials and moment two are minimal.
For two support monomials with unequal charges, fixed length and charge balance
determine the occupation vector uniquely.  If their charges agree, a balanced
return exists only when both charges are zero; the two distinct radial
monomials then have distinct heights, so fixed length and height again
determine the occupation vector.

Thus radial grading is a genuine positive lift but not a universal channel
separator.  Equal-height collisions require a grading not induced by source
exponents, or a nonlinear whole-face mechanism.

## 5. External coefficient characters always separate a finite fibre

For completeness, fix `m`, put `B=m+1`, and assign coefficient `c_i` the
mixed-radix weight

```text
g_i=B^(i-1),              N=B^k.                         (22)
```

The code

```text
gamma(r)=sum_i g_i r_i mod N                             (23)
```

is injective on every occupation vector with `0<=r_i<=m`.  If `xi` is a
primitive `N`-th root and one is supplied the externally twisted bank

```text
M_m^(u)=M_m(xi^(u g_1)c_1,...,xi^(u g_k)c_k),
                     0<=u<N,                             (24)
```

then root orthogonality gives the exact signed decoder

```text
(1/N)sum_(u=0)^(N-1) xi^(-u gamma(r0)) M_m^(u)
 = [m! A(r0)!/prod_i r0_i!]c^r0.                        (25)
```

This proves that the obstruction is forgotten grading, not an intrinsic
linear dependence among the coefficient monomials.  But (24) changes the
coefficient phases of the polynomial.  The hypothesis `M_m(P)=0` supplies
only the `u=0` equation and does not imply any other row in (24).  Therefore
(25) is an external tomography control, not a Gaussian-nullity argument.

THM-2020's finite-place carry filtration has the same logical boundary.  At
good primes an algebraic coefficient vector is unit-valued, so a unique
least-valuation channel cannot cancel.  That gives a nonvanishing
certificate when the unique-minimum condition holds; it does not turn all
raw scalar moments into an exact multichannel decoder.

## 6. Exact evidence and scope

Run

```text
python 04-computation/gmc2_radial_height_graded_decoder_thm2638.py
python -O 04-computation/gmc2_radial_height_graded_decoder_thm2638.py
```

The dependency-free companion enumerates the channels of (10) through
`j=12`, checks (11)--(15), verifies private radial rows at every enumerated
level, reconstructs the equal-height hostile (19)--(21), exhausts a bank of
source-torus characters against (18), and verifies mixed-radix injectivity
and exact root-orthogonality selectors.  All computations are integer or
rational; no floating point is used.  Normal and optimized runs LF-normalize
to the stored transcript.

The theorem proves a statement about the shell polynomial before radial
integration.  It does not infer shellwise vanishing from scalar moment
vanishing, manufacture the external twists (24), establish a finite-place
unique minimum on arbitrary support, improve an effective moment cutoff, or
prove NC2/GMC(2).  NC2 and GMC(2) are already proved by THM-2022.  This result
instead identifies the exact grading restored by the radial lift, the exact
forgetful map that destroys it, and the minimal collision beyond the reach of
all source characters.

`QED.`
