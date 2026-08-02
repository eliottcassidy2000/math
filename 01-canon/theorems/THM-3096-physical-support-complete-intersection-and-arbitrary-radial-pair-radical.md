---
id: THM-3096
title: "Physical-support complete intersection and arbitrary-radial pair radical"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  positive factorial support whose physical first-window resultant is
  nonzero, the first t factorial moment forms are a homogeneous complete
  intersection of degrees 1,...,t.  Consequently the two-charge family
  P=alpha W+Z C(ZW), with arbitrary complex radial coefficients on that
  support, has moment null variety exactly {alpha=0} union {C=0}; moments
  through 2t give literal pair-radical certificates with exponent
  t(t-1)/2+1.  This realizes the sharp HYP-8765 cutoff on this family, not on
  general multi-charge or neutral-channel supports.
source: root-gmc-radial-complete-intersection-2026-08-01
audit: >
  An independent hostile audit checked the determinant-one physical
  elimination, projective zero-locus equivalence, Cohen--Macaulay regular
  sequence and socle invoice, literal ideal membership, every Wick and
  factorial normalization, the radical intersection, and the neutral
  cancellation hostile.  Normal, optimized, and stored companions replayed
  identically; the declared LF hashes, evidence counts, and documentation
  checks all passed.
depends_on:
  - THM-1645-gmc2-angular-layer-is-dvdk-the-gap-is-purely-radial
related:
  - THM-2019-gmc2-affine-height-supports
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3091-arbitrary-gap-remote-pair-desuspension-and-exact-Jensen-contraction
  - THM-3093-arbitrary-gap-remote-cluster-monge-flag-compactification
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_radial_complete_intersection_thm3096.py
output: 05-knowledge/results/gmc_radial_complete_intersection_thm3096.out
script_sha256: 4cf0966c1ae6e6c0468d1ab318358504e72d107a59effe814bdd42884dca4826
output_sha256: ce8453d0dcc931c3e4b73bfc5a11498f0f5832fba8794926782403957a738ccb
hash_basis: LF-normalized bytes
---

# THM-3096 -- physical-support complete intersection and radial pair radical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The pure two-charge specialization of HYP-8765 is exactly the Strong
Factorial problem: for `H=sC(s)`, must one of `L(H),...,L(H^t)` detect every
nonzero `t`-slot coefficient vector?  The physical first-window resultant is
precisely the missing support-sensitive certificate.  Once it is nonzero,
the factorial moment forms are not merely set-theoretically faithful: they
are a complete intersection with an explicit socle bound.

## 1. The physical support resultant

Let

```text
Lambda={lambda_1<...<lambda_t} subset Z_(>0),
f_lambda=s^lambda/lambda!,
L(s^n)=n!.                                                (1)
```

For `x=(x_1,...,x_t)`, put

```text
H_x=sum_i x_i f_(lambda_i),
F_j(x)=L(H_x^j),                    1<=j<=t.              (2)
```

Thus `F_j` is homogeneous of degree `j`, and the normalization in `(1)`
gives

```text
F_1=x_1+...+x_t.                                         (3)
```

Use a determinant-one linear change with last coordinate `w=F_1`, restrict
`F_2,...,F_t` to `w=0`, and call their normalized resultant

```text
R_Lambda.                                                 (4)
```

This is the repository's physical first-window resultant for the factorial
support `Lambda`; multiplying its forms by their fixed positive normalization
constants does not change whether it vanishes.  The elementary linear-form
resultant quotient gives the exact equivalences

```text
R_Lambda!=0
 iff Res(F_1,...,F_t)!=0
 iff V(F_1,...,F_t)={0} in A^t.                          (5)
```

The theorem assumes `(5)`.  It is a property of the support, not of the
coefficient vector subsequently inserted.

## 2. Complete intersection and the socle invoice

Let

```text
I_Lambda=(F_1,...,F_t) subset C[x_1,...,x_t].            (6)
```

By `(5)`, this ideal has height `t`.  The polynomial ring is
Cohen--Macaulay, so the `t` homogeneous generators of degrees

```text
1,2,...,t                                                 (7)
```

form a regular sequence.  Hence the quotient is an Artin complete
intersection with Hilbert series

```text
Hilb_(C[x]/I)(z)=product_(j=1)^t(1-z^j)/(1-z)^t,         (8)
```

length `t!`, and socle degree

```text
D=sum_(j=1)^t(j-1)=t(t-1)/2.                             (9)
```

Put

```text
N_0=D+1.                                                  (10)
```

Every homogeneous class above degree `D` is zero.  In particular, for each
support coordinate there are homogeneous polynomials

```text
deg A_(i,j)=N_0-j
```

such that the literal ideal identity

```text
x_i^N_0=sum_(j=1)^t A_(i,j)(x)F_j(x)                    (11)
```

holds.  This is the finite algebraic defect/carry sidecar; no
Nullstellensatz exponent is left unspecified.

## 3. Arbitrary radial coefficients and exact Wick moments

Let `Z` be a standard complex Gaussian, write `W=conj(Z)` and `s=ZW`, and
choose arbitrary complex coefficients

```text
alpha,c_1,...,c_t in C,
C(s)=sum_i c_i s^(lambda_i-1),
x_i=lambda_i! c_i.                                      (12)
```

Consider the two-charge polynomial

```text
P=alpha W+Z C(s).                                        (13)
```

Every term in the first summand has charge `-1`; every radial term in the
second has charge `+1`.  Odd moments vanish.  In degree `2j`, balance forces
exactly `j` choices from each side, so

```text
M_(2j):=E[P^(2j)]
 =binom(2j,j)alpha^j L((sC(s))^j)
 =binom(2j,j)alpha^j F_j(x),            1<=j<=t.         (14)
```

This is one balanced channel, but its radial coefficients are completely
arbitrary and may have unrelated phases.

Let

```text
J_t=(M_2,M_4,...,M_(2t)) subset C[alpha,c_1,...,c_t].    (15)
```

If `alpha=0` or `c_1=...=c_t=0`, every generator vanishes.  Conversely, if
`alpha!=0`, equations `(14)` and `(5)` force `x=0` whenever all generators
vanish.  Therefore

```text
V(J_t)=V(alpha) union V(c_1,...,c_t),

sqrt(J_t)=(alpha c_1,...,alpha c_t).                     (16)
```

Thus this family satisfies NC2 exactly: the moment nullcone is the union of
the two strict charge-one-sided coordinate spaces.

The result is effective.  Multiply `(11)` by `alpha^N_0` and use `(14)`:

```text
(alpha c_i)^N_0
 =(lambda_i!)^(-N_0) sum_(j=1)^t
   binom(2j,j)^(-1) alpha^(N_0-j)
   A_(i,j)(lambda_1!c_1,...,lambda_t!c_t) M_(2j).       (17)
```

Hence every pair product has its **same explicit exponent** and is already
in the finite moment ideal, not merely its radical:

```text
(alpha c_i)^N_0 in J_t,               1<=i<=t.           (18)
```

The polynomial `(13)` has `K=t+1` possible monomials and primitive return
length `R=2`.  Its largest moment in `(15)` is

```text
2t=(K-1)R,                                               (19)
```

exactly the cutoff proposed by HYP-8765.  THM-2173 proves this cutoff is
dimension-sharp from below.

## 4. Support consequences

The conclusion applies to **every** positive support with `R_Lambda!=0`.
In particular:

1. one- and two-slot supports are elementary;
2. every arbitrary three-slot support is good by THM-2824, so arbitrary
   three-term radial `C` is detected by moment at most six;
3. every support supplied by the proved arbitrary-gap remote-pair theorem
   THM-3091 inherits the sharp radial cutoff and literal certificates; and
4. once candidate THM-3093 is independently promoted, the same holds above
   every one of its arbitrary-gap remote clusters, including its
   unconditional low-child tails.

This is the precise leverage from SFC geometry into NC2: a support resultant
is upgraded to an arbitrary-coefficient radial nullcone theorem.  It proves
the pure two-charge effective HYP-8765 obligation on every resultant-good
support, rather than only on monomial or sign-restricted radial data.

## 5. The neutral-channel boundary is real

The complete-intersection flag cannot be summed blindly with a neutral
channel.  Take

```text
P=Z+W+i sqrt(2)(s-1).                                    (20)
```

The visible pair channel at level two contributes `2L(s)=2`, but the neutral
square contributes `-2L((s-1)^2)=-2`.  Directly,

```text
M_1=0,
M_2=L([-2(s-1)^2]+2s)=0,
M_3=L([-2i sqrt(2)(s-1)^3]+6i sqrt(2)s(s-1))
    =2i sqrt(2)!=0.                                      (21)
```

Thus neither the Monge flag nor a nonzero pair channel preserves a scalar
sum across unrelated charge-zero contributions.  General multi-charge or
neutral-channel HYP-8765 still needs a joint channel filtration; the theorem
does not resurrect the false first-return isolation of MISTAKE-211.

## 6. Exact evidence and scope

The exact companion verifies:

1. seven one-, two-, and three-slot physical resultants after exact `F_1`
   elimination;
2. the complete-intersection Hilbert series, length, and socle cutoff through
   eight slots;
3. seventeen literal identities `(11)` using `241` exact multiplier columns;
4. `7,890` direct Wick expansions against `(14)`;
5. `1,415` exact finite-grid controls of the zero-locus equality `(16)`; and
6. the exact neutral hostile `(20)--(21)` over `Q(i sqrt(2))`.

Run

```text
python 04-computation/gmc_radial_complete_intersection_thm3096.py
python -O 04-computation/gmc_radial_complete_intersection_thm3096.py
```

Both modes must equal the stored transcript after LF normalization.

The theorem fixes a finite positive support with nonzero physical resultant.
It does not prove SFC for a support whose resultant vanishes, general
multi-charge or neutral-channel pair radicals, arbitrary-radial GMC(2) in
all charge patterns, LRC(14), JC(2), or DC(2).

**QED.**
