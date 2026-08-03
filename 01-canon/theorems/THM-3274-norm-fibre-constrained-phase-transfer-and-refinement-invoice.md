---
id: THM-3274
title: "Norm-fibre constrained phase transfer and refinement invoice"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the pinned
  F_169 norm-phase model, a varying norm fibre has an exact 12-state quotient
  with irreducible order-12 recurrence; a varying projective scalar line has
  coarsest phase refinement (phase,normalized transverse determinant squared),
  48 states and a three-mode closed form; and the abstract THM-3246 seam has
  no proper phase-refining equitable quotient, requiring all 168 points, but
  its punctured transfer has an exact Fourier/derivative characteristic
  polynomial and matrix recurrence order 168. Fixed-increment hostiles are
  separate. No physical LRC walk or LRC(14) decrement is asserted. A
  post-promotion VERIFIED-EXACT addendum proves that the seam's raw twelve
  target-phase counts already decode all 168 sources without being given the
  source phase, with seven coordinates necessary and sufficient; that delta
  is independently hostile-audited.
source: root/2026-08-03
audit: >
  The assertion-independent companion pins THM-3246, THM-3267 and THM-3268,
  reconstructs the finite field, all constrained alphabets, colour
  refinements, spectra, recurrences, Hankel controls, Fourier factors and
  fixed-increment hostiles. An independent scratch engine used integer-coded
  F_169 arithmetic, generic colour refinement, direct SymPy characteristic
  polynomials, resultant-based Fourier norms and separate Rabin powering. It
  reproduced every quotient, recurrence, factorization, digest and hostile
  without importing the target companion. Normal, optimized and stored
  outputs agree byte-for-byte. The post-promotion companion additionally
  exhausts all source pairs for the raw seam decoder, its exact separation
  metrics, every six- and seven-coordinate projection, degree sidecar, cyclic
  window and scalar-normalized straight-arc conjugacy; independent audit of
  this delta independently reconstructed every stated census and accepted the
  simple-module and scope conclusions without repair.
depends_on:
  - THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word
  - THM-3267-norm-phase-factorization-ladder-and-projective-determinant-blindness
  - THM-3268-nonzero-translation-norm-phase-walk-closed-form-and-rank-eleven-mixing-mode
related:
  - THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary
script: 04-computation/lrc_constrained_norm_phase_transfer_atlas_20260803.py
output: 05-knowledge/results/lrc_constrained_norm_phase_transfer_atlas_20260803.out
script_sha256: 24a77d063bc36b45c43d10427124a33b68e8972e63fa97ddfdeba305d0cdb523
output_sha256: bdd318c9cd54a7ac26aaacf6e7bd199288192aaf8d7c5413b7641deeb5915b99
hash_basis: LF-normalized bytes
---

# THM-3274 -- constrained norm-phase transfers and their exact state invoices

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Universe and two increment semantics

Work in the chosen model from
[THM-3268](THM-3268-nonzero-translation-norm-phase-walk-closed-form-and-rank-eleven-mixing-mode.md):

```text
F_169=F_13[u]/(u^2-2),        alpha=1+2u,
X=F_169^*,                    phi(alpha^j)=j mod 12.      (1)
```

For an increment alphabet `S subset X`, distinguish:

```text
varying-set walk: choose a fresh t in S at every step, with x+t!=0;
fixed-increment walk: choose one t once, repeat it, and stop on hitting 0. (2)
```

The distinction in `(2)` is load-bearing. Every transfer matrix and closed
form below uses varying-set semantics. Fixed translations supply hostile
controls and do not iterate those quotients.

## 2. One norm fibre: twelve phases and twelve irreducible modes

Let

```text
S_0={t in X:phi(t)=0},              |S_0|=14.             (3)
```

The partition of `X` by `phi(x)` is outgoing-equitable for the varying
`S_0` walk. In phase order `0,...,11`, its symmetric quotient is

```text
Q_0=
(0,0,2,2,2,0,2,2,0,2,1,0)
(0,0,2,0,2,0,0,2,2,2,2,2)
(2,2,1,2,0,2,1,2,2,0,0,0)
(2,0,2,2,2,0,0,2,0,0,2,2)
(2,2,0,2,0,2,2,0,1,2,1,0)
(0,0,2,0,2,2,2,2,2,0,2,0)
(2,0,1,0,2,2,2,0,1,0,2,2)
(2,2,2,2,0,2,0,0,0,2,0,2)
(0,2,2,0,1,2,1,0,2,2,0,2)
(2,2,0,0,2,0,0,2,2,2,0,2)
(1,2,0,2,1,2,2,0,0,0,2,2)
(0,2,0,2,0,0,2,2,2,2,2,0).                           (4)
```

Its row sums are 13 in source phase zero and 14 in the other phases. Its
characteristic polynomial is

```text
p(lambda)=lambda^12-13lambda^11-77lambda^10+884lambda^9
 +1508lambda^8-19811lambda^7-2325lambda^6+160510lambda^5
 -96880lambda^4-341392lambda^3+243616lambda^2
 +27968lambda-23104.                                    (5)
```

The exact Rabin--Frobenius certificate modulo 331 gives

```text
x^(331^12)=x mod p,
gcd(x^(331^6)-x,p)=gcd(x^(331^4)-x,p)=1 mod 331.         (6)
```

Hence `(5)` is irreducible over `Q`. Because `(4)` is symmetric, it has
twelve distinct real algebraic eigenvalues.

There is therefore no nonzero proper rational `Q_0`-invariant phase sector.
For every nonzero `v in Q^12`, the cyclic list

```text
v,Q_0v,...,Q_0^11v                                      (6a)
```

is a basis. Symmetry gives the dual observability statement as well. In this
precise rational sense, restricting to one norm fibre destroys THM-3268's
constant-plus-augmentation split maximally.

Let `C_n(d)` count paths aggregated over sources whose terminal/source phase
difference is `d`. Every `C_n(d)` obeys the order-12 recurrence defined by
`(5)`. For every `d in C_12`, the corresponding `12 by 12` Hankel determinant
is nonzero; their tuple digest is

```text
92dc06c80eebcc5059a1c83719b987b3b4223b5b0aa3fbc842216ec3b0f0c934. (7)
```

Thus every one of the twelve relative-phase sequences has minimal scalar
recurrence order exactly 12.

Multiplying source and increment by `alpha^(-r)` shows that every absolute
norm fibre `S_r={t:phi(t)=r}` has a phase-shift conjugate of `(4)`.

For the fixed increment `t=1`, the same-phase sources `alpha^0` and
`alpha^12` land in phases 10 and 9, so phase is not even pointwise equitable.
The fixed profile at length 13 is

```text
(156,0,0,0,0,0,0,0,0,0,0,0).                           (8)
```

This is the characteristic-13 hostile, not the varying-fibre recurrence.

## 3. One projective line: the exact squared-transverse repair

Take

```text
L^*=F_13^*(1,0)={(lambda,0):lambda!=0},        |L^*|=12. (9)
```

Norm phase alone is not outgoing-equitable: every source phase has four
distinct profiles. Writing `x=(a,b)`, varying increments in `(9)` preserve
`b`, while the norm `a^2-2b^2` sees the transverse coordinate through `b^2`.
Generic colour refinement proves that the coarsest outgoing-equitable
refinement of phase is exactly

```text
(phi(x),b^2).                                            (10)
```

It has 48 nonempty cells: 36 of size four and 12 of size two. One refinement
step produces `(10)`, and a second leaves it stable.

The `b^2=0` block is a cell split of loopless `K_12`; each of the six
nonzero-square blocks is a cell split of loopless `K_13`. Consequently the
48-state quotient spectrum is

```text
12 with multiplicity 6,
11 with multiplicity 1,
-1 with multiplicity 41.                                (11)
```

Its minimal polynomial is `(lambda-12)(lambda-11)(lambda+1)`. Every fixed
linear observable obeys

```text
C_(n+3)=22C_(n+2)-109C_(n+1)-132C_n.                    (12)
```

For relative phase difference the exact closed forms are

```text
C_n(0)=(300*12^n+26*11^n+1858*(-1)^n)/13,
C_n(d)=168*(12^n-(-1)^n)/13,                 d odd,
C_n(d)=(144*12^n+26*11^n-170*(-1)^n)/13,    d even, d!=0.
                                                                  (13)
```

Their total is `156*12^n+12*11^n`.

For an arbitrary projective direction `L=dF_13`, multiplication by `d^(-1)`
conjugates the walk to `(9)`. In original coordinates the invariant repair is

```text
(phi(x)-phi(d), (det(d,x)/Norm(d))^2).                   (14)
```

Thus the determinant which is globally phase-blind in
[THM-3267](THM-3267-norm-phase-factorization-ladder-and-projective-determinant-blindness.md)
becomes, after normalization and squaring, the exact dynamic sidecar for this
specific constrained operation. A fixed scalar increment again returns the
hostile `(8)` and not the 48-state transfer.

## 4. The abstract THM-3246 seam: no compression, exact spectrum

Reinterpret the negative-owner seam exponents of
[THM-3246](THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md)
only as the abstract varying increment alphabet

```text
N={alpha^j:j=0,...,5,162,...,167},              |N|=12. (15)
```

Multiplication by `alpha^6` normalizes `(15)` to the consecutive Singer arc
`{alpha^0,...,alpha^11}`. Although it contains one increment in each phase,
the phase partition is maximally non-equitable: within every source phase,
all fourteen points have distinct outgoing profiles. The first outgoing
refinement therefore has 168 singleton cells. Hence the coarsest
outgoing-equitable refinement which refines phase is the full point/Singer-
exponent partition. The setwise stabilizer of `N` in `GL_2(F_13)` is trivial.

The source phase used by that refinement is in fact redundant. Define the raw
target-phase histogram

```text
H_N(x)_r=#{t in N:x+t!=0, phi(x+t)=r},       r in C_12. (15a)
```

The 168 vectors `H_N(x)` are pairwise distinct. Their minimum pairwise `L1`
distance is 2, their minimum number of differing coordinates is 2, and their
affine span has dimension 12. Thus `(15a)` is a complete discrete tomography
chart `X -> Z_nonnegative^12`, without an externally supplied source phase.

Only seven count coordinates are needed, and seven is sharp. Exactly 21 of
the `binom(12,7)` coordinate sets remain injective; the lexicographically
first is

```text
{0,1,2,5,8,9,11}.                                      (15c)
```

No six-coordinate projection is injective. The best one,
`{0,1,3,5,6,8}`, has 166 signatures, consisting of 164 singleton fibres and
two double fibres. Appending the total degree `sum_r H_N(x)_r in {11,12}`
still leaves every six-coordinate projection noninjective. Among cyclic
consecutive windows the sharp length is eight, attained for starting phases
`4,7,8,9,10,11`.

Multiplication by `alpha^6` sends `N` to the straight arc
`{alpha^0,...,alpha^11}`. More generally, for `c in X`,

```text
H_(cN)(cx)=shift_(phi(c)) H_N(x).                       (15b)
```

Here `shift_s(h)_r=h_(r-s)`. Hence every scalar copy, including the straight
arc, is again a raw point decoder. Frobenius conjugacy preserves the same
conclusion. This is still an
abstract varying-increment query: it does not construct a physical seam
measurement or same-ancestry packet.

Restore zero and let `B_N` be the 169-vertex additive Cayley adjacency for
increment set `N`. Let `A_N` be its principal submatrix on `X`. Additive
characters diagonalize `B_N`; the 168 nontrivial characters form fourteen
projective Galois orbits of size 12. Therefore

```text
chi_(B_N)(lambda)=(lambda-12)*product_[c] F_[c](lambda), (16)
```

with fourteen explicit monic degree-12 integer factors. Thirteen are
distinct; directions `(1,7)` and `(0,1)` give the same factor.

Translation symmetry makes all 169 principal characteristic minors equal.
The derivative of a determinant is the sum of its principal minors, so

```text
chi_(A_N)(lambda)=chi_(B_N)'(lambda)/169.                (17)
```

A direct `168 by 168` characteristic polynomial and an independent
cyclotomic Fourier reconstruction agree. Their exact digests are

```text
chi_(B_N): b8025e76edb12075836f20791c13418030c09e73e87b2011c25dba1f0b8a27e3,
chi_(A_N): 653c5c1d8251e610ab744b6a72f5bba5e782cf2bb74738c0f1d33dbba1be1dd4.
                                                                  (18)
```

Over `Q`, `chi_(A_N)` factors in degrees 12 and 156. The degree-12 factor is

```text
lambda^12-lambda^11+27lambda^10-66lambda^9+66lambda^8
 +1039lambda^7-4029lambda^6+5784lambda^5+9062lambda^4
 -55849lambda^3+110475lambda^2-102233lambda+41003.       (19)
```

The degree-156 cofactor has digest

```text
ccee063f7dd21de43d7502f4d1db02153cd78bf0dac92a65307824238cbfbf98. (20)
```

An exact gcd with the derivative proves `chi_(A_N)` squarefree. Hence the
minimal polynomial of the full punctured transfer matrix has degree 168.
Every matrix entry and scalar observable satisfies the recurrence from
`(17)`; minimal order 168 is asserted for the matrix, not for every scalar
projection.

The one-per-phase alternate transversal with exponents

```text
(0,13,2,15,4,17,6,19,8,21,10,23)                       (21)
```

also refines immediately to singletons, but its punctured characteristic
digest is
`58aaacff1ebc624f20e95541ba3d214356389739137c0683feff88722c512609`.
One-per-phase balance does not determine the spectrum. Repeating the one seam
increment `1` gives `(8)`, not `A_N`.

## 5. Information-loss ladder and closed-form classification

The three alphabets give the exact state invoice

```text
14-element norm fibre       -> 12 states, irreducible order 12;
12-element projective line  -> 48 states, three spectral modes;
12-element Singer seam      -> 168 states, matrix order 168.      (22)
```

Alphabet size is not the controlling invariant. Compression is governed by
the interaction between the additive increment set, multiplicative norm
phase and additive stabilizers.

Equation `(14)` is the main structural repair: a coordinate globally blind to
phase can be exactly the missing dynamic sidecar for a specified consumer.
Equation `(17)` is the complementary spectral repair: deleting the translation
origin destroys Fourier diagonalization but the derivative identity recovers
the exact characteristic polynomial.

## 6. Scope and reproduction

The seam in `(15)` is used only as an abstract increment alphabet. Nothing
here says that THM-3246 supplies a physical ancestry/current transition. The
theorem constructs no LRC endpoint, safety certificate, row exclusion,
physical owner-to-phase map, or `LRC(14)` ledger decrement.

Run

```text
python3 04-computation/lrc_constrained_norm_phase_transfer_atlas_20260803.py
python3 -O 04-computation/lrc_constrained_norm_phase_transfer_atlas_20260803.py
```

and compare LF-normalized bytes with the declared output. The companion uses
exact integer, finite-field and polynomial arithmetic only, with no assertion
node, floating literal, fitted recurrence or randomness.

QED.
