---
id: THM-3208
title: "Positive sparse-unit bispectrum charged line and central-clutch obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Over an odd cyclic root group, every nonzero nonnegative one- or two-site
  fibre is automatically a group-algebra unit.  Its normalized bispectrum
  reconstructs the projective cyclic orbit, while its raw Hermitian
  bispectrum is the residual common-scalar charged line.  On the positive
  coefficient cone that line has a canonical positive orientation: equal
  normalized orbit and equal mass force equality, never a central sign
  reversal.  Fibrewise orbit-and-mass transport would therefore contradict a
  physical `-1` clutch, but no such LRC horn transport is constructed here.
audit: >
  The pure-integer cyclotomic companion pins THM-2312, THM-2802, and THM-3190;
  performs 5,605 automatic-unit Fourier checks and 5,605 positive
  scale/translation checks over p=5,7,13; rebuilds 541 positive whole-face
  orientations; and exhausts 10,416 normalized-bispectrum collision pairs
  over p=5,7.  It also checks the even-radix, broad-support, signed-mass,
  unequal-scale, and normalized-same/raw-opposite boundaries.  Normal and
  optimized replay agree with the stored transcript.  An independent hostile
  audit rederived the unit, equality-fibre, positive-orientation, varying-root
  translation, integrated nonvanishing, and conditional clutch arguments;
  replayed the immutable companion; and accepted every stated scope boundary.
source: root/2026-08-02
depends_on:
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-2802-normalized-unit-bispectrum-and-projective-cyclic-orbit-reconstruction
  - THM-3190-root-neutral-central-odd-bispectrum-clutch-criterion
related:
  - THM-2889-dicyclic-reverse-action-joint-carrier-and-skew-lift-separation
  - THM-2466-delayed-word-simultaneous-drift-service-retention
script: 04-computation/lrc14_positive_sparse_unit_bispectrum_clutch_thm3208.py
output: 05-knowledge/results/lrc14_positive_sparse_unit_bispectrum_clutch_thm3208.out
script_sha256: f4257d5b87b060f229c99857102465ff665aac453cf0b9159b0698c373717c18
output_sha256: 8d83920f6c4b2a9554b3b64a74277eaf9a2b8b6294708dc6fdcea434582729e8
hash_basis: LF-normalized bytes
---

# THM-3208 -- positive sparse-unit bispectrum charged line and central-clutch obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2312 proves that a complete raw bispectrum face is positive on every
active one/two-root word fibre.  THM-2802 proves, in a different part of the
repository, that the normalized bispectrum of a cyclic group-algebra unit
reconstructs its projective cyclic orbit.  These mechanisms compose exactly:
the sparse positive fibres are automatically units, and the raw bispectrum is
the scalar-charged line left above the reconstructed orbit.

The positive cone orients that line.  Consequently a physical pair of sparse
positive fibres with the same origin-free orbit and the same mass cannot
differ by THM-3190's central sign.  This replaces a chosen root-origin map by
a weaker orbit-and-mass transport obligation.  It does not construct that
transport or the missing THM-2889 horn edge.

## 1. Positive sparse fibres are automatic units

Let `p` be odd, let `zeta` be a primitive `p`th root of unity, and let

```text
f=sum_(r in C_p) f_r z^r in C[C_p],
fhat(k)=sum_r f_r zeta^(-kr).                              (1)
```

Assume that `f` is nonzero, every `f_r>=0`, and at most two coefficients are
nonzero.  Then

```text
fhat(k)!=0                         for every k in C_p.     (2)
```

Indeed the one-site case is immediate.  In the two-site case, after removing
one nonzero root phase, a vanishing Fourier coefficient would say

```text
a+b zeta^j=0,                     a,b>0.                  (3)
```

Taking absolute values forces `a=b`, and then `(3)` forces `zeta^j=-1`.
That is impossible at odd order.  By the splitting criterion in THM-2802,
`(2)` says exactly that `f` is a unit of `C[C_p]`.

The oddness and support hypotheses are sharp as automatic criteria.  At
`p=2`, `(1,1)` has a zero nontrivial Fourier mode.  At `p=13`, the positive
uniform vector has twelve zero modes.  The theorem does not say that every
positive unit is sparse; sparsity is a sufficient physical input inherited
from THM-2312.

## 2. The reconstructed base and the charged line

For a unit `f`, use THM-2802's normalized bispectrum

```text
Beta_f(k,l)=fhat(k)fhat(l)/(fhat(k+l)fhat(0)).             (4)
```

Also retain the raw Hermitian bispectrum

```text
H_f(k,l)=fhat(k)fhat(l) conjugate(fhat(k+l)).              (5)
```

If

```text
g=lambda z^h f,                  lambda in C^*, h in C_p, (6)
```

then direct Fourier covariance gives

```text
Beta_g=Beta_f,
H_g=|lambda|^2 lambda H_f.                                 (7)
```

Thus `Beta` is the origin-free projective base, while `H` retains one unit of
common scalar phase.  In particular `lambda=-1` fixes `(4)` and negates `(5)`.
This is the exact normalized-same/raw-opposite boundary required by
THM-3190's central-odd observable.

Conversely, THM-2802 says

```text
Beta_g=Beta_f   iff   g=lambda z^h f                    (8)
```

for units.  Therefore the fibre of the normalized-bispectrum map is precisely
the scalar-times-translation torsor, and `(5)` is its scalar-charged line.
This is a quotient reconstruction statement, not a nontrivial cohomology
class: the normalized array is a multiplicative coboundary and translations
are its exact kernel.

## 3. Positivity canonically orients the charged line

Suppose `f` and `g` are nonzero coefficientwise nonnegative real units and
`Beta_f=Beta_g`.  Equation `(8)` gives `g=lambda z^h f`.  Choose any positive
coefficient of `f`.  Its translated coefficient and the corresponding
coefficient of `g` show that

```text
lambda=c in R_(>0).                                        (9)
```

Consequently

```text
H_g=c^3 H_f.                                               (10)
```

If the masses agree, then their zero Fourier modes agree:

```text
sum_r g_r=sum_r f_r>0.                                    (11)
```

But `(6)` at frequency zero says `sum g=c sum f`, so `c=1`.  Hence

```text
Beta_g=Beta_f and equal positive mass
       imply g=z^h f and H_g=H_f.                          (12)
```

In particular `H_g=-H_f` is impossible because every Fourier coefficient is
nonzero, so every coordinate in `(5)` is nonzero.  Equivalently, the full
complex charged line has a `-1` point, but its coefficientwise positive locus
occupies only the canonically positive ray.

Without `(11)`, the strongest conclusion is the positive scaling `(10)`.
For example `g=2f` has the same normalized bispectrum and raw bispectrum
`8H_f`.  Normalization reconstructs projective scale; it does not invent a
mass anchor.

## 4. Origin-free packet criterion

Let `(Y,mu)` be a measured base and let `f_y,g_y` be measurable nonnegative
root fibres.  Assume for almost every `y` that either both vanish, or both are
nonzero and supported on at most two roots.  On the active fibres assume

```text
Beta_(f_y)=Beta_(g_y),
sum_r f_y(r)=sum_r g_y(r).                                (13)
```

The shift `h=h(y)` reconstructed by `(12)` may vary with `y`.  This causes no
loss: every coordinate of `(5)` is translation-invariant.  Thus

```text
H_(g_y)(k,l)=H_(f_y)(k,l)               a.e. in y,        (14)
```

and for every common measurable gate `Q subset Y`,

```text
integral_Q H_g(k,l)dmu=integral_Q H_f(k,l)dmu.             (15)
```

No root origin or measurable choice of `h(y)` is needed in `(14)--(15)`.
An orientation choice is still required: reversal sends `(k,l)` to
`(-k,-l)` rather than fixing the displayed array.

At `p=13`, if the active part of `Q` has positive measure, THM-2312's exact
whole-face identity gives

```text
sum_(k,l,k+l nonzero) integral_Q H_f(k,l)dmu>0.            (16)
```

Hence the integrated vector in `(15)` is nonzero.  This is stronger than
merely knowing one unspecified coordinate survives.

## 5. Conditional central-clutch obstruction

Suppose a proposed physical realization of the THM-2889/THM-3190 horn gives
direct and torus-normalized via root fields on one common measured word
packet.  If, after the known orientation/index action is removed, it proves
the fibrewise orbit-and-mass conditions `(13)`, then `(15)--(16)` force

```text
C_via=C_direct!=0.                                        (17)
```

If the same realization also implements the quaternionic central clutch,
THM-3190 instead requires

```text
C_via=-C_direct.                                          (18)
```

Equations `(17)--(18)` contradict characteristic zero.  Therefore no such
positive sparse physical horn can satisfy all of:

```text
common base and word;
fibrewise normalized-orbit transport;
fibrewise mass transport;
the claimed central -1 clutch.                            (19)
```

This is a proved conditional obstruction and a smaller connection contract
than a chosen amplitude isomorphism: cyclic origins may vary and are never
selected.  The load-bearing transport hypotheses in `(19)` are not supplied
by current canon.  An arbitrary linear action on the 132 integrated cubic
coordinates is not automatically a normalized-bispectrum index action.
THM-2466 supplies a conditional common delayed word, but not the fibrewise
orbit and mass identification required here.

## 6. Sharp boundaries and scope

1. **Fourier zeros.**  The normalized array is undefined off the unit locus.
   Broad positive support can have such zeros.
2. **Signed coefficients.**  At `p=13`, `(1,-1)` has zero mass and zero
   whole-face cubic sum.  Positivity is the orientation mechanism.
3. **Scalar sign.**  `f` and `-f` have equal normalized bispectra and opposite
   raw bispectra.  The negative packet lies outside the positive cone.
4. **Mass.**  Equal normalized orbit without `(11)` permits any positive
   cubic rescaling.
5. **Orientation.**  Reversal reindexes the array; it is not a cyclic origin
   change.
6. **Order of operations.**  A word-integrated raw current need not itself be
   a group-algebra unit, and there is no meaningful operation of normalizing
   the integral by pretending that integration commutes with `(4)`.  The
   criterion is fibrewise before integration.

The theorem proves neither the physical `e9=(-9,+9,QB)` edge, the common
oriented horn base, the torus action on root fibres, semantic owner transport,
a row exclusion, nor `LRC(14)`.  It gives the exact cheapest positive test:
transport only the origin-free normalized orbit and the zero-mode mass on one
common packet; if both survive, the central clutch is impossible.

## 7. Exact evidence

Run

```text
python 04-computation/lrc14_positive_sparse_unit_bispectrum_clutch_thm3208.py
python -O 04-computation/lrc14_positive_sparse_unit_bispectrum_clutch_thm3208.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
an explicit integer basis for `Q(zeta_p)`, pins all three dependencies, and
performs the counts stated in the audit field.  Its finite collision census is
an independent control of `(8)` on positive one/two-site banks; the proof of
the unrestricted unit statement remains THM-2802.  There is no floating point,
random sampling, imported executable, or assertion-sensitive test.

**QED.**
