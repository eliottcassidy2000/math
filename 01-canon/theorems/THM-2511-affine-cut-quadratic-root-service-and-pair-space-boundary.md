---
id: THM-2511
title: "Affine-cut quadratic root service and the pair-space boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Contracting the
  seven cut translates of one THM-2508 affine Radon bank quadratically gives
  a nonnegative C_13-valued energy E_(tau,a).  It is not an invariant scalar:
  the 72 energies transform equivariantly under the full physical CRT affine
  gauge.  Every nonconstant rational energy has all twelve nonzero C_13
  Fourier colours.  On the complete THM-2436 atlas, every one of 14,952
  nonflat defects has nonconstant energy for at least 64 of the 72 pairs
  (tau,a), sharply.  Hence, in any chosen measurable bundle trivialization,
  some fixed pair is nonconstant on at least 8/9 of the essential locus,
  giving inherited parent floors 16/63 and 8/21 and an unnormalized
  coefficient floor 2268^(-11).  The construction is quadratic
  signed pair-space data: it loses defect sign, cut character, and owner
  orientation, and varying root phase can still cancel after integration.
  It does not supply one Boolean owner/deep current, exclude a live row, or
  prove LRC(14).
source: codex-2026-07-27-affine-cut-quadratic-energy
depends_on:
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
script: 04-computation/lrc14_affine_cut_quadratic_root_service_thm2511.cpp
output: 05-knowledge/results/lrc14_affine_cut_quadratic_root_service_thm2511.out
script_sha256: cf5d6eb1976aa61a99667d5d403f45207c9952a0a6a9282e59767765045b3e0d
output_sha256: b6261657c874281919fd8cc8250d394eba7a2ce328060ddadf5009122a62cb45
hash_basis: working-tree bytes (LF)
---

# THM-2511 -- quadratic cut energy is a positive equivariant root bundle

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2508 closes the linear affine-carry seam by retaining a nontrivial cut
character.  Its full invariant scalar contraction vanishes, while its global
quadratic norm is positive but target-neutral.  There is an intermediate
object: contract only the cut-translation coordinate and retain the output
root.  This produces a nonnegative *equivariant vector*, not an invariant
scalar.

On the actual finite THM-2436 atlas that vector is almost always
nonconstant.  Thus quadratic contraction recovers positive root service with
an `8/9` component invoice.  One exact price is visible: the observable lives
on a signed pair space and forgets orientation data still needed by the
owner/deep-current argument.  No complete kernel classification is claimed.

## 1. The quadratic cut energy

Let

```text
d:F_13 x F_7 -> Q,

sum_(r in F_7)d(h,r)=0                    for every h,           (1)
```

and retain THM-2508's affine-cut Radon outputs

```text
R_(tau,a,c)(v)
 =sum_(r in F_7)d(v-tau rep(ar+c),r),                         (2)

tau in F_13^*,       a in F_7^*,       c in F_7.
```

Define the cut-contracted energy

```text
E_(tau,a)(v)=sum_(c in F_7)R_(tau,a,c)(v)^2.                  (3)
```

It is pointwise rational and nonnegative.  If `d` is integral, as in the
THM-2436 atlas, then `E_(tau,a)` is a nonnegative integer vector.

The square in (3) is important.  This is neither the zero linear contraction

```text
sum_c R_(tau,a,c)=0                                           (4)
```

from THM-2508 nor the root-contracted component scalar
`sum_v E_(tau,a)(v)`.  That scalar is translation-invariant but is permuted
across `(tau,a)` by the unit gauges; only the sum over all `72` labels is
fully gauge-invariant.  Equation (3) keeps the `F_13` output coordinate.

## 2. Exact affine covariance

For

```text
g=(A,H;B,C) in AGL_1(F_13) x AGL_1(F_7)
```

act by

```text
d^g(h,r)=d(A^(-1)(h-H),B^(-1)(r-C)).                         (5)
```

THM-2508 proves

```text
R^(d^g)_(tau,a,c)(v)
 =R^d_(A^(-1)tau,aB,aC+c)(A^(-1)(v-H)).                      (6)
```

Squaring, summing over `c`, and reindexing `c'=aC+c` gives

```text
E^(d^g)_(tau,a)(v)
 =E^d_(A^(-1)tau,aB)(A^(-1)(v-H)).                           (7)
```

Thus the `72=12*6` vectors `E_(tau,a)` form a finite equivariant bundle.
The translation carry permutes the contracted cut coordinate, the two unit
multipliers permute `(tau,a)`, and the root output transforms by the usual
affine action.  In particular, a root translation does not disappear: it
shifts the argument of `E`.  This is why (3) can retain target phase even
though the globally summed quadratic norm cannot.

## 3. Nonconstant is equivalent to all twelve root colours

For `alpha in F_13`, use the unnormalized transform

```text
Ehat_(tau,a)(alpha)
 =sum_(v in F_13)E_(tau,a)(v)zeta_13^(-alpha v).              (8)
```

For a rational vector `e:F_13 -> Q`, the following are equivalent:

```text
e is nonconstant;

ehat(alpha)!=0 for one alpha!=0;

ehat(alpha)!=0 for every alpha!=0.                            (9)
```

Indeed, if `ehat(alpha)=0`, then the degree-at-most-twelve polynomial

```text
Q(X)=sum_(v=0)^12 e(v)X^v
```

is divisible over `Q` by

```text
Phi_13(X)=1+X+...+X^12.
```

Its coefficients are therefore all equal.  The converse is immediate.
Equivalently, the twelve nontrivial transforms are the Galois conjugates of
one another.  Applying this to (3) proves

```text
E_(tau,a) nonconstant
  iff every one of its twelve nonzero C_13 colours survives. (10)
```

This conclusion is universal; the finite atlas is used only to force
nonconstancy for many `(tau,a)`.

## 4. The exact complete-atlas invoice

Retain the complete normalized THM-2436 universe used independently in the
THM-2507 finite atlas:

```text
guard [0,25] in Z/91Z;
five supports chosen from the 3,276 unoriented cyclic AP_13 supports
  with unit step modulo 91;
blocker-source columns;
all 28 source multisets s_0<=s_1 in F_7.                      (11)
```

The inherited exhaustive generator gives

```text
41,379 assignments;
 1,736 flat assignments;
39,643 nonflat assignments;
14,952 distinct nonflat defect matrices.                     (12)
```

For each distinct nonflat defect, the exact companion computes all `72`
vectors (3) and tests whether their thirteen integer entries are equal.  The
number of nonconstant pairs has the exact histogram

```text
64:      4 defects,
66:    107 defects,
68:  1,498 defects,
70:  1,977 defects,
72: 11,366 defects.                                          (13)
```

Consequently every nonflat atlas defect satisfies the sharp pointwise bound

```text
#{(tau,a):E_(tau,a) is nonconstant}>=64=8*8,                 (14)

64/72=8/9.                                                   (15)
```

Sharpness already occurs at the three-row defect

```text
d(0,-)=-e_0+e_5,
d(7,-)=-e_2+e_5,
d(8,-)=-e_4+e_5,                                             (16)
```

whose eight constant pairs are

```text
(1,1),(1,2),(2,2),(3,6),(10,1),(11,5),(12,5),(12,6).         (17)
```

This witness is an actual member of the audited atlas, not an arbitrary
row-zero hostile.  Conversely, `3,586` distinct defects have at least one
constant pair.  More sharply, over the `72` fixed labels, the number of
nonconstant defect shapes ranges from `14,632` to `14,903`, so every fixed
label has between `49` and `320` constant hostiles.  No fixed pair is a
pointwise universal detector on this atlas.

The finite computation proves (13)--(17) for the exact universe (11).  It
does not assert the universal `64/72` inequality for arbitrary rational
row-zero arrays.

## 5. Fixed-pair parent floors

Let `P minus mathcal E` be THM-2436's essential-parent locus.  At almost every
point of this locus, (14) holds for the local normalized defect.  Sum the
indicators over the `72` bundle labels and integrate.  Fubini gives one fixed
label `(tau,a)` in any chosen measurable bundle trivialization with

```text
mu{Y in P minus mathcal E:E_(tau,a;Y) nonconstant}
  >=(8/9)mu(P minus mathcal E).                               (18)
```

THM-2435/2436 give

```text
mu(P minus mathcal E)>=(1+k)/7,               k in {1,2}.     (19)
```

Therefore the fixed chart-label nonconstant-energy locus floors are

```text
k=1:       (8/9)(2/7)=16/63,

k=2:       (8/9)(3/7)= 8/21.                                (20)
```

On the selected locus, (10) supplies all twelve nonzero root colours.  These
are component-locus floors in a measurable affine-bundle chart.  Equation
(7) proves that changing the chart transports the label and phase lawfully;
it does not turn the selected component into a chart-free scalar current.

## 6. Exact coefficient floor

Every THM-2436 defect satisfies

```text
||d||_1<=18.                                                  (21)
```

For each cut `c`, the triangle inequality in (2) gives

```text
||R_(tau,a,c)||_1<=||d||_1<=18.
```

Since `R` is integral on the atlas,

```text
sum_v E_(tau,a)(v)
 =sum_(c,v)R_(tau,a,c)(v)^2
 <=sum_c ||R_(tau,a,c)||_1^2
 <=7*18^2
 =2268.                                                       (22)
```

If `E_(tau,a)` is nonconstant, each value (8) with `alpha!=0` is a nonzero
algebraic integer in `Q(zeta_13)`.  Its twelve conjugates are the twelve
nonzero root colours and, by nonnegativity and (22), each has modulus at most
`2268`.  The field norm is a nonzero integer.  Hence

```text
|Ehat_(tau,a)(alpha)|>=2268^(-11).                            (23)
```

Combining (18)--(23), the unnormalized fixed-colour `L^2` energies have floors

```text
k=1:       (16/63)2268^(-22),

k=2:        (8/21)2268^(-22).                                (24)
```

The exact atlas maximum of the left side of (22) is only `316`; `2268` is the
uniform theorem bound inherited directly from `||d||_1<=18`.

## 7. Pair-space interpretation and exact losses

Expanding (3) gives

```text
E_(tau,a)(v)
 =sum_(c,r,r')
   d(v-tau rep(ar+c),r)
   d(v-tau rep(ar'+c),r').                                   (25)
```

Thus the positive number `E` is assembled from products of a *signed defect*
on two copies of the stalk.  It is not yet the measure of one Boolean atom or
one positive owner-supported intersection.  The contraction has the
following exact information ledger:

- **preserved:** nonnegativity, the `F_13` output, full CRT affine covariance,
  nonconstancy, and hence all twelve nonzero root colours;
- **destroyed:** the sign `d -> -d`, linear response, the individual cut
  character `beta`, and the orientation of the two defect legs;
- **still required:** a lawful signed pair-space realization on one ancestry
  cospan, an absolute root mark through integration, owner/arrival support,
  and descent to the old deep sheet.

In cut Fourier coordinates, (25) pairs `beta` with `-beta`; that is precisely
why the energy is positive and why the primitive `F_7` orientation retained
by THM-2508 is gone.  Moreover, (7) shows that a normalization carry
`H(Y)` translates the root vector.  Pointwise all-colour nonvanishing does
not prevent the integrated coefficient from cancelling when this phase
varies with `Y`.

On every nonflat atlas defect each scalar `sum_v E_(tau,a)(v)` is positive
and target-neutral; the affine gauge permutes these `72` scalars, while their
full sum is invariant.  Keeping `v` repairs target colour, not temporal
ownership.  Therefore THM-2511 strengthens positive static root service but
does not supply THM-2471 owner-loop drift, transport THM-2478's old signed
charge, transplant the already-empty high-septimal packet to the live `165`
rows, exclude a scalar row, or prove LRC(14).

## 8. Exact companion

Compile and run

```text
c++ -std=c++20 -O2 \
  04-computation/lrc14_affine_cut_quadratic_root_service_thm2511.cpp \
  -o /tmp/thm2511 && /tmp/thm2511

c++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_affine_cut_quadratic_root_service_thm2511.cpp \
  -o /tmp/thm2511-opt && /tmp/thm2511-opt
```

Both builds reproduce

```text
05-knowledge/results/lrc14_affine_cut_quadratic_root_service_thm2511.out
```

byte-for-byte.  The companion reuses the audited THM-2436 exhaustive universe
but independently derives every defect and every quadratic cut energy.  It
checks the exact assignment counts, row-zero and `L^1` laws, positivity,
the `2268` bound, all `1,076,544=14,952*72` defect/pair incidences, histogram
(13), sharp witness (16)--(17), and the fixed-label hostile range.

An independent audit reran both optimization levels, reproduced the stored
transcript and hashes, and reimplemented the sharp witness directly in
Python.  It rederived covariance (7), the cyclotomic equivalence (9), the
`8/9` Fubini invoice, both parent floors, and the `2268^(-11)` norm exponent.
It also caught and repaired one scope ambiguity: a fixed
`sum_v E_(tau,a)(v)` is permuted by the unit gauges; only the sum over all
`72` labels is a fully gauge-invariant scalar.  The audit accepted the
explicit finite-atlas and signed-pair-space boundaries.  **QED.**
