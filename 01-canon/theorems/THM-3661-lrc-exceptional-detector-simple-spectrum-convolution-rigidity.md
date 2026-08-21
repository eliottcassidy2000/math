---
id: THM-3661
title: "LRC exceptional detector simple-spectrum convolution rigidity"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  THM-3660's signed eight-address detector has zero Fourier coefficient only
  at the trivial character, and its other 168 coefficients are pairwise
  distinct, for both the split F13^2 and assembled cyclic C169 group laws.
  Hence convolution by the detector has simple spectrum; its full linear
  centralizer is the polynomial algebra it generates, and its translates
  span the 168-dimensional augmentation ideal.  The vertical derivative
  spans exactly the 156-dimensional zero-vertical-sum space.  These are group-
  algebra and cyclotomic facts, not a physical translation law or LRC(14).
source: kps-s191 / THM-3660 Fourier continuation, 2026-08-21
depends_on:
  - THM-3660-lrc-exceptional-leakage-functional-and-fourteen-edge-boundary
related:
  - THM-3658-lrc-mod169-carry-fourier-block-intertwiner
  - THM-2334-relation-residue-current-and-character-twist-pushforward
script: 04-computation/lrc_exceptional_detector_simple_spectrum_thm3661.py
output: 05-knowledge/results/lrc_exceptional_detector_simple_spectrum_thm3661.out
script_sha256: 7c8ebacaff5318b58ed1588d2e38edf1bc3c14ea81beec069d781450db986251
output_sha256: 61893d6e628232e877c73d6102f042654327fddfe7d51d1cf79c1aa6bc65419f
semantic_sha256: c9ffab0b789d087cc094e137274c6a17a2ef4c06aaa165a6969d77e404bc8ec9
hash_basis: raw LF bytes
---

# THM-3661 -- one eight-point detector separates every nonconstant mode

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  The signed exceptional detector is simultaneously cyclic for both
natural convolution algebras on the 169-address set.  More strongly, its
convolution operator has no repeated eigenvalue.

## 1. The detector and the two group laws

Let

```text
G_sp=F13^2,
G_cy=C169,                                             (1)
```

identified as sets by `a=r0+13r1`.  From THM-3660 define

```text
g=1_(X_+)-1_(X_-),                                    (2)

X_+={(12,0),(0,11),(6,5),(9,3)},
X_-={(0,12),(12,1),(6,7),(3,9)}.
```

Thus `g` has support eight and mean zero.

Over the pinned field

```text
F=F_p, p=755373809845391722745761,                    (3)
```

use the exact roots

```text
omega=298763986285447441216949, order 169,
zeta =123453826432109539554819=omega^13, order 13.    (4)
```

The two Fourier transforms are

```text
g_hat_sp(u,v)=sum_(r0,r1) g(r0,r1) zeta^(-ur0-vr1),

g_hat_cy(k)=sum_(r0,r1) g(r0,r1)
                         omega^(-k(r0+13r1)).          (5)
```

## 2. Both spectra are complete and simple

Exact evaluation gives

```text
g_hat_sp(u,v)=0 iff (u,v)=(0,0),
g_hat_cy(k)=0   iff k=0.                              (6)
```

Moreover the 168 nonzero coefficients in each transform are pairwise
distinct.  Their exact products are

```text
product_((u,v)!=(0,0)) g_hat_sp(u,v)
  =669422013050837354847410 !=0 mod p,

product_(k=1)^168 g_hat_cy(k)
  =540112521324494664145058 !=0 mod p.                (7)
```

The complete Fourier-ledger digests are

```text
split:  9e2088a22f89c4e81ac10bd824c3e735f09f06a4e54f89a070485da7a831e61e,
cyclic: c579b898150678a32c36ec081f0adf5c599a770ffdd9caa4c2aaf41fa8f6c809.
                                                                    (8)
```

These finite-field inequalities also certify the corresponding cyclotomic
statements in characteristic zero.  Indeed the coefficients in (5) are
algebraic integers in `Z[zeta_13]` or `Z[omega_169]`, and (4) defines a
reduction homomorphism to `F_p`.  An equality or zero in the cyclotomic
integer ring would remain one after reduction.  Therefore every nontrivial
complex Fourier coefficient is nonzero, and the 168 values are pairwise
distinct, for both group laws.

## 3. Convolution centralizer rigidity

For either group law let `C_g` denote convolution by `g` on the
169-dimensional function space.  The characters diagonalize `C_g`, with
eigenvalues given by (5).  Equations (6)--(7) say that all 169 eigenvalues,
including the trivial eigenvalue zero, are distinct.  Hence

```text
minimal polynomial of C_g = characteristic polynomial of C_g,
degree 169, square-free.                               (9)
```

The characteristic-polynomial digests over `F_p` are

```text
split:  79074cae506283cdcdacfecd6b47f54e8464abdb415e410fc0aca89cfe296ee3,
cyclic: fc36447edf23c4027e422bd4ddf49a3dede5d78c8d143767517e5e5adae121b0.
                                                                    (10)
```

The standard simple-spectrum centralizer theorem now gives, over `F_p` and
over the cyclotomic characteristic-zero field,

```text
Cent_End(C_g)=F[C_g]
 ={P(C_g):deg P<169}.                                 (11)
```

For completeness, Lagrange interpolation gives every character projector as

```text
P_chi=
 product_(psi!=chi)
 (C_g-g_hat(psi)I)/(g_hat(chi)-g_hat(psi)).            (12)
```

Thus a linear endomorphism commuting with this single sparse convolution
observable is forced to be a polynomial in it.  No translation-equivariance
hypothesis is needed in (11).

## 4. Translate completeness and exact inverses

Equation (6) also implies that the 169 translates of `g` span precisely the
augmentation ideal

```text
V_0={f:sum_x f(x)=0},              dim V_0=168,         (13)
```

for either group law.  Their only linear relation is the sum of all
translates.  Equivalently, a mean-zero signal orthogonal to every translated
exceptional detector is zero.

On `V_0` convolution by `g` is invertible.  Defining the inverse multiplier
by

```text
h_hat(0)=0,
h_hat(chi)=1/g_hat(chi) for chi!=0,                   (14)
```

gives

```text
h*g=delta_0-(1/169)1.                                 (15)
```

Over the pinned field, each inverse kernel has support 168 and takes 169
distinct values.  The split inverse vanishes only at `(7,7)` and the cyclic
inverse only at assembled residue `85`.  Their digests are

```text
split:  57461b00c68e7debadec08403171b06bed70cf5c9b9dc9852e3e091b73b7e40c,
cyclic: fde5a03b6439bd4a06186adb313827fc48e97397c8b35a87c65cd8958069bb53.
                                                                    (16)
```

The distinct inverse values furnish an exact address-separating coordinate,
but it depends on the chosen group law and is not yet a physical address.

## 5. The derivative detects exactly high-digit variation

Let

```text
D_1g(r0,r1)=g(r0,r1+1)-g(r0,r1).                     (17)
```

In split Fourier coordinates,

```text
(D_1g)_hat(u,v)=(zeta^v-1)g_hat_sp(u,v).              (18)
```

By (6), its support is exactly

```text
{(u,v):v!=0},                         size 156.        (19)
```

Therefore the translates of the 14-edge boundary `D_1g` span exactly

```text
W={f:sum_(r1) f(r0,r1)=0 for every r0}, dim W=156.    (20)
```

This sharply identifies the information retained by THM-3660's carry
boundary: it detects all and only high-digit variation.  The Fourier ledger
digest for (18) is

```text
5eab6723e39aa7b96b67f803abff2108013f31bb1a0faa6b74e3eda9c0711551.
                                                                    (21)
```

## 6. Exact verification and boundary

Reproduce with

```bash
python3 -B 04-computation/lrc_exceptional_detector_simple_spectrum_thm3661.py
python3 -B -O 04-computation/lrc_exceptional_detector_simple_spectrum_thm3661.py
```

The assertion-free companion source-pins THM-3660, checks both root orders,
all 338 Fourier evaluations, distinctness, products, both square-free
characteristic polynomials, a direct rank-168 translate elimination, the
156 derivative modes, both inverse kernels, and both convolution identities.
Normal and optimized streams are byte-identical.

This theorem concerns two explicitly defined finite group algebras and their
cyclotomic lifts.  It does not identify THM-2334 target twists with either
address group, prove that physical chronology permits address translation,
transport a current, restore a translated detector to the fixed `X`, prove
exceptional entry, or imply LRC(14).  The centralizer resemblance to
horizontal-endomorphism rigidity in other problems is structural only.
**QED.**
