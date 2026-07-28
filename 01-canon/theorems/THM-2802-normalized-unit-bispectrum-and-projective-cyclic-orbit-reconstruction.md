---
id: THM-2802
title: "Normalized unit bispectrum and projective cyclic-orbit reconstruction"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  Over any splitting field for a finite cyclic group,
  the normalized Fourier bispectrum of a group-algebra unit is a complete
  invariant modulo nonzero scalar and cyclic translation.  Equality makes
  the Fourier ratio a character.  Every THM-2790 endpoint cycle therefore
  has a canonical projective-origin signature, and none equals the
  THM-2791 semantic two-point signature.  This reconstructs a coefficient
  orbit, not a positive carrier, physical frame map, semantic attachment,
  row exclusion, or LRC(14).  Until independent promotion, no proved result
  may depend on this candidate.
source: root/unit-bispectrum-holotopy-2026-07-28
depends_on:
  - THM-2790-universal-depth-two-central-response-and-carry-wall-spectrum
  - THM-2791-full-arm-orbit-transfer-and-lower-central-chord
related:
  - THM-2312-sparse-root-bispectrum-positive-word-current
  - THM-2439-cyclic-marker-replica-degree-and-homometric-gram-boundary
  - THM-2789-interval-gram-tomography-and-graceful-gap-tail-quadratic-detector
  - THM-2792-cyclic-unit-intertwiner-and-positive-naturality-boundary
script: 04-computation/lrc14_normalized_unit_bispectrum_thm2802.py
output: 05-knowledge/results/lrc14_normalized_unit_bispectrum_thm2802.out
script_sha256: c220b84cd9574d882640fb3fa9e186a1121722256b37530349c4bb4449713341
output_sha256: 95bc5935e0aaff968b37dc0c3a778c9232f15919fbc12dca2d7803c5e2b5df3c
hash_basis: LF-normalized bytes
---

# THM-2802 -- a unit cycle is reconstructed projectively by its bispectrum

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

THM-2790 proves that every central endpoint cycle has complete Fourier
support.  THM-2791 produces a very different full-spectrum object: one
positive two-point semantic chain.  THM-2792 gives the unique abstract
regular-module isomorphism between them, but correctly withholds a pointwise
physical identification.  The present theorem asks what can be recovered
without choosing a cyclic origin.

The answer is exact.  A normalized third-order Fourier ratio reconstructs
the entire projective cyclic orbit of any unit.  In cochain language it is a
multiplicative coboundary whose kernel is precisely the character group.
Thus the coefficient/origin ambiguity is solved, while the physical
allocation problem remains visible rather than being hidden inside an
arbitrary Fourier phase choice.

## 1. Splitting-field unit cycles

Let

```text
G=C_n=<z>,                    char(K) does not divide n,
mu_n subset K,
R=K[G]=K[z]/(z^n-1).                                    (1)
```

Fix `xi in K` of order `n`.  For

```text
f=sum_(j in C_n) f_j z^j
```

use the Fourier convention

```text
fhat(k)=sum_j f_j xi^(-kj),             k in C_n.         (2)
```

The splitting map

```text
R -> product_(k in C_n) K,              f |-> (fhat(k))_k (3)
```

is an algebra isomorphism.  Hence

```text
f is a unit of R  <==>  fhat(k)!=0 for every k.           (4)
```

For such a unit define its **normalized bispectrum**

```text
Beta_f(k,l)
 =fhat(k) fhat(l)/(fhat(k+l) fhat(0)),
                    k,l in C_n.                           (5)
```

Every denominator in `(5)` is nonzero by `(4)`.

## 2. Scalar and origin invariance

For `lambda in K^x` and `h in C_n`, put

```text
g=lambda z^h f.                                          (6)
```

Then

```text
ghat(k)=lambda xi^(-kh) fhat(k).                         (7)
```

Substitution in `(5)` cancels both copies of `lambda` and the total
translation phase:

```text
Beta_g(k,l)=Beta_f(k,l)             for all k,l.          (8)
```

Thus `(5)` is intrinsic to the projective cyclic orbit

```text
K^x <z> f.                                               (9)
```

Changing the displayed cycle origin, or removing one common nonzero
coefficient content, does not change it.

## 3. Completeness: equality forces a scalar translate

Let `f,g in R^x` and suppose

```text
Beta_f=Beta_g.                                           (10)
```

Define

```text
r(k)=ghat(k)/fhat(k),             q(k)=r(k)/r(0).         (11)
```

Taking the ratio of `(5)` for `g` and `f`, equation `(10)` becomes

```text
r(k)r(l)=r(k+l)r(0),
q(k+l)=q(k)q(l).                                         (12)
```

Therefore `q:C_n->K^x` is a character.  Since `K` contains `mu_n`, there is
a unique `h in C_n` with

```text
q(k)=xi^(-kh).                                           (13)
```

Put `lambda=r(0)`.  Equations `(11)--(13)` give

```text
ghat(k)=lambda xi^(-kh) fhat(k)
```

for every `k`, and Fourier inversion yields

```text
g=lambda z^h f.                                          (14)
```

Together with `(8)` this proves:

> **Normalized unit-bispectrum theorem.**  For units of a split cyclic group
> algebra, `Beta_f=Beta_g` if and only if `f` and `g` differ by one nonzero
> scalar and one cyclic translation.

No positivity, reality, or genericity hypothesis is used beyond the unit
condition.

## 4. The exact holotopy quotient

Normalize the Fourier one-cochain by

```text
a_f(k)=fhat(k)/fhat(0).                                  (15)
```

Then `(5)` is its multiplicative two-coboundary:

```text
Beta_f(k,l)=a_f(k)a_f(l)/a_f(k+l).                       (16)
```

Equation `(12)` says that two one-cochains have the same coboundary exactly
when their ratio is a one-cocycle, hence a character.  Under Fourier
duality those characters are precisely the physical translations `z^h`.
Consequently `(5)` gives an injection

```text
R^x/(K^x <z>)  ->  (K^x)^(C_n x C_n).                   (17)
```

This is a quotient reconstruction statement, not a nontrivial cohomology
class: every displayed bispectrum is already the coboundary `(16)`.  The
holotopy content is that the residual fibre of the coboundary map is exactly
the translation torsor, with no larger phase ambiguity.

If orientation is also forgotten, reversal

```text
f^vee(z)=f(z^(-1))
```

obeys

```text
Beta_(f^vee)(k,l)=Beta_f(-k,-l).                         (18)
```

Thus a dihedral quotient must additionally identify the two arrays in
`(18)`; cyclic-origin invariance alone does not forget orientation.

## 5. Application to the endpoint and semantic units

Let

```text
K=Q(zeta_N),          N=50,334,435,734,703,120,
```

the splitting field used by THM-2625/2790.  For every nonzero endpoint
direction `s`, determinant cycle `Delta`, and any chosen cycle origin,
THM-2790 gives

```text
J_(s,Delta)(z)=sum_(j=0)^12 P_(R_j+s)Q_(R_j)z^j
```

with all thirteen Fourier coordinates nonzero.  Hence each of the

```text
168*13=2,184                                             (19)
```

oriented endpoint cycles has a well-defined origin-independent signature

```text
Beta_(J_(s,Delta)).                                      (20)
```

This removes the arbitrary origin phase from comparisons inside the
endpoint bank.  It does not identify different physical directions:
rescaling `s` changes the ordered edge `P_(R+s)Q_R`, not merely the
coordinate `j`.  A formal reindexing by `u in (Z/13)^x` sends the array to

```text
(k,l) |-> Beta_f(u^(-1)k,u^(-1)l),                      (21)
```

but THM-2790 supplies no physical frame covariance asserting `(21)`.

THM-2791's normalized semantic chain is

```text
A=z^6(1+z).                                              (22)
```

It is a unit and has coefficient support two.  Every `J_(s,Delta)` has
coefficient support thirteen by THM-2790.  If their normalized bispectra
were equal, `(14)` would make `J_(s,Delta)` a scalar translate of `A`,
which preserves support size.  Therefore

```text
Beta_(J_(s,Delta)) != Beta_A
           for every one of the 2,184 endpoint cycles.   (23)
```

Equation `(23)` is an origin-free coefficient obstruction to a weighted
pointwise cyclic identification.  It does not contradict THM-2792: the
unique regular-module isomorphism there is convolution by an arbitrary unit,
not scalar translation.

## 6. Why third order is load-bearing

THM-2439 gives the non-dihedrally translated homometric pair

```text
U={0,1,3,9},                 V={1,2,5,7} in C_13.         (24)
```

Their complete cyclic autocorrelations and Fourier power spectra agree.
Both indicator polynomials are units over a characteristic-zero splitting
field: the zero mode is `4`, and every nonzero power is `3`.

They are not scalar translates, so the theorem forces

```text
Beta_(1_U) != Beta_(1_V).                                (25)
```

The exact companion checks `(25)` modulo `53`, where `C_13` splits.  Thus
the normalized bispectrum repairs precisely the projective phase/order loss
left by quadratic power data.  This parallels THM-2789:

```text
interval Gram data  -> incidence-column multiset, loses order;
unit bispectrum     -> projective cyclic orbit, loses physical semantics.
                                                                    (26)
```

The interval theorem uses convexity; the present theorem uses the
full-spectrum unit condition.  Neither hypothesis transfers automatically
to an arbitrary LRC packet.

## 7. Sharp boundaries

1. **Fourier zero.**  If one `fhat(k)` vanishes, `(5)` is undefined and the
   division proof fails.  Sparse-root bispectrum positivity alone does not
   imply reconstruction.
2. **Orientation.**  Reversal acts by `(18)` and is not silently quotiented.
3. **Physical frame.**  Equation `(21)` is a formal coefficient reindexing,
   not the physical replacement `s->us`.
4. **Positivity.**  Reconstruction of a projective coefficient orbit gives
   no nonnegative Markov kernel, Boolean carrier, or pointwise ancestry map.
5. **Depth.**  A `C_13` signature forgets THM-2791's higher
   `Z_2,...,Z_5` lower-central digits.
6. **Semantic attachment.**  The signature does not identify owner, word,
   root, endpoint chronology, or same physical support.

The cheapest next exact test is now finite and honest: compute `(20)` on the
canonical bank and classify its projective cyclic and dihedral collision
fibres.  A singleton census would prove actual cross-frame coefficient
separation; a collision would identify the precise residual frame symmetry.
Either outcome still requires a physical allocation sidecar.

## 8. Exact companion and scope

Run

```bash
python 04-computation/lrc14_normalized_unit_bispectrum_thm2802.py
python -O 04-computation/lrc14_normalized_unit_bispectrum_thm2802.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_normalized_unit_bispectrum_thm2802.out.
```

The companion exhausts every Fourier-unit vector on `C_3/F_7` and
`C_5/F_11`.  The normalized bispectrum fibres have respectively

```text
12 classes of size 18,        2,000 classes of size 50,  (27)
```

exactly the scalar-translation orbit sizes.  It also checks all `676`
scalar/origin transforms of `(22)` over `F_53`, and verifies that `(24)` has
equal power spectra but distinct normalized bispectra.  It uses explicit
exception gates and no truth-bearing Python assertions.

No same-ancestry physical allocation, positive intertwiner, global frame
map, endpoint chronology, row exclusion, or LRC(14) conclusion is proved.

**Awaiting independent audit; not QED.**
