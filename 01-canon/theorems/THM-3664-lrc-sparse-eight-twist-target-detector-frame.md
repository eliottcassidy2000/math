---
id: THM-3664
title: "LRC sparse eight-twist target detector frame"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  After choosing any basis of THM-2334's target-twist group, one transplanted
  copy of THM-3661's signed eight-point mask detects exactly the nonconstant
  twist profile.  Hence a nonzero target fibre exists if and only if one of
  169 signed eight-twist imbalances is nonzero, and the whole nonconstant
  profile is reconstructed from those imbalances.  Exact projective-line
  norms give a quantitative l2 frame bound.  This is a typed criterion on
  the target-twist group, not an identification with the two-current chart
  and not a proof that any covering-row imbalance is nonzero.
source: kps-s192 / THM-2334 x THM-3661 typed Fourier continuation, 2026-08-21
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-3661-lrc-exceptional-detector-simple-spectrum-convolution-rigidity
related:
  - THM-3660-lrc-exceptional-leakage-functional-and-fourteen-edge-boundary
  - THM-3662-lrc-eleven-cell-exceptional-flux-and-high-digit-variation-gate
  - THM-3663-lrc-minimal-carry-closed-exceptional-observer-bank
  - THM-3665-lrc-support-minimal-three-twist-target-detector
script: 04-computation/lrc_sparse_eight_twist_detector_frame_thm3664.py
output: 05-knowledge/results/lrc_sparse_eight_twist_detector_frame_thm3664.out
script_sha256: f60472d32b1ce92fba3ddedf82b3cef98dc6b45e7eb82e7ee49aecb5561fda89
output_sha256: 092490d1ac2903287af1208803834cc8fe79ade67ee77db04a588de7a98acc14
semantic_sha256: 24c09ef5a3681ccbf9878276e1e7fda2efabc947b0e86c06a6bf594e088561a6
hash_basis: raw LF bytes
---

# THM-3664 -- eight twists detect the complete nonzero target current

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  THM-2334 reduces target survival to nonconstancy of a function on
169 twists.  The present theorem replaces its dense variance test by a
sparse local witness involving four plus and four minus twists.

**Sharpening.**  THM-3665 proves that three support sites are optimal for an
abstract linear twist detector.  The continuing value of the present mask is
its inheritance from the physical exceptional-address pattern, its exact
cyclotomic norm atlas, and its comparison with THM-3660--3663; it is not
support-minimal on the freely based target group.

## 1. The typed transplant

Fix one positive strict shallow-owner word stratum in the scope of THM-2334.
Its target quotient and twist group are

```text
G isomorphic to F13^2,
T=G^ isomorphic to F13^2.                            (1)
```

Let `H:T->C` be the 169 boundary twist values in THM-2334 (41), and let
`A:G->C` be their inverse Fourier transform:

```text
A(q)=1/169 sum_(ell in T) e_13(-<ell,q>)H(ell).      (2)
```

Choose any group isomorphism

```text
theta:T -> F13^2.                                    (3)
```

This is only a basis of the already typed twist group.  It does **not**
identify a relation character with a two-current ancestry address.

On `F13^2`, use the signed mask of THM-3661:

```text
g=1_(X_+)-1_(X_-),

X_+={(12,0),(0,11),(6,5),(9,3)},
X_-={(0,12),(12,1),(6,7),(3,9)}.                    (4)
```

Pull it back to `T` by `gamma=g o theta` and define

```text
D_theta(s)=(gamma*H)(s)
 = sum_(x in X_+) H(s-theta^(-1)(x))
   -sum_(x in X_-) H(s-theta^(-1)(x)).              (5)
```

Thus every value of `D_theta` is a signed sum of exactly eight twists.

## 2. Exact sparse target equivalence

Use the Fourier convention

```text
f_hat(q)=sum_(ell in T)e_13(-<ell,q>)f(ell).         (6)
```

The basis (3) induces a dual basis on `G`.  The convolution theorem and (2)
give

```text
(D_theta)_hat(q)=gamma_hat(q)H_hat(q)
                =169 gamma_hat(q)A(q).              (7)
```

THM-3661 proves in characteristic zero that

```text
gamma_hat(q)=0 iff q=0,                              (8)
```

while the mean-zero identity for (4) gives `gamma_hat(0)=0`.  Therefore the
following statements are equivalent:

1. `A(q)!=0` for some nonzero target vector `q`;
2. `H` is nonconstant;
3. `D_theta` is not identically zero;
4. for some `s in T`, the signed eight-twist imbalance (5) is nonzero.

This is an if-and-only-if theorem, not merely a sufficient certificate.
In particular, one successful eight-value calculation proves nonzero target
survival without evaluating the dense sum of squares in THM-2334 (42).

## 3. Exact reconstruction

Define the cyclotomic inverse kernel on `T` by

```text
h_hat(0)=0,
h_hat(q)=1/gamma_hat(q), q!=0.                       (9)
```

Writing

```text
H_bar=1/169 sum_(ell in T)H(ell),                   (10)
```

Fourier inversion gives

```text
h*D_theta=H-H_bar.                                  (11)
```

Equivalently, every individual nonzero target fibre can be recovered as

```text
A(q)=(D_theta)_hat(q)/(169 gamma_hat(q)), q!=0.     (12)
```

Hence the 169 sparse imbalances retain precisely all nontrivial target
information and discard precisely the unavoidable constant mode.

## 4. Exact projective-line norms

The nonzero frequencies in `F13^2` split into fourteen Galois orbits, one
for each projective line.  Take representatives

```text
(1,t), 0<=t<=12, and (0,1).                         (13)
```

For a representative `q`, put

```text
P_q(z)=sum_x g(x)z^(-q.x) mod (z^13-1).             (14)
```

The exact norm `Norm_(Q(zeta_13)/Q)(P_q(zeta_13))`, computed as either the
resultant with `Phi_13` or the determinant of multiplication by `P_q`, is

```text
13 * (
  1, 64^2, 1, 131^2, 1, 27^2, 103^2,
  157^2, 1, 53^2, 27^2, 13^2, 1, 1
).                                                   (15)
```

In unreduced form the list is

```text
(13,53248,13,223093,13,9477,137917,
 320437,13,36517,9477,2197,13,13).                  (16)
```

Their product is the exact characteristic-zero determinant on the
augmentation ideal:

```text
18258988969361052598805681298174465892061184
 =4273053822427357700928^2.                         (17)
```

Reduction of (17) modulo THM-3661's pinned prime is

```text
669422013050837354847410,                           (18)
```

exactly the previously certified finite-field determinant.  Thus (15)--(18)
also give an independent exact lift of the nonvanishing part of THM-3661.

## 5. A quantitative frame, not only a detector

Every multiplier is a signed sum of eight roots of unity, so

```text
|gamma_hat(q)|<=8.                                  (19)
```

The other eleven members of a nonzero frequency's Galois orbit obey the
same bound.  Since its line norm is at least 13 by (15),

```text
|gamma_hat(q)|>=13/8^11, q!=0.                      (20)
```

Parseval and (7) now give the explicit stable-frame inequality

```text
(169/8^22) sum_ell |H(ell)-H_bar|^2
 <= sum_s |D_theta(s)|^2
 <=64 sum_ell |H(ell)-H_bar|^2.                     (21)
```

Using THM-2334 (42), this is equivalently

```text
(169^2/8^22) sum_(q!=0)|A(q)|^2
 <=sum_s |D_theta(s)|^2
 <=169*64 sum_(q!=0)|A(q)|^2.                       (22)
```

The lower constant is intentionally crude but fully explicit and uniform.

## 6. Basis covariance and the honest gauge boundary

Changing (3) by `M in GL(2,13)` replaces `g` by a linear pullback and merely
permutes its nontrivial Fourier multipliers.  Thus (7)--(22) hold for every
basis.  Exact enumeration gives

```text
|GL(2,13)|=26208,
number of distinct transformed signed masks=26208.  (23)
```

So the signed mask has trivial linear stabilizer.  There is no hidden
canonical basis: basis choice is a genuine gauge, even though every choice
gives an equivalent eight-twist test.  This is useful freedom for a future
proof, because one may choose a target basis that makes one proposed
imbalance easiest to analyze.  It is not permission to call the eight
selected twists the physical exceptional addresses of THM-3660.

## 7. Consequence and remaining frontier

For every THM-2334 stratum, the previously open unrestricted target step
can now be stated in one physically inherited sparse exact form:

```text
find one basis theta and one centre s for which (5) is nonzero.      (24)
```

Failure for all `s` in one basis is equivalent to failure in every basis
and to concentration of the full aggregate at `q=0`.  A promising next
move is therefore to choose `theta` so the eight coordinate translations
in (5) have the simplest interval-boundary geometry, then seek a unique
jump, sign-coherent sector, or first nonzero Taylor coefficient.

Nothing here proves (24) on a hypothetical covering row.  The theorem also
does not address THM-2334's all-`91`-unit aggregate `B(q)`, preserve visible
height or terminal phase, align relation twists with ancestry digits, or
close LRC(14).

## 8. Exact companion

Reproduce with

```bash
python3 -B 04-computation/lrc_sparse_eight_twist_detector_frame_thm3664.py
python3 -B -O 04-computation/lrc_sparse_eight_twist_detector_frame_thm3664.py
```

Both streams match the stored transcript.  The assertion-free companion
source-pins THM-2334 and THM-3661, verifies all `169^2` profile-basis transfer
identities, computes the fourteen integer cyclotomic norms by a
fraction-free multiplication determinant, lifts the augmentation
determinant, and exhausts all `26208` basis changes.  **QED.**
