---
id: THM-2375
title: "Gaussian angular complete-line charge tomography"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every
  finite abelian unitary action, demodulating by one labelled character
  and averaging the squared difference over the complete group orbit
  gives exactly twice the energy outside that isotypic component. The
  full labelled bank recovers every isotypic norm. Applied to circular
  Gaussian U(1)-charge and N>2d angular samples, it recovers every
  charge-sector norm of a degree-d polynomial and decides strict
  one-sided charge support. N>2d, character labels, and circular
  orthogonality are sharp: explicit aliasing, unlabelled, and
  anisotropic hostiles break the corresponding conclusions. These are
  Hermitian L2 measurements, not scalar Gaussian moments. Z and
  conjugate(Z) have identical zero scalar-moment sequences but opposite
  labelled tomography, so no new NC2 proof or LRC transfer follows.
source: codex-2026-07-25-gaussian-angular-charge-tomography
depends_on: []
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-1520-gmc2-one-sided-branch-closed-by-charge-telescoping
  - THM-2369-complete-line-target-dirichlet-and-balanced-observable-no-go
script: 04-computation/gaussian_angular_complete_line_thm2375.py
output: 05-knowledge/results/gaussian_angular_complete_line_thm2375.out
script_sha256: e6960ade066052a7a66afd169e21086b15419887133a38fd389b3d44a9357a43
output_sha256: 63e675ad6ae53fb4eb7599807d66443b930cdcaf96a5cf0477fea7ee6fa951ee
hash_basis: working-tree bytes (LF)
---

# THM-2375 -- demodulate, then average the complete angular orbit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2369 uses complete target-character lines to turn a covariant
difference into an exact support-energy multiplier. The same operation
has a reusable representation-theoretic form. Its Gaussian specialization
recovers the entire angular charge spectrum, but only because it is given
conjugated `L2` data which scalar Gaussian moments do not contain.

## 1. Finite-abelian complete-orbit tomography

Let a finite abelian group `H`, with `|H|>1`, act unitarily on a complex
Hilbert space by `U_h`. Write the orthogonal character decomposition

```text
f=sum_(chi in H^) f_chi,

U_h f_chi=chi(h)f_chi.                              (1)
```

For a labelled character `lambda`, define the complete-orbit
demodulated Dirichlet energy

```text
D_lambda(f)
 =1/|H| sum_(h in H)
   ||conjugate(lambda(h))U_h f-f||^2.               (2)
```

Character orthogonality gives the exact identity

```text
D_lambda(f)
 =2 sum_(chi!=lambda)||f_chi||^2
 =2(||f||^2-||f_lambda||^2).                       (3)
```

Indeed, the multiplier on `f_chi` is

```text
conjugate(lambda(h))chi(h)-1.
```

Its averaged squared magnitude is zero for `chi=lambda` and two
otherwise. Orthogonality of the isotypic components then proves (3).

Summing (3) over every `lambda in H^` gives

```text
sum_lambda D_lambda=2(|H|-1)||f||^2.               (4)
```

Consequently the labelled complete bank recovers

```text
||f||^2
 =1/(2(|H|-1))sum_lambda D_lambda,

||f_lambda||^2
 =||f||^2-D_lambda/2.                              (5)
```

The bank recovers isotypic norms and support. It does not recover the
vectors `f_lambda` or their relative phases.

## 2. Circular Gaussian charge sectors

Let `Z` be a standard circular complex Gaussian and put
`W=conjugate(Z)`. For a polynomial of total degree at most `d`, write

```text
P(Z,W)=sum_(q=-d)^d P_q(Z,W),                      (6)
```

where `P_q` is the sum of the monomials `Z^a W^b` with charge

```text
q=a-b.
```

The circle action

```text
(R_theta P)(Z,W)
 =P(exp(i theta)Z,exp(-i theta)W)
```

is unitary for the Hermitian Gaussian norm

```text
||P||_2^2=E[|P(Z,W)|^2].
```

The charge sectors in (6) are orthogonal. Choose an integer `N>=2` and
`zeta_N=exp(2 pi i/N)`. Put

```text
R_a P=sum_q zeta_N^(a q)P_q,

D_r
 =1/N sum_(a=0)^(N-1)
   ||zeta_N^(-a r)R_aP-P||_2^2,    r in Z/N.       (7)
```

Without any no-alias assumption, (3) gives

```text
W_r=sum_(q congruent r mod N)||P_q||_2^2,

D_r=2(||P||_2^2-W_r).                              (8)
```

Thus the exact object recovered at scale `N` is the residue-binned
charge spectrum.

## 3. The sharp no-alias bank

If

```text
N>2d,                                               (9)
```

then the charges `-d,...,d` have distinct residues. Equations (5) and
(8) become

```text
||P||_2^2
 =1/(2(N-1))sum_(r in Z/N)D_r,

||P_q||_2^2
 =||P||_2^2-D_(q mod N)/2.                         (10)
```

If only the `2d+1` possible charge labels are retained, then for `d>=1`

```text
||P||_2^2
 =1/(4d)sum_(r=-d)^d D_(r mod N).                  (11)
```

The band-only `{0}` sub-bank at `d=0` is identically zero and needs the
norm as an external datum. A full `Z/N` bank with `N>1` still recovers
the norm through (10).

Circular Gaussian measure has full support, so

```text
P_q!=0
 iff ||P_q||_2^2>0
 iff D_(q mod N)<2||P||_2^2.                       (12)
```

Therefore the labelled bank recovers the exact charge support. In
particular it decides whether that support is strictly positive,
strictly negative, nonnegative, nonpositive, or genuinely two-sided.
It still loses every radial coefficient and every phase inside or
between the nonzero sectors.

## 4. Three sharp boundaries

### 4a. Aliasing at `N=2d`

At the first excluded scale, the charges `d` and `-d` coincide. The two
polynomials

```text
P_mix=Z^d+W^d,

P_plus=sqrt(2)Z^d                                  (13)
```

have the same total norm and the same residue-binned norm at every
label modulo `N=2d`. Hence they have identical complete `D` banks, while
the first is two-sided and the second strictly positive-charge. The
strict inequality in (9) is sharp.

### 4b. Unlabelled banks

At `N=5`, the polynomials

```text
Z+Z^2/sqrt(2),

W+Z^2/sqrt(2)                                      (14)
```

have two charge sectors of norm one. Their unlabelled multisets of
`D_r` agree, but their supports are respectively `{1,2}` and `{-1,2}`.
The labels in (10) are load-bearing.

### 4c. Noncircular laws

If distinct charges are not orthogonal, the exact no-alias formula is
instead

```text
D_r
 =||sum_(q!=r)P_q||^2+sum_(q!=r)||P_q||^2.         (15)
```

Cross terms remain. For independent centered real Gaussians with

```text
E[X^2]=2,    E[Y^2]=1,

Z=X+iY,      W=X-iY,      P=Z+W,
```

one gets

```text
D_0=14,                   2||P||^2=16.             (16)
```

Thus circular unitarity, rather than the Gaussian factorial law, is the
essential hypothesis. The theorem extends verbatim to any faithful
circular `L2` law.

## 5. Why this is not scalar-moment progress

For every `m>=1`,

```text
E[Z^m]=E[W^m]=0.                                   (17)
```

The two scalar moment sequences are identical, even with the same norm.
Yet their labelled banks have

```text
D^Z_1=0,                     D^W_(-1)=0,           (18)
```

with the other relevant entries equal to two. Scaling `Z` shows that
the zero scalar-moment sequence does not even determine the magnitudes
in the Hermitian bank.

THM-2022 proves NC2 by preserving a whole balanced Wick face through
Kummer/Frobenius descent. It classifies the union of the two one-sided
loci; it does not recover their label or sector norms. The present
theorem assumes richer conjugated measurements and supplies no map from

```text
(E[P^m])_(m>=1)
```

to (7). It is therefore not a new proof of NC2 and does not shorten
THM-2022. It also gives no LRC relation-sum identity; MISTAKE-235 remains
fully operative.

## 6. Relation to complete-line target energy

The common reusable operation is

```text
demodulate by a prescribed character
  -> average a squared difference over a complete subgroup
  -> obtain an exact 0/2 isotypic multiplier.       (19)
```

THM-2369 applies (19) to the two target lines in `F_13^2`; this theorem
applies it to angular charge. The representation and measurements differ,
so neither application implies the other.

The operation is useful when the complete-orbit differences are lawful
observables. It is counterindicated when only balanced scalar moments,
unlabelled energies, or a nonunitary action are available.

## 7. Exact companion

The dependency-free `Fraction` companion:

- checks one nontrivial aliased `C_5` instance and full-bank recovery;
- verifies all-occupied and omitted-sector controls at `N=7>2d=6`;
- checks the exact `N=2d` aliasing hostile;
- checks the labelled/unlabelled support hostile;
- checks `200` nonzero monomial-charge instances underlying (17); and
- evaluates the anisotropic Gram hostile as `14` versus `16`.

Run

```bash
python3 04-computation/gaussian_angular_complete_line_thm2375.py
python3 -O 04-computation/gaussian_angular_complete_line_thm2375.py
```

Both transcripts must match

```text
05-knowledge/results/gaussian_angular_complete_line_thm2375.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent hostile audit rederived the finite-unitary identity, Gaussian
recovery formulas, and all three sharp boundaries; normal, optimized, and
stored transcripts agree. QED.
