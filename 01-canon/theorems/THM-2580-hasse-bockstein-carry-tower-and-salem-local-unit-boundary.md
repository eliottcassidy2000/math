---
id: THM-2580
title: "Hasse-Bockstein carry tower and Salem local-unit boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the augmentation lattice of a prime cycle, the first r binomial Hasse
  moments are the complete integral obstructions to r-fold cyclic-difference
  or Cayley filling for 1<=r<=p-1; the Smith form is
  1^(p-1-r) + p^r.  At p=13, C^12 Lambda=13 Lambda and the tower is
  12-periodic up to multiplication by 13.  On THM-2571's one globally
  primitive carrier, every target difference is once but not twice Cayley
  fillable: its secondary class is a faithful affine C13 displacement
  cocycle, nonzero on all 468 off-diagonal owner-target edges and supported
  on all 72 product coordinates.  HYP-9060's two cyclotomic signatures have
  exact carry depths 2 and 3; its Salem signature has 13-primary depth 1 and
  extra 5^4 torsion.  These are signed coefficient/filter statements, not a
  positive physical filling, semantic endpoint, scalar-row exclusion, or
  LRC(14) conclusion.
source: root-holotopy-2026-07-28-hasse-salem-tower
depends_on:
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2571-deep-colour-cayley-filling-bockstein-and-norm-curvature-split
  - THM-2579-socle-flat-target-torsor-and-integral-difference-filling
related:
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - HYP-9060
script: 04-computation/lrc14_hasse_bockstein_carry_tower_thm2580.py
output: 05-knowledge/results/lrc14_hasse_bockstein_carry_tower_thm2580.out
script_sha256: 2c9db264406bff002ae687f08f092258661100e18d06c59fd270ddf97225d75c
output_sha256: 5220628af6991239fb01d8997c2ce2982d7da9ae853237e3015fc525fdc2fa8a
hash_basis: LF-normalized bytes
---

# THM-2580 -- relation signatures are an integral carry spectrum

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

THM-2571 found one 13-primary obstruction after rational Cayley filling;
THM-2579 found that target differences kill it.  The killed class is not the
end of the story.  It is the first page of a finite integral filtration.
The next Hasse moment detects **every** nonzero target displacement.

For the canonical carrier the three relevant rungs are

```text
absolute profile c^(kappa,b):        first class nonzero;

pair difference c^(kappa,b)-c^(kappa,b'):
                                      first class zero, second class nonzero;

full unnormalized target Fourier numerator:
                                      literal factor 13, hence depth >=12.
                                                                    (1)
```

Thus the target family is not merely a flat torsor.  Its edges carry a
faithful affine displacement cocycle one filtration layer deeper, while the
complete Fourier numerator collapses far down the same tower.  This is a
literal relation-and-carry spectrum, not an analogy.

## 1. Prime-cycle Hasse exact sequence

Let `p` be an odd prime, let `P` act on `Z^p` by

```text
(Pa)_m=a_(m+1),                                             (2)
```

and put

```text
Lambda={a in Z^p: sum_m a_m=0},
D=P-I.                                                      (3)
```

For `1<=j<=p-1`, define the binomial Hasse moment

```text
beta_j(a)=sum_(m=0)^(p-1) binom(m,j)a_m       in F_p.       (4)
```

Then for every `1<=r<=p-1` there is an exact sequence

```text
0 -> Lambda --D^r--> Lambda
  --(beta_1,...,beta_r)--> F_p^r -> 0.                     (5)
```

Equivalently,

```text
a in D^r Lambda  <=>  beta_1(a)=...=beta_r(a)=0.           (6)
```

The Smith form is

```text
SNF(D^r|Lambda)=1^(p-1-r) + p^r.                           (7)
```

### Proof

Represent a colour vector by

```text
A(z)=sum_m a_m z^m in Z[C_p].                              (8)
```

The augmentation lattice is the ideal `I=(z-1)`.  Modulo `p`, put
`t=z-1`; then

```text
F_p[C_p]=F_p[t]/(t^p),       Lambda/pLambda=(t).            (9)
```

With convention (2), `P=z^(-1)`, so `D` is `-t/(1+t)`, namely `t` times a
unit.  The coefficient of `t^j` in `A(1+t)` is exactly `beta_j(a)`.
Therefore `D^r Lambda` is killed by the first `r` Hasse coordinates modulo
`p`.

The map in (5) is onto: on `e_m-e_0`, `1<=m<=r`, its matrix is
`(binom(m,j))`, lower triangular with diagonal one.  Finally,
`|det(D|Lambda)|=p`, so `D^r Lambda` has index `p^r`.  The kernel of the
surjective Hasse map has the same index; hence the two lattices agree.  The
same index and mod-`p` corank `r` give (7). QED.

The theorem is functorial under base extension from `Z` to any free
`Z`-module, in particular `Z[zeta_7,zeta_13]`.  Origin translation changes
the ordered Hasse coordinates by an invertible triangular transformation,
so the simultaneous vanishing condition is chart invariant.

## 2. Cayley filling and the twelve-step periodicity

Put

```text
C=(P-I)(P+I)^(-1) on Lambda.                               (10)
```

Since

```text
det(P+I|Lambda)=Phi_p(-1)=1,                               (11)
```

`P+I` is an integral unimodular automorphism commuting with `D`.  Hence

```text
C^r Lambda=D^r Lambda                                      (12)
```

for every `r`, and (5)--(7) classify iterated Cayley filling as well.

At `p=13`, the cyclotomic uniformizer `pi=zeta_13-1` satisfies
`13=unit*pi^12`.  Since `Lambda` corresponds to `(pi)` and `C` is `pi`
times a unit,

```text
C^12 Lambda=13 Lambda,                                     (13)

C^(12q+r)Lambda=13^q C^r Lambda,
q>=0, 0<=r<12.                                             (14)
```

The first two post-period Smith controls are

```text
SNF(C^13)=13^11 + 169,
SNF(C^24)=169^12.                                          (15)
```

More generally, an integral filter

```text
F(P)=D^r u(P),             u(1) not=0 mod p,               (16)
```

has the same `p`-primary carry depth `r`: in the local algebra (9), `u(P)`
is a unit.  Other primes can still occur through
`Res(u,Phi_p)`.  If that resultant is `+-1`, the equality is global rather
than merely `p`-primary.  This local/global distinction is load-bearing in
the Salem application below.

## 3. The canonical secondary target cocycle

Use THM-2571's globally primitive carrier

```text
x_(ell,s,r),       ell in F_7, s,r in F_13,                (17)
```

and its deep cycles

```text
c_m^(kappa,b)=sum_(ell,s,r) x_(ell,s,r)
 zeta_7^(kappa ell) zeta_13^(b s+m r),                     (18)
```

for `kappa in F_7^*`, `b in F_13`.  The single global primitive clearing is
essential; no profile is separately reprimitivized.

THM-2571 proves

```text
beta_1(c^(kappa,b))=Omega Y(zeta_7^kappa) !=0,             (19)
```

independently of `b`, where

```text
Omega=sum_m m zeta_13^m=-epsilon^11,
epsilon=zeta_13-1.                                         (20)
```

For the second Hasse moment put

```text
H_2(r)=sum_m binom(m,2)zeta_13^(mr).                       (21)
```

In `F_13[epsilon]/(epsilon^12)`, exact finite differentiation gives

```text
H_2(r)=zeta_13^(2r)(zeta_13^r-1)^10
      =r^(-2)epsilon^10+r^(-2)(7r-5)epsilon^11.            (22)
```

Consequently, for `b!=b'`,

```text
(zeta_13^(bs)-zeta_13^(b's))H_2(r)
  =-(b-b')s r^(-2) Omega.                                  (23)
```

Define the septimal secondary carrier

```text
Z_ell=sum_(s,r!=0) x_(ell,s,r) s r^(-2) mod 13.            (24)
```

The exact reconstruction gives

```text
Z=(0,5,3,9,4,10,8).                                       (25)
```

Modulo `Phi_7`,

```text
Z_red=(5,10,8,1,9,2),
Z_red^(-1)=(4,8,5,8,0,3).                                 (26)
```

Thus `Z` is a unit.  Equations (23)--(26) give the closed formula

```text
beta_2(c^(kappa,b)-c^(kappa,b'))
 =-(b-b') Omega Z(zeta_7^kappa) !=0.                       (27)
```

The six owner factors are

```text
(5,10,8,1,9,2),   (4,8,9,1,7,12),  (10,7,6,2,5,1),
(3,6,7,11,8,12),  (9,5,4,12,6,1),  (8,3,5,12,4,11).       (28)
```

Every coordinate in (28) is nonzero.  Hence all
`6*binom(13,2)=468` off-diagonal target differences have full
`6*12=72` product-basis support.  The owner factors span rank `3` over
`F_13`; target displacement only rescales them.

By THM-2579, each difference `delta` already belongs to `C Lambda`.  It does
not belong to `C^2 Lambda` by (5) and (27).  If `delta=C y`, its unique
integral first primitive satisfies the exact orientation law

```text
beta_1(y)=-2 beta_2(delta) !=0.                             (29)
```

Indeed the symbol of `C` at `z=1+t` is `-t/(2+t)`; because `y` is already
augmentation zero, the `t^2` coefficient of `Cy` is
`-beta_1(y)/2`.  The control `delta=(z-1)^2`, `y=1-z^2` realizes the sign.

Therefore every canonical target edge is **once fillable but not twice
fillable**.  The class

```text
sigma_kappa(b,b')=beta_2(c^(kappa,b)-c^(kappa,b'))         (30)
```

is antisymmetric, additive along target paths, and nonzero exactly off the
diagonal.  It faithfully records affine `C_13` displacement but does not
choose an absolute target origin.

## 4. Why the complete target transform still collapses

For the common carrier (17), target orthogonality gives the literal identity

```text
sum_b zeta_13^(qb)c^(kappa,b)
 =13 * (the s=-q integral slice).                          (31)
```

Therefore every unnormalized target Fourier numerator lies in

```text
13 Lambda=C^12 Lambda.                                     (32)
```

In particular its first and second Hasse classes vanish.  The new secondary
edge class is thus affine/pairwise; it is not an absolute selector surviving
the complete unnormalized transform.

Division by `13` changes the lattice question.  On the one common globally
primitive carrier, (31) makes that quotient coefficientwise integral and
recovers a target-shift slice.  On arbitrary independent cokernel
representatives, no such division is an integral group-ring operation.
THM-2585 is reserved for the exact saturated-slice census.  It is not a
dependency here, and no survival claim from that reserved theorem is used.
The singleton hostile in THM-2579 already proves that normalization can
restore a first class killed by its numerator.

## 5. HYP-9060: cyclotomic depth versus Salem torsion

HYP-9060 supplied three reciprocal recurrence signatures.  Regard each as a
filter on the 13-cycle:

```text
F_1(x)=(x-1)^2(x+1),

F_2(x)=(x-1)^3(x+1)^2(x^2+1),

F_S(x)=(x-1)S(x),
S(x)=x^6-x^4-2x^3-x^2+1.                                  (33)
```

Their exact Smith forms on `Lambda` are

```text
SNF(F_1(P))=1^10 + 13^2,

SNF(F_2(P))=1^9  + 13^3,

SNF(F_S(P))=1^8 + 5^3 + 65.                               (34)
```

For the first two filters,

```text
Res(x+1,Phi_13)=Res(x^2+1,Phi_13)=1,                      (35)
```

so their carry depths `2` and `3` are global integral statements.  For the
Salem sextic,

```text
S(1)=-2,
Res(S,Phi_13)=625=5^4.                                    (36)
```

Thus `S(P)` is a 13-local unit and `F_S` has 13-primary depth `1`, but it
also creates four powers of `5` in the global cokernel.  A Cayley-fillable
target difference need not be `F_S`-fillable without a separate mod-5
check.  Calling the Salem column simply “fillable” would be false globally.

Over `C`, all three filters have the same kernel on the full 13-cycle: the
constant line.  None of the remaining factors vanishes at a nontrivial
13th root.  Their complex drift detector is therefore identical even though
their integral carry depths differ.  The real Salem eigenvalue is not a
13-cycle eigenvalue; its growth cannot amplify or physicalize a `C_13` Abel
normal.  Its actual contribution here is the separate 5-primary torsion in
(34), not a positive dynamical gain.

The canonical target edge (27) is consequently visible to the depth-2 and
depth-3 cyclotomic signatures but invisible to a depth-1 13-primary test.
This comparison is 13-primary; the Salem column retains the global mod-5
qualification just stated.

## 6. Exact companion and stopping boundary

Run

```bash
python3 04-computation/lrc14_hasse_bockstein_carry_tower_thm2580.py
python3 -O 04-computation/lrc14_hasse_bockstein_carry_tower_thm2580.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_hasse_bockstein_carry_tower_thm2580.out
```

byte-for-byte.  The exact companion verifies:

- all twelve Hasse kernels, ranks, and Smith forms, plus Cayley depths
  `12,13,24` and the sign control (29);
- all three HYP-9060 Smith forms, resultants, and one-dimensional complex
  kernels;
- a fresh interval reconstruction of all `1183` canonical carrier cells,
  its `1092` positive cells, primitive content `13`, and digest;
- the closed `H_2` identity, vectors `Z`, `Z^(-1)`, all six owner factors,
  their rank `3`, all `468/468` secondary edges, and full `72`-coordinate
  support;
- all `169` target-orthogonality cells giving the literal factor in (31).

There are `3599` explicit checks, none implemented with `assert`.

The theorem proves a complete prime-cycle integral carry filtration, its
canonical secondary target-edge application, and the exact local/global
meaning of the three old reciprocal signatures.  It does not construct a
positive first or second filling, a lawful physical recurrence filter, an
absolute semantic target/owner/root reference, a same-carrier Abel-normal
composition, a scalar relation, a row exclusion, or LRC(14). **QED.**
