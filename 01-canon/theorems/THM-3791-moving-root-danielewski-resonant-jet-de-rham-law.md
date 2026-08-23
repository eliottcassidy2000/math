---
id: THM-3791
title: "Moving-root Danielewski resonant-jet de Rham law"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  For every
  smooth deformation c^n e=Sigma(c,b) with squarefree special fibre, the
  algebraic de Rham class of the physical symplectic form is exactly the
  resonant coefficient of the universal Hensel root in the finite etale
  special-fibre algebra, modulo base-field scalars.  It is exact iff that
  coefficient is scalar.  A nonscalar resonant class
  excludes every polynomial Darboux pair.  This unifies the exponent-one
  fixed-root obstruction, higher-exponent fixed-root exactness, and the
  moving-axis obstruction that closes THM-3789.
source: jc_quartic_c3_construct / THM-3789 moving-arm abstraction, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The proof reconstructs the canonical Hensel
  jets, polynomial-division arm charts, every overlap and transition, the
  Cech--de Rham hypercover rows, the truncated-simplex incidence quotient,
  finite-etale descent, and the local-primitive residue sign.  The exact
  companion verifies nonfactorized smooth Poisson packets and plane charts
  for n=1..5, an irreducible quadratic finite-etale control, incidence ranks
  through seven arms, all resonance controls, and the THM-3789 specialization.
  Normal and optimized executions byte-match the frozen transcript.
  Independent hostile audit remains due.
depends_on:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
related:
  - THM-3789-higher-pole-hermite-spectral-completion
script: 04-computation/jc2_moving_root_danielewski_resonant_jet_thm3791.py
output: 05-knowledge/results/jc2_moving_root_danielewski_resonant_jet_thm3791.out
script_sha256: 6ee083b05f540011bbed618adf241c708772709d25fc472e26c8e0c44b471eec
output_sha256: c6fd123e77cf8b297bb5ee252877f6a37f09aab090449ab9603d2858035a8c9a
semantic_sha256: 54c8a0a181045f8bc875a94e0b8d7b4800996ccfcb5891c28a5605727851c388
hash_basis: raw LF bytes
---

# THM-3791 -- only the resonant moving-root jet survives in de Rham cohomology

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  This is an
unconditional structure theorem for smooth moving-root Danielewski surfaces.
It turns three previously separate phenomena into one coefficient law:

- fixed roots at exponent one carry the classical logarithmic obstruction;
- fixed roots at exponent at least two have an exact physical form; and
- a moving root can recreate the logarithmic obstruction at any exponent.

Let `k` be any field of characteristic zero.  Fix an integer `n>=1` and a
polynomial `Sigma(c,b) in k[c,b]` whose special fibre

```text
Sigma_0(b)=Sigma(0,b)                                          (1)
```

is squarefree of degree `h>=1`.

Put

```text
A=k[c,b,e]/(c^n e-Sigma(c,b)),                 X=Spec A.       (2)
```

Give `X` the Jacobian Poisson bracket oriented by

```text
{c,b}=c^n,
{b,e}=n c^(n-1)e-partial_c Sigma,
{c,e}=partial_b Sigma.                                        (3)
```

On `c!=0`, its inverse form in the repository convention is

```text
omega=dc wedge db/c^n.                                        (4)
```

Let

```text
E=k[b_bar]=k[b]/(Sigma_0(b))                                   (5)
```

be the finite etale special-fibre algebra.  Since
`partial_b Sigma_0(b_bar)` is a unit in `E`, Hensel recursion gives a unique
universal root jet

```text
b_hat(c)=u_0+u_1c+...+u_(n-1)c^(n-1) in E[c]/(c^n),
u_0=b_bar,
Sigma(c,b_hat(c))=0 mod c^n.                                  (5a)
```

After the coefficients below order `m` are fixed, the coefficient of `c^m`
determines `u_m` by one linear equation with that unit coefficient.  Thus
`(5a)` is intrinsic data of the surface, not a choice of a splitting field or
global factorization.

Define the resonant Hensel element

```text
rho=u_(n-1) in E.                                              (5b)
```

Then `omega` extends to a regular symplectic form on all of `X`, and there
is a canonical identification

```text
H^1_dR(X)=0,
H^2_dR(X/k)=E/k.                                               (6)
```

under which

```text
[omega]=rho mod k.                                             (7)
```

In particular,

```text
omega is exact
iff rho in k.                                                   (8)
```

If the resonant Hensel element is not scalar, there are no
`P,Q in A` with `{P,Q}=1`.

The sign in `(7)` is fixed by the orientation `{c,b}=c^n` and the local
primitive chosen below.  Reversing the Poisson bracket reverses `[omega]`
but leaves `(8)` and every Darboux consequence unchanged.

## 1. Smoothness and the full moving-arm atlas

Let

```text
D=c^n e-Sigma(c,b).                                           (9)
```

If `c!=0`, then `D_e=c^n!=0`.  On `c=0`, equation `(9)` and the
squarefreeness of `Sigma_0` give

```text
Sigma_0(b)=0,                 D_b=-partial_b Sigma_0(b)!=0.    (10)
```

Thus `X` is smooth.  Formula `(3)` is the signed gradient packet of `(9)`, so
it is symplectic.

The arm calculation is most transparent after a faithfully flat scalar
extension.  Fix an algebraic closure `K` of `k`.  Then

```text
E tensor_k K ~= product_(i=1)^h K,                             (10a)
```

where the factors evaluate `b_bar` at the distinct roots `beta_i` of
`Sigma_0`.  The universal jet `(5a)` correspondingly splits into its unique
Hensel lifts

```text
beta_i^[n](c) in K[c],            deg_c beta_i^[n]<n.          (10b)
```

Write `X_K=X times_k K`.  Everything through `(25)` is first computed on
`X_K`; Section 3 then descends the result canonically to `k`.

The arm divisors are

```text
L_i=V(c,b-beta_i) in X_K,                 L_i ~= A1_(K),e.     (11)
```

Retain one arm and delete the others:

```text
U_i=X_K \ union_(j!=i)L_j,
a_i=[b-beta_i^[n](c)]/c^n.                                    (12)
```

Divide `Sigma` by the monic linear polynomial in `(12)`.  The Hensel
congruence `(5)` makes the remainder divisible by `c^n`, so uniquely

```text
Sigma(c,b)=[b-beta_i^[n](c)]S_i(c,b)+c^n R_i(c).               (13)
```

On `U_i`, the two expressions

```text
a_i=[b-beta_i^[n](c)]/c^n=[e-R_i(c)]/S_i(c,b)                 (14)
```

prove regularity: `D(c)` and `D(S_i)` cover `U_i`.  The mutually inverse
chart formulas are

```text
U_i ~= A2_(c,a_i),
b=beta_i^[n](c)+c^n a_i,
e=a_i S_i(c,b)+R_i(c).                                        (15)
```

At `c=0`, `S_i(0,beta_i)=partial_b Sigma(0,beta_i)!=0`, so `(15)` retains
the whole arm.
The `U_i` cover `X_K`.  For every set of at least two distinct indices,

```text
U_(i_0) intersection ... intersection U_(i_p)
 =D(c) ~= G_m,K,c x A1_K,b.                                   (16)
```

On a double overlap the exact transition is

```text
a_j=a_i+[beta_i^[n](c)-beta_j^[n](c)]c^(-n).                  (17)
```

Thus every squarefree deformation has precisely the same affine hypercover
as the fixed-root arm plane, but its singular shears carry the canonical
Hensel jets.  No global factorization of `Sigma(c,b)` is required.

The same cover also gives the exact split intersection

```text
A tensor_k K=intersection_(i=1)^h K[c,a_i]    inside K(c,b).   (18)
```

## 2. Local primitives and the unique resonant coefficient

In the `i`-chart, `(15)` gives

```text
omega=dc wedge da_i=d eta_i,                 eta_i=-a_i dc.    (19)
```

On an overlap,

```text
eta_j-eta_i=[beta_j^[n](c)-beta_i^[n](c)]c^(-n)dc.           (20)
```

For a monomial `c^m` in a root difference, the corresponding term in `(20)`
is `c^(m-n)dc`.  It is rationally exact unless

```text
m=n-1,                                                         (21)
```

when it is `dc/c`.  Therefore the class of `(20)` in

```text
H^1_dR(G_m,K x A1_K)=K[dc/c]                                  (22)
```

is exactly, with `r_i` the image of `rho` in the `i`-th factor of `(10a)`,

```text
(r_j-r_i)[dc/c].                                               (23)
```

This proves the core principle: lower and nonresonant Hensel jets are
Hamiltonian corrections; terms of order at least `n` never enter the chart
transition.  Only the `(n-1)`-jet survives.

## 3. The Cech--de Rham quotient

For clarity, `(7)` is not inferred merely from a pole in one rational
primitive.  Apply the algebraic Cech--de Rham hypercover of `X_K` from
`(12),(16)`.  Its
de Rham cohomology rows are

```text
H^q_dR(U_i/K)=K,0,0,...                       for q=0,1,2,...,
H^q_dR(D(c)/K)=K,K[dc/c],0,...                for q=0,1,2,... . (24)
```

The `q=0` row is the full simplex complex and is exact beyond degree zero.
The `q=1` row starts at edges because `H^1_dR(U_i)=0`.  Its edge-to-triangle
kernel is the image of the usual vertex-to-edge difference map, but there is
no vertex term in this truncated row to quotient it further.  Hence

```text
E_2^(1,1)=K^h/K(1,...,1).                                     (25)
```

There is no possible incoming or outgoing higher differential, so `(25)`
is `H^2_dR(X_K/K)`.  The local primitives `(19)` send `omega` to the edge cocycle
`(r_j-r_i)dc/c`; under `(25)` this is precisely the vertex class `[r]`.

This split residue identification is natural under further scalar extension:
the two pullbacks to `K tensor_k K` agree because uniqueness of the Hensel
jet and polynomial division identify the pulled-back roots, charts, and
residue entries.  In particular it carries the canonical faithfully flat
descent datum.  Algebraic de Rham cohomology of a smooth finite-type scheme
commutes with extension of the ground field, while

```text
(E/k) tensor_k K ~= K^h/K(1,...,1).                           (25a)
```

Faithfully flat descent therefore turns the split residue map into the
canonical `k`-linear isomorphism in `(6)`.  Under it, the base change of
`rho mod k` is exactly `[r]`; faithful flatness gives `(7)` over `k`.  The
same argument applied in degree one gives `H^1_dR(X/k)=0`.  This proves
`(6)--(8)` without choosing roots over the base field.

## 4. The Darboux obstruction

The bracket convention in `(3),(4)` gives, for all `P,Q in A`,

```text
dP wedge dQ={P,Q}omega.                                       (26)
```

If `{P,Q}=1`, then

```text
omega=dP wedge dQ=d(P dQ),                                    (27)
```

so `[omega]=0`.  A nonscalar `rho in E` contradicts `(7)`.  This proves the
all-degree Darboux nonentry.

## 5. Exact controls and sharp boundaries

### 5.1 Fixed roots

If `Sigma` is independent of `c` and `n=1`, its Hensel jets are the constant
roots after splitting and

```text
r=(beta_1,...,beta_h),                                        (28)
```

which is nonconstant when `h>=2`.  This is exactly the exponent-one
logarithmic obstruction of THM-3600, with the bracket orientation reversed.

For the same fixed-root surface with `n>=2`, one has `r=0`; the form is exact.  This
recovers THM-3600's higher-exponent primitive.  More generally, a common
`c^(n-1)` translation of the Hensel jets changes `r` by a constant vector and
has no effect on the class.

### 5.2 The THM-3789 moving axis

The surface of THM-3789 has `n=3` and Hensel root jets

```text
beta_0^[3](c)=c^2,                 beta_i^[3](c)=lambda_i.      (29)
```

Thus

```text
r=(1,0,...,0),                 [omega]!=0.                      (30)
```

Equation `(30)` is the exact all-support obstruction used there.

### 5.3 Boundaries

For one arm, `E=k` and `E/k=0`, so every physical form is exact.  If two special-
fibre roots collide, the surface can become singular and the arm cover no
longer has the form `(15),(16)`; squarefreeness is load-bearing.  Replacing
`b` by `b+u(c)` translates every Hensel jet by the same polynomial and changes
`rho` only by the scalar `[c^(n-1)]u(c)`, so `(7)` is coordinate-invariant.
The finite-etale formulation `(5)--(8)` includes split, nonsplit, and mixed
special fibres uniformly.  What remains outside the theorem is the ramified
boundary where `Sigma_0` is not squarefree; there `E` is not etale and the
arm hypercover degenerates.

### 5.4 A genuinely nonsplit control

Over `k=Q`, take `n=3` and

```text
Sigma(c,b)=b^2+1+cb+c^2,
E=Q[B]/(B^2+1).                                                (31)
```

Direct Hensel recursion gives

```text
b_hat(c)=B-c/2+(3B/8)c^2 mod c^3,
rho=3B/8 notin Q.                                             (32)
```

Thus `[omega]!=0` and this nonsplit surface has no polynomial Darboux pair.
After adjoining `i`, `(32)` becomes the two-entry residue vector
`(3i/8,-3i/8)` modulo constants, exactly as finite-etale descent predicts.

The exact companion verifies `(3),(9)--(17)` for nonfactorized deformations
at `n=1,...,5`, the residue extraction `(20)--(23)`, the truncated-simplex
ranks in `(25)`, the nonsplit recursion `(31),(32)`, and every split control
`(28)--(30)`.  Normal and optimized
executions byte-match the frozen transcript.  **QED, conditional only on
independent hostile audit.**
