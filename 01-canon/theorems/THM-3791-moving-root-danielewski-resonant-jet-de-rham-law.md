---
id: THM-3791
title: "Moving-root Danielewski resonant-jet de Rham law"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  For every
  smooth deformation c^n e=Sigma(c,b) with squarefree special fibre, the
  algebraic de Rham class of the physical symplectic form is exactly the
  vector of resonant coefficients in the canonical Hensel root jets modulo
  c^n, taken modulo constants.
  It is exact iff those coefficients agree.  A nonconstant resonant vector
  excludes every polynomial Darboux pair.  This unifies the exponent-one
  fixed-root obstruction, higher-exponent fixed-root exactness, and the
  moving-axis obstruction that closes THM-3789.
source: jc_quartic_c3_construct / THM-3789 moving-arm abstraction, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The proof reconstructs the canonical Hensel
  jets, polynomial-division arm charts, every overlap and transition, the
  Cech--de Rham hypercover rows, the truncated-simplex incidence quotient,
  and the local-primitive residue sign.  The exact companion verifies
  nonfactorized smooth Poisson packets and plane charts for n=1..5, incidence
  ranks through seven arms, all resonance controls, and the THM-3789
  specialization.  Normal and optimized executions byte-match the frozen
  transcript.  Independent hostile audit remains due.
depends_on:
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
related:
  - THM-3789-higher-pole-hermite-spectral-completion
script: 04-computation/jc2_moving_root_danielewski_resonant_jet_thm3791.py
output: 05-knowledge/results/jc2_moving_root_danielewski_resonant_jet_thm3791.out
script_sha256: d232d006065df8f65ae0785a3a1632ae629495805d32f7645c425628669753de
output_sha256: ed1b4fc68c37004acb3a1076bb069735d79dbd7a22011436a68b59023e16c1cb
semantic_sha256: 20a587137a4591a6b7b49d4a90d385f631d965522e35164c133a6372b1c916e4
hash_basis: raw LF bytes
---

# THM-3791 -- only the resonant moving-root jet survives in de Rham cohomology

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  This is an
unconditional structure theorem for smooth moving-root Danielewski surfaces.
It turns three previously separate phenomena into one coefficient law:

- fixed roots at exponent one carry the classical logarithmic obstruction;
- fixed roots at exponent at least two have an exact physical form; and
- a moving root can recreate the logarithmic obstruction at any exponent.

Let `k` be an algebraically closed field of characteristic zero.  Fix an
integer `n>=1` and a polynomial `Sigma(c,b) in k[c,b]` whose special fibre is

```text
Sigma(0,b)=gamma product_(i=1)^h (b-beta_i),
h>=1,                 gamma!=0,                 beta_i pairwise distinct. (1)
```

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

For each `beta_i`, Hensel lifting gives a unique polynomial jet

```text
beta_i^[n](c) in k[c],                 deg beta_i^[n]<n,
beta_i^[n](0)=beta_i,
Sigma(c,beta_i^[n](c))=0 mod c^n.                              (5)
```

Existence and uniqueness follow recursively because
`partial_b Sigma(0,beta_i)!=0`: after the coefficients below order `m` are
fixed, the coefficient of `c^m` determines the next root coefficient by one
nonzero linear equation.  Thus `(5)` is intrinsic data of the surface, not a
choice of a global factorization.

Define its resonant coefficient

```text
r_i=[c^(n-1)] beta_i^[n](c),                 r=(r_1,...,r_h).  (5a)
```

Then `omega` extends to a regular symplectic form on all of `X`, and there
is a canonical identification

```text
H^1_dR(X)=0,
H^2_dR(X)=k^h/k(1,...,1)                                    (6)
```

under which

```text
[omega]=[r].                                                    (7)
```

In particular,

```text
omega is exact
iff r_1=...=r_h.                                                (8)
```

If the resonant coefficients are not all equal, there are no
`P,Q in A` with `{P,Q}=1`.

The sign in `(7)` is fixed by the orientation `{c,b}=c^n` and the local
primitive chosen below.  Reversing the Poisson bracket reverses `[omega]`
but leaves `(8)` and every Darboux consequence unchanged.

## 1. Smoothness and the full moving-arm atlas

Let

```text
D=c^n e-Sigma(c,b).                                           (9)
```

If `c!=0`, then `D_e=c^n!=0`.  On `c=0`, equation `(9)` forces
`b=beta_i` for a unique `i`, and

```text
D_b=-partial_b Sigma(0,b)!=0                                  (10)
```

because the special fibre is squarefree.  Thus `X` is smooth.  Formula `(3)`
is the signed gradient packet of `(9)`, so it is symplectic.

The arm divisors are

```text
L_i=V(c,b-beta_i) ~= A1_e.                                    (11)
```

Retain one arm and delete the others:

```text
U_i=X \ union_(j!=i)L_j,
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
The `U_i` cover `X`.  For every set of at least two distinct indices,

```text
U_(i_0) intersection ... intersection U_(i_p)
 =D(c) ~= G_m,c x A1_b.                                       (16)
```

On a double overlap the exact transition is

```text
a_j=a_i+[beta_i^[n](c)-beta_j^[n](c)]c^(-n).                  (17)
```

Thus every squarefree deformation has precisely the same affine hypercover
as the fixed-root arm plane, but its singular shears carry the canonical
Hensel jets.  No global factorization of `Sigma(c,b)` is required.

The same cover also gives the exact intersection

```text
A=intersection_(i=1)^h k[c,a_i]              inside k(c,b).   (18)
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
H^1_dR(G_m x A1)=k[dc/c]                                      (22)
```

is exactly

```text
(r_j-r_i)[dc/c].                                               (23)
```

This proves the core principle: lower and nonresonant Hensel jets are
Hamiltonian corrections; terms of order at least `n` never enter the chart
transition.  Only the `(n-1)`-jet survives.

## 3. The Cech--de Rham quotient

For clarity, `(7)` is not inferred merely from a pole in one rational
primitive.  Apply the algebraic Cech--de Rham hypercover to `(12),(16)`.  Its
de Rham cohomology rows are

```text
H^q_dR(U_i)=k,0,0,...                         for q=0,1,2,...,
H^q_dR(D(c))=k,k[dc/c],0,...                  for q=0,1,2,... . (24)
```

The `q=0` row is the full simplex complex and is exact beyond degree zero.
The `q=1` row starts at edges because `H^1_dR(U_i)=0`.  Its edge-to-triangle
kernel is the image of the usual vertex-to-edge difference map, but there is
no vertex term in this truncated row to quotient it further.  Hence

```text
E_2^(1,1)=k^h/k(1,...,1).                                     (25)
```

There is no possible incoming or outgoing higher differential, so `(25)`
is `H^2_dR(X)`.  The local primitives `(19)` send `omega` to the edge cocycle
`(r_j-r_i)dc/c`; under `(25)` this is precisely the vertex class `[r]`.
This proves `(6),(7)` and therefore `(8)`.

## 4. The Darboux obstruction

The bracket convention in `(3),(4)` gives, for all `P,Q in A`,

```text
dP wedge dQ={P,Q}omega.                                       (26)
```

If `{P,Q}=1`, then

```text
omega=dP wedge dQ=d(P dQ),                                    (27)
```

so `[omega]=0`.  A nonconstant vector `(5a)` contradicts `(7)`.  This proves
the all-degree Darboux nonentry.

## 5. Exact controls and sharp boundaries

### 5.1 Fixed roots

If `Sigma` is independent of `c` and `n=1`, its Hensel jets are the constant
roots and

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

For one arm, `k^h/k(1)=0`, so every physical form is exact.  If two special-
fibre roots collide, the surface can become singular and the arm cover no
longer has the form `(15),(16)`; squarefreeness is load-bearing.  Replacing
`b` by `b+u(c)` translates every Hensel jet by the same polynomial and changes
`r` only by a constant vector, so `(7)` is coordinate-invariant.  Over a
non-algebraically-closed constant field the jets form a Galois set and need a
trace/isotypic statement; that descent is **OPEN**.

The exact companion verifies `(3),(9)--(17)` for nonfactorized deformations
at `n=1,...,5`, the residue extraction `(20)--(23)`, the truncated-simplex
ranks in `(25)`, and every control `(28)--(30)`.  Normal and optimized
executions byte-match the frozen transcript.  **QED, conditional only on
independent hostile audit.**
