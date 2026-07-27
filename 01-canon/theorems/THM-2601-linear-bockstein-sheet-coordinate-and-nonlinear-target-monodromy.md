---
id: THM-2601
title: "Linear Bockstein sheet coordinate and nonlinear target monodromy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  One explicit F_13-linear functional of THM-2585's six-dimensional
  septimal Bockstein factor is a bijective coordinate on all thirteen target
  sections.  Its six owner-conjugate rows recover the same coordinate in
  every labelled owner colour.  The inverse sheet map and target successor
  are explicit degree-eleven polynomials; the successor is one 13-cycle.
  A rank-eight affine system, certified by determinant 4, proves that no
  F_13-linear functional can make target translation affine.  Two CRT
  components separately distinguish all sheets, while the third has exactly
  the collisions {1,9} and {4,12}.  This scalarizes the coefficient atlas but
  does not identify it with a physical root deck, semantic endpoint, or
  common-ancestry transition.
source: root-holotopy-2026-07-28-bockstein-sheet-coordinate
depends_on:
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
related:
  - THM-2593-charged-target-section-atlas-and-minimal-c91-holonomy-trivialization
  - THM-2590-boolean-bockstein-and-theta-selector-incidence-spectrum
script: 04-computation/lrc14_linear_bockstein_sheet_coordinate_thm2601.py
output: 05-knowledge/results/lrc14_linear_bockstein_sheet_coordinate_thm2601.out
script_sha256: 99a993d03b649398c85c5698c1cd021322520ef0ceaa16485e4e5d5419a1c115
output_sha256: 988e100c989a334f0534bc5c9059de400d0270905b40ed87d080ca3356c51f0f
hash_basis: LF-normalized bytes
---

# THM-2601 -- a scalar Bockstein coordinate for the target-sheet atlas

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

THM-2585 constructs thirteen pairwise distinct unit factors

```text
Y_q in R_7=F_13[z]/(Phi_7(z)),                 q in F_13,       (1)
```

whose owner-coloured first Bocksteins are

```text
beta(D^(kappa,q))=Omega sigma_kappa(Y_q),
sigma_kappa:z |-> z^kappa.                                (2)
```

Pairwise distinctness says that the full six-coordinate factor remembers
the target section.  The stronger exact fact is that one scalar linear
measurement already remembers it.  The price is equally exact: target
translation becomes a nonlinear 13-cycle, and no linear measurement can
make that action affine.

This theorem lives on THM-2585's one globally primitive coefficient carrier.
It does not identify the target section with THM-2542's physical root deck or
construct an ordered transition kernel.

## 1. The six-dimensional factor table

Reduce a seven-row `(y_0,...,y_6)` modulo

```text
Phi_7(z)=1+z+...+z^6
```

in the power basis `1,z,...,z^5` by

```text
bar y=(y_0-y_6,...,y_5-y_6) in F_13^6.                    (3)
```

Let `bar Y_q` be the reductions of the thirteen exact rows printed in
THM-2585.  Define the linear functional

```text
lambda(x_0,...,x_5)
 =x_0+11x_1+8x_2+x_3+x_5.                                (4)
```

Direct substitution gives, in the order `q=0,...,12`,

```text
t_q=lambda(bar Y_q)
   =(9,3,10,1,8,0,7,2,12,11,4,6,5).                      (5)
```

This is a permutation of `F_13`.  Therefore

```text
q |-> t_q                                                     (6)
```

is a scalar coordinate on the entire target-section atlas.  No quotient,
division, or slice-by-slice reprimitivization is used: (4) is applied after
THM-2585's one global primitive reduction.

The coordinate is deliberately not called canonical.  The row (4) is one
explicit separating functional.  What is intrinsic is the verified
existence of a one-dimensional linear quotient that is lossless on these
thirteen distinguished factors.

## 2. Exact inverse and nonlinear monodromy

Since the values (5) exhaust `F_13`, the sheet label is recovered by the
unique polynomial function `P:F_13->F_13` of degree at most twelve.  In
low-to-high coefficient order,

```text
P(X)=5+4X+3X^2+7X^4+8X^5+12X^6
       +6X^7+11X^8+6X^9+6X^11,                           (7)
```

and exact evaluation gives

```text
P(t_q)=q                     for every q in F_13.          (8)
```

Transport the target successor `q->q+1` through (6).  Its scalar form is

```text
S(t_q)=t_(q+1),                                             (9)
```

with value table, indexed by `t=0,...,12`,

```text
(7,8,12,10,6,9,5,2,0,3,1,4,11).                          (10)
```

The unique reduced polynomial is

```text
S(X)=7+4X+X^3+6X^4+4X^5+6X^6
       +X^7+7X^8+6X^9+3X^10+2X^11.                       (11)
```

It has degree eleven.  Exhaustion of all `13^2` affine maps shows that none
equals (10).  Moreover (10) is one cycle of length thirteen, so

```text
S^13=identity,
S^n has no fixed point for 1<=n<13.                       (12)
```

Thus the scalar quotient loses no torsor information, but it linearizes
neither the torsor action nor the holonomy.  For a general target step `a`,
the action is `S^a`.  Equations (7)--(12) give a faithful nonlinear
permutation model of the complete `C_13` target atlas.

## 3. Every owner colour has the same scalar sheet coordinate

The owner action in (2) changes the power-basis coordinates.  It does not
destroy the scalarization.  Put

```text
lambda_kappa=lambda o sigma_(kappa^(-1)).                  (13)
```

The six rows of these conjugate functionals are

```text
kappa=1: (1,11,8,1,0,1),
kappa=2: (1,0,11,1,8,4),
kappa=3: (1,1,1,11,4,0),
kappa=4: (1,8,0,4,11,1),
kappa=5: (1,1,4,8,1,11),
kappa=6: (1,4,1,0,1,8).                                  (14)
```

By construction, and verified on all `78` pairs,

```text
lambda_kappa(sigma_kappa(bar Y_q))=t_q.                   (15)
```

THM-2590's socle calculation makes the Bockstein typing precise.  With
`u=zeta_13-1`,

```text
Omega=-u^11,
beta(D^(kappa,q))=-sigma_kappa(Y_q) tensor u^11.           (16)
```

Multiplication by `Omega` embeds `R_7` into that socle summand.  Hence, on
the image of the first Bockstein, extracting the `u^11` coefficient and
applying (14) recovers `t_q`, then (7) recovers `q`.  This is a linear
coefficient operation on the Bockstein image; it is not division by `Omega`
inside the nonreduced cyclotomic ring.

## 4. No linear coordinate can make target translation affine

The nonlinearity in (11) is not a poor choice of the separating row (4).
Suppose a linear functional `mu:F_13^6->F_13` and scalars `alpha,beta`
satisfied

```text
mu(bar Y_q)=alpha q+beta                 for all q.        (17)
```

The thirteen equations (17) have coefficient rows

```text
(bar Y_q,-q,-1) in F_13^8.                               (18)
```

Their rank is eight.  More concretely, the square minor using rows
`q=0,...,7` has determinant

```text
4 mod 13.                                                 (19)
```

Therefore (17) has only

```text
mu=0, alpha=0, beta=0.                                    (20)
```

In particular there is no nonconstant affine-covariant linear observable.
This also rules out an affine scalar successor for *any* linear separating
coordinate: an affine permutation of `F_13` having order thirteen must be a
nonzero translation, and its orbit would force (17) with `alpha!=0`.

The distinction is the exact holotopy boundary:

```text
linear scalar observation remembers every vertex;
nonlinear degree-11 permutation remembers the target action;
no linear gauge remembers both simultaneously.                         (21)
```

## 5. CRT boundary: one component can forget a sheet

The factorization

```text
Phi_7=(z^2+3z+1)(z^2+5z+1)(z^2+6z+1)                    (22)
```

gives three quadratic CRT components.  Reducing the thirteen `Y_q` in each
component gives the exact fibre census:

```text
z^2+3z+1: thirteen singleton fibres;
z^2+5z+1: thirteen singleton fibres;
z^2+6z+1: double fibres {1,9}, {4,12}, all others singleton. (23)
```

Thus two field components separately recover the target section, while the
third does not.  The full etale factor and the conjugate rows (14) repair
that ambiguity.  Equation (23) is the sharp hostile to the claim that an
arbitrary single CRT/owner component must be a complete sheet coordinate.

It is also why mere nonvanishing of one owner evaluation is weaker than
sheet recovery.  THM-2585 gives all six nonzero owner Bocksteins; (15)
uses the correctly conjugated functional on the full factor rather than
silently identifying one quadratic component with the whole atlas.

## 6. Map, loss, and next operation

The exact connection is

```text
source:    THM-2585's globally primitive factor Y_q;
map:       owner-normalize by sigma_(kappa^-1), then lambda;
target:    scalar t_q in F_13, with q=P(t_q);
preserved: target-section label, all six owner labels, C_13 orbit length;
destroyed/not supplied:
           positivity, physical x-support, semantic root/endpoint,
           relation charge, and an ordered common-ancestry edge;
action:    target translation becomes the nonlinear permutation S.       (24)
```

Consequently this theorem improves THM-2593's coefficient atlas: the full
six-coordinate unit factor can be compressed to one scalar without losing
the target sheet.  It does **not** remove THM-2593's chosen affine
root/target-torsor identification.  Indeed (19) proves that the available
linear coefficient observation cannot itself supply such an affine
identification.

The next physical object remains an ordered edge kernel

```text
M_k(q,q')                                                  (25)
```

on one common ancestry carrier.  Vertex factors determine `q` through this
theorem, but they do not determine which `q'` follows it.  A transition
matrix, not another vertex scalar, is needed to test the twisted seven-step
return from THM-2591/2593.

No row is excluded and LRC(14) remains open.

## 7. Exact verification

The dependency-free companion reconstructs every `Y_q` from THM-2585's
printed seven-row table and performs all arithmetic over `F_13`.  It checks:

- the separator (4), all thirteen values (5), and bijectivity;
- the inverse and successor polynomials at all field elements;
- all thirteen length-thirteen successor orbits and all `169` affine maps;
- all `78` owner-conjugate identities (15);
- rank eight and the determinant certificate (19); and
- every fibre in all three CRT factors (23).

Run

```text
python 04-computation/lrc14_linear_bockstein_sheet_coordinate_thm2601.py
python -O 04-computation/lrc14_linear_bockstein_sheet_coordinate_thm2601.py
```

Both modes must byte-match

```text
05-knowledge/results/lrc14_linear_bockstein_sheet_coordinate_thm2601.out
```

after LF normalization.  Every check raises explicitly under optimized
Python; no assertion disappears under `-O`.

QED pending independent hostile audit.
