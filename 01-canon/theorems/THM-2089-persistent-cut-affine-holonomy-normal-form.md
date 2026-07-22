---
id: THM-2089
title: "Persistent cut templates are flat affine holonomy systems"
status: >
  PROVED. Every cut edge transports one endpoint speed affinely in the other
  and the guard. A chord is persistent on the THM-2088 tree template exactly
  when its cycle has product-one multiplicative holonomy and zero signed
  offset holonomy. Hence every rank-six persistent cut admits the flat gauge
  q_i=u_i(z+v_i h), with explicit rational height bounds and a finite-index
  integrality lattice. A nonflat cycle is the rank-seven finite branch. This
  classifies the template geometry but does not discharge all flat templates.
source: codex-2026-07-22-LRC14-affine-cut-holonomy
depends_on:
  - THM-2088
related:
  - THM-1226
  - THM-1405
  - THM-2062
  - THM-2065
  - HYP-3832
script: 04-computation/lrc14_persistent_cut_affine_holonomy_referee_codex_20260722.py
output: 05-knowledge/results/lrc14_persistent_cut_affine_holonomy_referee_codex_20260722.out
script_sha256: dbc509b66e1d2e249bea94ebb01dcc5be38b3892d8f80bbb5d819b100b12758e
output_sha256: 90faa42a8ec6db5bfd5e0f0b1bb7c37f8663dc9a595d5a3c2b32169651311b05
hash_basis: normalized repository blobs (LF)
---

# THM-2089 -- persistent cuts are flat affine systems

Work in THM-2088's no-pair branch. Every chosen cut relation has the form

```text
a_e h+b_e q_i+c_e q_j=0,
b_e c_e!=0,
max(|a_e|,|b_e|,|c_e|)<=57.                            (1)
```

When the edge is traversed from `i` to `j`, put

```text
alpha_e=-b_e/c_e,             beta_e=-a_e/c_e.         (2)
```

Then (1) is the affine transport

```text
q_j=alpha_e q_i+beta_e h.                              (3)
```

For the reverse traversal use

```text
alpha_(e^-1)=1/alpha_e,
beta_(e^-1)=-beta_e/alpha_e,                           (4)
```

so the two transports are literal inverses.

## 1. Exact signed cycle holonomy

Consider an oriented cycle of length `r`, labelled so that

```text
a_i h+b_i s_i+c_i s_(i+1)=0,       s_(r+1)=s_1.       (5)
```

Define

```text
A_C=product_i(-b_i/c_i),
B_C=sum_(i=1)^r (-a_i/c_i)
                 product_(j=i+1)^r(-b_j/c_j).          (6)
```

Iterating (3) once around the cycle gives

```text
s_1=A_C s_1+B_C h.                                    (7)
```

For an integer form, put

```text
D_C=product_i c_i,
N_C=(-1)^r product_i b_i,
R_C=sum_(i=1)^r (-a_i)
       product_(j=i+1)^r(-b_j)
       product_(j=1)^(i-1)c_j.                        (8)
```

Then `A_C=N_C/D_C`, `B_C=R_C/D_C`, and (7) becomes

```text
(D_C-N_C)s_1=R_C h.                                   (9)
```

This is the signed-affine extension of THM-1226's multiplicative-holonomy
identity. No positivity of the edge coefficients is required.

## 2. Persistence is exactly two zero holonomies

Fix a cut spanning tree `T`. Its six rows are independent by THM-2088, so
their kernel `V_T` has dimension two. Choose a speed root `o`. Repeated tree
transport expresses every speed as

```text
q_i=A_i z+B_i h,                 z=q_o,                (10)
```

where `A_i!=0`. In particular `h,z` are free coordinates on `V_T`; the tree
imposes no relation between them.

Let a chord close a fundamental cycle `C`. Its row is persistent on `V_T`
if and only if the cycle transport is the identity for every `h,z`. By (7),
and because `h` and the base speed are independent on `V_T`, this is
equivalent to

```text
A_C=1 and B_C=0,
```

or, after clearing denominators,

```text
D_C=N_C and R_C=0.                                    (11)
```

Therefore THM-2088's rank-six branch is equivalent to (11) on every
fundamental cycle. One failure of (11) raises the cut matrix to rank seven
and enters THM-2088's finite maximal-minor branch.

Every simple cycle of `K_(A,B)` has length at most six. Thus

```text
|D_C|,|N_C|<=57^r,
|R_C|<=r*57^r,
r<=6.                                                  (12)
```

If a cycle is nonflat, (9) gives a nonzero guard/base-speed relation of
coefficient height at most

```text
6*57^6=205778683494.                                  (13)
```

This larger relation-height invoice is not needed for THM-2088's sharper
absolute finite bound; it is useful as a direct falsifier of persistence.

## 3. Flat gauge normal form

Assume (11) on all fundamental cycles. Then every cycle is flat, so the tree
formula (10) is independent of the chosen root-to-vertex path. Put

```text
u_i=A_i,                   v_i=B_i/A_i.                (14)
```

The whole terminal core has the normal form

```text
q_i=u_i(z+v_i h).                                      (15)
```

Equivalently, along an oriented edge,

```text
u_j=alpha_e u_i,
v_j-v_i=beta_e/u_j.                                   (16)
```

Thus `alpha` is an exact multiplicative cocycle and the normalized offset
`beta_e/u_j` is an exact additive cocycle. Conversely any rational data
`u_i!=0,v_i` satisfying (16) produce (3) and force every cycle holonomy to be
`(1,0)`. This proves that persistent cut templates are exactly flat affine
line bundles over the cut graph. QED.

This is the cut/cycle split of THM-1405 in a different coefficient group:
the spanning tree fixes a gauge, while chords measure holonomy in the affine
group `Q^x semidirect-product Q`.

## 4. Explicit height and integrality lattice

A simple tree path has length at most six. Clearing its transport gives

```text
q_i=(N_i z+R_i h)/D_i,                                 (17)
```

with

```text
0<|N_i|,|D_i|<=57^6,
|R_i|<=6*57^6.                                        (18)
```

After reduction the gauge in (15) may therefore be chosen with numerator and
denominator bounds

```text
height(u_i)<=57^6,
height(v_i)<=6*57^6.                                  (19)
```

The integer points of the terminal template are parametrized by the explicit
congruence lattice

```text
Lambda={ (h,z) in Z^2:
          D_i divides N_i z+R_i h for every i}.        (20)
```

It has finite index satisfying

```text
[Z^2:Lambda]<=product_i |D_i|<=57^36,                 (21)
```

because the root coordinate has `D_o=1` and there are only six other speed
coordinates. Saturating `Lambda` gives the integral coefficient-row template
to which the THM-2062 prime wheel can be attached.

After quotienting common positive scale, put `t=z/h`. Positivity is the
intersection of the seven rational half-line conditions

```text
u_i(t+v_i)>0,                                         (22)
```

For each pair, distinctness either deletes the rational wall

```text
u_i(t+v_i)=u_j(t+v_j).                                (23)
```

If two affine forms in (23) are identical, the template has no admissible
distinct point and is discarded. Otherwise there are at most 21 such walls.
Thus the real projective parameter space of an admissible persistent cut
template is an explicit interval, punctured by finitely many collision walls.
The unresolved difficulty is arithmetic and phase incidence on the lattice
points of that interval, not continuous template geometry.

## 5. Frontier effect and honest scope

THM-2088's opaque word "persistent" is now an exact pair of cycle equations.
For cut types `2+5` and `3+4`, only four or six fundamental-cycle pairs
`(D_C-N_C,R_C)` need be tested. Failure is finite; simultaneous vanishing
gives the gauge (15), bounds (18)--(21), and the one-parameter projective
interval (22).

This still concerns only the terminal seven-core and last guard. The earlier
three dyadic guards and two original seam tails remain outside (15). Nor do
flatness, positivity, or integrality imply a lonely time. The next target is
to splice those five outer coordinates into `Lambda`, then combine the
THM-2062/2069 deletion wheel, THM-2082 translated grids, THM-2079 owner
addresses, and THM-2081 restricted-edge margin on the resulting interval.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that a persistent circuit is an indivisible
exception. On a cut graph it is a flatness statement with two independently
checkable defects: multiplicative curvature `D_C-N_C` and affine-offset
curvature `R_C`.

Candidate vertices were speeds, cut edges, fundamental cycles, Fourier modes,
and proof obligations. Fundamental cycles are the smallest lossless falsifier
set. The binary relation "which cycle has the larger normalized curvature"
can orient a scheduling tournament, but its score histogram and SCCs discard
the ordered pair of signed rational holonomies. When every curvature is zero,
all tournament edges are ties and any chord order is a Hamiltonian path. The
faithful carrier is the affine cochain with its tree gauge and two-component
cycle holonomy. The complete-bipartite link analogy in HYP-3832 is structural
only; no coboundary-expansion theorem is invoked. QED.
