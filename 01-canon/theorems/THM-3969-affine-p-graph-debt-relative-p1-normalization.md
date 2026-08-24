---
id: THM-3969
title: "Collision-free affine-P graph debt has a finite-class relative-P1 normalization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For the first
  graph debt beyond quadratic P-depth,
  q=3rP-r^3+(P-r^2)^2(c+a(P-r^2)), the generic normalization is the graph
  of x=(3r+3v+cv^2)/(1-av^3) on a relative projective line. Its only
  basepoint polynomial is Xi=c^3+27a^2r^3-27acr+27a. If Xi is a nonzero
  scalar, the graph is already a finite smooth normalization: the complement
  of the cubic pole multisection in A1 times P1. Its class group is finite
  cyclic, so THM-3922 forbids a same-field Keller affine plane. This closes
  the collision-free family only; every zero of Xi creates an exceptional
  normalization curve and remains to be typed.
source: jc-cohn3709 / post-THM-3967 first cubic-P conductor seam, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_degree6_one_place, 2026-08-24). The
  audit independently checked relative-Proj affineness, smoothness, the
  fixed-degree resultant and separate W=0 basepoint row, degree three on
  every projective-line fibre, proper plus quasi-finite finiteness, the
  function-field recovery v=y/x and full-normalization bridge, and the
  Weil-localization computation Cl(Y)=Z/gcd(d_j). It also checked the
  nonreduced a=0 support and both moving controls. Normal and optimized
  32-gate runs byte-match the frozen output; hashes and documentation checks
  pass.
depends_on:
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
  - THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt
related:
  - THM-3964-polynomial-graph-hidden-double-root-normalization
  - THM-3967-quadratic-p-depth-natural-cubic-conductor-debt-closure
script: 04-computation/jc2_affine_p_graph_relative_p1_collision_thm3969.py
output: 05-knowledge/results/jc2_affine_p_graph_relative_p1_collision_thm3969.out
script_sha256: 947b7e94bb83accbfed9b7612dcfb16b4f2be112b58e386e94a1d6bc27457fb1
output_sha256: e1d137c24e6330b888cda5fbfaf968a6da404f0db2e1c664049d742338008a26
semantic_sha256: fa36713e47425b67ca0a12a38e82e3ab2c66c82911206acf322e750cac6c6e73
hash_basis: raw LF bytes
---

# THM-3969 -- three Kummer-conjugate punctures are a finite-class trap

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of
characteristic zero. Let `a,c,r in k[t]` and put

```text
y=P-r^2,
q=3rP-r^3+y^2(c+ay),
F=T^3-3PT-q,
A=k[P,t,T]/(F).                                         (1)
```

Define the collision polynomial

```text
Xi=c^3+27a^2r^3-27acr+27a.                              (2)
```

Assume `Xi in k^*`. Then `A` is a domain and its finite normalization is a
smooth affine surface `Y` whose class group is finite cyclic. Consequently
the same cubic function field admits no planar Keller chart with target
coordinates `(P,t)`. This proves the collision-free row only.

## 1. The relative-projective-line normalization

Set

```text
x=T+r.                                                   (3)
```

The equation becomes

```text
F=x^3-3rx^2-3xy-cy^2-ay^3.                              (4)
```

On a line `y=vx` through the repeated graph, the residual equation is

```text
D(v)=1-av^3,
N(v)=3r+3v+cv^2,
x=N/D,                 y=vN/D.                          (5)
```

The denominator cannot be discarded: its three geometric roots are the
three points of the normalized generic fibre above target infinity, while
`v=infinity` is generally a finite point with `(x,y)=(0,-c/a)`.

Homogenize on

```text
S=A1_t times P1_[V:W]
```

by putting

```text
D_h=W^3-aV^3,
N_h=3rW^2+3VW+cV^2.                                    (6)
```

Let

```text
E=V(D_h),                       Y=S minus E.             (7)
```

This is the standard relative-Proj open `D_+(D_h)`, so it is affine over
`A1_t`; it is smooth because it is an open of the smooth surface `S`. The
degree-zero functions

```text
x=W N_h/D_h,
y=V N_h/D_h,
P=r^2+y,                       T=x-r                    (8)
```

are regular on `Y` and satisfy `(4)`.

## 2. The resultant is exactly the finite-graph gate

Direct elimination gives

```text
Res_v(1-av^3,3r+3v+cv^2)=Xi.                            (9)
```

The homogeneous pair defining the target projective coordinate is

```text
[r^2D_h+V N_h : D_h].                                  (10)
```

A basepoint of `(10)` is exactly a common zero of `D_h,N_h`. The affine
resultant `(9)` also sees the point `W=0`: when `a=0`, a common point there
occurs exactly when `c=0`, and then `Xi=0`. Thus `Xi in k^*` makes `(10)`
basepoint-free on every fibre.

It follows that `(10)` defines a proper morphism

```text
S -> P1_P times A1_t                                   (11)
```

whose restriction to every projective-line fibre has degree three. It is
therefore quasi-finite and proper, hence finite. Its inverse image of the
affine target `A1_P times A1_t` is exactly `Y`, so base change makes

```text
pi:Y -> A2_(P,t)                                       (12)
```

finite of degree three.

On `W!=0`, equation `(5)` gives `v=y/x` in the function field. Conversely
`P,T` are the functions `(8)`. Hence

```text
k(Y)=k(P,T),                    [k(Y):k(P,t)]=3.         (13)
```

The monic polynomial `F` is therefore irreducible, and `A` is a domain.
The finite smooth surface `Y` is birational over `A`, so it is the full
finite normalization.

## 3. The pole multisection leaves only finite class group

Write the reduced divisor `E` as the union of its irreducible components

```text
E_1 union ... union E_s
```

and let `d_j>=1` be the degree of `E_j` over `A1_t`. The polynomial `D_h`
is primitive in `(V,W)`, so there is no vertical component. Since

```text
Cl(S)=Pic(S)=Z H,                                      (14)
```

where `H` is the class of a section of the projective-line bundle, one has

```text
[E_j]=d_j H.                                           (15)
```

The Weil localization sequence for `Y=S minus E` gives

```text
Cl(Y)=Z/gZ,                   g=gcd(d_1,...,d_s).        (16)
```

This includes the nonreduced endpoint `a=0`: if `c` is a unit, the support
of `D_h=W^3` is one degree-one section and `(16)` gives `Cl(Y)=0`.

If a same-field planar Keller source existed, normalization-form Zariski
Main would embed it as a dense `A2` open in `Y`. Because the finite degree is
three, this open would be proper. THM-3922 requires the class group of such
a finite normal completion to be free of positive rank. Equation `(16)` is
finite, a contradiction.

## 4. Why collisions are the exact residual

When `Xi` has a zero, `(10)` has a basepoint on the pole multisection. The
finite graph completion then adds an exceptional affine curve over that
fibre; its class can turn the finite cyclic group `(16)` into a free class.
Neither the model `Y` nor the collision-free class computation may be used
unchanged there.

The seam is nonempty. For example

```text
r=0, c=t, a=(1-t^3)/27
```

has `Xi=1`, so it is a genuinely moving collision-free row excluded by this
theorem. By contrast `a=t,c=1,r=0` has `Xi=1+27t` and one collision fibre;
it is deliberately outside the conclusion.

Thus the remaining affine-`P` graph problem is not an untyped bivariate
normalization. It is the explicit exceptional-curve packet over `V(Xi)`.
Subsequent THM-3972 constructs its blowup normalization whenever `Xi` is
nonconstant and squarefree, computes the irreducible-pole class/canonical/
ramification ledger, and closes every `a=t`, constant-`(c,r)` simple row.
General canonical-compatible squarefree rows, multiple `Xi` roots, reducible
pole multisections, higher `P`-depth, nongraph repeated factors, and `JC(2)`
remain open. **QED.**

## Reproduction

```bash
python3 04-computation/jc2_affine_p_graph_relative_p1_collision_thm3969.py
python3 -O 04-computation/jc2_affine_p_graph_relative_p1_collision_thm3969.py
sha256sum 04-computation/jc2_affine_p_graph_relative_p1_collision_thm3969.py \
  05-knowledge/results/jc2_affine_p_graph_relative_p1_collision_thm3969.out
python3 agents/check_docs.py
```
