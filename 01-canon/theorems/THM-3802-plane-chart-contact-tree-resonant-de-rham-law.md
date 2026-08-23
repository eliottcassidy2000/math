---
id: THM-3802
title: "Plane-chart contact-tree resonant de Rham law"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  Any smooth
  affine surface equipped with an actual finite plane-chart cover whose
  nontrivial intersections are the same punctured f-chart and whose leaf
  coordinates satisfy w=beta_i(f)+f^n a_i has H1_dR=0 and
  H2_dR=k^h/k(1,...,1).  The physical symplectic class is exactly the vector
  of f^(n-1) leaf coefficients modulo constants.  It is exact iff that vector
  is constant; otherwise no global polynomial Darboux pair exists.  The
  theorem assumes the affine atlas and makes no existence claim from a
  decorative contact tree alone.
source: jc_sparse_direct_search / THM-3791--3797 atlas abstraction, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The proof types the actual open cover,
  transition orientation, regular chart inverses, full Cech--de Rham rows,
  finite-etale descent variant, equality boundary, and Darboux consequence.
  The exact companion checks transitions and every nonresonant primitive for
  n=1..6, simplex incidence ranks through eight leaves, common-translation
  invariance, exact constant-resonance controls, and the THM-3791, THM-3789,
  and THM-3797 specializations.  Normal and optimized runs byte-match the
  frozen transcript.  Independent hostile audit remains due.
related:
  - THM-3791-moving-root-danielewski-resonant-jet-de-rham-law
  - THM-3797-confluent-quadratic-hermite-jet-completion-no-go
script: 04-computation/jc2_plane_chart_contact_tree_resonant_derham_thm3802.py
output: 05-knowledge/results/jc2_plane_chart_contact_tree_resonant_derham_thm3802.out
script_sha256: 7ada1a8e2ded5eec4925aa5abaca6bdb984d0ef037aea043c82f5523553700e0
output_sha256: 540660b79ca4e60501e6d3ce585f8d41113497497adcbbd782df0831a7aab6e2
semantic_sha256: 4a77818e9aea370b9661778aa9478919b40df16d976b75f5aa247e5b459911a7
hash_basis: raw LF bytes
---

# THM-3802 -- the resonant coefficient is an invariant of the actual leaf atlas

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  The
hypersurface in THM-3791 and the iterated affine modification in THM-3797 look
different globally, but their proofs use the same smaller object: a finite
cover by planes whose deleted-arm overlaps are one common punctured chart.
This theorem isolates exactly that object.  It deliberately does not infer a
surface from a list of formal roots.

Let `k` be a field of characteristic zero, let `n>=1`, and let `X` be a smooth,
affine, geometrically connected `k`-surface.  Suppose the following data are
given.

1. There are global functions `f,w in O(X)` and a finite Zariski open cover

   ```text
   X=U_1 union ... union U_h,                 h>=1.            (1)
   ```

2. Every chart comes with a specified isomorphism

   ```text
   U_i ~= A2_(f,a_i)                                           (2)
   ```

   carrying the global functions to

   ```text
   w=beta_i(f)+f^n a_i.                                       (3)
   ```

   Here `beta_i in k[f]` is the unique chosen representative of degree less
   than `n`.  For distinct leaves, `beta_i-beta_j` is nonzero modulo `f^n`.

3. For every set of at least two distinct indices, the actual scheme-theoretic
   intersection is

   ```text
   U_(i_0) intersection ... intersection U_(i_p)
    =D_X(f) ~= G_m,f x A1_w.                                  (4)
   ```

The isomorphisms `(2)`, the regular coordinates `a_i`, and the equalities
`(3),(4)` are hypotheses.  In particular, the theorem does not assert that
arbitrary jets `beta_i`, or an arbitrary labelled contact tree, glue to a
separated affine surface.

Define on `D_X(f)`

```text
omega=df wedge dw/f^n.                                        (5)
```

Then `(5)` extends uniquely to a regular symplectic form on `X`.  Put

```text
r_i=[f^(n-1)]beta_i(f),                 r=(r_1,...,r_h).       (6)
```

There are canonical identifications

```text
H^0_dR(X/k)=k,
H^1_dR(X/k)=0,
H^2_dR(X/k)=k^h/k(1,...,1),
H^q_dR(X/k)=0 for q>=3,                                      (7)
```

under which

```text
[omega]=[r].                                                   (8)
```

Consequently

```text
omega is algebraically exact
iff r_1=...=r_h.                                               (9)
```

Let `{ , }` be the Poisson bracket inverse to `omega`.  If `r` is
nonconstant, there are no `P,Q in O(X)` with `{P,Q}=1`.

## 1. The physical form and the exact transition sign

On `U_i`, equation `(3)` gives

```text
dw=beta_i'(f)df+n f^(n-1)a_i df+f^n da_i,
omega=df wedge da_i.                                         (10)
```

Thus `(5)` extends regularly, and `(10)` is nowhere zero.  The local
primitive compatible with the repository orientation is

```text
eta_i=-a_i df,                 d eta_i=omega.                 (11)
```

On an overlap, `(3)` gives the exact singular shear

```text
a_j-a_i=[beta_i(f)-beta_j(f)]f^(-n).                          (12)
```

Therefore

```text
eta_j-eta_i=[beta_j(f)-beta_i(f)]f^(-n)df.                   (13)
```

For a term `c_m f^m` in a root difference, the corresponding one-form in
`(13)` is `c_m f^(m-n)df`.  It has the rational primitive

```text
c_m f^(m-n+1)/(m-n+1)                                        (14)
```

unless `m=n-1`.  The exceptional term is `c_(n-1)df/f`.  Since

```text
H^1_dR(G_m x A1)=k[df/f],                                    (15)
```

the residue of `(13)` is exactly

```text
r_j-r_i.                                                       (16)
```

This also explains why only the truncated jets modulo `f^n` matter.  Adding
`f^n g_i(f)` to `beta_i` merely changes `a_i` by the regular function
`-g_i(f)` and does not affect `(16)`.

## 2. The complete Cech--de Rham calculation

Apply the algebraic Cech--de Rham spectral sequence to the actual affine cover
`(1)--(4)`.  The chart and overlap rows are

```text
H^q_dR(U_i)=k,0,0,...,
H^q_dR(D_X(f))=k,k[df/f],0,...                               (17)
```

for `q=0,1,2,...`.  The `q=0` row is the full simplex cochain complex on `h`
vertices.  It has cohomology `k` in degree zero and zero in every positive
degree.

The `q=1` row has no vertex term because `H^1_dR(U_i)=0`; it starts at the
edge space.  Its edge-to-triangle kernel is the image of the ordinary
vertex-difference map

```text
delta:k^h -> k^(edges),                 (s_i) |-> (s_j-s_i). (18)
```

Unlike the full simplex row, there is no vertex term in this row by which to
quotient that image.  Hence

```text
E_2^(1,1)=image(delta)=k^h/k(1,...,1).                        (19)
```

The higher part of the simplex complex is exact, and no higher differential
can enter or leave `(19)`.  Equations `(17)--(19)` prove all groups in `(7)`.
The local primitives `(11)` produce the edge cocycle `(r_j-r_i)df/f`, which
is the image under `(18)` of `r`.  This proves `(8)` with the asserted sign.

Because `X` is affine, vanishing of the class in `(8)` is equivalent to the
existence of a global algebraic one-form whose derivative is `omega`.  This
proves `(9)`, including sufficiency rather than only an obstruction.

## 3. Darboux consequence and coordinate invariance

For the Poisson bracket inverse to `omega`,

```text
dP wedge dQ={P,Q}omega                                       (20)
```

for all global functions `P,Q`.  If `{P,Q}=1`, then

```text
omega=dP wedge dQ=d(P dQ),                                   (21)
```

so `[omega]=0`.  Equations `(8),(9)` prove the Darboux no-go for nonconstant
`r`.

A common target translation

```text
w_new=w+u(f)                                                   (22)
```

replaces every `beta_i` by the same truncated polynomial `beta_i+u mod f^n`.
It adds the scalar `[f^(n-1)]u` to every entry of `r`, so its reduced class is
unchanged.  Relabelling the leaves only permutes the coordinates.  Thus `(8)`
is intrinsic to the oriented atlas, not to a chosen origin for `w`.

The lower coefficients and the pairwise contact orders

```text
ord_f(beta_i-beta_j)<n                                        (23)
```

record the affine-modification/contact tree and control how the leaves first
separate.  They affect the surface, but `(14)--(16)` show that the physical
de Rham class forgets all of them except the resonant coefficient at depth
`n-1`.

## 4. Finite-etale descent form

There is a basis-free version when the leaves are not individually defined
over `k`.  Let `E` be a finite etale `k`-algebra of rank `h`.  Suppose that
after a faithfully flat field extension `K/k` splitting `E`, the base change
`X_K` has an atlas satisfying `(1)--(4)`, indexed by the `K`-points of
`Spec E`, and that the atlas, `f,w,omega`, and leaf jets carry their canonical
descent data.  Suppose the split jets are the evaluations of

```text
beta_hat(f) in E[f]/(f^n),
rho=[f^(n-1)]beta_hat in E.                                  (24)
```

The split calculation is natural under both pullbacks to `K tensor_k K`.
Algebraic de Rham cohomology commutes with flat scalar extension, and

```text
(E/k) tensor_k K ~= K^h/K(1,...,1).                          (25)
```

Faithful descent therefore gives

```text
H^1_dR(X/k)=0,
H^2_dR(X/k)=E/k,
[omega]=rho mod k.                                            (26)
```

In this form, exactness is equivalent to `rho in k`.  The descent datum is an
explicit hypothesis here; a nonsplit formal leaf set by itself is not being
promoted to a geometric atlas.

## 5. Exact controls and hostile boundaries

### 5.1 THM-3791

For the smooth hypersurface

```text
c^n e=Sigma(c,b),                 Sigma(0,b) squarefree,      (27)
```

the Hensel-division charts of THM-3791 satisfy `(1)--(4)` after splitting the
finite-etale special-fibre algebra.  Their `beta_i` are the canonical Hensel
jets, so `(8)` is exactly THM-3791's resonant-root law.  Section 4 recovers its
root-free statement `H^2=E/k`.  For fixed roots, exponent one is obstructed,
whereas every exponent at least two has `r=0`.

### 5.2 THM-3789 and THM-3797

The higher-pole Hermite surface of THM-3789 has `n=3` and split jets

```text
(f^2,lambda_1,...,lambda_d),                 r=(1,0,...,0).   (28)
```

The confluent quadratic surface of THM-3797 is not a squarefree hypersurface.
Only after its two affine modifications does it have the actual three-chart
atlas

```text
beta=(f^2,0,ell),                 r=(1,0,0).                  (29)
```

Thus `(29)` is a genuine contact-tree control rather than another use of the
hypersurface theorem.

### 5.3 What is not covered

The singular one-equation intermediate in THM-3797 fails the smoothness and
plane-chart hypotheses.  Its normalization still has an `A_1` point.  Neither
object is eligible for this theorem.  More generally, the following data do
not suffice:

- a list of Laurent transitions without regular chart inverses;
- a rooted tree without a separated affine gluing;
- pairwise overlaps larger or smaller than the common `D(f)`; or
- a form that does not extend as `df wedge da_i` on every leaf chart.

These guards are load-bearing: changing the cover changes the Cech incidence
complex and can change both `H^2` and the equality criterion.  Subject to the
actual-atlas hypotheses, however, `(8)` is universal.  **QED, conditional only
on independent hostile audit.**
