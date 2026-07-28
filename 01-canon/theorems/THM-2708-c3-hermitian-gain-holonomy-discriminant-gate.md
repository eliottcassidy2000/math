---
id: THM-2708
title: "C3 Hermitian gain-holonomy discriminant gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every nonsingular integral symmetric matrix with a permutation
  C3 action, evaluate its free-orbit circulant blocks at omega in F4.  The
  resulting matrix B is Hermitian, independent of orbit representatives up
  to diagonal unitary gauge, and the standard-plane multiplicity in
  (Z^V/MZ^V)[2] is exactly nullity_F4(B).  Thus B is singular precisely when
  the two-torsion discriminant module contains the nontrivial irreducible
  F2[C3] plane.  THM-2703 is the balanced forest specialization.  On the
  minimal quotient triangle, balanced voltage lifts give 3K3, determinant
  8, and one standard plane; nontrivial holonomy gives C9, determinant 2,
  and none.  This removes the tree restriction from the full-rank rational-
  surface gate.  Independence of the boundary classes also forces every
  unit to be constant, so det(B)=1 directly excludes a quartic Kummer
  carrier in that scope.  Boundary relations, reflection completion,
  geometric realization, physical LRC homology, and all general
  A4/S4/JC2/DC2 conclusions remain open.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2700-danielewski-s3-resolvent-standard-plane-exclusion
  - THM-2703-c3-boundary-tree-arm-determinant-standard-plane-gate
script: 04-computation/c3_hermitian_gain_holonomy_discriminant_thm2708.py
output: 05-knowledge/results/c3_hermitian_gain_holonomy_discriminant_thm2708.out
script_sha256: f46c10dc19e207fe4ca95b9a541096799d031ae2109c84f31e86122aa3bdca66
output_sha256: 3df42b1eaf24f5505bc946261f771a0a0fcf8b35405fd5e0893c023dd044a356
hash_basis: LF-normalized bytes
---

# THM-2708 -- C3 Hermitian gain-holonomy discriminant gate

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2703 finds the quartic standard plane arm by arm when the boundary graph
is a tree.  A cycle prevents that direct sum, but it does not destroy the
calculation.  It introduces one extra coordinate: `C3` gain holonomy.  After
that coordinate is retained, the whole standard block is one Hermitian
matrix over `F4`.

This is the literal algebra behind the recurring two/three pairing:

```text
V4 characters form the additive plane F2^2;
C3 cycles its three nonzero vectors;
x^2+x+1 is irreducible over F2;
therefore the standard plane is the one-dimensional F4 module.       (1)
```

The theorem computes that module for an arbitrary invariant lattice matrix,
not only a tree or plumbing.

## 1. Free orbits and the `F4` block

Let `sigma` be a permutation of a finite set `V` with `sigma^3=1`, and let

```text
M:Z^V -> Z^V                                                (2)
```

be symmetric, integral, `sigma`-invariant, and nonsingular.  Let
`O_1,...,O_s` be the free three-element vertex orbits and choose a
representative `v_i` from each.  Fixed vertices are allowed.

Put

```text
F4=F2[omega]/(omega^2+omega+1),
bar(omega)=omega^2.                                        (3)
```

For free orbits define

```text
m_ij(a)=M_(v_i,sigma^a v_j) mod 2,
B_ij=sum_(a=0)^2 m_ij(a) omega^a in F4.                    (4)
```

Then

```text
B_ji=bar(B_ij),                                             (5)
```

so `B` is Hermitian.  Indeed, invariance and symmetry give

```text
m_ji(a)=m_ij(-a),                                          (6)
```

which is exactly `(5)`.  In particular,

```text
bar(det B)=det(bar B)=det(B^T)=det B,
det B belongs to F2.                                       (7)
```

Thus the determinant gate is a genuine zero/one invariant even though its
entries live in `F4`.

## 2. Exact standard-plane theorem

Let

```text
A_M=Z^V/MZ^V,
W=the nontrivial irreducible two-dimensional F2[C3]-module. (8)
```

Then

```text
(A_M[2])_std = underlying_F2 ker(B:F4^s -> F4^s),           (9)

multiplicity_W A_M[2]=nullity_F4(B).                       (10)
```

Consequently

```text
A_M[2] contains W  iff  det(B)=0.                          (11)
```

Here the empty Hermitian matrix has determinant one, so an action with no
free vertex orbit correctly contributes no standard plane.

### Proof

Modulo two, Maschke decomposition gives

```text
F2[C3]=F2 direct_sum W,       End_(F2[C3])(W)=F4.          (12)
```

Every fixed vertex contributes only to the trivial isotypic part.  Every
free vertex orbit contributes one copy of the permutation module and hence
one copy of `W` to the standard part.  Choose the `F4` coordinate on that
copy so that `sigma` acts by multiplication by `omega` (choosing
`omega^2` conjugates every displayed entry and changes no kernel).

The coupling from a fixed vertex to a free orbit is constant on its three
vertices.  Its standard restriction is therefore multiplied by

```text
1+omega+omega^2=0.                                        (13)
```

Hence all fixed/free blocks vanish on the standard isotypic part.  Between
two free orbits, the invariant `3 by 3` circulant block acts on `W` by its
evaluation `(4)`.  Therefore the restriction of `M mod 2` to the standard
part is precisely the `F4`-linear map `B`.

Finally, tensoring

```text
0 -> Z^V --M--> Z^V -> A_M -> 0                           (14)
```

with `F2` gives equivariantly

```text
A_M[2]=Tor_1^Z(A_M,F2)=ker(M mod 2).                       (15)
```

Taking standard parts proves `(9)`--`(11)`.

## 3. Switching and the exact gain coordinate

Replace the representative `v_i` by `sigma^(h_i)v_i`.  Equation `(4)`
changes to

```text
B'_ij=omega^(h_i-h_j) B_ij.                               (16)
```

If `D=diag(omega^(h_i))`, then

```text
B'=D B D^(-1).                                             (17)
```

This is a unitary diagonal gauge: `D^(-1)=bar(D)`.  Rank, nullity, determinant,
and `(9)` are independent of every representative choice.

When each nonzero off-diagonal block in `(4)` is a monomial `omega^h`, it is
the voltage or gain on the corresponding edge of the free-orbit quotient
graph.  Equation `(16)` is ordinary gain switching.  Products around cycles
are the surviving holonomies.

If the quotient component is a forest, choose a root and recursively choose
the `h_i` so every edge gain becomes one.  The Hermitian block then becomes
the ordinary mod-two weighted quotient matrix.  In the tree setting of
THM-2703, the free-orbit quotient splits into the representative moving
components and this reduction gives exactly

```text
B=direct_sum_C (M_C mod 2).                                (18)
```

Thus THM-2703 is the balanced-forest specialization of `(9)`, while cycle
holonomy is the precise coordinate its tree proof never needed.

## 4. Minimal sharp cycle: `3K3` versus `C9`

Take a quotient triangle with zero diagonal and one unit gain on each edge.
Orient its gains as `g_01,g_12,g_20 in {1,omega,omega^2}`.  Then

```text
B=[ 0          g_01       bar(g_20) ]
  [ bar(g_01)  0          g_12      ]
  [ g_20       bar(g_12)  0         ],                    (19)

det B=Tr_(F4/F2)(g_01 g_12 g_20).                         (20)
```

There are exactly nine balanced assignments, with gain product one, and
eighteen unbalanced assignments.  Since

```text
Tr(1)=0,             Tr(omega)=Tr(omega^2)=1,              (21)
```

the standard plane occurs exactly in the balanced cases.

The integral voltage lifts make the distinction visible:

```text
balanced holonomy:
  lift = 3 disjoint K3,
  abs(det M)=8,
  Smith(M)=(1^6,2^3),
  A_M[2]=F2^3=F2_triv direct_sum W;                         (22)

nontrivial holonomy:
  lift = C9,
  abs(det M)=2,
  Smith(M)=(1^8,2),
  A_M[2]=F2_triv.                                          (23)
```

Both lifts have the same unlabelled quotient triangle and the same local
degrees.  Forgetting the gain changes the answer.  This is the smallest
sharp obstruction to replacing `(4)` by the quotient graph alone.

## 5. Full-rank boundary-lattice corollary

Let `U` be a smooth complex surface with a `C3` action and a smooth
projective rational equivariant completion

```text
U subset Xbar,         D=Xbar\U.                           (24)
```

Assume the irreducible boundary component classes are independent and span
a finite-index `C3`-stable lattice `L` in `Pic(Xbar)`.  No tree hypothesis is
made.  Let `M` be their intersection matrix and form `B` by `(4)`.

The rational-surface Picard lattice `P=Pic(Xbar)` is unimodular.  Divisor
localization gives

```text
Pic(U)=P/L,
L subset P=P^* subset L^*,
Pic(U) injects C3-equivariantly into L^*/L=A_M.             (25)
```

The same localization row computes the units:

```text
O(U)^*/C^*
 =ker(Z^{components of D} -> Pic(Xbar)).                   (26)
```

The component classes are independent, so this kernel is zero.  Since every
nonzero complex constant is a square,

```text
O(U)^*/O(U)^(*)2=0.                                        (27)
```

Therefore

```text
Pic(U)[2] contains W  ==>  det(B)=0.                       (28)
```

For a quartic `A4/S4` Keller resolvent in the scope of THM-2655, its required
`V4` character plane lies either in unit squareclasses or in `Pic(U)[2]`.
Consequently

```text
det(B)=1
 ==> this completion cannot be the full quartic resolvent. (29)
```

This removes THM-2703's tree restriction and closes both Kummer alternatives
inside the stated completion box.  It does not remove the full-rank,
independence, or rational-completion hypotheses.  When boundary classes have
relations, `(26)` shows exactly how a unit-standard plane can return.

## 6. What the gain theorem does not identify

The word `holonomy` appears elsewhere in the repository, but the types must
remain separate.

THM-2658 uses **integer lift gains** of circular arcs to decide positive
common intersection.  Here `(4)` uses a **mod-three voltage** of boundary
component orbits to decide a two-torsion discriminant representation.  The
shared switching law is exact; positivity, mass, and physical components are
not transported.

Likewise THM-2672 shows that forgetting carry can turn thirteen disjoint
physical `Delta^11` strata into a virtual coarse `boundary Delta^12`.  The
pair `(22)`--`(23)` is the lattice analogue: forgetting voltage identifies
`3K3` with `C9`.  Neither result converts the other's virtual quotient into
physical homology.

The present theorem leaves open:

- unit-standard planes when boundary classes have relations;
- the reflection action needed for a prescribed full `S3` module;
- boundary component relations or non-full-rank lattices;
- geometric realization of a singular Hermitian block by a Keller
  resolvent;
- any physical LRC owner, phase, endpoint, or carry current;
- arbitrary `A4`, `S4`, planar Jacobian, or Dixmier candidates.

## 7. Exact companion

Run

```text
python 04-computation/c3_hermitian_gain_holonomy_discriminant_thm2708.py
python -O 04-computation/c3_hermitian_gain_holonomy_discriminant_thm2708.py
```

and compare with the stored output.  The exact script checks:

- every invariant symmetric binary circulant matrix on one, two, or three
  free orbits: `32,900` cases;
- `16,384` further cases with two fixed vertices and every fixed/free
  coupling, verifying that those couplings vanish on the standard block;
- all `10,588` labelled-tree gain assignments through five quotient
  vertices, explicitly switching every forest to unit gain;
- all `27` quotient-triangle gain assignments and all `729` representative
  gauges;
- the integral determinants, Smith forms, connected components, and standard
  kernels of the sharp `3K3/C9` pair.

Normal and optimized transcripts byte-match.  The script uses explicit
failure raises and exact `F2`, `F4`, integer determinant, and Smith arithmetic.

## 8. Connection contract and next test

```text
source object:
  an integral symmetric lattice with a permutation C3 action;

target object:
  a Hermitian matrix over F4 on the free vertex orbits;

map:
  evaluate each invariant circulant block at omega;

preserved predicate:
  exact multiplicity of the standard plane in discriminant two-torsion;

destroyed information:
  reflection, boundary-lattice saturation outside the independent case,
  geometric realization, and every physical LRC current;

needed sidecars:
  a full-rank independent equivariant completion for quartics;
  a physical component/owner/phase lift for LRC;

cheapest decisive next test:
  compute B for the boundary divisor or resolution graph of each surviving
  cyclic-cubic or S3 quartic-resolvent normal form, retaining all cycle
  voltages instead of only the unlabelled dual graph.
```

The determinant `det B in F2` is therefore the exact non-tree boundary gate,
not another scalar discriminant coincidence.
