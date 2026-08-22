---
id: THM-3269
title: "Scale-invariant clutch strength and canonical weighted bispanning polarization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The exact reset-link trap intervals canonically weight every edge of
  THM-3260's bispanning core by a positive row-scale- and orientation-invariant
  overlap ratio.  Among all 9,920 ordered complementary-tree charts, the
  product weight has one exact minimizer.  Hence the analytic clutch data
  selects a rigid ordered integral polarization with unique center 17.
  Combining its unique incident primitive direction with THM-3273's full
  cyclic Jacobian orders all twelve distinct Abel--Jacobi vertex classes and
  canonically removes both of THM-3260's abstract sidecars.  The resulting
  rank-compressed C12 label is not a group or physical phase intertwiner.
  The primitive direction separately identifies THM-3268's norm-phase
  augmentation with THM-3273's repaired rank-eleven edge sampler.  In fact,
  the genuine normalized J12 phases and the clutch order canonically select
  six undirected edges whose eleven directed phase increments form an
  integral unimodular sampler of the augmentation lattice.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-03
audit: >
  The assertion-independent exact companion pins promoted THM-3254,
  THM-3260, THM-3268 and THM-3273.  It recovers
  THM-3254's directly reconstructed response bank but
  trusts none of its chosen covers; rebuilds every endpoint trap and overlap
  strength; verifies orientation invariance; independently enumerates all
  complementary spanning trees and both orderings; proves the minimum unique
  by exact rational comparison; compares the intrinsic C2 image; and rebuilds
  the selected chart's two unimodular incidence matrices and GL11 transition.
  It also exhausts the selected tree's degree-compatible automorphisms and
  checks its unique center, radius, diameter and leaf set; then reconstructs
  the full critical coordinate, unique incident primitive generator, twelve
  distinct normalized exponents, circular order and rank-compressed C12
  bijection.  Finally it checks the generator-normalized J12 chart and the
  repaired 48-by-11 sampler's full rank; derives the six-edge selector from
  genuine J12 sign-orbits and exact clutch minima; and checks that its eleven
  directed evaluations have integral determinant minus one.
  Normal, optimized and stored replay plus LF hashes agree. An independent
  hostile audit rebuilt the response bank, all 74,748 spanning trees, all
  9,920 ordered bispanning charts, the unique weighted minimizer, incidence
  transition and rigid rooted geometry without importing this theorem. A
  second delta audit of the completed sidecar claims independently rebuilt
  the full critical Smith/adjugate coordinate, unique incident generator,
  twelve distinct exponents, rank compression, normalized edge sampler and
  augmentation embedding. A third delta audit independently rederived the
  five in-tree phase-orbit clutch minima, root-selected delayed edge, all
  eleven directed residues and integral determinant minus one. All audits
  passed with no repair.
depends_on:
  - THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go
  - THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary
  - THM-3268-nonzero-translation-norm-phase-walk-closed-form-and-rank-eleven-mixing-mode
  - THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary
related:
  - THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
script: 04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py
output: 05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out
script_sha256: 41ae9aeb01fea1384f59f3a2687b1a0482954bf202e5be2b6fc928ef579b116a
output_sha256: 65e6e1bb04f3d42d64c2a8e5322c0f3b37c9d05ed6b053670a2cd46293742e3c
hash_basis: LF-normalized bytes
---

# THM-3269 -- scale-invariant clutch strength and canonical weighted bispanning polarization

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3260 proves that the 22-edge reset-link core has 4,960 complementary
spanning-tree decompositions, connected by symmetric exchanges.  Any ordered
decomposition gives the integral polarization needed by its conditional
bridge to the eleven missing `C_12` modes, but the bare graph distinguishes no
chart.

The analytic clutch intervals do distinguish one.  The selector is intrinsic
under the harmless choices that a response atlas actually permits: positive
rescaling of either row and reversal of an undirected edge.

## 1. A scale-free weight on a response edge

For a first-link-blocked pair `e={i,j}`, use the orientation `(i,j)` and let

```text
I_sigma(i,j) subset [0,infinity]                       (1)
```

be THM-3254's closed interval of positive ratios `lambda` for which the reset
link state `sigma` traps `lambda f_i+f_j`.  Let

```text
A_e={U : I_sigma(i,j)=[0,U] for some sigma},
B_e={L : I_tau(i,j)=[L,infinity] for some tau}.        (2)
```

THM-3254 proves that every one of the 23 reset-link edges has a sharp
two-interval cover of the positive line and no one-state cover.  Thus the two
finite sets in `(2)` are nonempty and some `U>=L>0`.  Define the **clutch
strength**

```text
kappa_e=max {U/L : U in A_e, L in B_e}.                (3)
```

This is a positive rational number at least one.

The definition is independent of row normalization.  Replacing
`f_i,f_j` by `alpha f_i,beta f_j`, with `alpha,beta>0`, multiplies every ratio
endpoint by the common factor `beta/alpha`; every quotient `U/L` is unchanged.

It is also independent of edge orientation.  For the reversed pair,

```text
mu f_j+f_i = mu F_(1/mu).
```

Hence `[0,U]` becomes `[1/U,infinity]` and `[L,infinity]` becomes
`[0,1/L]`.  The corresponding overlap is

```text
(1/L)/(1/U)=U/L,                                      (4)
```

and inversion bijects all endpoint pairs.  Therefore `kappa_e` is an
invariant of the undirected, positively projectivized response edge.

On THM-3260's 22-edge leafless core, exact reconstruction gives

```text
12654135326210/11849988963879
  <= kappa_e <=
66479967673052544/306432210196733.                    (5)
```

The complete strength bank, including its maximizing state pairs, has digest

```text
c7b8441f54bd893f33f1767db913da0d40f392f20d74f6d177782a0490dc911a. (6)
```

## 2. The unique weighted chart

For an ordered complementary-tree chart `(T,T')`, put

```text
W(T)=product_(e in T) kappa_e.                         (7)
```

This avoids logarithms: comparing `W` is an exact rational version of
minimizing the sum of the logarithmic clutch strengths.  Exhaustion of all
4,960 unordered decompositions and both orderings gives a unique minimizer,

```text
T_*={(2,13),(2,17),(2,19),(2,21),(3,11),(3,16),
     (7,22),(10,11),(16,17),(18,19),(21,22)}.          (8)
```

Its canonical complement is

```text
T_*'={(2,10),(3,22),(7,18),(10,16),(10,22),(11,13),
      (13,18),(13,22),(16,21),(17,22),(19,22)}.        (9)
```

Both are spanning trees.  The ratio between the second-smallest ordered
weight and the minimum is exactly

```text
111771229679330735594490655
-------------------------------------------  > 1.     (10)
106891942612503954628585203
```

Thus uniqueness is an exact separation, not a floating-point tie break.
Since the product over all 22 edges is fixed, `(8)` simultaneously orders the
decomposition: `T_*` is the weak-clutch tree and `T_*'` the strong-clutch
cotree.

The nontrivial automorphism of the unweighted core swaps rows 17 and 21.  Its
image of `T_*` is not the minimizer; its exact weight ratio is

```text
W(C_2 T_*)/W(T_*)
 =719636186849558063/684759823142197440 > 1,           (11)
```

and it is the fourth ordered chart in strict weight order.  Hence the analytic
weights break the whole intrinsic `C_2`; the weighted response graph has no
residual chart ambiguity.

The selected weak tree is itself rigid:

```text
Aut(T_*)=1,     center(T_*)={17},
radius(T_*)=4,  diameter(T_*)=8,
leaves(T_*)={7,10,13,18}.                              (11a)
```

Thus the analytic selector canonically supplies row 17 as a root as well as
the ordered polarization.  Rigidity does **not** by itself provide a cyclic
ordering: imposing a generic graph-canonicalization algorithm and numbering
its output would be a convention, not an intertwiner with the physical phase
module.

## 3. Canonical integral polarization

Root at the canonical center row 17 and orient every edge by increasing
label.  The two reduced
incidence matrices for `(8)--(9)` satisfy

```text
det B_(T_*)=1,       det B_(T_*')=-1.                 (12)
```

Consequently

```text
U_*=B_(T_*)^(-1) B_(T_*') in GL_11(Z),
det U_*=-1.                                            (13)
```

The exact matrix digest is

```text
93f9ae1643ad9606f15826d8f1d3809fd0b04ce1d74c11d92d354edbdd737d4f. (14)
```

As in THM-3260, the columns `(-U_*z,z)` give an integral cycle frame and the
two incidence matrices give dual cut coordinates.  Equations `(3)` and `(8)`
make this polarization canonical from the normalized response problem itself;
no external tree-pair choice remains.

## 4. The full Jacobian orders the vertices

The proved critical-group portion of THM-3273 gives

```text
Jac(G_0)=C_74748.                                      (15)
```

Use the canonical tree center 17 as Abel--Jacobi root.  Its two neighbors in
`T_*` are 2 and 16.  In THM-3273's exact primitive cyclic coordinate their
rooted classes satisfy

```text
[2-17]=53470,       gcd(53470,74748)=2,
[16-17]=55507,      gcd(55507,74748)=1.               (16)
```

Thus `16` is the unique neighbor whose difference from the canonical root
generates the full Jacobian.  Put

```text
g=[16-17],       g |--> 1 in Z/74748,
g^(-1)=4891 in the displayed cyclic coordinate.       (17)
```

This normalization is independent of THM-3273's initially chosen cyclic
coordinate: multiplying every class and `g` by the same unit cancels in the
quotient.  The twelve normalized Abel--Jacobi exponents are distinct.  In
increasing cyclic order they are

```text
vertex:   17 16  18   13   11    22    19    10    7     21    2     3
exponent:  0  1 289 2344 9088 21481 25012 39646 48259 49832 53266 60022.
                                                               (18)
```

Let `ell(v)` be the rank `0,...,11` of `v` in `(18)`.  Equivalently,

```text
ell={17:0,16:1,18:2,13:3,11:4,22:5,
     19:6,10:7,7:8,21:9,2:10,3:11}.                  (19)
```

Every ingredient--the weighted tree, its center, its unique incident
primitive direction, Abel--Jacobi differences and the generator order--is
functorial under isomorphism of the weighted response graph.  Hence `(19)` is
a canonical bijection

```text
V(G_0) --> C_12,       17 |--> 0.                     (20)
```

This is rank compression of a twelve-point subset of `C_74748`, not a group
homomorphism.  Indeed every homomorphism `C_74748 -> C_12` factors through the
order-twelve quotient of THM-3273, whose vertex landing has only six classes;
`(19)` has twelve.  Therefore `(20)` is an abstract combinatorial phase label,
not the missing physical owner-phase intertwiner.

## 5. The abstract `C_12` bridge is now canonical

THM-3260's integral bridge

```text
Aug(Z[C_12]) --> H_1(G_0;Z)                           (21)
```

was conditional on a cyclic vertex bijection and an ordered tree-pair chart.
Equations `(8)--(9)` provide the chart, while `(20)` provides the bijection and
root.  Thus both formal sidecars are now internal to the exact response graph.
The load-bearing cyclic-difference map `Z[C_12]/Z N -> Aug(Z[C_12])` from
THM-3260 remains unchanged.

This is a genuine positive holotopy move: the unweighted exchange atlas has
no preferred origin, while the exact analytic clutch supplies a discrete
potential with one global minimum, and the full Jacobian separates the
otherwise colliding order-twelve phases.  The result is a canonical abstract
integral bridge, not a physical one.

## 6. Norm-phase mixing embeds in the repaired edge sampler

THM-3268 identifies the nonconstant norm-phase sector of the freely varying
translation quotient with

```text
A_norm=Aug(Q[C_12]),       Q|_(A_norm)=-I.             (22)
```

The primitive direction `g` in `(17)` also fixes the generator of
THM-3273's genuine critical quotient `J_12`.  It therefore gives a canonical
group isomorphism `C_12 -> J_12`, distinct from the nonlinear vertex-rank map
`(19)`.  THM-3273's directed core-edge sampler has one-dimensional kernel,
and the two delayed edges `(11,17),(11,21)` fill it.  Consequently the
composition

```text
A_norm  --> Aug(Q[J_12]) --> Q^(directed core + two delayed)       (23)
```

is a canonical injective `48 by 11` coefficient map.  Its image carries the
transported scalar action `-I` from `(22)`.

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

There is a smaller integral form of this sampler.  For an edge `e={u,v}`, let

```text
d(e)=min((j(v)-j(u)) mod 12,(j(u)-j(v)) mod 12),       (24)
```

where `j` is the genuine generator-normalized `J_12` phase, not the nonlinear
rank label `ell` from `(19)`.  The sign-orbit types occurring in `T_*` are
exactly `1,2,3,5,6`.  Within each type choose the unique edge of least exact
clutch strength.  Type 4 is absent from `T_*`; of THM-3273's two delayed
repair edges, exactly `(11,17)` meets the canonical root and so supplies that
type.  The resulting selection is

```text
d:       1        2        3         4         5         6
edge: (16,17)  (2,21)  (18,19)  (11,17)  (21,22)   (7,22). (25)
```

Use both orientations of each selected edge, except that the self-opposite
type 6 uses only increasing-label orientation.  Ordered by the resulting
nonzero phase increment, the eleven directed edges are

```text
1:(17,16), 2:(21,2), 3:(18,19), 4:(17,11), 5:(21,22),
6:(7,22),  7:(22,21), 8:(11,17), 9:(19,18), 10:(2,21),
11:(16,17).                                             (26)
```

They evaluate an augmentation coefficient at each nonzero residue exactly
once.  Relative to the standard integral basis of `Aug(Z[C_12])`, the
resulting `11 by 11` evaluation matrix has

```text
det=-1.                                                  (27)
```

Thus the canonical sampler is already an integral unimodular six-edge
subsampler; the 48-coordinate rational rank calculation is not needed for
injectivity once the delayed type-4 edge is admitted.

This is an abstract representation embedding only.  The edge sampler
evaluates coefficient functions at phase increments; it is not a transition
operator on the response graph, does not realize THM-3268's freely varying
walk, and does not make the rank label `(19)` into a physical phase observable.
In particular, `(25)--(27)` are coefficient evaluations, not six lawful
physical transitions or an owner-labelled current.

## 7. Failure boundary and scope

The selector uses only THM-3254's eleven reset-link states and the 22
first-link-blocked core edges.  It is special to the complete support-(1,3),
bank-I2 response normalization up to positive projective row scaling.

The argument would fail if the minimum in `(7)` were tied.  Merely knowing
the graph, its Betti number, its critical group, or the unordered interval
cover does not imply uniqueness; the exact rational endpoint values are
load-bearing.

The theorem does not turn the rank label `(19)` into an owner/phase observable,
does not identify it with THM-3255's physical phase marker, and does not
provide a `C_7` carrier on the delayed relative layer, a positive physical
current, or a scalar Gaussian-moment response.  It proves no row exclusion,
no arbitrary-radial NC2 theorem, no Gaussian Moment Conjecture consequence,
and no `LRC(14)` decrement.

## 8. Exact verification

Run

```text
python 04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py
python -O 04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py
```

and compare LF-normalized bytes with

```text
05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out.
```

The companion uses exact integer and rational arithmetic only.  It has no
floating point, randomness, discovery cache, graph-library dependency, or
optimization-sensitive assertion.  It independently reconstructs all
complementary tree pairs and the selected integral transition.

QED.
