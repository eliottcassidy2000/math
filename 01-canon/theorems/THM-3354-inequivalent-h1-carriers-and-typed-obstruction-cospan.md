---
id: THM-3354
title: "Inequivalent H1 carriers and the typed-obstruction cospan"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT. The proved LRC
  seven-chart class is an F13-valued graph-Cech class; its C91 mapping torus
  is an orbit carrier, not a mixed Z/91 coefficient class. Berggren-tree
  graph H1 and every pulled Gaussian ancestry charge are exact. All
  coefficient homomorphisms between the odd LRC groups and the sporadic S3,
  quartic V4, or characteristic-zero Hamiltonian-response carriers are
  trivial in both directions. The lawful D5 object is therefore a
  type-retaining comparison cospan, not a cohomology-class transport. The
  proof and exact finite controls below remain outside the proof graph pending
  independent hostile audit.
source: codex-2026-08-12-d5-typed-obstruction-cospan
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2703-c3-boundary-tree-arm-determinant-standard-plane-gate
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2708-c3-hermitian-gain-holonomy-discriminant-gate
  - THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z
  - THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion
  - THM-3353-split-prime-parabolic-branch-transplant-and-unary-transducer-compiler
  - HYP-9031-d5-h1-dictionary-lrc-word-current-vs-jc-flux
script: 04-computation/d5_typed_h1_no_go_thm3354.py
output: 05-knowledge/results/d5_typed_h1_no_go_thm3354.out
script_sha256: 6bc5389a51bea7eac74f586c6af3c471a591e8f18ca0f930c33bb3382c93d71c
output_sha256: 907d91f08e005e8e63411a3a7929c98ce46a780372a89155f64ad90a5e1f6fa8
hash_basis: LF-normalized bytes
---

# THM-3354 -- inequivalent H1 carriers and the typed-obstruction cospan

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**

This theorem audits the proposed D5 dictionary at the level of its actual
sites, coefficient objects, classes, target predicates, and missing sidecars.
The useful parallel survives, but a direct LRC-to-JC cohomology map does not.

No literature-priority claim is made. All statements are deductions from the
listed canon plus elementary graph cohomology and finite-group arithmetic.

## 1. The five carriers are not one H1

The current repo contains five relevant but inequivalent objects.

1. **LRC chart holonomy.** THM-2542 constructs

   ~~~text
   [g] in H^1(C_7^graph;F_13),             [g]=7a != 0.       (1)
   ~~~

   This is simplicial/Cech cohomology of a seven-cycle nerve.

2. **Gaussian multiplication charges.** THM-3336 constructs

   ~~~text
   tau in H^1(G;F_2),       lambda_p in H^1(G;Z),             (2)
   ~~~

   on the primitive Gaussian associate multiplication group G.

3. **Sporadic cubic monodromy.** On the finite-etale locus of the fixed
   three-dimensional sporadic Keller map, the Galois-closure torsor of
   THM-2473 has class

   ~~~text
   [rho_F] in H^1_et(U;S_3),             Mon(F)=S_3.          (3)
   ~~~

4. **Quartic Kummer data.** THM-2655/2685 give the standard plane

   ~~~text
   W=Hom(V_4,C_2) -> H^1_et(R_reg;mu_2),                     (4)
   ~~~

   together with its divisor-parity boundary and unit/class-group completion.

5. **Planar Hamiltonian response.** For R=K[x,z],

   ~~~text
   C_P=R/D_P(R),                    D_P(q)=Jac(P,q),          (5)
   ~~~

   is a K[P]-module. In the linear-in-z scope of THM-3348,

   ~~~text
   C_P tensor_(K[P]) K(P)
      ~= H^1_dR(Spec K(P)[x,g^(-1)]/K(P)).                  (6)
   ~~~

Equations (1)--(6) use the glyph H1, but their sites, variance, coefficients,
algebraic structures, and target predicates differ. They cannot be identified
by notation.

## 2. The Berggren ancestry carrier is exact

Let T be any tree, orient its edges arbitrarily, and let A be an abelian
coefficient group. Fix a root o. Given an A-valued one-cochain c, define
f(o)=0 and integrate c along the unique path from o to each vertex. Then on
every oriented edge

~~~text
c(x,y)=f(y)-f(x).                                          (7)
~~~

Thus

~~~text
H^1_graph(T;A)=0                                           (8)
~~~

for every constant abelian coefficient group A.

For the rooted Berggren tree, write w_x for the root word of vertex x. The
unique reduced path is

~~~text
P(x,y)=w_x^(-1)w_y.                                       (9)
~~~

If a character chi is evaluated on endpoint labels, put
f_chi(x)=chi(w_x), or f_chi(x)=chi(z_x) for Gaussian labels z_x. Its path
pullback is

~~~text
c_chi(x,y)=f_chi(y)-f_chi(x),                             (10)
~~~

so it is a coboundary. In particular, the genuine multiplication-group
classes tau,lambda_p from (2) do not become Berggren-tree classes after
endpoint pullback. THM-3345's source-dependent XOR path functor and
THM-3353's fixed-prime unary compilers therefore add exact ancestry arrows,
not ancestry holonomy.

## 3. Ninety-one is the carrier length, not the coefficient order

For a in F_13^*, THM-2542 has the constant overlap cochain

~~~text
g_k=a,                    k in F_7,                         (11)
~~~

and the sum is

~~~text
sum_k g_k=7a != 0 in F_13.                                (12)
~~~

Its minimal trivializing base cover has degree thirteen. The corresponding
skew successor is

~~~text
(k,r) |-> (k+1,r+a)             on F_7 x F_13.             (13)
~~~

Because (1,a) has order lcm(7,13)=91, (13) is one cycle of length ninety-one.
This proves a C91 mapping-torus **carrier**. It does not change the coefficient
group in (1).

More sharply, every coefficient homomorphism

~~~text
iota:C_13 -> C_91                                         (14)
~~~

has image in the unique order-thirteen subgroup 7C_91. Under

~~~text
C_91 ~= C_7 x C_13,                                       (15)
~~~

the C7 projection of every element in that image is zero. A nonzero iota
transports (1) only to an order-thirteen class. It cannot create a primitive
mixed order-ninety-one class or jointly nonzero 7-by-13 character. Conversely,
projection C91 -> C13 can recover the thirteen-primary coordinate while
discarding the seven-coordinate; it cannot prove that the latter existed.

Hence the statement

~~~text
rho_LRC in H^1(chart;Z/91) with both CRT components nonzero              (16)
~~~

is not a consequence of THM-2542. The proved statement is exactly (1), plus
the carrier (13).

## 4. The coefficient no-go is exact

Write C_n for the cyclic group of order n. The following coefficient
homomorphism sets contain only the trivial map:

~~~text
Hom(C_13,S_3),     Hom(C_91,S_3),
Hom(S_3,C_13),     Hom(S_3,C_91),                          (17)

Hom(C_13,V_4),     Hom(C_91,V_4),
Hom(V_4,C_13),     Hom(V_4,C_91).                          (18)
~~~

For maps from a cyclic group, the image order divides both the odd number 13
or 91 and the target order 6 or 4; it is therefore one. Every map from S3 to
an abelian group factors through S3^ab=C2, and every map from V4 has
two-primary image. An odd cyclic target has no such nonzero subgroup. This
proves (17)--(18).

There is an equally sharp characteristic-zero boundary. Since D_P(R) is a
K-linear subspace, C_P in (5) is a K-vector space. Its additive group is
torsion-free and divisible. Therefore

~~~text
Hom(C_13,C_P)=Hom(C_91,C_P)=0,                            (19)
Hom(C_P,C_13)=Hom(C_P,C_91)=0.                            (20)
~~~

Indeed, (19) sends a finite-order generator to torsion, while the image in
(20) would be a divisible subgroup of a finite group.

Cohomology transport also needs a map of bases/sites with the proper
variance. No canonical map is supplied between the seven-cycle nerve, the
Gaussian multiplication group, the sporadic etale locus, the quartic
resolvent normalization, and the generic punctured line. Even if one
fabricated a base map, every coefficient-only route from the LRC class to
(3)--(6), or back, is zero by (17)--(20).

This is a no-go for the proposed **direct coefficient-induced D5 map**. It
does not deny that a future construction could introduce a new correspondence
space and a new coefficient-changing sidecar; such data would be the substance
of a new theorem, not a consequence of the present carriers.

## 5. The Hamiltonian class has a different target predicate

The relevant planar mate observer is

~~~text
theta=[1] in C_P.                                         (21)
~~~

A polynomial Q with Jac(P,Q)=1 exists exactly when

~~~text
theta=0.                                                  (22)
~~~

In THM-3348's generic identification, theta maps to [-dx/g]. The full generic
de Rham group can be nonzero while this one observer vanishes; for a one-root
coefficient, the missing obstruction can instead survive as vertical
K[P]-torsion. Conversely, a Keller pair already has D_P(Q)=1, so (21)
vanishes tautologically even though other cover or monodromy data may be
nontrivial.

Thus neither nonzero generic de Rham H1 nor nonzero etale monodromy is
synonymous with failure of the planar mate equation. The observer, integral
lattice, and localization loss must be retained.

## 6. The same tree can carry opposite Kummer answers

Graph H1 also cannot replace THM-2655's weighted boundary lattice. Consider
the same three-arm star with one fixed centre and three leaves, cycled by C3.

For the negative D4 Cartan weights, its matrix is

~~~text
M_2=[-2  1  1  1;
      1 -2  0  0;
      1  0 -2  0;
      1  0  0 -2].                                       (23)
~~~

THM-2703 gives

~~~text
|det M_2|=4,          A_(M_2)[2] ~= W.                    (24)
~~~

Keep the same unweighted star and centre weight, but change every leaf weight
to -3:

~~~text
M_3=[-2  1  1  1;
      1 -3  0  0;
      1  0 -3  0;
      1  0  0 -3].                                       (25)
~~~

Then

~~~text
|det M_3|=27,         A_(M_3)[2]=0.                       (26)
~~~

Both underlying graphs are trees and have zero graph H1. One weighted lattice
has exactly the quartic standard plane and the other has none. Therefore the
intersection matrix, divisor parity, saturation, and gain sidecars are
load-bearing. An unweighted ancestry or boundary graph cannot stand in for
the Kummer carrier.

## 7. The corrected D5 map is a typed comparison cospan

Define TypedObstruction to be the set of records

~~~text
(site,
 coefficient object,
 class or distinguished observer,
 reversible/directed variance,
 proved status of the target implication,
 missing realization or boundary sidecar,
 quantifier tag).                                         (27)
~~~

There are explicit encoding maps

~~~text
LRC_record  --Phi_L-->  TypedObstruction  <--Phi_J--  JC_records,            (28)
~~~

where JC_records is the disjoint union of the sporadic monodromy, quartic
Kummer, and Hamiltonian-response records in (3)--(6). Their images are:

~~~text
Phi_L:
  (C_7^graph,F_13,[7a],horizontal reversible,
   semantic arrival OPEN,
   positive vertical 2-cell or twisted lift missing,
   local proved / global LRC universal target);

Phi_J,sporadic:
  (U_et,S_3,[rho_F],reversible etale,
   fixed three-dimensional noninjectivity REALIZED,
   nonproper escape boundary retained,
   existential fixed-map counterexample);

Phi_J,quartic:
  (R_reg,mu_2,W,reversible etale,
   A4/S4 Keller exclusion OPEN,
   weighted saturated boundary and Kummer realization missing,
   universal quartic target);

Phi_J,response:
  (Spec K(P)[x,g^(-1)]_dR,K(P),theta,
   additive/localized,
   mate iff theta=0,
   integral torsion/observer sidecar retained,
   fixed P before universal JC(2)).                        (29)
~~~

The maps in (28) are type-retaining encodings, **not** maps on cohomology
classes. A coarser projection of (27) preserves only this grammar:

~~~text
a quotient/localization carries an obstruction;
the target needs a named realization sidecar;
discarding that sidecar invalidates the target implication.                (30)
~~~

It destroys the site, coefficient arithmetic, group law, class order,
orientation variance, and quantifier. Consequently equality after the coarse
projection in (30) supports analogy and experiment design, but no theorem
transfer.

## 8. Exact sidecars and cheapest lawful next tests

The cospan identifies different next tests rather than one shared map.

1. **LRC flat fill.** A root-trivializing semantic correction must provide

   ~~~text
   c_k in F_13,              sum_k c_k=-7a,                (31)
   ~~~

   together with THM-2542's target-role intertwiner killing neutral roles and
   detecting a positive target-active edge.

2. **LRC twisted lift.** Instead of cancelling (1), construct a positive
   semantic vertical path on the nontrivial carrier (13), preserving the old
   common ancestry/deep/source face.

3. **Quartic JC.** Construct or exclude the Q-equivariant Kummer standard
   plane using the weighted saturated boundary lattice and its divisor
   residues. The unweighted graph is insufficient by (23)--(26).

4. **Planar response.** Test the distinguished integral class theta, not
   arbitrary generic de Rham rank, and retain the vertical torsion discarded
   by localization.

5. **Gaussian ancestry.** Any nonzero current must add a quotient, boundary,
   or physical incidence sidecar not present in the tree, because (7)--(10)
   make every endpoint charge exact.

These are parallel realization problems in inequivalent cohomology theories.
Solving one does not solve another.

## 9. Exact finite referee and stopping boundary

The companion independently enumerates every finite coefficient homomorphism
in (17)--(18), all thirteen maps C13 -> C91 and their twelve embeddings, all
twelve nonzero mapping-torus steps, the graph Betti controls, and the two
weighted-star determinants and mod-two nullities. Normal and optimized
transcripts must agree with the stored output:

~~~bash
python3 04-computation/d5_typed_h1_no_go_thm3354.py
python3 -O 04-computation/d5_typed_h1_no_go_thm3354.py
~~~

The universal conclusions rest on (7), group-order/abelianization arguments,
and the characteristic-zero vector-space argument, not on bounded inference.

No LRC current, semantic 2-cell, scalar-row exclusion, quartic Kummer
realization, planar mate, JC(2), DC(2), or cross-domain reduction follows.
The direct-map form of HYP-9031 must be read as superseded if this candidate
passes audit; its useful survivor is precisely the typed grammar (27)--(30).
