---
id: THM-3354
title: "Inequivalent H1 carriers, the direct-map no-go, and a typed comparison table"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The proved LRC
  seven-chart class is an F13-valued graph-Cech class; its C91 mapping torus
  is an orbit carrier, not a mixed Z/91 coefficient class. Berggren-tree
  graph H1 vanishes and every endpoint-pulled Gaussian ancestry charge is
  exact. All direct coefficient homomorphisms between the odd LRC groups and
  the sporadic S3, quartic V4/mu2, or characteristic-zero integral
  Hamiltonian-response carriers are trivial in both directions. The generic
  Hamiltonian localization loses vertical torsion, so generic vanishing is
  only necessary, not sufficient, for a polynomial mate. The typed comparison
  cospan below is DEFINITIONAL bookkeeping, not theorem content, and no new
  correspondence space or coefficient-changing sidecar is excluded.
source: codex-2026-08-12-d5-typed-obstruction-cospan
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2703-c3-boundary-tree-arm-determinant-standard-plane-gate
  - THM-2708-c3-hermitian-gain-holonomy-discriminant-gate
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
  - THM-3353-split-prime-parabolic-branch-transplant-and-unary-transducer-compiler
related:
  - THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z
  - THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion
  - HYP-9031-d5-h1-dictionary-lrc-word-current-vs-jc-flux
script: 04-computation/d5_typed_h1_no_go_thm3354.py
output: 05-knowledge/results/d5_typed_h1_no_go_thm3354.out
script_sha256: b0b86426cedd58ca3e23ff46e9c8af7fe8504d0dcffb0a16c16dc9cf3fa02a48
output_sha256: a6d340d8870ecccc0cdca5cf20dbfdbfcea997be7b85243ca7b5dce761035cb9
hash_basis: LF-normalized bytes
---

# THM-3354 -- inequivalent H1 carriers, the direct-map no-go, and a typed comparison table

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

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

3. **Sporadic cubic monodromy.** In the target coordinates of THM-2473, put

   ~~~text
   L=27a^2c^2-18abc+16a+b^3c-b^2,
   U=A^3_C \ Z(L).
   ~~~

   The Galois closure of the degree-three finite-etale cover over `U` is a
   finite-etale `S3`-torsor. It is classified by the nonabelian pointed set

   ~~~text
   [rho_F] in H^1_et(U;underline S_3),
   geometric monodromy S_3.                                  (3)
   ~~~

4. **Quartic Kummer data.** Here the coefficient sheaf is `mu_2`, while
   the parameter plane is

   ~~~text
   W=Hom(V_4,C_2) ~= F_2^2.                                  (4)
   ~~~

   THM-2655/2685 place the three nonzero characters in
   `H^1_et(R_reg;mu_2)`, subject to the Kummer unit and class-group rows and
   the divisor-parity completion. The additive group of `W` happens to be
   isomorphic to `V4`; neither is the coefficient sheaf.

5. **Planar Hamiltonian response.** For `R=K[x,z]`, the integral response
   record is

   ~~~text
   C_P=R/D_P(R),                    D_P(q)=Jac(P,q),          (5)
   theta_int=[1] in C_P.
   ~~~

   This is a `K[P]`-module. In the linear-in-`z` scope of THM-3348, its
   generic localization is a different record:

   ~~~text
   C_P tensor_(K[P]) K(P)
      ~= H^1_dR(Spec K(P)[x,g^(-1)]/K(P)),
   theta_gen=image(theta_int).                              (6)
   ~~~

These records are sometimes compressed under one `H1` or flux vocabulary,
but their sites, variance, coefficients, algebraic structures, and target
predicates differ. They cannot be identified by notation.

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

For the rooted Berggren tree, write `w_x` for the root word of vertex `x`. The
unique reduced path is

~~~text
P(x,y)=w_x^(-1)w_y.                                       (9)
~~~

If `chi` is a character of the free Berggren word group, put
`f_chi(x)=chi(w_x)`. Separately, after choosing endpoint Gaussian labels
`z_x`, evaluate the genuine multiplication-group characters as
`f_tau(x)=tau(z_x)` and `f_p(x)=lambda_p(z_x)`. Every such endpoint function
has path pullback

~~~text
c_chi(x,y)=f_chi(y)-f_chi(x),                             (10)
~~~

so it is a coboundary. In particular, the genuine multiplication-group
classes `tau,lambda_p` from (2) do not become Berggren-tree classes after
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

Although `C91` is not the coefficient object in (1), suppose hypothetically
that its abstract carrier-cycle group is treated as coefficients. Every
coefficient homomorphism

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

Hom(C_13,C_2),     Hom(C_91,C_2),
Hom(C_2,C_13),     Hom(C_2,C_91).                          (18a)
~~~

For maps from a cyclic group, the image order divides both the odd number 13
or 91 and the target order 6 or 4; it is therefore one. Every map from S3 to
an abelian group factors through S3^ab=C2, and every map from V4 has
two-primary image. An odd cyclic target has no such nonzero subgroup. The
same order argument proves (18a); over `C`, the abstract geometric fibre of
`mu_2` is `C2`. This proves (17)--(18a). It also keeps separate the quartic
parameter plane `W`, its abstract additive group `V4`, and its coefficient
sheaf `mu_2`.

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

## 5. Integral and generic Hamiltonian observers have different predicates

The relevant integral planar mate observer is

~~~text
theta_int=[1] in C_P.                                     (21)
~~~

A polynomial Q with Jac(P,Q)=1 exists exactly when

~~~text
theta_int=0.                                              (22)
~~~

Under THM-3348's generic identification, let

~~~text
theta_gen=image(theta_int)=[-dx/g].                        (22a)
~~~

A polynomial mate implies `theta_gen=0`, but the converse is false because
localization discards vertical `K[P]`-torsion. The minimal linear-in-`z`
hostile is

~~~text
P=x+x^2 z,                 D_P=(1+2xz)partial_z-x^2 partial_x.
~~~

In the localized ring, `Q=1/x` satisfies `D_P(Q)=1`, so `theta_gen=0`.
Integrally, THM-3348 gives

~~~text
Ann_(K[P])(theta_int)=(P),
~~~

so `theta_int!=0` and no polynomial mate exists. Thus the full generic de
Rham group can be nonzero while this observer vanishes; the missing
obstruction survives as vertical torsion. Conversely, a Keller pair already
has `D_P(Q)=1`, so (21) vanishes tautologically even though other cover or
monodromy data may be nontrivial.

Thus neither nonzero generic de Rham H1, generic observer vanishing, nor
nonzero etale monodromy is equivalent to the integral planar mate predicate.
The observer, integral module, and localization loss must be retained.

## 6. The same tree can carry opposite Kummer answers

Graph H1 also cannot replace THM-2655's weighted boundary lattice. Consider
the same three-arm star with one fixed centre and three leaves, cycled by C3.
The lattice computations below are unconditional; their geometric Kummer
application is only under the full-rank, independent rational-surface
completion hypotheses of THM-2703.

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
intersection matrix, divisor parity, and saturation sidecars are load-bearing
in the stated completion scope. For quotient cycles, THM-2708 adds gain
holonomy as a further independent coordinate. An unweighted ancestry or
boundary graph cannot stand in for the Kummer carrier.

## 7. Definition: the typed comparison table

**DEFINITIONAL ONLY.** This section supplies bookkeeping, not an additional
theorem. It asserts no naturality, universal property, or cross-domain
morphism.

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

There are tautological set-theoretic encoding maps

~~~text
LRC_record  --Phi_L-->  TypedObstruction  <--Phi_J--  JC_records,            (28)
~~~

where `JC_records` is the disjoint union of the sporadic monodromy, quartic
Kummer, integral-response, and generic-response records in (3)--(6). Their
images are:

~~~text
Phi_L:
  (C_7^graph,F_13,[g]=7a,horizontal reversible,
   semantic arrival OPEN,
   positive vertical 2-cell or twisted lift missing,
   local proved / global LRC universal target);

Phi_J,sporadic:
  (U,underline S_3,[rho_F],reversible finite etale,
   fixed three-dimensional noninjectivity REALIZED,
   nonproper escape boundary retained,
   existential fixed-map counterexample);

Phi_J,quartic:
  (R_reg,mu_2,image of W,reversible etale,
   A4/S4 Keller exclusion OPEN,
   unit row plus Cl(R_reg)[2], weighted saturated boundary,
   and Kummer realization missing,
   universal quartic target);

Phi_J,response-integral:
  (K[P] subset R,C_P,theta_int,
   additive/integral,
   mate iff theta_int=0,
   polynomial realization and vertical torsion retained,
   fixed P before universal JC(2));

Phi_J,response-generic:
  (Spec K(P)[x,g^(-1)]_dR,K(P)-de-Rham module,theta_gen,
   additive/localized,
   mate implies theta_gen=0 only,
   integral module and vertical torsion destroyed,
   fixed linear-in-z P before universal JC(2)).            (29)
~~~

The maps in (28) merely package already-defined records; their existence is
tautological. They are **not** maps on cohomology classes. A coarser
projection of (27) preserves only this grammar:

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

   together with an intertwiner of the type required by THM-2542: it must kill
   neutral roles and detect a positive target-active edge. No such
   intertwiner is constructed there.

2. **LRC twisted lift.** Instead of cancelling (1), construct a positive
   semantic vertical path on the nontrivial carrier (13), preserving the old
   common ancestry/deep/source face.

3. **Quartic JC.** Construct or exclude the `Q`-equivariant Kummer standard
   plane in the full universal Kummer row: units together with
   `Cl(R_reg)[2]`. When the completion hypotheses of THM-2703/2708 hold, use
   the weighted saturated boundary lattice, divisor residues, and any cyclic
   gain coordinate. The unweighted graph is insufficient by (23)--(26).

4. **Planar response.** Test the distinguished integral class `theta_int`,
   not arbitrary generic de Rham rank or `theta_gen` alone, and retain the
   vertical torsion discarded by localization.

5. **Gaussian ancestry.** Any non-exact ancestry class or holonomy must add a
   quotient, boundary, or physical incidence sidecar not present in the tree,
   because (7)--(10) make every endpoint charge exact. Nonzero exact
   one-cochains may of course exist on a tree.

These are parallel realization problems in inequivalent cohomology theories.
Solving one does not solve another.

## 9. Exact finite referee and stopping boundary

The companion independently enumerates every finite coefficient homomorphism
in (17)--(18a), all thirteen maps `C13 -> C91` and their twelve embeddings,
all twelve nonzero mapping-torus steps, the graph Betti controls, and the two
weighted-star determinants and mod-two nullities. It also checks the exact
localized hostile `D_(x+x^2z)(1/x)=1` and the one-root annihilator exponent.
Normal and optimized transcripts agree with the stored output:

~~~bash
python3 04-computation/d5_typed_h1_no_go_thm3354.py
python3 -O 04-computation/d5_typed_h1_no_go_thm3354.py
~~~

The universal conclusions rest on (7), group-order/abelianization arguments,
the characteristic-zero vector-space argument, and the integral definition
of `C_P`, not on bounded inference.

No LRC current, semantic 2-cell, scalar-row exclusion, quartic Kummer
realization, planar mate, JC(2), DC(2), or cross-domain reduction follows.
The direct-map form of HYP-9031 is superseded. Its useful survivor is the
typed grammar (27)--(30), which remains definitional and transfers no theorem.
