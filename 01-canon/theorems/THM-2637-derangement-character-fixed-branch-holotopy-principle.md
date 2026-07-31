---
id: THM-2637
title: "Derangement-character fixed-branch holotopy principle"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  Let a finite
  group G act on a finite sheet set Omega and let chi:G->C_l be a
  prime-cyclic character whose nonzero support consists of derangements.
  Every holonomy element with a fixed sheet then lies in ker chi.  Hence a
  G-local system on a connected graph whose fundamental-cycle holonomies
  each have a fixed sheet has character-trivial holonomy.  If chi is
  injective on the holonomy image, the local system is trivial and has
  |Omega| parallel sections.  For the regular C13 carry torsor, a
  common-carry fixed branch on each fundamental-cycle holonomy forces zero
  holonomy and thirteen sections.  The common carry gauge and the fact that
  one is testing one selected holonomy component are load-bearing: separate endpoint
  marginals, full determinant-sector support, a private row at only one
  endpoint, or a diagonal contribution inside a mixture of translations do
  not imply trivial holonomy.  This abstracts the abelian
  derangement-character hinge inside the stronger D4 Keller exclusion and
  identifies the exact conditional consumer of an LRC endpoint pair twist;
  it proves no common-carry branch, row exclusion, or LRC(14) conclusion.
source: root-2026-07-28-derangement-holotopy-bridge
depends_on: []
related:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
  - THM-2635-half-tooth-opposite-graph-unit-section-and-reversed-digit-closure
---

# THM-2637 -- a charged derangement cannot carry a fixed branch

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**

The same elementary group-theoretic hinge appears in two apparently remote
frontiers.  In the Jacobian problem, a boundary character must occur on an
inertia element while openness supplies a finite fixed inverse branch.  In
the LRC carry problem, chronology asks whether a principal carry transition
has a branch with the same physical carry at both ends.  The first situation
is closed for `D_4`; the second is still open.  This theorem isolates the
common implication and, equally importantly, the exact extra data needed to
apply it.

## 1. The fixed-branch character lemma

Let a finite group `G` act on a finite set `Omega`, and let

```text
chi:G -> C_ell                                             (1)
```

be a homomorphism to the additive cyclic group of prime order `ell`.  Assume
the **derangement-support condition**

```text
chi(g)!=0  ==>  Fix_Omega(g)=empty.                        (2)
```

Then

```text
Fix_Omega(g)!=empty  ==>  chi(g)=0.                        (3)
```

This is the contrapositive of (2), but its typing matters: the fixed point
must be a branch of the *same action element* `g`, in one declared sheet
gauge.  A fixed point of another element, or two separately nonempty sheet
marginals, is not a witness for (3).

Condition (2) is conjugacy invariant.  Indeed `chi(aga^(-1))=chi(g)`, and
conjugation carries `Fix(g)` bijectively to `Fix(aga^(-1))`.  Thus (3) is
intrinsic for a holonomy element, whose representative is defined only up to
conjugacy.

## 2. Graph holotopy and a cycle-basis criterion

Let `Gamma` be a connected finite graph carrying a `G`-local system with
sheet fibre `Omega`.  Choose a base vertex.  Parallel transport gives a
holonomy representation

```text
rho:pi_1(Gamma) -> G.                                     (4)
```

Choose a spanning tree.  The fundamental loops associated to the non-tree
edges freely generate `pi_1(Gamma)`; write their holonomies as

```text
H_1,...,H_b.                                              (5)
```

If each fundamental holonomy has a fixed branch,

```text
Fix_Omega(H_i)!=empty                 for i=1,...,b,       (6)
```

then (3) gives `chi(H_i)=0` for every generator.  Hence

```text
chi o rho =0 on pi_1(Gamma).                              (7)
```

If `chi` is injective on `rho(pi_1(Gamma))`, equation (7) makes every
holonomy element the identity.  Choosing an arbitrary sheet over the base
and transporting along the spanning tree is then path independent.  Thus
the local system has exactly

```text
|Omega| global parallel sections.                         (8)
```

The argument needs fixed branches only for one cycle basis, not for every
closed path.  Conversely, one fixed fundamental loop says nothing about the
other generators.

## 3. The regular `C_13` carry specialization

Take

```text
G=Omega=C_13,
a.c=c+a,
chi(a)=a.                                                 (9)
```

Every `a!=0` is a derangement and `chi` is injective.  Therefore a connected
principal carry local system whose fundamental holonomies each have a
same-carry fixed branch has zero holonomy and exactly thirteen parallel
sections.

After trivializing transport on a spanning tree, a selected non-tree edge
represents its fundamental-cycle holonomy.  Equivalently, a common physical
cospan may identify both endpoint fibres in the same gauge.  If that selected
holonomy is translation by `a`, a physical joint array `K_a(c_L,c_R)` is
supported on

```text
c_R=c_L+a.                                                (10)
```

A common-carry entry in this common gauge

```text
K_a(c,c)>0                                                (11)
```

is possible only for `a=0`.  Thus a lawful pair twist or cross-correlation
which recovers (11) for the *selected translation component*, after the
spanning-tree trivialization or equivalent cospan, supplies the fixed branch
required by Section 2.  An arbitrary open-edge diagonal without this gauge
identification is not itself a cycle-holonomy witness.

This is the exact consumer of a two-carry endpoint refinement.  It is not a
claim that such a branch is already present in the LRC canon.  In particular,
the following implications are false:

```text
both endpoint carry marginals nonempty   ==> common carry;
all endpoint/determinant sectors nonzero ==> common carry;
one endpoint has a private root row      ==> fixed transition;
some diagonal mass in a mixed bank       ==> every component is trivial.
                                                               (12)
```

The first two failures are compatibility-fibre loss.  The third forgets the
target endpoint.  The fourth confuses existence of an identity component
with purity of the transition label.

## 4. Sharp finite hostiles

All four boundaries in (12) occur already on `C_13`.

1. Let the transition be `c -> c+1` and give every source and target carry
   positive marginal weight.  Both marginals have full support, while the
   graph (10) is disjoint from the diagonal.
2. A private source atom at `c=0` is transported to target carry `1`; it is
   private but not fixed.
3. Mix the identity graph and the `+1` graph.  The mixture has positive
   diagonal mass, but its charged component remains a derangement.
4. Relabel the target fibre by `c_R' = c_R-1`.  The `+1` graph becomes a
   diagonal.  This is not a contradiction: the relabeling changed the
   endpoint cospan.  A common physical carry gauge is therefore part of the
   hypothesis, not a cosmetic coordinate choice.

The two-array hostile motivating THM-2634 is even stronger and can be checked
directly here.  Let `l,r` range over `F_13^2`.  In model A put

```text
P_0(l)=1,       Q_1(r)=1,                                 (15)
```

for every endpoint and set all other carry bins to zero.  In model B keep
`P_0(l)=1` and replace `Q_1(r)` by `Q_0(r)=1`.  Both models have the same
carry-forgotten marginals and the same full endpoint array

```text
L(l)=R(r)=J(l,r)=1.                                      (16)
```

Consequently every target/determinant pushforward of `J` is identical in the
two models.  In particular, for `q!=0` every `(q,Delta)` sector has value
`13`, while `(q,Delta)=(0,0)` has value `169`.  Yet the common-carry
contraction is respectively

```text
sum_c P_c(l)Q_c(r)=0       and       1.                   (17)
```

Thus full carry-forgotten support, even with its exact sector weights, does
not determine a fixed branch.  Section 3 explains why the missing cross
term, rather than another marginal, is the relevant datum.

The derangement-support hypothesis itself is sharp.  In the natural action
of `S_4`, the sign character is nonzero on a transposition, and a
transposition fixes two sheets.  Hence (3) would be false without (2).

## 5. The two applications and their different missing hypotheses

For a polynomial Keller map, THM-2633 supplies both nontrivial inputs:

```text
Kummer localization: a nonzero character occurs on some inertia element;
open affine image:    that same inertia element fixes a finite sheet.    (13)
```

For the `D_4` source-deck character, (2) holds, so (13) is contradictory.
This is the abelian derangement-character hinge of the `D_4` exclusion,
without replacing THM-2633's geometric inputs or its stronger conclusion
that affine Jelonek inertia normally generates the full monodromy group.

For LRC, THM-2622 classifies the consequence once a carry holonomy is
supplied, but current results do not provide the analogue of both lines in
(13).  Proved and independently audited THM-2623 records carry-forgetting
rows.  Proved and independently audited THM-2635 records imposed
half-tooth/seven-clock data, but lacks both paired endpoint carries and a
selected cycle-holonomy component.  Neither is used as a dependency here.
THM-2625 fills allocated endpoint sectors after separate endpoint
aggregation.  A theorem applying Section 3 must
instead retain

```text
(common word/clock component, c_L, c_R, transition label)                (14)
```

before summing, prove positive diagonal service `c_L=c_R` on a cycle basis,
and show that the serviced entries belong to the selected principal
transition rather than to an unlabeled mixture.  Pair-twist recovery of the
two-carry cross term is a candidate measurement; positivity, semantic-root
typing, and adjacent-clock transport remain separate obligations.

Therefore this theorem is a conditional holotopy bridge, not an LRC row
exclusion.  It changes the residual from the vague request for “more root
support” to the typed fixed-branch service (14).  The scalar ledger remains
`165`.

## 6. Independent audit

An independent hostile audit rederived the fixed-branch contrapositive and
its conjugacy invariance, the free fundamental-cycle/cycle-basis argument,
the injective-holonomy conclusion, and the count of global sections.  It
also checked the regular-`C_13` specialization, all four gauge/mixture/private
row hostiles, the complete endpoint-plane models (15)--(17), and the
`S_4` sign-character boundary.  In particular, the audit checked the
load-bearing quantifier that different fundamental loops may fix different
sheets: this proves only character-triviality until injectivity on the whole
holonomy image is invoked.

QED.
