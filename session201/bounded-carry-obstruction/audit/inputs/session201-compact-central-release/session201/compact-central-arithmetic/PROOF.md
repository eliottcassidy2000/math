# Twelve explicit arithmetic templates suffice, conditionally

Status: **PROVED, independently audited in the stated scope** for the finite-menu reduction
and conditional consequences. Every uniform occupancy statement below is
**OPEN**. This companion does not change the separately reviewed central
arithmetic proof. Session 201.

Use its g, tau, Q, chronological matrices M_W, physical graph compiler, and
complete inverse. In particular tau(u)=(g(u),u mod2,u,1), u>=2, and

    Psi_W(u,v)=tau(u)^t Q M_W tau(v).

There are no free profile coordinates. Both integer indices may be in either
order in each reduced menu below.

## 1. Dagger quotient of the finite normal menu

The reviewed boundary-absorption reduction retains235 central words through
length8, with counts1,1,1,1,3,8,21,55,144 by length. Its value-set equality is
over all actual unit index pairs, not just ordered pairs.

The dagger maps a JZ=OJ rewrite component to another such component. It also
preserves whether the component exposes a binary boundary, since it reverses
the two boundaries and interchanges O,Z. Thus it permutes the retained
components and is an involution on them. The exact identity

    Psi_W(u,v)=Psi_dagger(W)(v,u)

allows one representative from each pair of retained components. This
establishes equality of the union of values, because both index orders are
allowed. It does not license restricting u<=v afterward.

Exact finite generation in the companion verifier gives these sizes:

    Maximum central length    Raw words    Retained before dagger    After
              5                  42                 15                12
              6                 123                 36                26
              8                1095                235               137

The twelve length-five representatives can be chosen explicitly as

    empty, J, JJ, JJJ, JOZJ, JJZJ, JJJJ,
    JOOZJ, JJZZJ, JOZJJ, JJZJJ, JJJJJ.                       (1)

The nontrivial pairings of retained representatives at this cutoff are
JOOZJ with JOZZJ, JOZJJ with JJOZJ, and JJZJJ with JJJZJ. All other retained
components are fixed. `MENUS.json` records the exact expanded matrices and
the corresponding26- and137-template lists. No minimality of these lists
among all possible arithmetic identities is claimed.

## 2. A particularly small remaining arithmetic lemma

**OPEN twelve-template occupancy lemma.** Every odd n>2134117 has

    n=Psi_W(u,v) for some u,v>=2 and W in the twelve-word list (1). (2)

Equation (2) quantifies only over two ordinary integers and twelve fixed
integer bilinear forms in g(u),u,u mod2,g(v),v,v mod2. The recurrence for g
is explicit, and its exact physical compiler is already proved. This is a
sufficient condition for the desired tournament coverage, not an asserted
equivalence to the unrestricted problem.

For searching (2), one may restore the full transpose-closed35-matrix
length-five family and then use the reviewed ordered-index product bound
and congruence inverse. Alternatively keep the twelve representatives and
allow both orders, with 2<=u<=floor(2n/b_W) and
2<=v<=floor(4n/(b_W u)). The two complete searches have the same value set.

The exact inventory through2134117 finds65 gaps in this family. Of these,
29 are the previously established small C1 gaps. Every other one of the
36 gaps has an explicit longer-central-word C1 certificate in `REPAIRS.json`.
The existing27 C2/C3 adjacency-list repairs cover the29 small C1 gaps other
than7 and21. A transitive tournament covers1. Consequently (2) would imply
positive existence for every positive odd count except7 and21.

This finite repair list has no hidden occupancy assumption: the36 listed
source words are replayed coordinate by coordinate and compiled into actual
graphs using the already proved construction. They are not independently
chosen cone points. The27 inherited repair graphs are counted literally by
the reviewed core checker.

## 3. Hostile endpoints and optional equivalent-size alternatives

The complete finite inverse proves2134117 has no representation at central
length at most5, at any index height. It is repaired by

    2134117=Psi_JOOZZJ(5,14256).

Its chronological graph word has24 vertices. Thus the threshold in (2)
cannot be lowered past that target without expanding the template menu or
allowing an explicit exception. No claim that there are no larger gaps is
proved by this finite obstruction.

The26-template length-six family similarly suggests the separate OPEN
lemma that every odd n>17383 is represented. Its six additional finite C1
repairs beyond the29 small gaps are1811,1877,1897,6151,6721,17383, all included
in the36-word repair list. The137-template length-eight statement is exactly
the core's OPEN length-eight lemma after this valid dagger quotient.

These are three different sufficient all-target hypotheses, not proved
equivalent to one another. In particular a long-central witness does not
automatically produce a short-central witness for the same target.

## 4. Conditional size bound and evidence limits

For represented targets in the length-five or length-six family, the core's
exact product constants give graph order at most floor(log2 n)+6. Each of
the36 displayed C1 repairs and27 inherited graph repairs also satisfies
that bound, by direct integer checks on their orders. Hence either the
twelve-template lemma or the26-template lemma would imply a construction
of order at most floor(log2 n)+6 for every allowed positive odd count, using
at most three Hamilton cycles. Above the respective threshold the generic
construction has exactly one Hamilton cycle.

The core's independently replayed finite scan through10000001 has65 gaps
at length5,35 at length6, and29 at length8, with respective last gaps2134117,
17383,2647. That is finite evidence only. The separate author exploration
through100000001 is not needed for this implication or for any finite repair.
Universal occupancy remains OPEN in every formulation.
