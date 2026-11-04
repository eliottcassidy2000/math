# Original-module collisions require highest prime quotients of order at least eight

## 1. Action, response, and reviewed inputs

Fix a tournament T on a fixed finite set of physical names. A family consists of modules of the **original T**. Its simultaneous action toggles an old edge exactly when an odd number of selected original modules contain both endpoints. Selections may overlap, nest or cross, and repetitions cancel modulo two. Write U for the resulting tournament on those same names.

The full response F_T(a,c) counts Hamilton paths after adjoining actual ports x->y, with a_v=[v->x] and c_v=[v->y]. Full equality means equality at every physical attachment pair. It does not mean equality of selected endpoint marginals, a few coefficients, or independent unlabeled isomorphism classes.

This integration uses the following **separately reviewed** inputs:

1. The Session188 highest-prime frontier at immutable commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7`, `prime_seven_frontier/PROOF.md`, with its independent review. It excludes every highest changed linear quotient and every highest changed prime quotient of order3,5,6. Thus every still-possible highest changed quotient is prime of order at least seven.
2. The Session189 sibling `profiles/PROOF.md`: arbitrary nonempty children, matching corresponding complete profiles, and every fixed ambient are covered for every prime-seven quotient except the specified carousel with seven even children of order at least four. Its complete census checks all2^21 labeled quotients and all197 prime isomorphism types. Its singleton-position and odd-total branches are proved using complete augmented-profile contrasts, not scalar nonvanishing alone.
3. The Session189 sibling `prime7/PROOF.md`: the missing all-even carousel branch has a fixed finite actual-query certificate in every fixed ambient. This is the separate modulo-sixteen theorem, not an extrapolation from a modulo-eight zero bank.
4. The joint independent review in `prime7_audit/REVIEW.md` approves inputs2 and3, including completeness, physical labels, odd-total safety and the all-even certificate. All three source files and the review are pinned in this packet's verification receipt.
5. The reviewed hereditary original-module theorem and highest-changed-node reduction already used in input1: the original strong-module tree remains a common recursive partition; scalar Hamilton equality forces equality of every complete ordered-cover profile on every fixed induced subset; highest changed quotient nodes are disjoint; all changed old edges lie inside them; and a changed originally prime quotient is its converse. The reviewed core converse-complement duality and profile monotonicity are used only for the separate root-family corollary below.

A strong module is a module that overlaps no other module. This is not the same as a strongly connected tournament. The original strong-module tree need not be the canonical strong tree of U. The inherited proofs and independent review are linked and pinned in `INHERITED.json`; their original text is retained without modification.

Together inputs2 and3 give the following complete baseline statement: **every prime-seven quotient converse with arbitrary nonempty fixed physical children and arbitrary profile-preserving target child changes is separated by at least one actual child-uniform query in every fixed ambient.** For exterior vertices, both ports have the source core's original exterior incidences. No seven-child profile or singleton-position assignment is omitted.

## 2. Excluding one highest changed seven-child node

Suppose T differs from U and F_T=F_U. Query a=c=0 forces x,y before all old names, so H(T)=F_T(0,0)=F_U(0,0)=H(U). The hereditary theorem therefore supplies

    q_h(T[S])=q_h(U[S])

for every fixed old-name subset S and every h. Here q_h counts ordered spanning covers by h nonempty directed paths.

By input1, all highest changed original strong quotients are prime of order at least seven. Suppose one such node M has exactly seven original children A_0,...,A_6. These exact physical sets remain modules on both sides. The target quotient is the converse of the source quotient Q, and every corresponding child has the same complete q profile by hereditary equality.

Form an intermediate tournament V by replacing only T[M] by U[M], while leaving every exterior name and arc as in T. The exterior incidences of M have not changed: every old edge that changes lies inside a highest changed node. Apply the complete prime-seven baseline theorem from Section1 in this fixed ambient. It supplies an actual query z=(a,c) such that

    F_T(z) != F_V(z).                                    (1)

Both ports imitate M at every exterior vertex. The query retains the actual names of the A_i and both full physical attachment masks. Descendant changes inside the A_i are allowed because the theorem expressly permits arbitrary profile-preserving target interiors. No literal equality of child interiors is assumed.

Now let N be any other highest changed node. It is disjoint from M and from all the other highest changed nodes. Since both M and N are old modules, their old cross arcs are uniform. Both physical ports have M's incidences to N and therefore see N uniformly. Thus N remains a module in the actual augmented tournament for this **same query z**.

The complete profiles of T[N] and U[N] agree. Exact occupied-run substitution consequently allows replacing T[N] by U[N] without changing the queried Hamilton count. The equality is exact in the integers, not merely in a selected residue. Restore every other highest changed node successively. Their boundaries and pairwise disjointness are unchanged, so each replacement is valid. Every changed old edge lies inside one of these highest nodes; the final tournament is therefore exactly U. Hence

    F_V(z)=F_U(z),

contradicting (1) and the assumed F_T=F_U.

The other highest node N may itself have an unresolved larger prime quotient. Its own detectability is irrelevant: only its complete equal profile and uniform exposure to the same ports are used. Changes outside M are not silently included in a theorem that restricts all selected operations to M.

## 3. Sharpened necessary frontier

Combining Section2 with the previously reviewed exclusions proves:

**Theorem.** If a nonidentity simultaneous XOR of original modules takes a named tournament T to U with the same full named two-port response, then **every highest changed quotient in T's original strong-module tree is prime of order at least eight**.

In particular, if every prime quotient in the original strong-module tree has order at most seven, every original-module XOR target different from its source is detected by the full named response. Linear nodes of arbitrary order were already excluded. The two-child convention is linear, and there is no prime quotient of order four; these inherited facts do not leave a missing lower-order case.

This is a necessary condition on a still-possible collision. It does not assert that a collision exists at order eight, close arbitrary prime orders, or prove unrestricted tournament reconstruction. It is source-relative. It does not assert pairwise injectivity among arbitrary outputs of the original source's module family, because those original modules need not be modules of an intermediate target. It does not cover arbitrary sequences whose later/current modules were not original modules.

## 4. Separate source-relative root-family consequence

Fix a source core Q[A_0,...,A_6] with prime-seven quotient and any fixed rooted ambient. Allow any family of original modules contained inside this core and selecting the whole core oddly. Do not select any module involving an exterior name.

Every proper original module is contained in one child. The target quotient is Q^op, while hereditary original-module monotonicity makes every target child profile coordinatewise at least the corresponding source profile. A whole-child converse caused by the root selection does not change its profile.

Choose a baseline separator z supplied by the complete prime-seven theorem, and pair it with Jz=(1-c,1-a) on **core bits only**. In both queries, keep the same exterior incidences for both ports. Reviewed whole-augmented-profile duality gives opposite exact baseline ambient contrasts. Fix the member whose source-minus-baseline contrast is negative. At that same physical query, uniform child incidence and nonnegative occupied-run substitution make every child-profile increase unable to cancel the strict negative contrast.

Thus the fixed pair {z,Jz} detects every root-containing family of this kind relative to the source; after choosing its baseline-negative member, that one query works for all those targets in the specified ambient. The sign selection can depend on the ambient. No equality between distinct source children is needed, and crossing descendants are permitted by hereditary monotonicity.

This restricted root-family statement is not used to absorb arbitrary exterior changes in Section2. The global collision proof instead uses scalar equality and exact restoration of the other highest nodes.

## 5. Exact physical controls and scope of verification

The checker retains two named examples with two highest changed original strong nodes:

- An ordered ambient containing a prime-seven core and a triangle, both selected as original modules.
- A strong cyclic ambient containing a prime-seven core with a transitive-triple child, a second triangle, and one outer singleton. The original core, a proper two-name interval inside that child, and the second triangle are selected simultaneously. The child interval operation cancels the root's reversal on just its internal edge; the triple child still changes to a different transitive order with the same complete profile. This demonstrates an actual profile-preserving descendant change under a proper selected original module.

For each example it computes every original module and the complete original strong-node structure, verifies that every selected set is an original module, executes simultaneous edge XOR, identifies the exact highest-node antichain, and checks equality of complete ordered-cover profiles on **every fixed induced subset**. It retains all source, intermediate and target rows, all physical child sets, the single lifted actual query and its three exact Hamilton values. Exact restoration of the unrelated highest node is checked at that same query.

Every named bidegree-(2,2) coefficient and basis value of all six physical source/intermediate/target responses is retained, including all zero, cubic and quartic coordinates. Physical labels are never independently relabeled. The focused examples illustrate the hypotheses and restoration mechanism; the theorem follows from Sections1–3 and the separately reviewed all-order inputs, not from extrapolating the examples.

The standard-library verifier and installed C++17 compiler perform the same explicit checks in normal and optimized Python. No software is installed. At most15 vertices occur in a counted augmented control, so exact Hamilton counts fit safely in64 bits. Source pins bind the new baseline proofs, their independent review, the inherited frontier/review and the executable certificate. This integration requires its own independent review before the highest-prime frontier is announced as updated.
