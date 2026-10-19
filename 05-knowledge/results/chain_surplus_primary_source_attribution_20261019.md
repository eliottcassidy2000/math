# Chain-surplus Gold: primary-source attribution

Date: 2026-10-19. Scope: comparison of packet88 with the sources below; no computation re-audit and no literature-priority claim.

## Source statements actually checked

**Peczarski (2017), pp.107–109, Section2.** The convention is e(P)≥t1+t2, allowing comparisons already determined by the poset. His structural class has a pair with identical lower covers and both full principal upper sets chains. The p.109 relaxation retains the chain above x but replaces the chain above y by P(x<y)≥1/2. The proof chooses adjacent chain vertices around a half-probability crossing and produces a balanced triple by a swap injection. Fact1 converts that certificate to Gold; Fact2 separately covers an oriented pair with at most one slave and probability at least one half. [Publisher PDF](https://cdm.ucalgary.ca/article/download/62349/46785/178291).

**Zaguia, arXiv:1610.00809v3, Definitions1,5 and Section2.** A good pair (a,b), possibly in the dual, satisfies D(a)⊆D(b), U(b)\U(a) a chain, and P(a<b)≤1/2. A very good pair has equal downsets and both upset differences chains. Theorem2 gives a balanced pair from a good pair. Its proof partitions extensions by insertion of a into {b}∪(U(b)\U(a)), proves decreasing bin masses by swapping a with the preceding chain element, and uses a half-mass crossing. [Author-deposited full text](https://arxiv.org/html/1610.00809v3), [PDF](https://arxiv.org/pdf/1610.00809v3).

The earlier insertion-bin argument already appears in **Zaguia (2012), Section2, pp.4–5**, for the selected chain in an N-free poset. [Publisher PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v19i2p29/pdf/).

## Exact match with packet88

Our assumptions D(u)=D(v), U(u)=U(v)∪T, with T a chain, make {u,v} a very good pair: the two differences are T and the empty set. The ordered pair (a,b)=(v,u) is good. Indeed (u,v) is critical, so the ordinary swap injection gives P(u<v)≥1/2 and hence P(v<u)≤1/2.

Therefore the ambient neighborhood class and its balanced-pair consequence are already covered by Zaguia. The monotone insertion law is also already present. The additional equality of our first two bins follows from our stronger equal-downset/nested-upset assumptions.

Peczarski's Section2 is not literally the same marked-pair hypothesis: our common successors may branch, or be incomparable with some extras, so U(u) need not be a chain. Conversely, his two upper chains need not be nested. This is a comparison of hypotheses, not a proof that the collections of entire posets certified by all results in his paper are disjoint or strictly ordered. For one extra successor, his Fact2 already supplies the unweighted Gold conclusion directly.

## Independent deduction from the checked arguments

The following is our deduction, rather than a quotation or a claim of priority:

**Every good pair in Zaguia's sense supplies a weak Gold strategy.** His decreasing insertion bins q0≥q1≥...≥q_(m−1) have q0≤E/2. If q0=E/2, the first threshold is a half split and any legal second comparison works. Otherwise take the first prefix reaching E/2. Write its preceding mass as a, its newly added bin as x, and its remaining suffix as b. Then a<E/2, b≤E/2, and x≤q0≤a. The crossing is neither the first nor the last bin, because the prefix before the last has mass at least E−q0>E/2.

If b≤a, compare at prefix a and then prefix a+x: the heavy suffix splits into x and b, each at most a. If a≤b, reverse the comparisons: the heavy prefix splits into a and x, each at most b. The other first branch is already at most E/2. Thus every path satisfies E≥t1+t2. These are comparisons between a and adjacent vertices of {b}∪(U(b)\U(a)). If the difference is empty, the critical-pair injection and the good-pair probability hypothesis force a half split, which handles the two-bin boundary.

**Corollary.** Every poset with a very good pair has Gold: orient that pair so its forward probability is at most one half. Equal downsets give the required inclusion in either orientation, and both upset differences are chains. Apply the deduction above (or its dual).

This shows that packet88's first-two-equal requirement is stronger than necessary for the numerical Gold step. It upgrades the conclusion explicitly stated in Zaguia's Theorem2 by combining his insertion law with the balanced-triple/Gold mechanism used by Peczarski. It should be described as a short synthesis or corollary of those arguments, without claiming that this synthesis has never appeared elsewhere.

## What our presentation adds, with appropriate limits

Packet88 makes the two-threshold ordering explicit and records all weak-boundary cases. Its monotone-bin proof is an explicit realization of the older half-crossing/balanced-triple mechanism, not a new injection principle.

Our proof also states a nonnegative outside-skeleton-weight version. All swaps retain outside labels and marked occupied positions, so their inequalities survive such weights. For the general good-pair version the weighted hypothesis q0≤E/2 must still be imposed; for packet88 it follows automatically from q0=q1. Neither of the checked papers states this weighted formulation. That supports identifying it as an additional formulation in our presentation, but does not establish literature novelty.

Recommended attribution: “We specialize Zaguia's good-pair insertion-bin argument and make its Gold consequence explicit using Peczarski's half-crossing/balanced-triple mechanism. We also record the outside-skeleton invariant-weight formulation.”

The official University of Warsaw author-homepage endpoint returned HTTP403 during this check. The publisher-hosted Peczarski paper and Zaguia's author-deposited full text were accessible and read directly. No absence claim rests on search-result snippets.
