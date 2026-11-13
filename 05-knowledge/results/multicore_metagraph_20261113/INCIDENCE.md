# What the metagraph retains and what its nodes forget

**PROVED:** the boundary-window obstruction and simultaneous-pair recovery
below. The tiling, reflection and line-count identities are inherited from
[THM-796, three-sorted recursive incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)
and [THM-1430, the tiling-class metagraph dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md).

## Four levels with different purposes

|Object|Information retained|What it controls|
|---|---|---|
|Ordinary node|T up to isomorphism|A000568 class census|
|Merged node|T up to isomorphism and converse|The user's metagraph nodes|
|End state|(T,last vertex) up to rooted isomorphism|Exact forward growth of presentation weights|
|Path presentation|(T,ordered Hamiltonian path)|Complement partner, line colour, every face deletion|

The **node weight** is its number of path presentations, not H(T) alone.
For a self-converse node it is H(T)/|Aut(T)|; otherwise it is twice that.
The **blue degree** counts incident reflection-fixed presentations. The
**black degree** is the remaining weight. A loop contributes two incidences.
This terminology separates the number of nodes, their weights, and the
number of distinct neighbours; these are different counts.

For a presentation with path 0→1→...→n−1, tau reverses every nonconsecutive
arc and keeps the path arcs. A **complement match** is {t,tau(t)}. Reflection
sigma reverses path order and converses the tournament. A match is blue
exactly when sigma(t)=t, and otherwise black. Sigma and tau commute.
This graph is not the graph of single-arc flips.

The inherited parity mechanism explains the colours without guessing a
drawing. Reflection pairs nonfixed presentations within each merged node.
A self-converse node has odd weight and hence an odd positive blue degree.
A non-self-converse node has zero blue degree. Each black line has a
distinct reflected companion with the same two merged endpoints, so every
black edge multiplicity is even, including loops. Consequently, if M is
the number of merged nodes and B_active the number incident to a blue line,

    A000568(n) = 2 M − B_active.

These constraints do not specify which nodes meet. The exact incidence
still comes from the complement matches of presentations.

## Even bounded information at both ends cannot recover the colour

At order four, use the bit convention of `rooted_growth.py`, where a 1
orients a free arc forward. Bits 2 and 4 give two path presentations of
the same rooted tournament, with last-root canonical rows (2,12,9,1).
The first is blue and the second black. Thus the end state of G1 does not
determine line colour or, as the receipt checks, the merged partner.

There is an all-size strengthening. Fix any q≥0. Let A and B be these two
four-vertex presentations and place each between a source chain and a sink
chain of length q:

    TT_q → A → TT_q,                 TT_q → B → TT_q.       (I1)

They are isomorphic by a map fixing every vertex in both outer chains.
In particular they agree as tournaments with their first q and last q
path vertices individually marked. Reflection swaps the two chains and
restricts to the old four-vertex reflection in the middle. Thus the first
presentation is blue and the second is black. No fixed amount of endpoint
window information determines colour at every size.

The same obstruction holds **among strong tournaments**. For q≥1, reverse
all comparisons between the outer chains, giving the cyclic substitution
C3[TT_q,A,TT_q] and C3[TT_q,B,TT_q]. The quotient cycle makes both
tournaments strong. Reflection plus converse still swaps the outer chains
and acts on the middle as before. The identical boundary-fixed isomorphism
works. At q=0 the original four-vertex tournament is already strong. Thus
the obstruction is not caused by detachable source or sink padding.

This is a precise obstruction to expanding the end state by a fixed number
of extra end markers. The complete rooted graph still controls forward
attachment; the missing information is the internal Hamiltonian ordering.

## The coupled partner remembers the path exactly

Keep T and its complementary partner Q on the **same** vertex set, with
their relative alignment. Their common directed arcs are exactly the
Hamiltonian spine: every other arc was reversed. A directed spanning path
has a unique order, so the simultaneous isomorphism class of the ordered
pair (T,Q) recovers the entire presentation. Its simultaneous automorphism
group is trivial. Forgetting which member comes first gives a complement
match, with its two sheets identified.

Thus the full aligned pair is faithful but does not compress the tiling
space. Merely keeping the two endpoint node names loses their alignment;
THM-796 already records concrete failures of such independently sorted
couplings. This is the boundary between a useful quotient and a complete
description. The all-order two-face pullback in that theorem remains the
exact full-metagraph growth rule; G1 supplies a smaller quotient for the
specific task of forward node-weight enumeration.

Reproduction: `python -B 05-knowledge/results/multicore_metagraph_20261113/blue_incidence.py 7`.
The receipt retains every node weight and coloured edge multiplicity through
seven, and explicit boundary-fixed isomorphisms for q=0,...,8. The proof of
I1, rather than that finite list, gives its all-q quantifier.
