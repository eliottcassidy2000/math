# Two births cannot repair one missing binary digit

**PROVED by explicit labelled injections.** This is a theorem about the
source-preserving operation B_x of BINARY_PORTS.md, not a restriction on
unrestricted tournament Hamilton-path values. H always means the raw count.

## Endpoint concentration lemma

Let e_x be the number of Hamilton paths ending at a non-sink vertex x of
T, and F=H(T). Then

    e_x <= (F-1)/2.                                      (1)

To see this, take a path p_1,...,p_(n-1),x. Choose its first vertex p_i
beaten by x; one exists because x is not a sink. Move x from the end to
immediately before p_i. If i>1, the preceding p_(i-1) beats x by minimality
of i, so this is still a Hamilton path. Its last vertex differs from x.
The map is injective: delete x and append it to recover the input path.
Thus e_x<=F-e_x; the inherited oddness of F strengthens this to(1).
This elementary lemma is proved here for use, with no external novelty claim.

## The two-birth theorem

Suppose T has **no sink**, put F=H(T), and make a birth B_x adding z.
For any next port y, old or new, let K=H(B_y(B_x(T))). Then

    y old:  K >= (5F+3)/2,
    y=z:    K >= (5F+1)/2.                               (2)

In particular, since F>=3 and K is odd,

    K >= 2F+3.                                          (3)

The hypothesis includes every nontrivial strong source, but does not require
strongness. Both births retain the complete old adjacency and labels.

**Proof.** A Hamilton path P of T has the following two possible certified
insertions of z: immediately before x, always; or at the end, provided P
ends somewhere other than x. A second birth at y similarly allows insertion
of its new vertex w immediately before y, or at the end if the first
extended path does not end at y. Every resulting path is legal. All these
paths are distinct: deletion of z,w recovers P, and their positions then
identify the insertion choices.

Counting this explicit subfamily, according to P's endpoint, gives

    y old, y!=x:   K >= 4F-2e_x-e_y,
    y=x:           K >= 4F-3e_x,
    y=z:           K >= 3F-e_x.                          (4)

For example, when y!=x is old, a P ending at x supplies two paths, one ending
at y supplies three, and every other endpoint supplies four. Since T has no
sink, (1) applies to every old port. Substitution proves(2). Also a tournament
without a sink contains a directed cycle and has F>=3. The smaller bound
in(2) is at least 2F+2; oddness proves(3). QED.

This argument counts a labelled injection, not disjoint scalar reserves.
It is therefore unaffected by automorphisms, re-marking, or repeated ports.
The value in(3) is sharp: a directed triangle admits two births giving9.

## Consequence: arbitrarily long words add no binary transition

For any initial tournament with F>=3, if any word of these births reaches
2F-1 or 2F+1, a single birth already does so. No assumption is made that
intermediate counts themselves follow binary updates.

First, induced-extension path counts never decrease. One direct proof is
that every old path admits insertion of a new tournament vertex; choose its
first beaten old vertex, or put it last if it beats none. Deleting the new
vertex recovers the old path, giving an injection.

If the initial source has no sink, births preserve that property, and(3)
excludes every word of length at least two. If it has a sink, a birth at that
sink merely inserts a new singleton immediately before it, preserving F and
the sink. A birth anywhere else has a nonempty SCC suffix. The suffix lemma
in BINARY_PORTS.md gives I_x>=F, hence the new count is at least3F>2F+1.
Consequently neutral sink insertions cannot prepare a useful later birth.
This proves the assertion in all cases.

Combined with the SCC reduction, the full collection of arbitrary-length
B words realizing a binary target is classified by the single-port tests on
strong cores. For example **no B word at all**, from any source or order with
H=31, can end at61 or63. Likewise13 cannot reach27 and19 cannot reach37 by
this language. The target counts themselves are attained by other operations.

## Audit and route switch

`audit_two_birth.cpp` independently enumerates the Hamilton paths of all687
strong sources below125, explicitly constructs every path in(4), checks its
arcs and uniqueness, and compares with a direct subset-DP count of the final
adjacency. It also builds and checks each endpoint injection in(1).
`word_search.cpp` was the earlier bounded search on five prime source counts;
the proof above supersedes its finite depth limitation.

The positive reroute is not a longer one-port word. Permit the fresh vertex
to beat two old vertices. Exact examples are55->111,61->123,73->145, recorded
in AUDIT_PORTS.json. All three one-port transitions are globally impossible,
so two is the minimum number of outgoing old ports for these one-vertex
constructions. A three-port example gives31->61, and a four-port example
gives63->127; no global minimality of three or four is asserted.

The remaining task is to find a family of these multi-port sources closed
under enough target-dependent attachments, with a genuine inverse-descent
certificate. The finite examples do not supply that closure.
