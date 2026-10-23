# Complementary graphs and tournaments: the observer is indispensable

23 October 2026. Elementary all-size proofs, plus exhaustive controls through order five. The no-common-invariant and four-order witness statements below are proved here; no claim of priority in the literature is made. The older ordered-chart and marked-path dictionaries are explicitly inherited.

## 1. Inherited dictionary and the precise question

The prior packet `05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`, sections 1–3, already proves the ordered graph chart, the affine relabeling correction, and the spanning-tree switching trijection. The separate `05-knowledge/results/trijection_reset_20261020/graph_audit.md` proves the marked-path compression onto graphs on one fewer vertex. These are different carriers.

Fix a finite named set V and an observer order pi. Let R_pi point from earlier to later. For a simple graph G on V, define T(G,pi) by reversing precisely the edges of G in R_pi. Thus (G,pi) determines a tournament exactly, and every tournament has one graph coordinate in every observer order.

The question is whether forgetting the observer and replacing G by the unordered complementary pair {G,G^c} leaves a unique tournament isomorphism class. It does not. More strongly, except for the empty/complete pair, it does not even leave a unique class up to converse.

This does not forbid an arbitrary convention that selects a canonical graph adjacency word and then decodes it. Such a rule defines a function, but it is a chosen encoding, not an observer-independent identification or an automatic transfer of ordinary graph theory.

## 2. The two reflections and their diagonal kernel

The exact identities are

    T(G^c,pi) = T(G,pi)^op,
    T(G,reverse(pi)) = T(G,pi)^op,
    T(G^c,reverse(pi)) = T(G,pi).

Graph complementation and observer reversal are two commuting involutions. Each induces converse; their simultaneous product is invisible to the decoded tournament. Therefore retaining the correlated pair

    {(G,pi), (G^c,reverse(pi))}

preserves the decoded tournament, whereas forgetting which graph is paired with which observer loses this distinction. With a fixed observer, an unordered complementary graph pair corresponds exactly to a labeled tournament/converse pair, not necessarily one tournament isomorphism class.

For example, at order four the backward graph with the single edge 13 encodes C3 -> sink. Its complement encodes source -> C3. Their score multisets are respectively (0,2,2,2) and (1,1,1,3), so they are not isomorphic.

This full-n-vertex chart must not be confused with the path-normalized tiling chart. In the latter, complementing the compressed graph reverses only non-path arcs while retaining the directed Hamiltonian path. That operation is not full converse of the original tournament. The inherited path dictionary retains both the actual path and its positions.

## 3. Minimal witness and a complete order-stability classification

At order three, a graph with one edge has only one ordinary graph type. If its edge joins adjacent observer positions, T is transitive; if it joins the first and third positions, T is a directed triangle. The complementary two-edge graphs have the same converse types. Thus the ambiguity remains after pairing G with its complement, and remains after allowing tournament converse. No smaller order can distinguish tournament types.

### Theorem 1: only the trivial graph pair is observer-stable

For n >= 2, the following are equivalent for a simple graph G:

1. T(G,pi) has the same tournament isomorphism type for every observer order pi.
2. T(G,pi) has the same type up to converse for every pi.
3. The number of cyclic triangles in T(G,pi) is independent of pi.
4. G is empty or complete.

The implication 1 -> 2 -> 3 is immediate, because converse preserves cyclic triangles. Empty and complete graphs always encode a transitive tournament, so 4 -> 1.

To prove 3 -> 4, suppose G is neither empty nor complete. There exist distinct a,b,c with g(ac) != g(bc). Otherwise every vertex sees all other vertices identically, forcing all edges of the complete graph to have the same bit. Fix any order of the remaining vertices, denoted rest, and consider

    pi1  = (a,b,c,rest),    pi1' = (b,a,c,rest),
    pi2  = (c,a,b,rest),    pi2' = (c,b,a,rest).

Swapping adjacent a,b changes only the tournament arc ab. Hence the change in cyclic-triangle count depends only on triples containing ab. Every vertex in rest contributes the same change to the two swaps, because it follows both a and b in both orders. Vertex c contributes opposite signs: g(ac) != g(bc) makes exactly one of the two arc-ab orientations close its directed two-step path, and moving c from after to before reverses that path. Consequently

    (cyc(T_pi1') - cyc(T_pi1))
      - (cyc(T_pi2') - cyc(T_pi2)) = +2 or -2.

At least one of the two swaps changes the cyclic-triangle count. This contradicts 3, proving the theorem.

The proof is constructive: four observer orders always exhibit the obstruction. It does not require enumerating all n! orders. Graph complementation preserves the cyclic-triangle count of the decoded tournament, so it cannot repair this witness.

## 4. No nonconstant invariant survives both ordinary quotient rules

Let X = F2^E be the backward-edge cube for one fixed reference order. Ordinary graph relabeling by g acts by the edge-coordinate permutation P_g. Tournament relabeling acts instead by

    A_g(x) = P_g(x) + epsilon_g,

where epsilon_g is the inversion graph of the relabeled reference order. All sums here are modulo two.

### Theorem 2: the common-quotient collapse

For n >= 2, any function f on X that is invariant under every ordinary graph relabeling P_g and every tournament relabeling A_g is constant. The conclusion remains true if either invariance is strengthened to include complement/converse.

Proof: take the transposition exchanging the first two reference positions. Its inversion mask is the single edge 12. Thus A_g P_g^{-1} translates x by that single coordinate. Invariance under both actions makes f invariant under this edge flip. Conjugating by ordinary permutations gives translation by every individual edge. These translations generate all of X, so f is constant. QED.

Equivalently, there is no nonconstant statistic that can be simultaneously an ordinary unlabeled graph invariant and an unlabeled tournament invariant through this fixed-order decoding. This is stronger than the mismatch of object counts. It does not rule out a graph predicate explicitly using the observer; such predicates are exactly how useful transfers work.

The positive common universe is the full ordered chart, not a shared nontrivial quotient obtained by forgetting the orders on both sides.

## 5. The connected incidence universe between the two quotients

There is still a meaningful, many-to-many relation. Its left vertices are ordinary graph-complement isomorphism classes, its right vertices are tournament-converse isomorphism classes, and a backward-edge word x gives an incidence between its two classes. Give each incidence weight equal to the number of words producing it.

Theorem 2 implies this bipartite incidence graph is connected for every n >= 2: components would otherwise define a nonconstant invariant shared by both actions. Theorem 1 says its only left vertex of degree one is the empty/complete pair. Every other complementary graph class reaches at least two tournament-converse classes.

Complete finite controls give:

| n | graph types | complementary graph types | tournament types | converse tournament types | incidence edges |
|---|---:|---:|---:|---:|---:|
| 2 | 2 | 1 | 1 | 1 | 1 |
| 3 | 4 | 2 | 2 | 2 | 3 |
| 4 | 11 | 6 | 4 | 3 | 14 |
| 5 | 34 | 18 | 12 | 10 | 120 |

At order five, five of the eighteen complementary graph types reach all ten tournament-converse types. This is a finite census, not an asserted asymptotic rule. All weights and degree histograms are in receipt.json.

The counts at order four already rule out a bijection of the entire proposed unlabeled object classes: six graph-complement types are neither four tournament types nor three tournament-converse types. The connected incidence construction preserves the discrepancy rather than hiding it.

## 6. What actually transfers

Every theorem about the full ordered graph object can be conjugated through the bijection, but its predicate must be translated. For example, for i<j<k with backward bits a=x_ij, b=x_jk, c=x_ik, the tournament triangle is cyclic exactly when a=b!=c. Its integer indicator is

    c + ab - ac - bc.

Ordinary graph triangles do not represent this predicate. Ordinary graph connectivity is likewise not tournament strong connectivity. In the marked directed-path chart, strong connectivity translates to coverage of every path cut by a backward chord interval, as proved in the inherited graph audit.

A valid map can therefore retain a genuinely useful theorem about interval coverage, joint attachment consistency, or a transported algebraic operator. It cannot infer preservation merely because both objects have the same number of pair coordinates.

The inherited spanning-tree gauge trijection is another valid positive alternative: labeled switching classes correspond bijectively to free graph coordinates relative to each chosen spanning tree. Its lost datum is the cut representative. Hamiltonian path count is not invariant under switching: a transitive triangle and a directed triangle lie in one switching class but have path counts 1 and 3.

## 7. Fibonacci addresses: a valid order, not a color-only repair

Use standard Zeckendorf addresses on positive integers with weights F2=1,F3=2,..., no adjacent selected indices. Compare finite digit strings from highest index down, padding leading zeros. This lexicographic comparison agrees with numerical order: the maximum legal sum using indices 2,...,k-1 is F_k-1, so the highest differing digit dominates every lower digit.

Thus a distinct numerical Zeckendorf address for each vertex supplies an exact observer order. Retaining those addresses makes the ordered graph/tournament dictionary lossless and supplies actual positions for shared ports. This is a genuine structural connection with the Fibonacci address system.

Replacing addresses by only their three colors does not supply that order. More generally, no finite tag alphabet on an infinite address set can repair the unrestricted correspondence: some tag class contains three distinct addresses. Put a one-edge graph on those three vertices and use the minimal witness from section 3. The ordinary colored graph is unchanged up to a tag-preserving permutation while the tournament type changes. The same obstruction applies to any fixed combination of colors and residue tags, including mod 8, unless actual addresses/order are retained as additional data.

This does not invalidate finite tags as pruning invariants or the repository's electrical guards. It proves they are not a replacement for the full observer in this particular identification.

## 8. Reproduction and scope

Run `python3 check.py` and `python3 -O check.py`; both produce identical JSON receipts. Checks use explicit RuntimeError failures, so optimization cannot disable them. They enumerate all backward words and all permutations at n=2,3,4,5; compute both ordinary and affine canonical forms independently; tabulate every weighted quotient incidence; check connectedness and the degree-one classification; and verify the constructive four-order identity for all 1,090 nontrivial labeled graphs in that range.

No assertion here resolves unrestricted tournament metagraph wiring, the Hamiltonian-count spectrum, or a natural all-simple-graph/tournament-type bijection. The payoff is a complete obstruction, a constructive witness, and an exact larger relational universe that retains both sides of the discrepancy.
