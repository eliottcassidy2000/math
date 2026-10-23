# Actual-path flags transport every rooted motif at every order

23 October 2026. **PROVED:** restriction and deletion commute with the
decorated complement rule below, at every tournament order. **FINITE-EXACT:**
the complete five-vertex flag carrier and its realized hostile controls.
**OPEN:** whether any fixed small motif observer determines all metagraph
wiring. The transfer rule is exact; its small aggregated observations can
still lose information.

## Inheritance and the missing coordinate

The closest mechanism is the actual-path complement in
[the signed-curl packet](../metagraph_port_rigidity_20261023/PROOF.md),
and the finite observer is the
[rooted five-vertex packet](../metagraph_minimal_signatures_20261023/README.md).
The canonical hostile is a same-margin port switch. A corrected near miss
is that a rooted motif which distinguishes two large regular tournaments
need not repair a different small port collision. The least-used sidecar
here is the subset of **actual consecutive path edges** retained inside a
selected motif. The restricted ordering alone does not record those edges.

Board: source path; local induced tournament; retained path forest; rooted
observer; reversible transport; incidence aggregation. The anchor is the
actual-port rule, the niche is its smallest finite carrier, and the
wildcard is deletion-compatible reconstruction.

## Exact local rule

Let T have a specified directed Hamiltonian path
P=(v_0,...,v_(n-1)). Define K_P(T) by keeping those n-1 path arcs and reversing
every other arc. For any vertex subset S let

    F_S = {v_i -> v_(i+1) : both endpoints lie in S}.

This is a directed linear forest, possibly with isolated vertices. It is
not the compressed path obtained by listing S in the order inherited from P.
Write a local flag as (T[S], F_S), with an optional marked root r in S.
Define U on such flags by reversing every arc outside F_S and fixing F_S.
Then

    (K_P(T)[S], F_S, r) = U(T[S], F_S, r),
    U squared = identity.

**Proof.** Each pair in S is a consecutive pair of P precisely when it is
in F_S. That is exactly the condition for its arc to be retained in both
the global and local constructions. The remaining pairs are reversed in
both. Applying the same operation twice restores every arc. A relabelling
of the flagged rooted object commutes with this operation, so U descends
to isomorphism classes of flags. QED.

Thus an integer histogram of source flags transforms by a permutation
matrix. A target rooted-motif histogram is the projection that forgets F
after this permutation. Keep the actual vertex as root while doing this;
an unordered global histogram does not generally reconstruct which
vertices own which records. This is a finite reversible classical carrier,
not an assertion about a physical quantum system.

### Deletion and recursion

For R contained in S, the correct smaller forest is F_S restricted to R.
Deleting a vertex **does not insert a shortcut** between its former
neighbors. With that convention U commutes with every deletion map.
Consequently it commutes with the alternating sum of deletion maps on
ordered flags as well. This gives an all-size compatible hierarchy: the
k-vertex observer is the restriction of the same pairwise operation used
at every larger size. It does not imply that the k-vertex marginal
determines a larger source.

### Why the sidecar is genuinely needed

The frozen receipt contains pairs of actual local sources with identical
rooted tournament type but different target rooted type after transport.
Even retaining the number of path edges does not repair this: a second
witness has the same source rooted type and the same forest cardinality,
yet different targets. The embedding of the marked forest in the source
tournament is the missing information. These are realized examples, not
formal arrays outside the domain of tournaments.

## Complete five-vertex carrier

The same labelled carrier has a closed count at every local size k:

    sum_F 2^(binom(k,2)-|F|)
      = 2^(binom(k,2)-(k-1)) * 3^(k-1).

There are k-1 possible consecutive gaps; each contributes either a retained
forced-forward arc or an unretained arc with two orientations. All other
pairs are free. Every state embeds in a global path presentation on at most
2k-1 vertices by the same gap-insertion construction. Marking a root
multiplies the labelled count by k. This general count does not count
isomorphism classes or remove the need for the forest embedding.

List a selected five-set in path order. There are four successive gaps.
Each gap records whether the corresponding pair really is consecutive in
the full path. A retained pair is forced forward; all other pairs are free.
Therefore the number of labelled admissible (forest, tournament) states is

    sum_F 2^(10-|F|) = 2^6 * 3^4 = 5184.

Each state has a global realization on at most nine vertices: insert one
extra path vertex in every nonconsecutive gap, prescribe the selected
arcs, and choose all remaining non-path arcs arbitrarily. This proves
realizability of the entire local universe. It does not claim every such
state occurs at every smaller global order.

Exhaustive permutation canonicalization and a separate literal global
construction give:

| Carrier | Number |
|---|---:|
| Ordinary rooted tournament types on five vertices | 48 |
| Labelled forest/tournament states | 5184 |
| The same states with one marked root | 25920 |
| Unlabelled forest/tournament flags | 676 |
| Rooted unlabelled forest/tournament flags | 3360 |
| Fixed points of U on unlabelled flags | 24 |
| Fixed points of U on rooted unlabelled flags | 24 |

The quotient tables are stored in the receipt. They retain the forest
**embedding**, not just its abstract component sizes. These numbers are
finite carrier sizes, not counts of tournament metagraph vertices.

## Acceptance tests and reproduction

```text
python -B 05-knowledge/results/metagraph_flag_transport_20261023/check.py
python -B -O 05-knowledge/results/metagraph_flag_transport_20261023/check.py
```

The checker constructs all 5184 admissible local states, embeds each into
a literal global Hamiltonian-path presentation, and independently
restricts its globally complemented adjacency matrix. It canonicalizes
the flags over all permutations, separately over root-fixing permutations,
checks that the quotient maps are well-defined involutions and commute
with all single-vertex deletions, and retains the two sidecar-forgetting
witnesses. Normal replay has 72,583 checks including equality with the
saved receipt (72,582 generating checks).

The labelled convention is a little-endian upper-triangle tournament word;
the four-bit `path_edges` records which consecutive local positions were
also consecutive in the original path. The canonical flag code is
`(forest_pair_mask << 10) | tournament_word`; in rooted canonicalization
the root occupies position zero. Forest directions are supplied by the
tournament, so the undirected pair mask plus the tournament loses no
marked-arc orientation. Every marked arc is retained by U.

## Consequence for actual ports

For any native source/path port, the rule computes each rooted five-vertex
record of its actual target without guessing a target class from margins.
If those records together with other retained features separate a finite
target class universe, the target class of every such port is forced.
The companion five-port experiment tests that separation and authenticates
the affected ports. The all-size theorem here supplies transport; it does
not supply an all-size separation theorem. Full source/path information
remains an exact fallback when a small observer collides.
