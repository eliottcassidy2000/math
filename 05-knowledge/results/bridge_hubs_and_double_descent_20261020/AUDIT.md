# Independent recovery audit: boundary hubs and double descent

Date: 2026-10-20. Status: **APPROVED**, with the explicit `ell=1` clarification recorded below. This is a recovery review of the saved proof, not a new research extension. The research review was confined to its isolated audit stage. Publication packaging subsequently changed only paths, links and read-only result handling in a separate stage; the candidate, frozen packets, shared repository and remote were not modified.

## 1. Verdict and exact scope

The following mathematical statements are correct:

1. The common-black-hub sufficient criterion, including its general shared-edge form, certifies cycles of length at most four in the **paired-black multigraph Q**.
2. Every incident Q edge at a rigidly guarded nonstrong node satisfying the stated degree hypothesis lies on a cycle of length at most four.
3. Every strong non-self-converse tournament of order at least three has at least two distinct strictly descending Q edges.
4. A terminal-free side of a Q bridge has a unique triangle-potential minimum at its boundary endpoint; that endpoint is self-converse and has only one strictly descending Q edge. The other side contains a nonstrong node.

There is no substantive gap in these conclusions. The double-descent proof should explicitly exclude `ell=1` after choosing an odd prime divisor of the rotation subgroup. The omitted explanation is immediate: rotation by `n/p=1` would act transitively, forcing all scores equal and contradicting nonregularity. A reviewed copy adds exactly this case explanation. It is a proof-completion sentence, not an additional hypothesis.

This approval does **not** establish all-order black-core connectivity, classify every Q bridge, prove all nonleaf Q edges have short cycles, or exclude every mixed-to-mixed blue bridge in the literal full graph.

## 2. Provenance and source audit

Reviewed original:

- `session15-bridge-continuation/BOUNDARY_HUB_AND_DOUBLE_DESCENT.md`
- SHA-256: `cb1246527eb2cecc663c4d4125f2febe9ca8210a75c1f7684302aebe12b208a3`

Read exact file contents with `git show`, without checkout operations, from commit `7d21b9fdcfaa30e723b92a51ebe02d35ff57b250`:

- `05-knowledge/results/black_descent_20261020c/cycles.md`
- `05-knowledge/results/black_descent_20261020c/scc.md`
- `05-knowledge/results/defect_transitions_20261020/blue.md`
- `05-knowledge/results/blue_quotient_20261020b/structure.md`
- `05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`

The prime-orbit block lemma was checked from its statement and proof, not accepted from its label. The all-order descent argument and the needed black-closed-path input were also checked, as detailed below. The order-nine reconciliation was checked against the precise census text at commit `126830b14ba947b188f5fa6f1c7a080ffe636a82` and the retained local bridge receipts; this audit does not rerun the order-nine census.

A separate freshness comparison uses exact files at `0fccdd7bc00ab288a7dd182b0c8fc50fd72bd4ae`. It is not retroactive provenance for the saved candidate and is not a broad audit of that newer packet.

## 3. Graph conventions and edge identification

There are three different graphs:

- K: literal black-line multigraph, retaining both reflected twins.
- Q: one edge for each reflected pair of black complement lines, retaining loops and all remaining multiplicities. K doubles every Q edge.
- F: K plus the actual blue lines. The actual blue-only subgraph is yet another graph.

The blue quotient obtained by contracting black connected components is not this Q. No inference in the approved proof silently substitutes one graph for another.

In fixed-path coordinates, a Q edge is the orbit of a black tiling `t` under `{1, sigma, tau, sigma tau}`, where sigma reflects the path and tau complements free bits. The underlying source tournaments of `sigma t`, `tau t`, and `sigma tau t` are respectively the converse source, complement partner, and converse complement partner.

For a strictly descending edge from a non-self-converse T:

- sigma cannot identify two path presentations inside the ordinary T fibre, because it lies in the distinct ordinary converse fibre;
- tau and sigma tau cannot identify them inside T, because their triangle potential is strictly smaller, while converse preserves that potential.

Therefore distinct ordinary Hamiltonian-path orbits in T give distinct descending Q edges. Conversely, two words give the same ordinary presentation precisely when the positional map is an automorphism. For cuts of one fixed Hamiltonian cycle, that positional map is exactly the corresponding cycle rotation. This proves the candidate's if-and-only-if identification, including possible equal destinations: parallel Q edges are not collapsed.

For completeness, black V4 edge orbits have four elements. A fixed reflected free coordinate, for example the nonpath endpoint pair, excludes `sigma t=tau t`; blackness excludes `sigma t=t`; tau has no fixed tiling.

At order nine, the retained receipts agree on 97,144 merged vertices, 97,138 black-active vertices, seven Q bridges, five full F bridges, and 24 actual blue-only bridges. All 24 actual blue-only bridges have a blue-degree-one endpoint. The incoming 19 mixed-to-mixed bridges lack a **full F-degree-one** endpoint, not a blue-degree-one endpoint. Its phrase “non-leaf bridge of the blue subgraph” is therefore using an ambient leaf convention and should not be repeated without qualification. The saved candidate correctly disambiguates it.

## 4. Common hub and rigid guards

### Literal common hub

A Hamiltonian path of a nonstrong tournament starts in the first SCC and ends in the last, so its first vertex a beats its last vertex z. For n at least three this is a nonpath arc. Complementing nonpath arcs therefore closes the retained path into a directed Hamiltonian cycle; in particular its partner is strong.

Rotating after a shared directed path edge `u -> v` removes that edge and adds the closing endpoint edge. Applying the second complement consequently reverses exactly `az` and `uv` in the original U. This is a literal labelled equality, not merely matching scores or isomorphism types.

The rotated endpoints are v,u. Their score sum in the common hub is

    d(u)+d(v)-[u=a]+[v=z].

Its inequality with n-1 excludes positional anti-automorphism. Blackness is preserved by the complement involution, so the same rotation line is black at both ends. The score-charge vector is

    q=-e_a+e_z-e_u+e_v,

with sum zero; expanding the cyclic-triangle score identity gives exactly

    c3(V)-c3(U)=-sum d(i)q(i)-(1/2)sum q(i)^2.

This remains valid when u=a or v=z; the charge multiplicities correctly combine. The two reversed unordered pairs remain distinct.

Nonzero triangle change implies `[V] != [U]`, even after converse merging. The partner nodes are strong and hence also different from nonstrong U. The alternate walk uses the explicitly different source edge and two black rotation edges having no endpoint U, so none can be the designated source edge. Removing repeated vertices gives an alternate path of length at most three, and hence a Q cycle of length at most four, including allowed digons. No claim of four distinct nodes is needed.

For `u=a, v=b`, the displayed formulas reduce to the candidate's two tests and to delta=`2d(a)-d(b)-d(z)-3`.

### Rigid guards

The automorphism group acts transitively on the Hamiltonian paths inside either a singleton or a directed triangle. Because every Hamiltonian path traverses SCCs consecutively in their fixed order, independent outer-component automorphisms normalize every source path to fixed literal paths on A and D. They do not change the Q edge. With `|A|>=2`, all Q edges consequently have representatives with the same first, second and last vertices.

If all SCCs were rigid, there would be one ordinary path orbit. A converse-palindromic rigid word is pure blue; a nonpalindromic one has Q degree one. Thus Q degree at least two forces a nonrigid strong SCC of order at least four in B, and hence n at least seven.

The degree bounds in the candidate are valid in every boundary case: first C3 versus first singleton, second singleton versus second C3, and last singleton versus last C3. They yield

    d(a)+d(b) >= 2n-5 > n,
    delta >= n-6 > 0.

Converse can interchange the prefix and suffix to arrange `|A|>=2`. All incident Q edges are nonloops because their other ends are strong. Hence an arbitrary second edge is available under the degree hypothesis and the common-hub criterion applies to every incident edge.

## 5. Double descent and its incoming dependency

### Prime-orbit block lemma

The exact incoming lemma assumes a directed Hamiltonian cycle, n=`p ell`, odd prime p, ell at least three, and automorphism given by rotation through ell positions. These are precisely the hypotheses supplied by the candidate when ell at least three.

Each prime orbit is circulant and has a Hamiltonian word using step +1 or -1. Consecutive orbit cross sets are nonempty. If any is full, their sum is all of F_p. Otherwise reversing the cyclic block order replaces the cross sets by their negated complements; the two total cardinalities add to `p ell`. Choosing the larger direction gives total size at least `ceil(p ell/2) >= p+ell-1`, including the borderline pair `(p,ell)=(3,3)`. The supplied Lagrange-coefficient proof of Cauchy-Davenport is valid and gives the complete sumset. Choosing cross differences summing to the sum of the within-block steps closes all boundary equations. The construction requires neither a prescribed edge on a Hamiltonian cycle nor a false arbitrary-edge Hamiltonicity claim.

The subgroup of cycle rotations preserving T has odd order, because an even-order tournament automorphism group would contain an involution, which reverses the arc on any transposed pair. Descending cut positions are invariant under this subgroup. If there is only one descending Q edge, the Q-edge identification above makes those positions exactly one coset. A trivial subgroup then gives n distinct increasing scores after cutting at the sole descent, impossible in a strong tournament with scores only 1 through n-2.

For a nontrivial subgroup choose an odd prime p dividing its order. The case ell=1 is impossible by transitivity and nonregularity. For ell at least three, the audited block lemma creates a cycle with each prime orbit consecutive. Nonregularity supplies two orbit scores that differ. Cutting an internal edge of each yields two closed paths with equal endpoint scores, each decreasing c3 by one; different endpoint scores forbid their identification by an automorphism. Since the source is non-self-converse, the Q-edge argument applies.

### Exactly two prime orbits

The candidate's alternating cycle implies `0 in S` and `-1 not in S` with its stated indexing. If `alpha in A`, `beta in B`, and `alpha+beta != 0`, a nonempty proper subset of F_p cannot be invariant under that nonzero translation. Thus some x satisfies `x in S`, `x-alpha-beta not in S`. With `t=x-alpha`, the first block boundary requires `t+alpha in S` and the closing boundary requires `t-beta not in S`; both hold. The induced step paths visit all p orbit vertices because p is prime. Internal cuts in the two blocks are inequivalent equal-score descents.

If every alpha+beta vanishes, both nonempty skew sets are singletons, forcing p=3 and B=-A. Every nonempty proper subset of F_3 is invariant under a reflection `d -> 2c-d`. The proposed map swapping U and W has the exact cross test `S=2c-S`; on the induced triangles B=-A gives the required arc reversal. It is therefore an anti-automorphism, contradicting the source hypothesis. This exhausts the two-orbit case.

For a regular source all cycle cuts descend. If only one Q edge were present among them, all rotations would be automorphisms; negation of cyclic labels is then an anti-automorphism of the circulant, again a contradiction. Thus the regular case is also complete.

### All-order one-descent input for bridge cuts

The exact incoming `cycles.md` proves the needed theorem for all strong black-active tournaments, including self-converse ones. Its odd case is the blue-cut-coset reduction plus the audited block lemma. Its even case reduces to exactly two opposite blue cuts with strictly increasing score runs, then performs a terminal-block insertion. The bound `2 <= t <= delta+1 < m`, the prefix-predecessor count, and the exclusion of the old predecessor of the maximum-score endpoint give the required pair of insertion arcs. The resulting closed path has minimum-score start and a last score strictly below maximum, so it descends and cannot be blue. The separate n=4 argument is valid.

The regular-source and all-blue-cycle inputs in `structure.md` were also checked. In even order its explicit two-orbit construction gives equal-score black cuts. In odd order all-blue cuts force a circulant; the supplied prefix obstruction produces black closed paths for n at least seven. The three forbidden prefix families close with increments respectively 3, k and 2, positive under their hypotheses; unit reindexing preserves closure. The four n=7 certificates close, and C3/R5 are correctly excluded as pure blue. Thus the bridge-cut argument does not rely on an unchecked assertion that descent already implies connectivity.

## 6. Bridge-cut consequence

At a minimum of c3 on a terminal-free bridge side, the all-order descent theorem supplies a smaller black neighbor. The only edge leaving the side is the bridge, so the minimum must be its boundary endpoint u, and its opposite endpoint v has smaller potential. A second minimum would have no available boundary edge; it is therefore impossible. Every descending edge at u must leave the side, so it is the bridge alone. Double descent forces u to be self-converse.

Starting at v, any descending path cannot cross the bridge upwards to u. Finite strict descent must end at a nonstrong black-active vertex on the opposite side. The equivalent restrictions for equal-potential bridges and for a non-self-converse higher-potential endpoint follow. This is a restriction, not a complete bridge classification: both-terminal-side bridges and self-converse single-descent bottlenecks remain possible.

## 7. Finite controls and reproducibility

All new checkers use explicit raising conditions, not removable Python assertions. Normal and `python -O` executions produced identical outputs.

### Independent controls

[independent_controls.py](independent_controls.py) imports no producer implementation. It reads the previously saved class representatives through order seven, reconstructs literal paths and quotient-edge keys, and independently checks:

- All 150 strong non-self-converse merged classes through n=7: 12,506 literal paths, 4,771 descending ordinary path orbits, no identification of distinct descending ordinary path orbits in Q, and minimum four descending Q edges at both n=6 and n=7.
- 18,448 ordered descending cycle-cut pair comparisons: equality of Q-edge keys if and only if the corresponding cycle rotation is an automorphism.
- Every skew internal connection-set pair and every cross set with `0 in S, -1 not in S` for p=3,5,7: 2,184 systems, comprising 2,180 verified literal block cycles with two distinct one-triangle descents and four p=3 exceptions with the stated literal anti-automorphism.
- Every nonempty cross-set system for `(p,ell)=(3,3),(3,4),(5,3)`, and every signed within-orbit step vector: 279,488 boundary certificates. An independent finite sumset dynamic program checks coverage, every boundary equation and cyclic closure.

[independent_common_hubs.py](independent_common_hubs.py) uses only this audit's helper code. For every nonstrong source class through n=7, it generates all literal Hamiltonian paths, groups black edge representatives by shared endpoints and shared directed edge, applies the stated sufficient tests, and verifies the literal hub equality, blackness, strong partner, triangle shift and avoidance of each designated Q edge. It checks 1,339 qualifying shared-edge groups and 3,777 alternate routes, including 723 first/second/last routes.

These controls are finite exact checks on the declared source tables, not new independent class enumerations or all-order proofs.

### Producer replay

A copied producer checker was run inside the isolated audit directory, leaving the original output files untouched. Its path-only publication adaptation is [replay_full_atlas.py](replay_full_atlas.py), which requires the authenticated full atlas directory. It replays the rigid-guard certificates through order eight: one eligible n=7 node with five incident edges, and eight eligible n=8 nodes with 60 incident edges. All 65 literal alternate walks pass. Its strong non-self-converse descent checks through n=8 also pass. This is explicitly a replay, not a third implementation.

Portable reproduction from this published directory:

    python -B independent_controls.py
    python -B -O independent_controls.py
    python -B independent_common_hubs.py
    python -B -O independent_common_hubs.py
    python -B replay_guard_certificates.py
    python -B -O replay_guard_certificates.py

These commands verify bundled expected results without changing this packet. Minimal representative projections through order seven are included under inputs, with exact source hashes. The 65 saved rigid-guard walks have a standalone literal replay. The optional historical full-atlas replay is:

    python -B replay_full_atlas.py --atlas-dir /path/to/full/atlas
    python -B -O replay_full_atlas.py --atlas-dir /path/to/full/atlas

The full atlas must include n8/census_n8.json with the hash in [DEPENDENCIES.json](DEPENDENCIES.json); its regeneration source is the [prior atlas packet](../merged_metagraph_black_quotient_bridges_20261020/README.md). Standalone certificate replay does not independently re-enumerate the set of eligible n8 nodes. Results are retained under results and controls. Packaging adaptations are documented in [README.md](README.md). No order-nine census was rerun in this recovery review.

## 8. Narrow freshness, overlap and credit note

Read exact `scc_basins_20261020d/README.md` and `moves.md` at commit `0fccdd7bc00ab288a7dd182b0c8fc50fd72bd4ae` after the original proof audit. The new packet reports paired-R4 blue bypasses improved to three black path reversals, end dissolution improved to four moves, arbitrary-context R4 normalization, and a singleton-padded mixed-source blue bypass of at most five black moves. The all-order singleton statement is explicitly relative to CITED Thomassen input whose original proof the source could not retrieve; this review does not independently certify that new dependency.

Compatibility and credit:

- The literal common-hub mechanism was already in `defect_transitions_20261020/blue.md` at the candidate's original 7d21 commit. The candidate's added role is the nonzero triangle change preventing folding back onto a **black target edge**, and its paired-black Q applications. Credit for the hub operation must remain with that earlier source.
- A black bypass of a blue edge certifies that edge is not a bridge in F. It does not remove its possible bridge status in the actual blue-only graph, and does not by itself certify a Q cycle containing any chosen black edge of the bypass.
- The new R4 outer-SCC and singleton-padding blue statements do not subsume the rigidly guarded Q theorem. Their colors, edge targets and guard hypotheses differ. In particular the singleton-padded **Q short-cycle problem** remains outside the candidate's theorem even where the newer blue bypass theorem applies.
- The new finite coverage is every **nonstrong mixed source through total order eleven**, not a full order-eleven metagraph census; blue lines with two strong endpoints are outside that stated finite coverage.
- The new packet itself leaves arbitrary black-core connectivity and genuine mixed-to-mixed blue bridges open. It is compatible with the candidate's explicit open all-order Q bridge classification.

The candidate's original provenance is unchanged. Newer results belong in a separately dated update or navigation note, not as retroactively claimed inputs to its saved proof.
