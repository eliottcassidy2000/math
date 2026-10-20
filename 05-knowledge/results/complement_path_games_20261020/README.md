# Four complementary path games, with the starting-vertex quantifier exposed

**Status:** PROVED by the elementary arguments below; FINITE-EXACT validation in `verification.json`. Independent review pending. No novelty claim is made for classical graph facts or extremal bounds.

## 1. The game is a complementary-graph problem

Let `n >= 1`, let `R` be the red graph on all vertices of `K_n`, and let `B = complement(R)` be the graph of remaining blue edges. Player 1 selects `R` in one move; Player 2 then seeks a route in `B`. Coloring a subset of blue edges cannot create a route absent from `B`, so the second move is exactly a graph-existence question, not an additional strategic coloring choice.

- `H(G)`: an undirected Hamilton path exists, visiting each vertex exactly once.
- `W(G)`: a spanning walk exists, with repeated vertices **and repeated edges allowed**.
- `H_s(G)`: a Hamilton path exists with specified first vertex `s`.
- `E(G)`: `H_s(G)` holds for every vertex `s`. This is sometimes called homogeneously traceable; it is weaker than requiring a path between every ordered pair of prescribed endpoints.

A one-vertex path/walk of length zero is allowed. `n=0` is excluded. If repeated edges are forbidden, the object is a trail, not the walk analyzed here.

Write `A/D` for the variant in which red must **not** have property `A`, and blue must **not** have property `D` for Player 1 to win. Thus the four winning conditions are `not H(R) and not H(B)`, `not H(R) and not W(B)`, `not W(R) and not H(B)`, and `not W(R) and not W(B)`.

“Player 2 starts on any node” has two different quantifier readings:

1. **Player 2 chooses a favorable start.** Success is `exists s: H_s(B)`, namely `H(B)`. Blocking Player 2 means *every* start fails.
2. **Player 2 must work from every prescribed start.** Success is `forall s: H_s(B)`, namely `E(B)`. Defeating that guarantee means *some* start fails. A fixed start known before Player 1 chooses edges has the same existence thresholds and minimum costs, by relabeling the constructions below.

If Player 2 is free to choose the start *after* seeing the coloring, reading 1 applies. The red prohibition in the user's wording is unambiguously “there must not be a way,” so it remains `not H(R)` even when blue is tested with `E`.

## 2. Complete classification by order

### Freely chosen start

| Variant | Player 1 can win exactly when |
|---|---|
| H/H | `n >= 4` |
| H/W | `n >= 4` |
| W/H | `n >= 4` |
| W/W | Never |

For `n >= 4`, color a star `K_{1,n-1}` red. It has at least three leaves, so cannot have a Hamilton path: every degree-one vertex would have to be one of only two endpoints. Blue is a clique `K_{n-1}` plus an isolated vertex, so has no spanning walk. This simultaneously witnesses H/H and H/W. Swap the colors to witness W/H.

For `n=2`, a graph without a Hamilton path is empty, whose complement is complete. For `n=3`, it has at most one edge, whose complement contains a three-vertex path. Replacing a prohibition on `H` by a prohibition on `W` cannot make the red legal class larger or the blue failure condition easier, so these small-order obstructions apply to all four variants. At `n=1`, red already has a length-zero Hamilton path and spanning walk, so no legal first move exists.

**Why W/W is impossible.** In a finite undirected graph, `W(G)` is equivalent to connectivity: a spanning walk connects every pair of visited vertices, and a depth-first traversal of a spanning tree gives a spanning walk from any prescribed start. If `R` is disconnected, vertices in different red components are blue-adjacent. Two vertices in the same red component are joined in blue through any vertex in another red component. Thus `B` is connected, in fact of diameter at most two. This also explains why walk-start quantifiers do not matter.

### Every prescribed blue start

Replace the blue `H` column by `E`, not the red prohibition:

| Red prohibition / blue target | Player 1 can defeat the guarantee exactly when |
|---|---|
| not H / E | `n >= 3` |
| not H / W | `n >= 4` |
| not W / E | `n >= 3` |
| not W / W | Never |

At `n=3`, let blue be `P_3`, and prescribe its middle vertex as the start. A Hamilton path exists, but cannot begin there. Red is one edge plus one isolated vertex, so has neither a Hamilton path nor a spanning walk. For all `n >= 3`, use a blue `K_{n-1}` with one extra pendant vertex `a` attached at `s`; there is no Hamilton path starting at `s`. Its red complement is `K_{1,n-2}` plus the isolated vertex `s`, so both red prohibitions hold. The `n<=2` boundary is immediate.

**Important direction:** replacing “no Hamilton path” by “no spanning walk” on Player 1's red graph makes the restriction *stronger*, because it forces disconnection.

## 3. Exact minimum number of red edges

These minima count red edges in a winning legal first move, not blue edges and not the length of a route. They are infinite when the relevant variant is impossible.

For freely chosen starts and `n >= 4`:

| Variant | Minimum red edges |
|---|---|
| H/H | `n-1` |
| H/W | `n-1` |
| W/H | `f(n)` |
| W/W | Impossible |

Here

`f(n) = 2n-5`, except `f(6)=6` and `f(8)=10`.

For every prescribed blue start, the corresponding first and third minima are `n-2` for all `n>=3`; H/W remains `n-1` for `n>=4`.

### A classical degree obstruction, with a self-contained proof

If an `n`-vertex graph `G` has no Hamilton path, then for some integer `1 <= k <= floor(n/2)` at least `k` vertices have degree at most `k-1`. If `G` is connected, then necessarily `k>=2`.

To prove this, adjoin a universal vertex to `G`, producing a non-Hamiltonian graph on `N=n+1` vertices, and add edges until it is maximal non-Hamiltonian. Call the result `J`. For any missing edge `uv`, adding it creates a Hamilton cycle; removing that edge gives a Hamilton path from `u` to `v` in `J`. The familiar path-reversal splice proves `d_J(u)+d_J(v)<=N-1`: along this path, the positions immediately before a neighbor of `u` and the positions occupied by a neighbor of `v` must be disjoint, or they close a Hamilton cycle.

Among all missing edges, choose `uv` with the sorted pair `(min(d_J(u),d_J(v)), max(d_J(u),d_J(v)))` lexicographically maximal, and put `k=d_J(u)<=d_J(v)`. Every nonneighbor of `v` has degree at most `k`, or its missing edge to `v` improves that pair. There are `N-1-d_J(v)>=k` such vertices, and `k<=floor((N-1)/2)`. None is the adjoined universal vertex. Each had degree at most `k-1` in the original `G`, proving the claim.

This is the usual Pósa-type degree obstruction. Historical primary source: [Pósa, *A theorem concerning Hamilton lines* (1962)](https://real.mtak.hu/201261/1/cut_MATKUTINT_7_1_-_2_1962_pp225_-_226.pdf). The present proof and specialization are included so that the results do not depend on an unchecked citation.

### Counting the complement edges

If a set `S` of `k` blue vertices has blue degree at most `k-1`, its total red degree is at least `k(n-k)`. Red edges inside `S` are counted twice, so the total number `m` of red edges satisfies

`m >= k(n-k) - binom(k,2) = g_n(k) := kn - (3k^2-k)/2`.

The right side is concave in `k`; its minimum on an integer interval occurs at an endpoint.

- On `1<=k<=floor(n/2)`, the minimum is `n-1` for `n>=2`. Hence at least `n-1` red edges are needed to destroy all blue Hamilton paths. The red star attains this and is legal for H/H when `n>=4`.
- If red is disconnected, blue is connected, so `2<=k<=q=floor(n/2)`. The lower bound is `min(2n-5,g_n(q))`. For even `n=2q`, `g_n(q)=q(q+1)/2`; for odd `n=2q+1`, `g_n(q)=q(q+3)/2`. Comparing these expressions yields precisely `f(n)` above.
- To disconnect blue, every edge across some nontrivial cut must be red. The minimum cut size in `K_n` is `n-1`, attained by a red star.

For the W/H upper bound, take blue to be a clique on `n-2` vertices with two extra leaves attached to the same clique vertex. The two leaves would force both Hamilton endpoints, leaving no way to visit the rest of the clique, so blue has no Hamilton path. Its complement has an isolated vertex and exactly `2n-5` edges. At `n=6,8`, instead take red to be `K_{n/2+1}` plus `n/2-1` isolated vertices, with respectively `6,10` edges. Blue has an independent set of size `n/2+1`, too large to separate with only `n/2-1` other vertices.

The connected-nontraceable extremal problem is classical, not a new discovery. Independently, substituting `k=n-1` in Theorem 1.3 of [Balister, Győri, Lehel and Schelp, *Connected graphs without long paths*](https://www.memphis.edu/msci/people/pbalistr/nopk.pdf) gives the same bound and extremal constructions.

For the every-start cost, a Hamilton cycle gives a Hamilton path from every vertex. If blue lacks a Hamilton cycle, the same maximal-supergraph splice above, now on `n` vertices, gives a blue nonedge `uv` with blue degree sum at most `n-1`. Its red degree sum is at least `n-1`, so there are at least `n-2` red edges. The pendant-clique construction in Section 2 attains this with legal disconnected red. This argument is used only for `n>=3`.

## 4. Exact separator/path-cover theorem for joins

This gives structural information about **every** disconnected-red coloring, rather than only a winning example.

Let `G=G_1 vee ... vee G_k` be the join of nonempty graphs: all cross-block edges are present. Write `n_i=|V(G_i)|`, `n=sum n_i`, and `p_i=pc(G_i)`, the minimum number of vertex-disjoint undirected paths whose vertices partition `V(G_i)`, with singleton paths allowed. For a vertex `s` of `G_i`, let `p_i(s)` be the corresponding minimum under the additional requirement that `s` is an endpoint of one of the paths.

**Join theorem.** For `k>=1`,

1. `H(G)` iff `p_i <= n-n_i+1` for all `i`.
2. For `s in G_a`, `H_s(G)` iff `p_a(s) <= n-n_a+1` and `p_i <= n-n_i` for every `i != a`.
3. For `k>=2`, `E(G)` iff `p_i <= n-n_i` for all `i`.

The `k>=2` hypothesis in statement 3 is essential. At `k=1`, use statement 2: `E(G_1)` iff `p_1(s)=1` for every `s`. For example `C_5` satisfies `E`, although its path-cover number is not zero.

**Necessity.** Cut any spanning path at its changes of block. The runs inside block `i` form a spanning path cover of `G_i`, separated by outside vertices. There are at most `n-n_i+1` such runs, or at most `n-n_i` if the path starts outside that block. In the start block, the initial run makes the specified vertex an endpoint of its cover. For statement 3, choose a start outside each block, possible exactly because `k>=2`.

**Sufficiency.** If every block has size at most `n/2`, the complete multipartite graph using only cross-block edges has a Hamilton cycle for `n>=3`; an elementary circular arrangement of the block labels proves this. For `n=2`, the one cross-edge suffices directly. Thus every root works in this balanced case.

Otherwise there is a unique dominant block `G_a` with `n_a>n/2`. Put `m=n-n_a`. All inequalities for other blocks are automatic. If no root is prescribed, a cover of `G_a` with at most `m+1` paths can be split into exactly `m+1` paths, then alternated with the `m` outside vertices. If the root lies inside `G_a`, use a rooted cover and place its root path first. If the root lies outside, a cover with at most `m` paths can be split into exactly `m` paths and alternated beginning with the chosen outside vertex. For `k=1`, this is simply a single rooted/unrooted spanning path. Finally, `p_i(s)<=p_i+1`, because a path containing `s` internally can be cut once to make `s` an endpoint. Therefore the inequalities in statement 3 imply every rooted inequality in statement 2. This proves all statements.

The circular arrangement fact can also be obtained directly from the classical Ore degree-sum splice already proved above: vertices in one multipartite part of size `n_i` have degree `n-n_i`, so every missing edge has endpoint degree sum `2(n-n_i)>=n`. A maximal non-Hamiltonian supergraph would contain a missing edge with degree sum at most `n-1`, a contradiction.

**Apply to the game.** If `R_1,...,R_k` are the red connected components, then `B` is the join of the induced graphs `complement(R_i)`. For disconnected red (`k>=2`), the theorem is an exact characterization of blue H, every blue start, and each individual prescribed blue start.

In particular, when red is a disjoint union of cliques, blue is complete multipartite and `p_i=n_i`. Thus blue has a Hamilton path iff `max n_i<=ceil(n/2)`, and blue works from every start iff `max n_i<=floor(n/2)`. In the odd, tight bipartite case, starts in the smaller part fail while starts in the larger part work.

## 5. Not just a disconnected-color trick

For every `n>=7`, both colors can be connected and still neither has a Hamilton path. The following sharper construction was supplied during the independent review and checked separately.

Take a blue clique `C` on `n-2` vertices, choose two distinct vertices `h,v`, delete its edge `hv`, and add two blue leaves `a,b`, both attached only to `h`. Blue is connected for `n>=5`. Removing `h` leaves three blue components, whereas removing one vertex from a Hamilton path leaves at most two, so blue has no Hamilton path.

In red, `h` is joined only to `v`; `a,b` are adjacent and each is joined to every vertex of `C` except `h`. This red graph is connected. Removing `{a,b}` leaves the edge `hv` and `n-4` isolated vertices, a total of `n-3` components. Deleting two vertices from a Hamilton path leaves at most three components. Since `n-3>3` for `n>=7`, red has no Hamilton path either.

Thus the separator obstruction `components(G-S)<=|S|+1`, necessary for traceability, simultaneously defeats both complementary connected graphs. The exact verifier rules out connected-both counterexamples through `n=6`, and independently checks this family through `n=16`. Consequently the first order is **seven**, with the small-order exclusion explicitly **FINITE-EXACT**, rather than an asserted all-order structural classification of those exclusions.

## 6. Honest connections to the tournament/path-count program

### The actual transport: split into maximal runs

The proof of Section 4 uses exactly the maximal-block-run operation behind the repository's **THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`**. The source here is an undirected graph join; the quotient is a complete **bidirected** graph, not a tournament. Replace every undirected edge by both arcs. This preserves undirected Hamilton-path existence, and ordered directed paths in the replacement are exactly oriented traversals of undirected paths. It does not preserve the tournament category, or the number of unoriented paths without reversal factors.

A precise counting lift retains `q_i(c)`, the number of unordered covers of block `i` by `c` **oriented** undirected paths. Each nonsingleton component has two orientations, each singleton one. If `w(c_1,...,c_k)` counts words with exactly `c_i` occurrences of each label and no equal adjacent labels, then

`h_ordered(G_1 vee ... vee G_k) = sum_c w(c) product_i [c_i! q_i(c_i)]`,

where each `c_i` ranges from `1` to `n_i`.

Cutting a path into maximal runs proves the bijection; assigning its distinct cover paths to the occurrences of each label gives the factorial. The sidecar is the entire oriented path-cover profile, not merely the scalar Hamilton-path count. Section 4 is the support/minimum-path-cover shadow of this operation. This is a proved transport, not an identification of directed and undirected path counts.

### A coloring polynomial and a one-edge update

Define `F_R(t)=sum_pi t^(number of red consecutive pairs in pi)`, summing over all `n!` vertex permutations. Then `F_R(0)` counts ordered blue Hamilton paths, and `[t^(n-1)]F_R(t)` counts ordered red Hamilton paths. Thus H/H asks for the simultaneous vanishing of these two extreme coefficients. Complementation gives the exact reciprocity

`F_complement(R)(t) = t^(n-1) F_R(1/t)`.

When a single blue edge `e` becomes red, let `D_e(t)` sum over permutations containing `e`, weighted by the red count on their other consecutive pairs. Then `F_new(t)-F_old(t)=(t-1)D_e(t)`. This is the relevant edge-color update. It resembles the edge-use partition in **THM-082, `01-canon/theorems/THM-082-deletion-contraction-ham-paths.md`**, but is not an arc-reversal identity; use counts and endpoint data, not a tournament degree-only formula.

### The hostile against a cosmetic tournament

Fix an order `1<2<3`, orient blue edges forward and red edges backward. Red edges `12,13` form an undirected Hamilton path `2-1-3`, but become arcs `2->1` and `3->1`, with no directed red Hamilton path. The full tournament still has the mixed-color path `2->3->1`. Thus a fixed-order red/backward gauge destroys the target predicate. The retained order plus color data can reconstruct the graph, but a scalar tournament path count cannot answer the original game.

Every complete tournament has a directed Hamilton path, whereas the permitted monochromatic subgraphs in this game are generally incomplete digraphs. The repository's **HYP-2202, `05-knowledge/hypotheses/HYP-2202-unit-distance-tournament-hampath-flop.md`**, already records the relevant warning: graph-level path existence, complement path existence, and path existence in a fixed tiling orientation are separate questions. No LRC or tournament-spectrum conclusion follows from these game results.

## 7. Reproducibility, inheritance, and boundaries

Run `python3 verify_path_games.py` in this directory. The standard-library-only verifier:

- enumerates all labeled simple undirected red graphs on `1..6` vertices;
- computes all induced-subgraph Hamilton endpoint sets by subset DP;
- checks both quantifier variants, winning counts, and minimum red costs;
- checks all three join criteria on every disconnected-red graph, including every prescribed root and rooted path-cover number;
- independently checks endpoint sets against permutation enumeration for all graphs on `1..6` vertices;
- verifies the all-order witness templates numerically through `n=16` and the connected-both family for `n=7..16`.

These are validation, not a substitute for the all-order proofs. The exact universe, counts, and comparison results are in `verification.json`.

Inheritance pass: closest proved mechanism is THM-3121's run decomposition; canonical hostile is the star/complement-clique pair; corrected near miss is the fixed-gauge orientation in HYP-2202; least-used sidecar here is the rooted path-cover number. The live concept board is complement, connectivity, traceability, root endpoints, path covers, and edge-color response. The decisive lost coordinate is endpoint eligibility, not total edge count.
