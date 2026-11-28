# Trijection carriers for the moving root and a single-chord repair

2026-11-28. This combines the [original path/cycle-plus-solo/star trijection](../trijection_four_nodes_three_edges_20261020.md), [tree-edge compression](../session90_even_graphs_odd_descent_20261029/trijection/TREE_COMPRESSION.md), the [functional disconnected charts and their obstruction](../session333_smooth_sums_and_trijection_20261114/FUNCTIONAL_TRIJECTION_AND_ROOT_RETURNS.md), and the [moving-root grammar](../chanpak_moving_root_and_module_interfaces_20261128). Distinguish a **marked carrier** (it can contain forward arcs) from a **disagreement graph** (every selected pair denotes a backward arc). Neither is the tournament itself.

## 1. The moving-root grammar lies over a transitive compressed tournament

In chronological order, the letters P, K and M append a vertex b to the current rooted tournament (T,a). All earlier vertices other than a point to b. The marked pair ab is forward for P and backward for K or M. M alone moves the root to b. Mark **every** such pair, including the forward P pairs. The mark is a spanning tree: each step adds one leaf. Its spine is the path followed by the successive roots, with all other vertices pendant at their parent on that spine.

For a fixed tree F, the older endpoint chart maps a free pair uv to the first and last edges f,g on its unique tree path, and transports u→v to f→g. It is bijective: given f,g, choose the outer endpoints of f and g by deleting each edge and taking the side away from the other. Here list the edges in their birth order. At a leaf append, the new marked edge b–a is last. For any older edge f, the free pair represented by (f,b–a) has outer endpoints x,b, with x distinct from a. Hence x→b regardless of the letter. Old free pairs are unchanged. **The entire compressed tournament on the n−1 edge slots is transitive, in birth order**, including for P, K and M words with large Hamiltonian-path counts.

The proof holds for any increasing parent tree in which each new vertex can reverse only its designated parent contact, whether or not the parent pointer is monotone. Conversely, give such an ordered tree and its marked orientations, and take the forward transitive tournament in edge birth order; expansion restores exactly that ordered one-parent construction. One can therefore vary the Hamiltonian count while keeping the smaller tournament completely fixed. The smaller tournament alone is not a scalar inverse.

The [rooted identity MK≃PM](../chanpak_dyadic_tribonacci_and_rooted_rewrite_20261128) genuinely forgets part of this mark. On a singleton, MK marks a three-vertex path with the final root at its middle; PM marks a three-vertex path with the final root at an end. The tournaments are rooted-isomorphic, but the rooted carriers are not. So the rewrite is safe for height search and not a rewrite of marked tree histories.

## 2. What changes for arbitrary n−1-edge graphs

Let D be any simple marked graph on n vertices with n−1 edges, with prescribed orientations on its edges and c connected components. Its [cycle rank](https://www.cip.ifi.lmu.de/~grinberg/t/22s/graphs.pdf) is c−1. It is a tree only if c=1. It is a pseudoforest precisely when all but one components are unicyclic and the remaining one is a tree. In the general case, if t is the number of tree components and \(\beta_i\) is the cycle rank of component i, then
\[
t-1=\sum_{\beta_i\ge2}(\beta_i-1).
\]
This uses \(\sum_i(\beta_i-1)=-1\); it is not a claim that there are only c−1 simple cycles.

**A marked chart for every D.** Pick a spanning forest A inside D, then join its components by c−1 virtual bridges B to get a spanning tree F. Let S be the c−1 edges of D outside A. Use the tree endpoint chart with F. In a completion of D the orientations on B are free bits, while the \(|S|\) corresponding slots in the compressed tournament are prescribed. Match B bijectively to S and move the free bridge bits into those prescribed slots, with any fixed sign convention retained in the mark. This gives a literal bijection
\[
\operatorname{Comp}(D)\longleftrightarrow\operatorname{Tour}(E(F)).
\]
The inverse takes those bits back out as orientations of B and reinstates the prescribed S orientations before expanding along F. The map is natural under isomorphisms preserving the forest, bridges, pairing and sign marks. It does **not** contradict the earlier K4-plus-isolates obstruction to forgetting such choices.

**Exactly how many assignments have to cross components?** Mark a vertex r and try to assign each edge of D to a different vertex of \(V-\{r\}\). If r is in a tree component, the least possible number of edges assigned to a vertex outside their component is \(t-1\). If r is in a nontree component it is t. For the lower bound every tree component not containing r has one fewer edge than available vertices, so needs at least one external assignment. For equality, assign edges of each tree not containing r to the nonroot endpoints of its rooted orientation, choosing one vacant vertex in each such tree. In a cyclic component not containing r choose a spanning unicycle, direct its cycle and direct all attached trees toward it; assign each chosen edge to its tail, covering every vertex. In the component containing r use a spanning tree rooted at r, including the star of all edges incident to r (an acyclic subset). Assign the remaining surplus edges to the vacant sites in the other tree components. The surplus count equals the vacancy count because the total edge count is n−1. Every local assignment is genuinely incident.

This gives a second marked compressed chart. All root edges were assigned to their nonroot endpoints. For each free pair rv, shift its bit into the prescribed nonroot edge assigned to v; keep free pairs not meeting r. This is the earlier tail-slot construction once no cross-component assignment is needed and the edges are oriented by their tails. In general use the prescribed orientation of each receiving slot as its sign convention. Thus \(t-1\), with a root in a tree component, is the sharp number of nonlocal transfers for this particular assignment model, not the number c−1 of forest-to-tree substitutions, and not a claim about all equivariant charts. For example the six-site diamond plus two isolates needs just one nonlocal assignment, while a unicyclic component plus an isolate needs none.

## 3. One chord repairs all 27 finite grammar misses other than 7 and 21

Start with a moving-root word and reverse one arc uv **outside** its marked tree. In the tree-edge chart that flips exactly one arc of the smaller, initially transitive tournament. On the original vertex set add uv to the carrier, and remove one carrier bridge outside the unique cycle. The marked graph now has n−1 edges and two components, one tree and one unicyclic. In the saved witnesses we can always remove a pendant edge, so the tree is a single isolated vertex and nothing is padded. The unicyclic component may have pendants of its own. Its unique undirected cycle is not required to be a directed cycle in the tournament.

A direct single-arc count avoids building a new complete subset table for every chosen chord. For a fixed ordered pair u,v and a subset A of the other vertices, count paths through A that can feed into u (ending at a vertex beating u), and paths through the complementary subset that can follow v (starting at a vertex beaten by v). Empty prefix/suffix counts are one. The product counts precisely the Hamiltonian paths using the block u→v with that split, independently of the arc uv on the underlying two-vertex block. Sum these products for both orders and take the difference. The saved witnesses are rechecked by direct Hamiltonian-path DP *after* the reversal, not only by the difference formula.

**Result.** Every one of the **27** finite targets missed by the original grammar besides 7 and 21 is obtained by this one-chord construction, on at most **10 vertices**. If the cycle is required to have length three, all 27 are still obtained, using at most **11 vertices**. The unrestricted saved carriers have 24 triangular and 3 quadrilateral cycles; the triangle-only bank is separate. The smaller tournament differs from its transitive order at exactly one pair in each witness. For example PPKP gives height 5 on five vertices; reversing the off-tree arc 1→4 (birth labels 0 through 4) gives height **13**. Mark 1–4 and release 0–2, leaving the triangle 0,1,4, one pendant vertex and the isolate 2.

Combined with the prior exact interval, every allowed odd count through **300,001** now has a uniform construction: a moving-root word followed by at most one off-tree flip, instead of a grammar plus unrelated fallback cores. Neither 7 nor 21 is filled, and the construction is not proved cofinite. Subsequent grammar operations on a repaired rooted core are legal as well; the new edge stays identifiable. The next exact question is whether the one-chord seed language admits a positive integer descent, not whether every graph in the n−1 family is reachable by one exchange.

## 4. Why the distinction between graph and carrier matters

If n−1 is only an edge-count restriction on an ordered *disagreement graph*, it cannot narrow the Chan–Pak spectrum at all, even if the graph must be disconnected. Start with any ordered tournament on n vertices with b backward arcs. Append an ordinal, transitive block of s sites; the Hamiltonian-path count is unchanged. If b>n−1, use s=b−n+1 forward sites. Otherwise choose s large enough that \(s+n-1-b\le\binom{s}{2}\), and order that transitive block with exactly \(s+n-1-b\) inversions. Every integer from zero through the maximum occurs as an inversion count, by the elementary process of inserting the next largest letter with any contribution from zero to the number of existing letters. No pair crossing the block boundary is backward, and the combined disagreement graph is disconnected with exactly n+s−1 edges. So it is **the low number of mark exchanges and the controlled smaller tournament**, not bare edge balance, that gives a useful restriction.

## Replay

~~~sh
g++ -O3 -std=c++17 one_chord.cpp -o /tmp/one-chord
/tmp/one-chord | sort -n
/tmp/one-chord --triangle | sort -n
python3 verify.py
python3 verify.py triangle_witnesses.txt
~~~

The independent verifier checks tree-edge transitivity on all 1,093 words through six appends, visibly checks the single off-tree bit and carrier swap for both banks, and recomputes every before/after count. The search omits MK by the rooted rewrite and only asks for witnesses; it does not claim minimum order over all marked or unmarked representations.
