# Session 208: paper scope and exact bridge reduction

This is an independent elementary reduction motivated by the attached preprint, not a claim to have certified its 45-page geometric proof or a literature-priority claim.

## Source and exact scope

Mohammad Hajiaghayi, *Logarithmic Embeddings into Probabilistic Spanning Trees for Distances and Capacities*, arXiv:2611.05189v1, manuscript dated November 2, 2026. PDF supplied locally as 2611.05189v1.pdf.

- Theorem 1.1 (pp.2–3) claims deterministic polynomial construction, for rational inputs, of separate distributions on actual spanning trees of a connected undirected N-vertex, M-edge graph. Expected all-pairs distance stretch and expected per-original-edge fundamental-cut congestion are O(log N), respectively. Each support has at most M trees.
- The central deterministic oracle is Theorem 6.5 (pp.26–27): arbitrary nonnegative demands a_e independent of positive lengths ell_e admit a tree with sum a_e d_T(e) <= rho_N sum a_e ell_e. The geometry and entropy proof is in sections 4–6. It is not independently verified here.
- Theorem 8.1 and Corollary 8.2 (pp.30–33) convert weighted response oracles to rational distributions, then compress support to M without increasing any coordinate.
- Lemma 9.1 (p.35) is the exact finite-sum identity sum_e c_e d_(T,w/c)(e)=sum_f w_f L_T(f)/c_f, where L_T(f)=sum_e c_e 1[f in P_T(e)]. For f in T this is c_G(delta(S_(T,f))); for f outside T it is zero.
- Propositions 9.3 and 9.5 (pp.36–37) establish all-cut domination and the coordinatewise minimality of fundamental-cut capacities. These are elementary, independently checkable reductions. Theorem 9.6 (pp.37–38) is fractional oblivious routing, not a guarantee for the congestion of one sampled tree.
- Theorem 1.2/section10.2 (pp.3,41–43) separates support requirements from cut approximation: the two-heavy-stars family needs at least ceil(k²/rho) supported trees despite near-unit optimal expected congestion.
- Proposition10.4 (pp.43–44) gives a unit-length/unit-capacity graph on k(k-1)+2 vertices where every one distribution has either stretch or congestion >=(k+1)/2. The two distributions must not be silently identified.
- Section11 (pp.44–45) permits C0=10^15, so this is an asymptotic claim, not an optimized practical constant.

All uses of the new logarithmic theorem should be labelled applications of the preprint's stated theorem, unless its proof is separately certified. No global tournament spectrum assertion follows.

## Exact bridge reduction, independent of the preprint

Let G be any finite connected loopless undirected multigraph with at least two vertices, with positive edge lengths ell and capacities c. Delete all bridges, and denote by G_i its nontrivial connected components (isolated vertices have no optimization coordinates). Every G_i is bridgeless; parallel edges are retained as different physical edges. Write B for the set of deleted bridges.

**Theorem.** Spanning trees of G are in bijection with tuples of spanning trees T_i of the G_i, by

    T = B union (union_i T_i).

For an original edge e in G_i its tree path is exactly its T_i path. For f in G_i its global fundamental-cut load is exactly its local load in G_i. For b in B, d_T(b)=ell_b and L_T(b)=c_b.

**Proof.** Every bridge belongs to every connected spanning subgraph. Contracting every component of G minus its bridges, including isolated-vertex components omitted from the indexed G_i, gives a tree. A simple path between two vertices in one G_i cannot exit that component and return, since the contracted path would traverse some bridge twice. Consequently the unique global tree path between endpoints of any internal original edge stays in T_i. An original bridge's tree path consists of that edge alone. The exact load claim now follows by summing c_e 1[f in P_T(e)]: an external edge cannot use an internal f, and an internal edge cannot use a bridge. The spanning-tree bijection follows by connectedness and the edge count, or by absence of cycles between components. QED.

For a distribution q define S(q)=max_e E[d_T(e)/ell_e] and C(q)=max_f E[L_T(f)/c_f]. Write S_*(G),C_*(G) for their separate minima. Then

    S_*(G)=max(1,max_i S_*(G_i)),
    C_*(G)=max(1,max_i C_*(G_i)).

For a graph with no cyclic component both equal1. The lower inequalities follow by taking marginals. Conversely couple independently any optimal component distributions and insert all bridges; every coordinate has its local expectation or is identically1. Finite families ensure the optima exist. The stretch optimum also equals optimal maximum all-pairs distortion, since edge bounds sum along graph shortest paths and each edge distance is at most its given length.

More strongly, any prescribed bounds on both coordinate families are globally feasible if and only if locally feasible in every component, together with the bridge floor1 when bridges exist. This does NOT assert that a distribution optimal for stretch is also optimal for congestion.

**Support coupling.** If the chosen component distributions have s_i atoms, they can be coupled with at most

    1+sum_i(s_i-1)

global atoms while preserving every local marginal exactly. Lay each distribution's successive probability intervals on [0,1), take their common refinement, and on each interval use the tuple of locally constant trees. There are at most sum_i(s_i-1) interior endpoints. Rational inputs give rational probabilities. This avoids multiplying support sizes. In particular, applying the paper componentwise gives support at most 1+sum_i(|E(G_i)|-1), and distortion O(log max_i |V(G_i)|), with an all-tree graph requiring just one tree. These improve over using total graph size when the cyclic core is small.

## Connection to the actual tournament metagraph

The operation is K_P(T): preserve a directed Hamilton path P and reverse every nonpath arc. Do not replace this by single-arc reversal. Work separately on connected components, without claiming all-order metagraph connectivity. For an ordinary isomorphism-class line multigraph the fibre degree is H(T)/|Aut(T)|, not H(T). Loops carry zero cut load and may be deleted; parallel lines must remain if their capacities are distinct or each line is being charged.

Session207's proved source classification gives T_ab=Tr_a -> C3 -> Tr_b with H=3, |Aut|=3 and ordinary line degree1. Its sole target U_ab has H=g(a,b), is rigid, and g grows tribonacci in either padding parameter. For n>=4 the two endpoints are distinct (source nonstrong, target strong). Therefore the sole line at T_ab is a bridge. Under EVERY spanning tree of its connected component, that edge has stretch exactly1 and congestion exactly1, independently of g(a,b). The transitive class likewise has H=1 and automorphism order1 and hence a degree1 edge (for n>=3).

Thus the unbounded neighbour-degree/count phenomenon from session207 is compatible with exact preservation of those incident edges. It cannot by itself force tree distortion or congestion. The embedded geometry's difficult part lies in cyclic core edges, not in these pendant source edges. This statement concerns ordinary class line multigraphs; labelled H=3 vertices have degree3, so they are not automatically leaves. Converse merging also needs a separate check before the degree1 conclusion is used.

At tournament order n the full labelled metagraph has N=2^(n choose2) vertices before splitting into components. Applying an O(log N) theorem gives O(n²), not O(log n). Passing to a connected component or bridge core can reduce N; no general polynomial-in-n bound on its number of vertices has been proved here.

## Suggested verification

Exhaustively enumerate actual path-complement class line multigraphs for n<=5, verify the H=3 source leaf edges and their exact local loads on each spanning tree for tractable components. Separately test the bridge theorem on small multigraphs with unequal rational capacities and lengths. The proof above is all-order and does not depend on those tests.
