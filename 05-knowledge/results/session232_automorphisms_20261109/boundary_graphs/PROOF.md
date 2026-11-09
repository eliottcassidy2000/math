# Endpoint boundaries and faithful automorphism graphs

Status: PROVED elementary reductions and constructions; FINITE-EXACT witness counts and minimum-order certificates, independently checked by subset dynamic programming and permutation enumeration. No claim of literature novelty for the standard group-action or graph-encoding mechanisms. No new proof of the all-order Chan–Pak or presentation gaps is claimed.

## 1. Inheritance and the exact new distinction

The starting facts, not new results, are:

- THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`: tournament Hamiltonian-path counts H never equal 7 or 21.
- THM-4812, `01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md`: the full path-orbit count tc(T)=H(T)/|Aut(T)| also never equals 7 or 21, at every order, with the finite-certificate dependencies and audit scope stated there.
- LEM-003, `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`, and THM-4747, `01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md`: the automorphism action on directed Hamiltonian paths is free.
- THM-4870, `01-canon/theorems/THM-4870-thomassen-merge-reducible-or-obstruction-dichotomy-asymmetric-path-spectrum-3-7-21-square-sum-vertex-deletions.md`: asymmetric spectrum and the distinction between group orders and path counts.
- THM-4986, `01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md`: an exact common multiplicative lemma for two different spectra, explicitly not a graph–tournament object map.

This note adds an explicit object map and a sharp boundary counterexample: the full path-action graph has neither 7 nor 21 connected components, but an intrinsic fixed-endpoint part can have exactly 7 or exactly 21 components. The smallest possible tournament orders for those boundary counts are 6 and 7, respectively. There are infinite nontrivial-automorphism towers of both examples.

## 2. Endpoint quotient, with the necessary stabilizer coordinate

Let P(T) be the set of directed Hamiltonian paths, G=Aut(T), and P_v(T) the paths ending at v. Write e_T(v)=|P_v(T)| and G_v={g in G:g(v)=v}.

**Theorem A (boundary decomposition).** The endpoint map P(T) -> V(T) is equivariant. It induces P(T)/G -> V(T)/G, whose fibre over the orbit of v is naturally bijective to P_v(T)/G_v. In particular

    q_T(v) := |P_v(T)/G_v| = e_T(v)/|G_v|,
    tc(T) = sum over v in V(T)/G of q_T(v).

The sum takes one representative of each vertex orbit. Individual q_T(v) can be zero or even; neither Rédei parity nor the global 7/21 exclusion applies to every summand.

*Proof.* Every path whose endpoint lies in Gv can be transported to one ending at v. Two paths ending at v lie in the same G-orbit exactly when the transporting automorphism fixes v. Freeness on ordered Hamiltonian paths restricts to freeness of G_v on P_v(T). This proves both statements. QED.

The stabilizer is load-bearing. For the Paley tournament P_7, e(v)=27, |G|=21, and |G_v|=3. Its genuine endpoint quotient has size 9. Dividing a single endpoint count by the whole group produces 27/21=9/7, which is not a count. This is an explicit hostile to the naive quotient.

## 3. Sharp small boundary realizations

Specify a labelled tournament by ordering unordered pairs (i,j), i<j, lexicographically, with bit 1 meaning i -> j and bit 0 meaning j -> i.

**Witness B_7.** Six vertices, mask 26. Equivalently, start with every larger label pointing to every smaller label, and reverse only (0,2), (0,4), (0,5). Then

    Aut(B_7)=1,
    (e(0),...,e(5))=(1,5,7,3,2,1),
    H(B_7)=19, q_B7(2)=7.

**Witness B_21.** Seven vertices, mask 290. Start with the same decreasing orientation and reverse only (0,2), (0,6), (1,4). Then

    Aut(B_21)=1,
    (e(0),...,e(6))=(3,9,21,9,7,3,3),
    H(B_21)=55, q_B21(2)=21.

The names B_7 and B_21 refer to the boundary counts, not their orders or total path counts. Both claims are FINITE-EXACT: `verify.py` enumerates all vertex permutations for these two witnesses, directly enumerates paths, and independently applies the subset recurrence

    D({v},v)=1,
    D(S,v)=sum_{u in S\{v}, u->v} D(S\{v},u).

It enumerates every automorphism by degree-class permutations and checks preservation of all arcs.

**Theorem B (minimum orders).** The minimum order of a tournament possessing an endpoint quotient q_T(v)=7 is 6; for q_T(v)=21 it is 7. The same minima hold for the raw endpoint counts e_T(v).

*Proof.* Deleting the final vertex injects P_v(T) into P(T-v), so e_T(v)<=H(T-v). Complete labelled enumeration gives maximum H on orders 1,...,5 equal to 1,1,3,5,15. Thus for n<=5, q_T(v)<=e_T(v)<=5; for n<=6, q_T(v)<=e_T(v)<=15. The displayed witnesses attain 7 and 21 on orders 6 and 7. The only finite input is the exhaustive maxima through order 5, reproduced by the companion program on all 1+2+8+64+1024 labelled tournaments. QED.

## 4. Nontrivial action graphs and infinite towers

For a group G acting on a set X and an inverse-closed generating set S, the undirected action graph has vertices X and edges x--sx, s in S. When the action is free, each connected component is a Cayley graph of G; the number of components is |X|/|G|. This graph is an action graph, not a graph whose full graph-automorphism group is necessarily G. Generator labels are needed if one wants to remember the chosen group action, rather than just the underlying abstract graph.

Let A=>B denote ordinal sum, with every vertex of A dominating every vertex of B. Every Hamiltonian path is exactly a concatenation of a Hamiltonian path in A and one in B. Strong components form a linearly ordered list, so every automorphism fixes each component setwise. Hence H and Aut are multiplicative over the ordinal sum.

For k>=0, define

    T_{k,7} = C3 => ... => C3 => B_7,
    T_{k,21} = C3 => ... => C3 => B_21,

with k copies of C3. Keep the root r corresponding to vertex 2 in the final block. Then

    Aut(T_{k,h}) = (C3)^k,
    |G_r|=|G|=3^k,
    H(T_{k,7})=19*3^k,    e(r)=7*3^k,   q(r)=7,
    H(T_{k,21})=55*3^k,   e(r)=21*3^k,  q(r)=21.

With the standard generators rotating one C3 component at a time, the complete path-action graph is respectively 19 or 55 disjoint copies of the Hamming graph H(k,3), equivalently the k-fold Cartesian power of a triangle. Its endpoint-r induced subgraph consists of exactly 7 or 21 copies of the same Hamming graph. At k=0 interpret H(0,3) as one isolated point. At k=1 the examples are especially concrete:

- n=9: whole action graph 19 triangles, endpoint-r subgraph 7 triangles;
- n=10: whole action graph 55 triangles, endpoint-r subgraph 21 triangles.

The program independently enumerates all paths by adjacency DFS for those two cases and partitions them by the explicitly computed automorphism group: component-size histograms are respectively {3:19}/{3:7} and {3:55}/{3:21}.

Prepending any number of singleton source components preserves all these boundary counts and groups, producing every larger order above each displayed example. This is a boundary-preserving stabilization, not an assertion about a paper's technical stabilization theorem.

## 5. A faithful ordinary simple-graph carrier

The preceding action graph uses paths as vertices. A different graph construction preserves the original tournament itself and its entire automorphism group.

For every tournament T on n>=3 vertices, construct an uncoloured simple graph Gamma(T):

1. Replace every arc u->v by a fresh three-edge path u--L_uv--R_uv--v.
2. Attach a pendant path of length 1 to each L_uv.
3. Attach a pendant path of length 2 to each R_uv.
4. Attach no pendant path to an original vertex.

There are n+5*binom(n,2) vertices, 6*binom(n,2) edges, and 2*binom(n,2) leaves.

**Theorem C (faithful graphification).** For n>=3, T can be recovered canonically from Gamma(T), up to isomorphism, and restriction gives

    Aut(Gamma(T)) ~= Aut(T).

*Proof.* The 2-core of Gamma(T) is exactly the subdivision of K_n obtained in step 1: every vertex in that core has degree at least two within the core, while iteratively stripping leaves removes all attached paths. At every core vertex, the length of its attached pendant path is an intrinsic invariant: 0 for original vertices, 1 for L vertices, 2 for R vertices. Thus every automorphism preserves these three roles. A three-edge core path u--L--R--v recovers exactly the arc u->v. Restriction therefore gives a tournament automorphism. Conversely, every tournament automorphism extends uniquely to arc gadgets and their rigid pendant paths. These operations are inverse and respect composition. QED.

Define an admissible spanning trace in Gamma(T) to be a core walk that visits every original vertex exactly once, moving only in steps original--L--R--original. Such traces biject equivariantly with P(T). Endpoint fibres and stabilizers are also preserved. Thus H, tc, the boundary q values, and the entire path-action constructions transport to an ordinary simple graph with an intrinsic trace predicate.

The predicate is indispensable: Gamma(T) has n(n-1)>=6 leaves, so it has no ordinary undirected Hamiltonian path. This directly rules out the incorrect statement that H(T) becomes the ordinary Hamiltonian-path count of Gamma(T). The graph encoding itself is fully intrinsic, but the counted object is a specially specified trace.

**Concrete order-7 and order-21 graph automorphism groups.** Let R_7=Cay(Z_7,{1,2,3}) and P_7=Cay(Z_7,{1,2,4}). Exact enumeration gives

    (H,|Aut|,tc)(R_7)=(175,7,25),
    (H,|Aut|,tc)(P_7)=(189,21,9).

Their carriers each have 112 vertices, 126 edges, and 42 leaves. Their graph automorphism groups have orders 7 and 21, respectively, while their admissible-trace counts and full trace-orbit counts retain the tournament gaps. These do not claim new abstract realizability of finite groups as graph automorphism groups; they exhibit a concrete faithful connection with the required path-count sidecar.

## 6. Transfer and failure ledger

- Source: tournament plus ordered Hamiltonian paths and endpoint map.
- Target 1: generator-labelled path-action graph. Map: each path becomes a graph vertex and group generators become edges. Preserved: orbits, stabilizers, endpoints, component counts. Lost if labels are removed: which group generators defined the action and identification of the original tournament.
- Target 2: intrinsic ordinary simple graph Gamma(T). Map: each arc becomes an L/R gadget. Preserved: tournament isomorphism type, full automorphism group, admissible traces, endpoint quotient data. Lost by dropping the trace predicate: the intended meaning of path count; ordinary graph Hamiltonian paths are absent.
- Whole-space exclusion: H and tc are not 7 or 21, by inherited theorems.
- Boundary failure: q(v) is 7 or 21 in the explicit minimum-order examples and in the infinite towers.
- Graph group-order failure: 7 and 21 occur as |Aut(Gamma(T))|. This does not contradict a count-spectrum exclusion.
- Naive quotient hostile: P_7 has e(v)/|G|=9/7 but e(v)/|G_v|=9.
- Open: no endpoint-spectrum completeness theorem is asserted. This construction does not settle the open Chan–Pak completeness direction or the corresponding full presentation-spectrum question.

## 7. Reproduction and validation scope

Run `python3 verify.py > results.json` from this directory. Also run `python3 -O verify.py > results-optimized.json` and compare the outputs: all validation uses explicit `require` guards and survives optimized Python. The program uses only the Python standard library.

- Every labelled tournament of orders 1 through 5 is included in the maximum-H certificate.
- For B_7 and B_21, Held–Karp-type subset DP is checked against direct permutation enumeration of paths.
- The two nontrivial-action examples have a third DFS path enumeration and explicit orbit partition.
- A generic uncoloured graph-automorphism enumerator, using degree/color refinement and bijection backtracking, verifies all graph carrier automorphism counts independently of the tournament encoding proof.
- This graph verification covers all 72 labelled tournaments of orders 3 and 4, both boundary witnesses, both C3-stabilized witnesses, and R_7/P_7. It receives only the uncoloured adjacency graph, not the original/port roles.
- An independent worker checked R_7/P_7 with its own paired-refinement graph-automorphism enumerator and intrinsic 2-core decoder: graph orders 112, edges 126, leaves 42, automorphism groups 7 and 21, restriction kernel 1, and recovered tournament adjacency all agree. Its separate code is the sibling `session232-aut-spectrum/verify.py`, functions `encode_graph`, `decode_graph`, `graph_automorphisms`, and `graph_checks`; this script does not depend on it.

No arithmetic estimate, random sample, graph-isomorphism heuristic, or unverified external theorem is used for the new finite claims. The polynomial-size graph encoding theorem is valid at every n>=3; the graph-group computations are checks, not its proof.
