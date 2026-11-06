# Size padding amplifies the opposite endpoint of an actual metagraph edge

Status: PROVED elementary construction; finite controls supplied. No literature priority claim.

## Objects and inherited definitions

For a tournament T on n vertices and a directed Hamilton path P, let K_P(T) preserve the n-1 path arcs and reverse every other arc. This is the repository's labelled path-complement metagraph move, not a single-arc reversal or a vertex-extension edge. For n>=3 its distinct neighbours are in bijection with Hamilton paths, so its degree is the ordinary path count H(T). On unlabelled classes the corresponding line degrees are H(T)/|Aut(T)|; after converse merging multiply by one or two. These are different statistics.

Inherited sources: THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md; THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md; session151_spectrum_metagraph_20261103/metagraph/STRUCTURE.md and bridges/OVERLAP_AND_PATH_WALK.md. Inspected via local cached origin/main on 2026-11-06; no claim that the shared old worktree was current.

The first two sources already prove nonstrong vertices form an independent set and identify the transitive vertex's tribonacci partner. Those facts are inherited, not new results here.

## The padding theorem

Let n>=1, T contain P=(v1,...,vn), and U=K_P(T). For k>=2 define

T' = Tr_k -> T -> Tr_k,

where the left block has forward path a1,...,ak and the right block b1,...,bk. Let P' be their concatenation with P and U'=K_{P'}(T'). Then

1. |T'|=|U'|=n+2k and H(T')=H(T).
2. T' is nonstrong, U' is strong.
3. Write R_k=K_{standard path}(Tr_k), and let alpha_k be its number of Hamilton paths not starting at its last vertex. Then

   H(U') >= alpha_k^2 H(U)+1.

In particular alpha_2=1 and alpha_3=2, giving

   k=2: H(U') >= H(U)+2,
   k=3: H(U') >= 4 H(U)+1.

The first rounding uses the oddness of tournament Hamilton-path counts.

### Proof

Every Hamilton path of T' must traverse the left block in its unique order, then one Hamilton path of T, then the right block in its unique order. This gives a bijection and proves the count preservation.

The induced left and right blocks in U' are R_k; the middle block is U. Between blocks, all arcs point from right toward left, except the two retained path junctions ak->v1 and vn->b1. (All left-right arcs point right-to-left.) The reversal of the vertex order is an anti-automorphism of R_k, hence alpha_k also counts its Hamilton paths not ending at its first vertex.

Choose a right-block Hamilton path B that does not end at b1, an arbitrary middle Hamilton path Q in U, and a left-block Hamilton path A that does not start at ak. The concatenation (B,Q,A) is a directed Hamilton path of U': its two junctions avoid precisely the exceptional arcs. Distinct triples give distinct paths, yielding alpha_k^2 H(U) paths. None is P', because these paths start in the right block whereas P' starts in the left. P' supplies one more path. This proves the inequality without needing to classify paths with repeated visits to a block.

For k=2 there is one permitted path in either exterior block; for k=3 the block is a directed triangle with two permitted starts or ends. Finally the retained P' closes through bk->a1, which is a reversed nonpath arc. Hence U' contains a directed spanning cycle and is strong.

## Iteration and the joint size-count consequence

Apply k=3 repeatedly, always retaining the extended distinguished path. Starting from the pair (T,U), after r>=1 steps there is an actual metagraph edge (T_r,U_r) with

   |T_r|=|U_r|=n+6r,
   H(T_r)=h:=H(T),
   H(U_r) >= 4^r H(U)+(4^r-1)/3.

For every attained ordinary count h, choose a realizing T and any Hamilton path P. Therefore the stratum H=h contains vertices, at arbitrarily large orders, adjacent to vertices whose ordinary counts grow at least exponentially in order along this family. The low-count vertices are nonstrong after the first step, and the opposite vertices are strong.

Consequently there is no order-independent finite function F(h) such that every actual metagraph edge satisfies H(U)<=F(H(T)). Nor can one bound adjacent-count difference or ratio using the low endpoint's count alone. This is an explicit relation between size spectrum, ordinary count spectrum and the real adjacency operation, stronger than count-preserving padding by itself. It does not say all odd counts occur, prove all-order connectivity, or provide arbitrary prescribed adjacent-count pairs.

The four-vertex version similarly gives H(U_r)>=H(U)+2r at order n+4r. Two-vertex padding was experimentally tested but no such general amplification claim is made for it here: the naive path (new sink,Q,new source) can hit both exceptional junctions and is not always legal.

## Controls

Run `python3 check.py > results.json`. It enumerates every fixed-base-path tiling of orders 1..4 and pads each by k=1,2,3. A direct subset dynamic program checks counts and the claimed inequalities for k=2,3; k=1 is an explicitly finite extra check. The output has 36 cases. These finite tests supplement, not replace, the injection proof. The checker does not claim an unlabelled class census.

## The amplifier also controls automorphisms and tiling-fibre degrees

For k=3, write N=n+6. In U', the only outdegree-one vertices are a1,a2, and the only indegree-one vertices are b2,b3. Indeed a3 has outdegree two; each middle vertex v has outdegree

   d_U(v)+3-1[v=v1]+1[v=vn],

which lies between two and N-3 (also valid when n=1, with cancelling indicators); b1 has indegree two. Every automorphism fixes the oriented pair a1,a2 pointwise, then fixes a3 as the unique outneighbour of a2. Dually it fixes b2,b3 and b1. The exceptional cross arcs then fix v1 and vn. Conversely every automorphism of U fixing v1,vn extends by fixing the six exterior vertices. Hence

   Aut(U') is isomorphic to Aut(U)_{v1,vn}.

In particular its order is at most |Aut(U)|. Iterating does not enlarge this group, so the normalized fibre H(U_r)/|Aut(U_r)| is also exponentially unbounded.

For T', the totally ordered strong-component decomposition forces every added singleton to be fixed and preserves each old component. Thus Aut(T') is isomorphic to Aut(T), and the low endpoint's fibre H(T')/|Aut(T')| is exactly constant. Padding by the same transitive blocks on both sides also preserves whether T is self-converse: the ordered strong-component word is reversed with each component conversed, and equal prefixes/suffixes cancel. Therefore the low endpoint's converse-merged fibre (one or two times H/|Aut|) is constant too, whereas the adjacent merged fibre is at least its unmerged fibre and grows without bound.

This rules out order-independent degree-based bounds both for the labelled graph and for the actual line multigraph of class fibres (with loops counted twice), including its converse merge. It does not concern the degrees of the graph obtained by discarding parallel lines.

The checker additionally verifies the two automorphism identities for all fixed-path tilings at orders 1..5 (76 cases), by independent degree-cell permutations, not a graph-isomorphism package. The JSON output is emitted only after all checks pass and reports both universes.
