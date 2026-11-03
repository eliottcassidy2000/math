# Square-sum interval chains and a Fibonacci carry defect

Status: proved, independently reviewed. All H below count directed Hamilton paths. This extends ZECKENDORF.md from nonconsecutive to arbitrary finite digit sets; it does not claim spectrum coverage or novelty.

## Attachment lemma

Let T be a tournament on [1,u] in which u's only outgoing edge to a smaller label is u→p. Assume no Hamilton path of T−u ends at p. Add k≥1 new interior labels u+1,...,u+k and a new extreme v=u+k+1. Orient all new edges forward except v→u. Call the new tournament T'. Then

    H(T')=(1+2^k) H(T).                                          (1)

Proof: Without v→u, all original vertices must precede all added vertices; append the new labels in increasing order to each path of T. The transition is valid even when its final core vertex is u, since k≥1. Thus this case gives H(T).

For paths using v→u, delete v and the k interior vertices. All surviving adjacencies form a directed path in T except possibly the adjacency just before u. Any invalid such adjacency must be p→u since u→p is its only backward outgoing edge. But then u must terminate the projected order, since its only possible outgoing edge to a lower vertex is  u→p and p is already used. An invalid projection would therefore give a Hamilton path of T−u ending at p, contrary to assumption. Each projection is consequently a path P of T. Conversely P expands by putting each interior vertex either after the part of P before u and before v, or after the part of P after u, both in increasing order. No other placements are possible: the interior vertices are beaten by all original labels, and v→u is the unique backward edge incident to the new labels. This gives 2^k expansions for each P, proving (1).

After extension, the endpoint v has exactly one backward outgoing edge v→u. The tournament T'−v has no Hamilton path ending at u, because each of the k vertices above u has no backward outgoing edges at all, so once such a vertex occurs, no smaller vertex can occur afterward. Hence the lemma iterates.

## Arbitrary square-sum indices, including contiguous strong runs

Let T_n=n(n+1)/2. Begin with core [1,10], all edges forward except 8→1, 10→6. Its count is 265 by NOTE.md. Its top endpoint 10 has the sole backward outgoing edge 10→6, and the core with 10 removed has no Hamilton path ending at 6 (the only remaining reversed edge is 8→1, and vertex 9 cannot appear before 8).

Choose any finite S⊂{4,5,...}, not necessarily nonconsecutive. Work on labels [1,max(10,max_(n∈S) T_(n+1))] if S is nonempty, and [1,10] if S is empty. Reverse only 8→1, 10→6, and each T_(n+1)→T_n for n∈S. All reversed-edge endpoint sums are squares, since T_n+T_(n+1)=(n+1)². Then

    H = 265 ∏_(n∈S) (1+2^n).                                   (2)

Proof: decompose S into consecutive runs. Each run comprises edge intervals which successively share one endpoint. The attachment lemma iterates within each run because every triangular gap is at least two, and the gap-interior vertices ensure its hypothesis for the next step. A run beginning at n=4 attaches to the 265 core using its special zero-endpoint property; a run beginning later starts with the ordinary one-backedge count 1+2^n. Distinct runs and the core if not connected are separated by forward-oriented intervals, so their Hamilton-path counts multiply. In particular S={4,5,...,m} produces a strong tournament: each edge interval contains a spanning directed cycle and successive cycles, including the intersecting core cycles, share an endpoint, covering every vertex.

Finite check: `g++ -O3 strong21.cpp -o /tmp/strong21 && /tmp/strong21` directly counts the 21-vertex strong case S={4,5} by endpoint-subset DP and reports 148665=265·17·33. The 15-vertex case S={4} is independently checked in zeckendorf_check.py. These checks do not replace the proof.

## Zeckendorf carry is a count defect

Index Fibonacci numbers F₂=1,F₃=2,... and map index i to the square-sum interval index n=i+2. Let i≥2, put x=2^(i+2), and let R be an arbitrary finite set of other indices disjoint from {i,i+1,i+2}. The two subsets R∪{i,i+1} and R∪{i+2} have the same Fibonacci sum by F_i+F_(i+1)=F_(i+2), but (2) says their distinct tournament path counts differ by

    265 [∏_(j∈R)(1+2^(j+2))] ((1+x)(1+2x)−(1+4x))
      = 265 [∏_(j∈R)(1+2^(j+2))] x(2x−1).

Thus the binary-weighted square-sum tournament code does not preserve Fibonacci carries: it measures their defect exactly. Choosing the nonconsecutive Zeckendorf normal form removes this representation dependence and gives the injective code in ZECKENDORF.md.
