# Reading a count from the metagraph: a three-coordinate tail and an exact inverse

**All-order:** a literal metaline with a transitive pad reads three marked counts, and three successive pad lengths recover them. Iterated point marking recovers absolute counts from normalized incidence in the **framed, cross-order** metagraph. A fixed positive profile has a unique inverse tail length; sufficiently long tails give exact algebraic rounding, but that regime is sparse. The remaining Chan–Pak question is supply of actual profiles at a prescribed target, not a rounding error.

## 1. Three coordinates, not one

A framed tournament is written in the order of a directed Hamiltonian path. Its line partner reverses all arcs not in that path. Let U be either endpoint, with carried frame ending at v. Take the opposite endpoint and append a transitive tail of m vertices by ordinal sum, with its natural frame. Cross the padded frame line. The new endpoint W has old induced tournament U, tail Z_m (upward consecutive arcs, every longer arc points down), and all cross arcs tail to old except v→s₁, where s₁ is the first tail vertex. This is a literal line of the metagraph; an abstract unweighted edge does not retain its frame or the cross-order lift.

Use the data

    a = H(U)
    b = number of ordered covers of U by one or two paths with the first ending at v
        (the second path may be absent)
    c = number of Hamiltonian paths of U starting at v.

Put h_m=H(Z_m), t_m=1+sum_{1≤i<m} h_i, and e_m equal to the number of full tail paths ending at s₁. The exact law is

    H(W)= a h_m + b t_m − c e_m.                                (1)

Without the exceptional arc a path consists of tail then old: subtract the forbidden adjacency s₁→v, giving a h_m−c e_m. With the exceptional arc, the old cover consists of a run ending at v and an optional subsequent run. Independently the tail has an optional preceding run and then a run starting at s₁. Starting at s₁ forces that run to be a consecutive initial interval, leaving an unrestricted Hamilton path on the remaining tail; this gives t_m. This argument also accounts for seam use followed by returns to both blocks. The carried frame and last-tail-to-first-old arc form a spanning directed cycle for m≥2. The two-vertex case m=n=1 uses (1) but is not strong.

The three series, with indices starting at 1, have a common denominator D=1−x−x²−x³ and numerators x+x³, x+x², x−x² respectively. The tribonacci law for h and e_m=h_(m−1)−e_(m−1), e₁=1, are inherited from [THM-4731](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md); t is a prefix sum. For any three consecutive m the rows (h_m,t_m,−e_m) have determinant +2: check m=1,2,3, then shift using the tribonacci companion matrix of determinant 1. For the first triple, putting s=(H₃−H₁)/2 gives b=H₂−s, a=2s−H₂, c=s−H₁. In particular a whole fixed-port sequence determines exactly the three readings. The sequence generating function is

    [x(a+b−c)+x²(b+c)+x³a]/D.

Splitting a full path just after v proves b≥a. Also c≤a and b≡c mod 2, the latter by (1) at m=1 and Rédei parity. Since h is odd and t−e is even, that parity is precisely what the denominator-2 inverse detects. For the natural last port of Z_n, (a,b,c)=(h_n,t_n,e_n) by reflecting the label order. Thus there is an exact addition law

    h_(n+m)=h_n h_m+t_n t_m−e_n e_m.

## 2. Recovering the symmetry lost by projection

Use the **raw** multigraph of all tile lines, before identifying black twins. Count loops twice for endpoint incidence. Its incidence count at a converse-merged vertex is H/|Aut| for self-converse classes, and 2H/|Aut| for chiral pairs. Every class has odd H/|Aut| since the automorphism action on frames is free and H is odd. Reversal acts on frame tilings within a self-converse class, so at least one is fixed. Hence a blue incidence identifies self-converseness and removes the factor 2. The corresponding formula for the already twin-quotiented multigraph is obtained by doubling its black incidences first.

For m≥4, n=|U|≥2, all old vertices have outdegree below n: v gains one arc but has an old predecessor. Each tail vertex has outdegree at least n (s₁ loses one arc but beats s₂). Therefore every automorphism preserves the two sets, fixes the tail pointwise by rigidity of Z_m, and fixes v, the only old vertex beating any tail vertex. Conversely every old automorphism fixing v extends. For n=1 everything is rigid. Thus three normalized incidence readings at successive m≥4 invert to

    (a,b,c)/|Aut(U)_v|.

In particular the ratio of the first decoded coordinate to the original normalized incidence a/|Aut(U)| gives the orbit size of v. There is an all-order **absolute** read by iterating. Choose a terminal port, recover this stabilizer index, then pass to its decorated W with m=4. W is strong; every vertex is an eligible terminal port by rotating a Hamiltonian cycle. Mark the remaining original vertices in succession. At every step the new symmetry group is exactly a point stabilizer of the previous one, in its action on the old vertices; all padding vertices remain fixed. After marking all original vertices the group is trivial. Multiply the decoded indices and the original normalized incidence to get H(U). This uses the marked lifts between metagraph orders, not a claim that one unlabelled fixed-order shadow already recovers its frame fibre.

Another exact restriction on a simple shadow: represent each merged step by reversing its path arcs instead of its nonpath arcs (take the converse of the latter). Each step transfers one outdegree unit from the first to the last vertex, changing no other score. Carry the same labels across a walk. After r steps from the transitive node the L¹ score displacement is at most 2r, hence

    r ≥ (1/2) sum_i |s_i−i|,

with sorted outdegrees s_i and i=0,…,n−1. Sorting gives the least matching cost. A regular tournament on 2k+1 vertices therefore lies at distance at least k(k+1)/2 from the blue gate's transitive origin, if reachable. Connectivity does not make the walk a unit count selector.

## 3. Invert the absolute target, and when not to use the tail

The first two increments of H_m are b+c and 2a+b−c, both strictly positive, and all later increments are positive by the recurrence. For a fixed profile binary search therefore gives the unique tail length or rejects the target. The exact physical example

    U outgoing bit rows (20,13,16,5,10,31), v=4 (zero-based)
    (a,b,c)=(11,26,0), m=2, H(W)=11+2·26=63

puts the first surviving point of the 7-chain inside the grammar. Here H=63 is not extracted by dividing another witness, and the 11-core is not used with its endpoint data erased.

For arbitrary target N, at each candidate m solve the integer plane (1) and require an actual terminal profile. Always N≥a(h_m+t_m−e_m)=a h_(m+1) and N≥b t_m, using b≥a, c≤a and e_m≤h_m. In particular any m≥2 gives a strictly smaller a and a bounded candidate plane, but **a smaller H witness alone is not a smaller marked-profile witness**.

There is an exact real approximation, with a sharp change of regime. Let α>1 be the real root of X³−X²−X−1; the other roots are complex conjugates of modulus α^(−1/2). The cubic is irreducible by the rational root test and has one real root (its two critical values are negative). Solving the recurrence expresses H_m=Cα^m plus the two conjugate terms, with C in Q(α) and linear in the profile. For B=max(a,b,c), their sum is O(Bα^(−m/2)) with an absolute constant. Thus H_m is the nearest integer to Cα^m for sufficiently large m, and uniformly so when α^m is larger than a fixed constant times N^(2/3), since B=b≤N/t_m.

The same scale controls collisions. For the difference of two distinct integer profiles, the dominant coefficient C cannot vanish: Galois conjugacy would kill every coefficient, contradicting the determinant 2. After multiplication by a fixed denominator it is an algebraic integer; its conjugates are O(B) for profile differences bounded in magnitude by B. The nonzero algebraic norm yields |C|≫1/B². Equality of the two readings at height m would then force α^(3m/2)≪B³. In particular the linear read is injective on integer profiles in a sufficiently small fixed multiple of the α^(m/2) box.

The long-tail sector cannot itself fill the spectrum. For 0<ε<1/3, targets N in [X,2X] with α^m≥N^(2/3+ε) have at most O(X^(1−3ε) log X) possible values even allowing every integer profile in the bounding box; there are O(log X) heights and O(X^(1/3−ε)) choices per coordinate. Exact rounding and injectivity do not imply coverage. Product closure and the known holes reduce coverage to odd 3-free targets other than 7, together with one witness 63. They do not imply that every such target has a representation by this particular terminal-tail grammar.

## 4. The first 32 as two different kinds of boundary

The first 32 odds have holes at positions 4 and 11. This yields triangular run sizes T₂,T₁,T₃,T₁,T₆ = 3,1,6,1,21. Also exactly 21 of the first 32 odds are 3-units, since the divisible positions are 2,5,…,32. These are equal cardinalities, not the same subset or a new periodic exclusion. The precise infinite decomposition is into chains b·3^j with b odd and 3∤b. Prepending a directed triangle multiplies the count by 3. The 7-chain has exactly the initials absent and all j≥2 present, by the 63 witness; every realized unit b supplies its entire chain. Thus the global unresolved problem is the other 3-unit roots, not additional exceptional sites on the 7-chain.

The square-sum graph gives a more useful analogy than the matching number 32. Its first triangle, at 30, is (6,19,30), from sums 25,36,49. This inserts 30 between an already adjacent pair of old vertices. Vertex 32 has precisely two old neighbors, 4 and 17 (sums 36 and 49), which are not adjacent. A spanning cycle at 32 is therefore exactly a spanning path in G₃₁ with prescribed endpoints 4 and 17, sealed by 32. The independent replay finds one such cycle and it uses the triangle subdivision 6–30–19. This is the distinction between merely having a local return and a return closing the required absolute boundary. See the repository's prior square-sum arrival anatomy and endpoint-aware kernel (Session 32, and THM-4505).

There is even a literal target-selected local construction on the square side. Fix x>0. Pick odd a large enough that b=(a²−2x−1)/2>a, put c=b+1, y=a²−x, z=b²−x. Then (x,y,z) is a square-sum triangle with roots (a,b,c), because y+z=(b+1)². It is primitive since a prime dividing all vertices would divide the consecutive roots b,c. Taking the first sufficient odd a gives z=O(x). The identity forces the factor 1 in c²−b²=(c−b)(c+b); this freely selectable factor is the feature missing on the tournament side. Also in the graph with sums 3 times a square, the ±1 mod 3 residues form a bipartite block, so every odd cycle lies on multiples of 3 and divides by 3 to a square-sum odd cycle. A second tripling puts the scaling into the square graph. This is an exact channel switch, not a proof about the tournament hole values.

A related tournament obstruction explains why local insertion is harder. For a quotient G at a marked vertex r, let h₀=H(G), g₀=H(G−r) (empty convention 1), and W₂ the number of Hamiltonian quotient words with exactly two visits to r. Replacing r by TT₂ yields h₀+2W₂ paths: one visit has only one internal path, and two visits give the two singleton orders. By the [odd-cycle packing identity](https://arxiv.org/abs/2307.05569), every old pack touching r lifts to at least two packs, with one or the other twin on that cycle. Untouched packs lift at least once. Thus 2W₂≥h₀−g₀. Inserting any B instead gives

    H(G[B]) ≥ h₀ H(B) + (h₀−g₀)/2 · (P_B(2)−2H(B)),

where P_B(2) counts two ordered slots allowing empties. If G is nontrivial strong, r lies on a directed triangle (otherwise the arcs between its in- and out-neighbors form a one-way separation), so W₂≥1. The (|B|−1)H(B) distinct cuts of full B paths then give H(G[B])≥(h₀+|B|−1)H(B). Consequently a universal strong constant multiplier cannot be obtained simply by widening one homogeneous module. In contrast ordinal ×3 needs no return at all.

## Replay

Run `python check.py` and `python -O check.py` (same output). This verifies literal padding/DP counts, pointed automorphism orders, determinant controls, the order-eight 63 witness, and the square-sum cycle. It is a finite check of the stated examples, not a certificate of complete coverage.
