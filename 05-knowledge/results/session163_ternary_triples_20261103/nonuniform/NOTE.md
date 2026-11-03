# Resilience and escape of run-divisible modules

Status: the statements below are proved, except the finite controls, which are explicitly identified. No claim is made that every positive odd integer is a tournament path count.

## General exceptional-interface theorem

Let a finite tournament V contain t pairwise disjoint, nonempty vertex sets B_1,...,B_t and a disjoint set S of s exceptional vertices. Assume that, for each i, every outside vertex not in S sees all of B_i in the same direction. (In particular, the edges between distinct B_i are uniform.) Suppose m is a positive integer for which m divides C_k(B_i) for every i and k, where C_k counts ordered covers by k nonempty directed paths. Then

    m^max(t−2s,0) divides H(V).

Proof. Replace each vertex in a Hamilton path by its block symbol B_i, or by its own individual symbol if it lies outside the union, then collapse each maximal consecutive run of the same B_i symbol to one symbol. Partition by this reduced word, retaining no information on run lengths; let k_i be the number of occurrences of B_i. The reduced word fixes all nonexceptional cross-block adjacencies, which are uniform. Only B_i symbols neighboring a vertex of S in the word can be affected by nonuniform edges; at most 2s distinct blocks are touched in this way.

For every untouched block, internal choices for its k_i ordered runs contribute exactly C_{k_i}(B_i), independently of all other blocks. For a touched block, they contribute some integer depending on the word and the relevant boundary conditions. The word class size is therefore the product of these factors (or zero if a fixed cross-edge is wrong). Each untouched factor is divisible by m, so each word class and their sum have the claimed divisibility. The argument includes composite m, arbitrary interactions inside S, and arbitrary other outside vertices. The conclusion for exponent zero is tautological.

A weighted variant allows each block its own modulus m_i dividing all its C_k. With d=max(t−2s,0), the universal divisor is gcd over all d-element U of the products ∏_{i∈U}m_i, with empty product 1. Indeed each word has at least d untouched blocks, and its product is divisible by one of these products. Equivalently, at each prime p the guaranteed valuation is the sum of the d smallest values v_p(m_i).

A sharper bound uses the actual interface graph, with an edge from x∈S to B_i only when x sees B_i nonuniformly. Let ν be the maximum number of blocks covered by a bipartite matching in which each block has capacity 1 and each exceptional vertex has capacity 2. Then m^(t−ν) divides H. For the weighted case set w_i=v_p(m_i). The guaranteed p-adic exponent is at least

    ∑_i w_i − max_matching ∑_{i matched} w_i.

Indeed, only blocks adjacent in the reduced word to an exception that actually sees them nonuniformly need be excluded. Assign each excluded block one such adjacent exception. Each exception has at most two neighbors in the word, so these assignments form a feasible capacity-2 matching. Applying this separately for each prime gives the stated divisor; it is a lower bound, not an assertion that every matching is realizable by a Hamilton path.

S162 supplies a sufficient hypothesis on each B_i: F_B(r) vanishes mod m for 1≤r<κ(m), where κ(m) is the least k with m|k!. Integer binomial inversion then gives m|C_k for k<κ(m), and k!|C_k gives the rest.

## Sharpness for 3: two triangles per exception

Let A and B be directed triangles, put A→B, and add x with x→0 and 1→x, 2→x in each triangle (local cyclic edges 0→1→2→0). This seven-vertex tournament has exactly 41 Hamilton paths. Directly, partition by how many of the six old vertices occur before x: the seven class sizes for k=0,...,6 are 3,5,7,8,7,5,6, summing to 41. The independent endpoint-subset DP in `check.py` verifies the full sequence and all 64 interface choices.

For any t≥2s, take the ordinal sum of s copies of that 41-path gadget and t−2s directed triangles. Then the exact count is

    H=41^s 3^(t−2s),

and the t triangle blocks satisfy the theorem with the s designated interface vertices. This realizes exact 3-adic valuation t−2s, so the bound and the factor 2 are sharp in general. For t<2s, pair as many triangles as possible into copies of the gadget, use a five-path triangle-plus-x gadget for an odd leftover, and add transitive singleton exceptional vertices as needed; the count is not divisible by 3. For strong sharpness, the seven-vertex gadget has F(2)=712, by the independent subset convolution in `check.py`. Given t≥2s and t≥1, put d=t−2s and let Y be the ordinal sum of s gadgets and d triangles. The strong wrapper of Y⊕T_n has exact count

    H = 41^s 3^d + 2^n 712^s 12^d.

Take n≡s (mod 2). The bracket after factoring out 3^d is nonzero mod 3, so this strong tournament has v₃(H)=d. The two wrapper vertices and all added transitive vertices see each triangle uniformly, and the same s exceptional vertices witness the claimed interface. The direct nine-vertex wrapper count 753 in `check.py` independently controls the formula.

For t<2s, the count of the ordinal-sum construction is not divisible by 3 and its F(2) is also not divisible by 3: the 41 gadget has F(2)=712, the single-triangle escape X has F(2)=32, and transitive singletons have F(2)=2. Choose n parity so that H(Y) and 2^n F_Y(2) agree mod 3. The same wrapper then gives a strong count indivisible by 3. Thus sharpness holds for all t≥1, s≥0, including t<2s.

## Prior all-modulus universality retains the specific escape interface

This is a containment corollary of Session 157's strong one-apex modular universality, not a new proof of modular coverage. Write M=2^e m, m odd; S157 proves that every residue mod M, odd when e>0, occurs infinitely often as the count of a strong one-apex tournament. Its construction has orientation word

    1 (0^(k−1)1)^n (10)^a,

with k≥2, n≥1 and a≥1. Choose the initial 1, any two later zeros, and the apex. The apex beats the initial 1 and loses to both zeros; among the three base vertices all edges point forward. These four vertices form precisely X from S162: the apex and first two selected base vertices form a directed triangle, while the other selected vertex has a nonuniform interface to that triangle. Hence all witnesses in this specific S157 construction retain X as an induced subgraph. The same exact count formula given in S157 applies. This does not embed an arbitrary prescribed S162 block B, and the preceding resilience theorem shows that one should not infer arbitrary-block escape from this example.

## A distinct exact family

For b≥1, n≥0, let Y be the ordinal sum of b copies of X and n transitive vertices; apply S162's two-vertex strong wrapper G. Since H(X)=5 and F_X(2)=32,

    H(G(Y)) = 5^b + 2^(5b+n).

In particular, b=1 and n=0,1 give strong counts 37 and 69, whose gcd is 1. For any odd prime power q=p^e for which 5 is primitive modulo q, this single family realizes every residue mod q: select t∈{1,2} so r−t is a unit mod p, choose b≥1 with 5^b=r−t mod q, and then choose n≥0 with 2^(5b+n)=t mod q. This includes q=3^e and 7^e. This corollary is weaker than the all-modulus coverage inherited from S157 but gives an elementary, compact exact-count family.

Source: S157 universal/NOTE.md and S162 square_transfer/MODULE_BARRIER.md and NOTE.md, at repository revision 3b532f92e3170b15771d9ad38182688efe1f526e.
