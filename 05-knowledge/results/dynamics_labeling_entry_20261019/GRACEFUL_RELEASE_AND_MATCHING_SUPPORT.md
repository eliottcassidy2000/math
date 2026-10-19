# Bounded release certificates and matching obstructions for graceful fibers

19 October 2026. Status: PROVED structural criteria, with a complete finite-exact audit of the earlier 13-label obstruction. These results reconfigure already graceful trees. They do not preserve abstract tree shape and do not prove the Graceful Tree Conjecture. The matching-switch mechanism is standard; the point here is its exact transfer to fixed-edge graceful fibers and a local certificate for releasing a small number of fixed edges.

## 1. Setup and the inherited rank-staircase input

Use vertices 0,...,m. Color d consists of all edges {i,i+d}, 0<=i<=m-d. Let F be an acyclic fixed selection containing one edge of each fixed color. Let C be the r remaining colors. Since |F|=m-r, F has r+1 components. Contract F and write

    q_F(d)=rank(F union E_d)-rank(F),

where E_d is the complete difference-d edge set and rank is graphic-matroid rank. Thus q_F(d) is the rank of the color-d contact graph on the contracted components. All actual edges, including parallel realizations, remain distinct choices.

The existing rank-staircase theorem says: if the free colors can be ordered c_1,...,c_r with

    q_F(c_i)>=r-i+1,

then the entire completion fiber is nonempty and connected by exchanges changing at most two colors. Its diameter is at most 2^(r-1) for r>=1. The first color has full contracted rank r and supplies the connected buffer; the remaining colors satisfy the partial staircase. This criterion is sufficient, not necessary.

The present result tells exactly what releasing a fixed edge does to these ranks, so a small release can be certified without enumerating the completion graph.

## 2. The exact rank-gain formula

### Lemma 1: one-edge gain

Let e={u,v} be an edge of F, put G=F minus {e}, and let d be an originally free color. Then

    q_G(d)-q_F(d) is either 0 or 1.

It equals 1 precisely when u and v are connected in the graph G union E_d.

Proof. Both F and G are forests, with rank(F)=rank(G)+1. Adding e to G union E_d raises its rank by 0 if its endpoints are already connected, and by 1 otherwise. Subtracting the two definitions of q gives the formula. QED.

The test is not whether a single d-edge directly joins u to v. The permitted reconnection may alternate many remaining fixed edges and d-edges. It is a connected-component computation on one explicit graph.

### Theorem 2: one released edge repairs a one-defect staircase

Suppose the original free colors can be ordered c_1,...,c_r so that

    q_F(c_i)>=r-i     for every i.

Call row i deficient when q_F(c_i)=r-i. Suppose a fixed edge e of color b has both properties:

1. The complete color-b graph becomes connected after contracting G=F minus {e}, equivalently q_G(b)=r+1.
2. For every deficient row i, the two endpoints of e are connected in G union E_(c_i).

Then every two completions containing F can be joined by two-color exchanges that keep every edge of F except e unchanged. More strongly, the entire completion fiber with only G fixed is connected, of diameter at most 2^r.

Proof. Removing e creates r+2 components and r+1 free colors. By Lemma 1, every deficient old row gains one rank; all other old ranks do not decrease. Thus

    q_G(c_i)>=r-i+1.

Put b first, then c_1,...,c_r. Property 1 gives the first rank r+1, and the displayed inequalities are exactly the remaining staircase. Apply the inherited theorem. Both original completions are vertices of this enlarged fiber. Their endpoints still contain e, so a path between them restores it automatically. QED.

The theorem bounds the number of originally fixed edges that may ever be disturbed: exactly one allowed buffer. It does not say that the buffer edge changes on only one step, nor that target agreement is monotone.

### Corollary 2.1: a general b-edge certificate

Let B be any set of b fixed edges, G=F minus B, and C' the r+b free colors. If the q_G values admit a staircase ordering of length r+b, the completion fiber is connected with diameter at most 2^(r+b-1). At most b originally fixed edges are ever released.

The exact gain of an originally free color d is

    q_G(d)-q_F(d)
      = b - [rank(G union E_d union B)-rank(G union E_d)].

It lies between 0 and b. This formula permits a parameterized certificate search over small subsets B using rank tests alone. It is a sufficient-condition search, not a claimed polynomial-time algorithm for finding the smallest successful release.

For arbitrary F a coarse certificate always exists: release every fixed color at most h, where h is the largest originally free color. The remaining fixed colors form an upper tail, so the prior upper-tail theorem applies. This uses at most h-r released edges. The local rank-gain test can be much smaller; the example below needs only one.

## 3. A complementary exact obstruction: tree-supported contacts

Suppose the union of all nonloop free-color contacts, after contracting F and forgetting parallel realizations, is a tree J on r+1 vertices. Construct a bipartite graph H:

- its left vertices are the r free colors;
- its right vertices are the r edges of J;
- color d is adjacent to a tree edge a when at least one actual d-edge has contact type a.

### Theorem 3: exact matching reduction

Every completion projects to a perfect matching of H, and every perfect matching lifts to completions by independently choosing an actual edge realizing each matched contact. Under at-most-two-color exchanges, the connected components of the completion fiber project exactly to the connected components of the perfect-matching graph generated by alternating four-cycle flips. Each matching fiber itself is connected by one-color exchanges.

Proof. A completion uses r contact edges and must connect all r+1 contracted vertices. As every contact belongs to J, it must use every edge of J exactly once. Assigning its color to each tree edge is exactly a perfect matching. Conversely, any actual realization of each assigned type gives a spanning tree on the contracted components, hence a completion together with F.

Changing realizations while retaining a type cannot affect the contracted tree and is a legal single-color move. Two different perfect matchings differ on disjoint alternating cycles. If at most two colors change, the only possible nontrivial matching change is one alternating four-cycle. Conversely such a flip changes just its two colors, and any chosen actual realizations still give J. Therefore the projection and lifting preserve precisely the stated connected components. QED.

This is stronger than the unsigned statement that the free contact union is connected. The color-to-contact incidence graph is the information that controls exchange.

### Corollary 3.1: chordal-bipartite support is a second sufficient condition

If every cycle of H of length at least six has a chord, the completion fiber is two-color connected.

More generally, if every induced cycle of H has length at most 2k, it is connected by at-most-k-color moves, with single-realization changes allowed as well.

Proof. The symmetric difference of two perfect matchings splits into alternating cycles. A chord in an alternating even cycle splits it into two smaller even cycles. Exactly one of the two chord-plus-path cycles is alternating in the current matching: its path begins and ends with matching edges. Flip it; then the other chord-plus-path cycle is alternating. Flipping both has the same net effect as the original cycle. Induct until all flips are on induced cycles. Their number of colors is at most k. For the chordal-bipartite case k=2. QED.

The alternating-cycle argument is classical. Related primary literature is Martin Dyer, Mark Jerrum, and Haiko Müller, “On the switch Markov chain for perfect matchings,” https://arxiv.org/abs/1501.07725 . The same chord mechanism appears in Hidefumi Ohsugi and Takayuki Hibi, “Toric Ideals Generated by Quadratic Binomials,” Journal of Algebra 218 (1999), 509-527, https://page.mi.fu-berlin.de/~chaase/lehre/ohsugiHibiQuadratic.pdf . No rapid-mixing claim is made here.

If H itself is a chordless cycle of length 2k, it has exactly two perfect matchings. Changing between them requires all k colors at once. Thus arbitrary exchange width is genuinely necessary for general tree-supported colored graphic fibers. This last observation is not an arithmetic realization theorem for integer-difference colors; the arithmetic constraints may prevent a proposed support graph from occurring.

Together, Theorem 2 and Theorem 3 give two useful routes. Rank gain certifies that a small release creates enough flexibility. Matching support identifies a concrete circuit obstruction when the frozen components force a tree-shaped contact union.

## 4. Complete one-edge release audit of the arithmetic witness

The existing counterexample has labels 0,...,12, free colors 6,9,11, and fixed edges

    color 1: 5--6       color 2: 1--3       color 3: 8--11
    color 4: 3--7       color 5: 0--5       color 7: 4--11
    color 8: 1--9       color 10: 0--10     color 12: 0--12.

Its fixed components are

    A={0,5,6,10,12}, B={1,3,7,9}, C={2}, D={4,8,11}.

The free supports are

    color 6: AD,CD; color 9: AB,CD; color 11: AB,AD.

Their union is the tree B--A--D--C, and H is a chordless six-cycle. The two perfect matchings lift respectively to three and two actual completions. Theorem 3 explains the disconnected fiber without a state-space search.

A complete enumeration after releasing each possible one fixed edge gives:

    release 1: 54 completions, one component
    release 2: 56 completions, one component
    release 3: 24 completions, one component
    release 4: 20 completions, one component
    release 5: 54 completions, one component
    release 7: 14 completions, one component
    release 8: 20 completions, one component
    release 10: 16 completions, one component
    release 12: 5 completions, components of sizes 3 and 2.

Releasing color 12 cannot help because 0--12 is its only possible edge, so it is still forced.

The one-release rank certificate proves the cases 1,2,4,5,7,8 without enumeration. In all six, use the old ordering (6,9,11). Its ranks before release are (2,2,2), so only the first row is deficient. The released color becomes connected with rank 4, and color 6 gains rank to 3. Colors 9 and 11 already meet their required ranks 2 and 1.

Releases 3 and 10 are connected despite failing every full rank staircase. For release 3 the free-color ranks are (4,2,2,2); for release 10 they are (3,3,3,2). Thus staircase failure is not an obstruction theorem. These two cases are certified by the complete finite enumeration, not by Theorem 2.

### An explicit shortest three-move repair

Retain every fixed edge except 5--6. List states in color order (1,6,9,11):

    S0: (5--6, 2--8, 0--9, 0--11)
    S1: (1--2, 0--6, 0--9, 0--11)
    S2: (1--2, 0--6, 2--11, 1--12)
    S3: (5--6, 4--10, 2--11, 1--12).

Every listed state, together with the eight unlisted fixed edges, is a tree with all differences 1,...,12. The moves change color pairs (1,6), (9,11), (1,6). The endpoints belong to opposite components of the original fiber, and the released edge is restored.

No two-move repair can join these endpoints while keeping the other fixed edges. Without changing color 1 the original obstruction applies. If color 1 changes and returns in two moves, it consumes two of the four possible coordinate changes, leaving only two changes for the three differing original free colors. Therefore three moves are optimal.

## 5. Remaining gap and falsifiable targets

The rank-gain theorem is a concrete sufficient condition for bounded release. It does not prove that every arbitrary fixed-color graceful fiber becomes connected after releasing one edge. The example's success with eight choices is not evidence sufficient to assert that universal statement.

A targeted next test is the following bounded-release question: does every disconnected nonempty integer-difference completion fiber admit a release of one nonforced fixed edge whose enlarged fiber is connected? A negative search should retain all actual integer-difference edges, rather than substituting an abstract support graph. Theorem 3 supplies a disciplined way to search for larger alternating-cycle obstructions, while Theorem 2 rapidly certifies many positive releases.

The square-sum connection is structural: both path lifts and these exchanges succeed by releasing selected constraints while retaining the exact incidence data needed by later moves. Signed ports in the Pythagorean fan construction and color-contact incidence here cannot be dropped. Neither representation transfer preserves all original labels, all unused interval vertices, or a prescribed abstract tree shape.
