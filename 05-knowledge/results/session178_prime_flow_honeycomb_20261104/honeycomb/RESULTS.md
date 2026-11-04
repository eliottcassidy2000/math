# Two constructive honeycomb models for the ordinary Collatz map

Provenance: the 704-versus-694 obstruction is inherited from prior session work (session 167) and is only rechecked here. The dyadic-ray chart and explicit finite wire layouts are the constructive contributions developed in this session; the dilation inequality is a direct extension of the inherited counting argument. No historical novelty claim is made.

Use the ordinary map T(n)=n/2 for even n and T(n)=3n+1 for odd n. These constructions do not assume convergence. They concern the brick-wall realization of the honeycomb graph H: vertices are integer pairs (x,y), all horizontal unit edges are present, and the vertical edge out of (x,y) goes up when x+y is even and down when it is odd. This is combinatorially the honeycomb lattice, not the triangular lattice of hexagon centers.

## 1. An exact global arithmetic chart: local halving, explicit odd jumps

Write n=2^a q with q odd. Define

    Phi(n) = (a,q-1).

This is a bijection from the positive integers onto the vertices with x>=0 and y>=0 even. Its inverse is L(x,y)=2^x(y+1). The exact outgoing arrow detector is

    x>0: (x,y) -> (x-1,y),
    x=0: (0,y) -> (a,(3y+4)/2^a-1), a=v_2(3y+4).

Proof: the first rule removes one factor of two; the second factors 3(y+1)+1 into its odd and power-of-two parts. The destination row is again nonnegative and even. Thus Phi(T(n)) is exactly the indicated destination, without orbit computation or a convergence assumption.

Every even arrow is a nearest-neighbor honeycomb edge. Odd arrows are straight arithmetic jumps, not necessarily honeycomb edges or disjoint routed paths. This is an injective labeling with long edges, not a full nearest-neighbor embedding. It makes all halving chains literal horizontal rays. An odd step followed by all its halvings travels from one ray's endpoint to another ray and returns horizontally to that ray's endpoint. The resulting return map on the boundary x=0 is exactly the accelerated odd Collatz map.

The odd-jump test can also be written without an explicit valuation call: from (0,y) to (a,z), require a>=1, z>=0 even, and 2^a(z+1)=3y+4. These conditions uniquely determine the target. This is an arithmetic arrow detector, not a bounded-radius finite-state detector on unlabeled lattice geometry.

For y=0 mod 4, the jump has a>=2 and destination row z<=3y/4. For y=2 mod 4, a=1 and z=3y/2+1. In other words, the odd return dynamics consist of explicitly visible contracting and expanding boundary moves. This reformulation alone proves no net drift.

## 2. Explicit finite nearest-neighbor realizations after subdividing edges

Every finite rooted binary tree of height d has a vertex-injective topological embedding in H. Original tree vertices are distinguished branch vertices; each tree edge is replaced by a finite nearest-neighbor path. Paths meet only at their common tree endpoints. Extra host edges can be omitted: this is a subgraph realization, not an induced-subgraph claim.

Here is an explicit construction, sufficient to implement without a planar-embedding oracle. Embed the complete binary tree of height d and delete unused branches. Put the root at (1,0). For a vertex at depth t<d and position (x,2t), set D=2^(d-t). Give its left and right children positions (x-D,2t+2) and (x+D,2t+2). All branch x-coordinates are odd. For a child at (X,2t+2), route its edge from the parent horizontally at height 2t to (X-1,2t), then use

    (X-1,2t) -> (X-1,2t+1) -> (X,2t+1) -> (X,2t+2).

Both vertical edges exist by parity. The horizontal run has positive length because D>=2. The left and right paths from one parent have no common vertex except that parent. Distinct parents at the same depth have x-separation at least 4D, while their outgoing paths are contained in [x-D-1,x+D] at that level. Their intervals are disjoint. Different levels use disjoint open horizontal strips, meeting only at prescribed child/parent endpoints. This proves the required disjointness and injectivity.

A left edge has length D+4 and a right edge length D+2. The total number of inserted wire vertices for the complete height-d tree is

    (d+2) 2^(d+1) - 4.

Indeed, at level t there are 2^t parents, each contributing (D+3)+(D+1)=2D+4 interior vertices. Summing for t=0,...,d-1 gives the formula. The drawing uses d+1 branch levels and O(2^d) width, with O(d 2^d) total vertices.

Apply this to a finite inverse Collatz tree rooted at a nonperiodic point, including the inverse tree of 5. The two possible incoming children of n are 2n and, only when n=4 mod 6, (n-1)/3. Label the corresponding branches consistently and orient every wire toward the root. All original Collatz arrows are now detectable as directed paths between arithmetic-labeled branch vertices. Intermediate vertices are auxiliary wire sites and are not additional positive integers in the Collatz dynamics. If the chosen root belongs to a cycle, the cyclic edge must be removed or separately accounted for before calling the result a tree.

The distinction matters: this construction is fully nearest-neighbor along wires but does not make one Collatz step equal to one lattice edge. Also, the height-d layouts are not nested as d increases. Consequently, this construction does not by itself furnish a global infinite subdivision embedding by taking a limit.

## 3. The prior volume obstruction becomes a quantitative dilation constraint

A vertex-injective drawing of a depth-d rooted tree in H whose edge paths have length at most L places every tree vertex within lattice distance Ld of the root. Hence necessarily

    |B_tree(d)| <= 1 + 3(Ld)(Ld+1)/2.

This uses the honeycomb ball count |B_H(r)|=1+3r(r+1)/2. In particular, the 704 vertices in the depth-21 inverse ball of 5 do not fit inside the 694 vertices of a radius-21 honeycomb ball: L=1 fails, while this count does not obstruct L=2. The inverse-ball count is separately recomputed in the accompanying test rather than assumed from prior notes.

For complete binary trees the same inequality forces the maximum edge length somewhere in a depth-d realization to grow at least on the order of 2^(d/2)/d. Thus subdividing is a genuine geometric cost, not just a cosmetic correction. This binary-tree statement should not be transferred to Collatz inverse trees without proving an appropriate growth lower bound for those trees.

## Verification

`check_honeycomb.py` uses only Python's standard library. It checks 100,000 chart/conjugacy identities and even-edge adjacencies; complete binary embeddings through depth 10, including injectivity and pairwise interior disjointness; and independently computes both 704 and 694 counts. All checks passed. The proofs above, rather than the finite tests, support the general construction.

## What this adds and what it leaves open

The arithmetic chart is global and exact but contains long odd arrows. The wire construction is finite, injective, and nearest-neighbor after subdivision, with explicit coordinates and resource bounds. Neither is a solution to Collatz, a global one-step honeycomb embedding, or evidence of positive drift. A useful next question is how much the specific Collatz inverse trees permit reducing total wire length relative to the complete-binary construction, with branch labels retained as a locally checkable arithmetic certificate.
