# Odd chains, square shells, and the tournament grid

Status: **PROVED elementary constructions and obstructions below; FINITE-EXACT tests.** No novelty claim, no theorem-ID reservation. This is a precise realization of the user's picture, not a claim that the user's intended curve was uniquely specified. No LRC or Collatz closure follows.

## Recovered definition, verified against current remote

The exact prior definition is in `05-knowledge/results/decoder_halving_20260925.md`, sections 1–3. GitHub fetch_file verified lines 1–124 at commit `aa977c095d7651e1f27c3b626fa61bed75cf4ae9` (blob `cb6891af3d0aa67bbc22bacba8c1b262e7bcb048`). It defines

    n = 2^k q = 2^k(2j+1), j,k >= 0.

The odd chain is q→q+2, and doubling raises k. This is NOT a new definition invented here. The same source already gives the regular odd cores H_q and the intrinsically halvable family E(2^k q)=H_q[TT_(2^k)]. The current note supplies an explicit shell traversal and obstructions around that existing grid. Other historical operation-grid prose was inspected only for provenance and is not a dependency.

Inheritance: closest proved mechanism is that recovered decoder; hostile is literal +2 off the odd row; corrected near miss is confusing address coverage with a continuous plane-filling curve; least-used sidecar is the scale k. Portfolio: anchor = exact odd-chain geometry; niche = tournament-family commuting square; wildcard = a full-plane square-shell spiral. Concept board: dyadic address / adjacency; odd square-shell count / coverage; tournament substitution / induced-core preservation; literal arithmetic / noncommutation; continuum limit / area.

## 1. The canonical grid is exact, but its horizontal increments are scaled

Define L(j,k)=2^k(2j+1). Unique factorization by the largest power of two dividing n gives a bijection N_0²→N_+. The inverse is

    k=v2(n), j=(n/2^k−1)/2.

Vertical neighbors satisfy L(j,k+1)=2L(j,k). Horizontal neighbors satisfy

    L(j+1,k)=L(j,k)+2^(k+1).

Thus horizontal movement is +2 **on the odd core**, and literal +2 on integer labels only on row k=0. Define A(n)=n+2^(v2(n)+1) and D(n)=2n. Then A and D commute, since A(Dn)=D(An)=2n+2^(k+2). Their address graph is exactly the quadrant square grid.

By contrast, S(n)=n+2 obeys DS=S²D, not DS=SD. If an injective coordinate map made S and D fixed translations by u and v, this identity would force u+v=2u+v, hence u=0, contradicting injectivity. More strongly, the undirected graph of literal +2 and doubling contains the 5-cycle

    1 — 3 — 6 — 4 — 2 — 1.

Every edge is a literal +2 or ×2 edge in one direction. Since the square grid is bipartite, there is no injective assignment of ALL positive integers to square-grid vertices that sends EVERY literal +2 and doubling edge to a unit grid edge. This obstruction is not just a flaw in the chosen dyadic coordinates. Directed-only positive A,D moves also cannot exhaust the quadrant in a single path: after leaving a coordinate axis it cannot return. A covering path must allow inverse moves or jumps.

## 2. Odd shell lengths produce a nearest-neighbor path through every address

For visit time t≥0 put r=floor(sqrt(t)), s=t−r², so 0≤s≤2r. First define

    (x,y)=(r,s)                 if 0≤s≤r,
    (x,y)=(2r−s,r)              if r≤s≤2r.

Set P(t)=(x,y) when r is even, and P(t)=(y,x) when r is odd. The repeated corner formula agrees. Shell r consists of all (j,k) with max(j,k)=r and contains exactly 2r+1 cells. Alternating shell direction joins their endpoints by one edge:

- even r ends at (0,r), next shell starts (0,r+1);
- odd r ends at (r,0), next shell starts (r+1,0).

Within a shell each move is a unit horizontal/vertical edge. The shells partition N_0², so P is a bijection, and its first (r+1)² values are exactly [0,r]². This proves the grid-filling path, not merely a finite sample.

Inverse: r=max(j,k); put (x,y)=(j,k) for even r and (k,j) for odd r; set s=y if x=r and s=2r−x otherwise. Then t=r²+s.

Now B(t)=L(P(t)) is a permutation of ALL positive integers. It begins

    1,2,6,3,5,10,20,12,4,8,24,40,56,28,14,7,9,18,36,72,144,...

Every transition is ×2, ÷2, or ±2^(k+1) in a fixed row. Crucially, these labels are not consecutive integers. Natural numerical order already fails unit-grid adjacency at 2→3: (0,1)→(1,0) has Manhattan distance 2.

The inverse visit time of integer n is rank(address(n)), explicitly given above. The first (r+1)² visited integers are exactly

    {2^k(2j+1): 0≤j,k≤r},

not the interval 1..(r+1)². Their maximum is 2^r(2r+1). Every integer appears, but size and visit time are very different scales. The odd shell count 1+3+...+(2r+1)=(r+1)² is the exact meeting of +2 shell-length growth and square-area growth. An odd shell is not itself an odd-core chain; these are two distinct uses of “odd.”

## 3. The recovered tournament family makes both axes substantive

Recover H_1 as one vertex. For odd q, choose S⊂V(H_q), |S|=(q−1)/2, add a,b, and orient a→b, a→S→b, b→(V\S)→a. Every old vertex gains one out-neighbor, and each new vertex has (q+1)/2 out-neighbors. Hence H_(q+2) is regular and contains H_q induced. A deterministic choice of S produces a fixed nested chain.

At cell (j,k) place

    F(j,k)=H_(2j+1)[TT_(2^k)].

Horizontal growth adds TWO whole transitive fibers, each of size 2^k, retaining all old arcs. Vertical growth substitutes TT_2 into every vertex. Associativity and TT_b[TT_2]=TT_(2b) give the commuting square, with compatible core/fiber labels:

    grow odd core, then double = double, then grow odd core.

Numerically it is n→n+2^(k+1) and n→2n, exactly the same A,D grid. This is a genuine family of tournaments with all positive orders and compatible graph operations, not a scalar analogy. Traversal P visits exactly one chosen tournament E(n) at each positive order n, not every tournament of that order.

For k≥1 the graph of two-vertex modules is the disjoint union of (2j+1) paths on 2^k vertices. Its unique perfect matching yields the canonical halving quotient F(j,k−1) up to isomorphism (prior source proof). Horizontal deletion, however, requires identifying the two newest core fibers. The old core labels or growth history are an essential sidecar; intrinsic vertical halving does not automatically make horizontal deletion intrinsic.

## 4. Quadrant, whole-plane tiling, and continuous curves are different claims

The half-open unit cells [j,j+1)×[k,k+1) tile the quadrant. The above path visits one center of every such cell exactly once. The path itself does not fill cell interiors.

For an actual whole-plane grid traversal, here is a separate explicit spiral with odd shell lengths. Put again r=floor(sqrt(t)), s=t−r² and m=floor(r/2). Define Z(t) by

    r even: (-m,m−s) for s≤r; (-m+s−r,−m) for s≥r;
    r odd:  (m+1,−m+s) for s≤r; (m+1−s+r,m+1) for s≥r.

Its first (r+1)² values are the nested square

    [-floor(r/2), ceil(r/2)]² ∩ Z².

Each new shell adds a connected L of 2r+1 cells on alternating sides; direct endpoints show the joins have unit length. These nested squares exhaust Z². An inverse is

    r=max(2 max(−x,−y,0), 2 max(x,y,0)−1), m=floor(r/2),
    r even: s=m−y if x=−m, otherwise s=r+x+m;
    r odd: s=y+m if x=m+1, otherwise s=r+m+1−x;
    t=r²+s.

Placing cell labels B(t) at position Z(t) gives a bijection between positive integers and full-plane unit cells whose scheduled successive labels follow the same arithmetic path B. This transport preserves the chosen path and bijectivity, but NOT all horizontal/vertical adjacency of the original dyadic grid. For an explicit hostile, the dyadic neighbors (1,1) and (2,1) occur at times 2 and 5; their plane images Z(2)=(1,1) and Z(5)=(−1,0) have Manhattan distance 3. The corner of the quadrant has degree 2 whereas every whole-plane grid vertex has degree 4, so no full grid-graph isomorphism exists either. The full transported adjacency requires retaining P/rank as a sidecar.

Finally, linearly joining a countable sequence of unit grid edges gives a countable union of line segments, of planar area zero; it cannot surject onto an open planar square. A continuous surjection from a compact interval also cannot cover the unbounded whole plane. A bounded-square Hilbert/Peano limit would require a separate compatible refinement construction, not merely this infinite grid traversal. Therefore the rigorous conclusion is: **a cell-filling nearest-neighbor grid path and a plane tiling indexed by integers**, not a continuous surjection from [0,1] onto R².

## 5. Loss ledger and reproducibility

- n↔(j,k): lossless arithmetic; natural numerical successor is not grid adjacency.
- (j,k)→odd core: loses valuation and therefore actual label increments; restore k.
- dyadic grid→B traversal: retains every vertex and chosen path, omits most grid edges; restore (j,k) or inverse rank.
- E(n)→n: loses tournament orientation/core choice; restore the chosen H_q family.
- E(n)→halving quotient: retains lower-level structure up to isomorphism, loses original vertex labels; restore fiber labels.
- quadrant path→whole-plane spiral: preserves visit time, coverage, and path adjacency; not a grid isomorphism.
- discrete path→drawn segments: does not preserve cell-area coverage; retain cells, and do not infer a continuous space-filling limit.

Run `python3 oddplane.py` and `python3 -O oddplane.py`. Both pass 560,200 explicit checks (checks remain active under -O), covering integer addresses through 100,000, 40,000 traversal positions, complete squares through shell 99, and the literal-operation odd-cycle hostile. Definitions, inverse formulas, and the all-scale proofs above supply what the finite tests cannot. `checks.json` records exact bounds and first labels; `oddplane.svg` is a 7×7 chart labelled by integers rather than visit times.

## 6. Positive continuous realization: a half-line surjects onto the plane

**PROVED, self-contained.** The preceding area-zero obstruction applies to the polygonal path, not to every continuous realization of the cell schedule. There exists a continuous surjection G:[0,∞)→R² whose image on [t,t+1] is exactly the closed unit cell centered at Z(t). Thus the discrete construction can be extended to a genuinely space-filling continuous map after adding an explicit cell-filling operation. This extra operation is substantive and is not implied by simply joining successive integer labels.

### 6.1 An explicit square-filling map, without a literature dependency

Let C⊂[0,1] be the ternary Cantor set. Every x∈C has a unique expansion

    x = sum_(i≥1) 2 b_i / 3^i, b_i∈{0,1},

where the all-0/2 expansion is used, including at endpoints. Define

    f(x) = (sum_(i≥1) b_(2i−1)/2^i, sum_(i≥1) b_(2i)/2^i).

This is onto [0,1]²: choose a binary expansion of each desired coordinate (using 0.111... for 1), interleave their digits, and use them as b_i. It is continuous. Indeed, if two Cantor expansions first differ at digit l, their numerical values differ by at least

    2/3^l − sum_(i>l) 2/3^i = 3^(−l).

Therefore |x−y|<3^(−m) forces their first m Cantor digits to agree, and each coordinate of f(x)−f(y) has absolute value at most 2^(−floor(m/2)). Also f(0)=(0,0) and f(1)=(1,1).

Extend f to F:[0,1]→[0,1]² by affine interpolation on each component (a,b) of [0,1]\C:

    F(x) = ((b−x)/(b−a)) f(a) + ((x−a)/(b−a)) f(b).

Keep F=f on C. Convexity keeps its image in the square, and F(C) already covers the square. For completeness, continuity at a point x∈C follows by a sequential argument. If y_n→x are in complement intervals (a_n,b_n) whose lengths tend to zero, both endpoints tend to x, so continuity of f and convexity imply F(y_n)→f(x). Otherwise, pass to a subsequence whose intervals have lengths bounded below. Only finitely many disjoint complement intervals have that length, so a further subsequence lies in one fixed interval; x must then be its endpoint, where the affine formula converges to f(x). Points in C and in interval interiors are already covered. Any hypothetical failure of continuity would have one of these subsequences, a contradiction. Thus F is continuous and onto the square, with the stated corner endpoints.

This is an elementary Cantor-digit/interpolation construction, not a novelty claim or an imported theorem. It supplies the square-filling existence proof directly, so no external reference is needed.

### 6.2 Fill the spiral cells, with connectors that glue

Write z_t=Z(t) for t∈N_0 and let

    Q_t = z_t + [−1/2,1/2]²,
    c_t = z_t − (1/2,1/2), d_t = z_t + (1/2,1/2),
    b_0 = z_0,
    b_t = (z_(t−1)+z_t)/2, t≥1.

Successive z_t are unit-grid neighbors, so b_t lies on the common edge of Q_(t−1),Q_t. For s∈[0,1], define G(t+s) by three stages:

    0≤s≤1/3:     (1−3s)b_t + 3s c_t;
    1/3≤s≤2/3:   c_t + F(3s−1);
    2/3≤s≤1:     (3−3s)d_t + (3s−2)b_(t+1).

At s=1/3 and 2/3 the corner endpoints of F give agreement. At s=1 the value is b_(t+1), exactly the value at s=0 for the next piece. The intervals [t,t+1] are locally finite on [0,∞), hence these continuous pieces glue to a continuous global map. No additional continuity condition “at infinity” is part of this domain.

The middle third covers all of Q_t; both connectors stay inside Q_t. Thus

    G([t,t+1]) = Q_t.

Since Z is a bijection N_0→Z², the union of the Q_t is R², proving surjectivity. Closed cells have disjoint interiors and shared boundaries, the usual geometric tiling convention; half-open cells give a disjoint ownership partition, but closed cells are needed for this continuous construction.

The construction is necessarily noninjective, and here this is visible directly: pick a point on the shared edge of Q_0 and Q_1 different from b_1. Surjectivity of each square-filling piece supplies distinct preimages inside [1/3,2/3] and [4/3,5/3]. This does not contradict the earlier obstruction: the domain is the noncompact half-line, and the map contains the nonpolygonal cell-filling stages. Reparameterization gives a continuous surjection [0,1)→R², but it has no continuous finite-valued extension at 1, since a continuous image of compact [0,1] is bounded.

### 6.3 Exactly what the continuous extension preserves

It preserves the discrete visit schedule, the integer label B(t) assigned to cell Q_t, and the consecutive-cell shared-edge transition. It does not preserve a one-to-one correspondence between integer labels and curve points: each time interval now contains uncountably many parameters, cells share boundaries, and the filler has repeated values. It also does not make ×2 a global geometric vertical translation in the plane. The continuous plane-surjection is therefore a positive completion of the user's geometric interpretation with an explicit added rule, rather than a silent strengthening of the original discrete-grid theorem.
