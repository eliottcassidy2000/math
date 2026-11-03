# Protected sources and sparse interface escape

Status: proved structural claims and exact examples below, with independent finite controls. The exact state obtained by adding a vertex is generally not the original marked state, and no ternary-state closure or spectrum result is asserted. Source states and q2/b transfer matrices are from [the strong-return note](../strong_returns/NOTE.md).

## Protected-prefix and protected-module lemmas

For a marked tournament (X,p), append a nonempty atom by making old vertices beat all new vertices except that the atom's first vertex beats p; then mark any vertex of that new atom. The a-query is the one-vertex special case; the S163 q2 and b atoms use their final new vertex as port.

If X=B⊕Z with p∈Z, every such operation leaves the literal prefix B unchanged and yields B⊕(Z with the operation), including all subsequent a-queries and suffixes. Hence v(B⊕Z)=H(B)v(Z) for the observer v=(h,x,y), and this relation persists through every suffix. In particular C3⊕Z is permanently nonstrong and its full observer is divisible by 3. No q2/b/a sequence can restore strong connectivity or escape that 3-divisibility without changing the permitted attachment or the port policy.

More generally, take run-divisible blocks B_i with some set of preexisting exceptions S, as in `NOTE.md`, and any marked port p. Every appended vertex sees every B_i uniformly except possibly the FIRST appended vertex, which may see nonuniformly the unique B_i containing p. Since all subsequent ports are newly appended and outside the fixed blocks, no later vertex has a nonuniform interface to any B_i. Thus the interface graph acquires at most one new exception of degree one. If p is outside the blocks, it acquires none; if p is inside a block, the matching-based guaranteed exponent drops by at most one. The conclusion covers arbitrarily long suffixes and the two observer queries, not just counts at the original source.

For odd m, if the interface bound gives m^d | H of the source and its a- and aa-queries, then m^d divides each observer coordinate because 2 is invertible mod m^d. All moduli >1 dividing the full path count of a nonempty tournament are odd. Even-modulus care: generally only m^d/gcd(m^d,2) follows by division from the three count congruences, although even m cannot meet the nonempty tournament block hypothesis.

## Quadratic interface criterion for arbitrary modulus

Fix an old tournament A, a designated subset B of b vertices, an old physical port p, and every edge between a new vertex z and A\B. For i∈B put d_i=1 if z→i and d_i=0 if i→z. For each D⊆B let A_D denote the resulting tournament, retaining port p.

The exact count h(D)=H(A_D) is an integer multilinear polynomial of degree at most 2 in the Boolean variables d_i. Indeed, partition a Hamilton path at z into its left and right old-vertex paths. Their internal edges are fixed and only the incoming and outgoing edge at z depend on d, each affinely. More explicitly, with s_v and e_u the old Hamilton path counts starting and ending at the indicated vertex, and E_uv the number of ordered two-path covers of A whose first path ends at u and second starts at v,

 h(d)=∑_v d_v s_v + ∑_u (1−d_u)e_u + ∑_{u≠v}(1−d_u)d_v E_uv,

where d_v is also used for the fixed directions on A\B. Thus h is completely determined by the 1+b+binom(b,2) values h(∅), h({i}), h({i,j}). The coefficient of d_i is h({i})−h(∅); the coefficient of d_i d_j is h({i,j})−h({i})−h({j})+h(∅). This identity is over the integers and hence valid modulo every modulus, including composites.

If h(∅)≡0 mod m, SOME interface makes h nonzero mod m iff a singleton or pair does. The forward implication follows because if all singleton and pair values vanish, all polynomial coefficients vanish modulo m; the reverse is immediate. In particular, if B is an m-run-divisible module, the uniform baseline A_∅ satisfies the hypothesis.

This polynomial gives a finite complete decision for simultaneous escape of all prime divisors of any chosen m: compute the coefficients using O(b²) exact-count evaluations, then test the Boolean assignments for gcd(h(D),m)=1. A singleton or pair with unit count is a quick sufficient certificate; a unit need not occur among them even if a larger D gives one. Nor does escaping each prime separately imply that one D escapes them all. No general bound or uniform success for this search is asserted.

The whole marked observer has the same degree bound, not just h. For the fixed old port p, in each a- or aa-query, orientations of z to the query vertices are fixed independently of D. The same two-neighbor argument makes h, h_a, h_aa quadratic. Their observer differences are always even, so coefficientwise division by 2 gives integer quadratic polynomials for x and y as well. Consequently O(b²) literal observer evaluations determine every new physical state, and arbitrary target or unit conditions can be decided from these polynomials by a finite Boolean search. This changes the state; it does not say the original state is retained.

## The literal R4 return escapes with one edge and q2

Use the S163 strong ten-vertex returned source with rows

    R4=(996,993,994,999,975,912,800,576,128,256), port=9,
    v(R4)=(291,120,225).

Vertices B={0,1,2} form a directed triangle and a literal module; the port is outside. Hence the matching barrier with no exceptions applies to R4 and persists through every ordinary q2/b/a suffix. Every such descendant's count is divisible by 3, so none is either 3·291+2 or 3·291+4. This is only a statement about append descendants of this physical R4, not alternative seed replacements.

Add z, keeping original port 9, with z→i for i∈D⊆B and every other old vertex beating z. Retain R4 as an induced subgraph. Direct subset DP gives the following observer, depending only on |D|:

    |D|=0: (291,345,345)
    |D|=1: (1965,1103,1564)
    |D|=2: (3445,1781,2633)
    |D|=3: (4731,2379,3552).

Equivalently, with k=|D| and c=binom(k,2), the exact state is

    (291+1674k−194c, 345+758k−80c, 345+1219k−150c).

For |D|=1, 2 or 3, the new tournament is strong: z has edges both ways to the old strong tournament. Two reversed incidences escape immediately, with H=3445≡1 mod3; one gives H=1965≡0 but x=1103≡2 and y=1564≡1. Append literal q2, whose first output coordinate is 5h+8x+6y. The resulting strong physical tournament has exactly

    H=5·1965+8·1103+6·1564=28033≡1 (mod 3).

Thus one nonuniform edge into the protected triangle suffices to escape after one marked-state transfer, although it does not suffice immediately. The four possible counts also give a concrete failure of simultaneous escape modulo 15: 291,1965,3445,4731 are each divisible by 3 or 5. Yet the strong |D|=2 and |D|=3 cases have coprime counts 3445 and 4731; the former escapes 3 and the latter escapes 5. Thus gcd 1 across a family, and escape from each prime separately, do not force a single interface that escapes every prime at once. This is a restriction on changing exactly those three triangle edges: an unrestricted interface escapes 15. If z beats exactly old vertices 0 and 3 and loses to the other eight, retaining original port 9, literal DP gives the strong marked state (6011,2980,4531), and gcd(6011,15)=1. Vertex 3 is outside B, so this choice is not in the restricted family. Note that adding a uniform sink |D|=0 preserves H=291 but changes the marked state from (291,120,225) to (291,345,345); the original and new states must not be conflated.

`marked_escape_check.py` independently checks R4's tournament module, all eight literal interface states, strong connectivity, the appended q2 tournament, and its exact count. Run under both Python and Python -O.

## All one-vertex extensions of R4

The immediate ternary-offset target counts from R4 are 875 and 877. No tournament formed by adding a single vertex to this literal R4 has either count, regardless of which port is chosen afterwards. `search_r4_onevertex.py` reconstructs the count polynomial from 56 values and evaluates every one of the 1024 choices of edges between the new vertex and R4. It also independently runs literal endpoint-subset DP on all 1024 graphs and checks equality with the polynomial, in both normal and optimized Python. This is a finite exhaustion, not a claim about two or more new vertices or different source replacements.

The same finite exhaustion extends to all five marked sources U,U2,U4,R2,R4 in the return note. For each source, define its two exact targets 3v+t and 3v+2t, where t=(2,−1,0). No single-vertex extension achieves either target at any port. Target-COUNT hits do occur: for U the new vertex's out-neighbor mask is 1,33,117,125, giving counts 97,95,95,97 respectively; for U2 the only hit is mask 4, giving 289. None of those count hits has a target observer at any port. U4, R2, and R4 have no target-count hit at all. Here a mask bit i means the new vertex beats old vertex i. The script `five_source_branches_check.py` independently runs literal Hamilton endpoint DP on every orientation and cross-checks the count polynomial; it tests the full observer at every possible port for each target-count hit. Run it both normally and with `python -O`. This finite result is limited to induced extensions by one vertex, and does not rule out replacing a seed, adding two or more vertices, or changing the branch family.

An unrestricted but finite one-vertex search also cannot guarantee a unit for every modulus: the 1024 extensions have 454 distinct counts, all at least 291. Their product is a modulus for which none is a unit. More generally, for any finite family of counts all greater than 1, the product of those counts (or its radical) obstructs simultaneous coprimality. The m=15 obstruction above is much smaller but applies only to its explicitly restricted three-edge family.
