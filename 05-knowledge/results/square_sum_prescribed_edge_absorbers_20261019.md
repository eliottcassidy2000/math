# Square-sum absorbers with prescribed host edges

Date: 2026-10-19. Status: exact local identities and conditional packing statements; no claim of a new unrestricted interval Hamiltonicity theorem.

## 0. Existing interval theorem and the present scope

Robert Gerbicz's 2018 construction already supplies square-sum Hamiltonian paths on every interval `[1,n]` with `n >= 25`, and cycles for `n >= 32`. His own [OEIS comment](https://oeis.org/A090461) records these thresholds. His [original discussion](https://www.mersenneforum.org/node/17331/page2), especially the later explanations in posts 30 and 40–46, uses alternating affine lifts of parity-compatible sequence pairs, with modulus 49 for cycles and modulus 25 sufficient for paths. We do not claim novelty for these moduli or for unrestricted interval coverage.

This note records independently checked local insertion identities, their exact host-edge matching condition, and a matching-preserving block construction. Its possible use is in extensions which must preserve prescribed seed edges or an existing backbone. No theorem here asserts that the required host edges occur in an arbitrary backbone, or that Gerbicz's construction preserves an arbitrary prescribed seed.

The square-sum graph on a finite set V of positive integers joins distinct x,y when x+y is a square.

## 1. Two-vertex insertion

Suppose x,y are distinct positive integers with x+y=s². Let A,B >= 2 and C be positive integers satisfying

    A²+B²=C²+1.

Put

    u=A²(x+y)−x,       v=B²(x+y)−y.

Then

    u+x=(As)²,  x+y=s²,  y+v=(Bs)²,  u+v=(Cs)².

Consequently an existing edge u–v may be replaced by u–x–y–v. Both host vertices exceed max(x,y). They are distinct: if A=B their difference is y−x; if A≠B, the nonzero term (A²−B²)(x+y) has greater magnitude than x−y. If x,y <= L, then max(u,v) <= (2 max(A²,B²)−1)L.

The symmetric choice (A,B,C)=(5,5,7) gives

    u=24x+25y,       v=25x+24y,

with max(u,v) <= 49L. More generally, every positive solution C²−2A²=−1 with A>1 gives a symmetric insertion. Starting from (A,C)=(1,1), the recurrence

    (A,C) -> (3A+2C, 4A+3C)

produces (5,7),(29,41),.... In edge coordinates (square root s, oriented difference x−y), the host edge has coordinates (Cs,−(x−y)); its gap is unchanged.

## 2. A one-vertex insertion for every positive integer

For any positive integer x, choose p to be the least odd integer at least ceil(sqrt(2x))+2, and put

    q=(p²−2x−1)/2,       u=p²−x,       v=q²−x.

Then q,u,v are integers and

    u+x=p²,       v+x=q²,       u+v=(q+1)².

Moreover x<u<v<196x. Indeed p>=sqrt(2x)+2 and p<sqrt(2x)+4; hence q>p and q<4sqrt(2x)+15/2<14sqrt(x). Thus an existing edge u–v can absorb x by becoming u–x–v, with a linear numerical budget. Other sufficiently large odd p give further candidate hosts. Existence of the numbers u,v is not existence of their edge in a prescribed backbone.

## 3. Conditional simultaneous absorption

Let P be a square-sum path or cycle whose vertex set is V\O. Partition O into disjoint blocks, each either a singleton or a distinct square-sum pair. For each block specify candidate host edges of P whose replacement is furnished by Section 1 or 2. All generated host vertices must already be vertices of P and must belong to the prescribed ambient set V.

Form the bipartite graph with blocks on the left and eligible edges of P on the right. If every family of blocks has at least as many neighboring host edges as blocks, Hall's theorem supplies distinct host edges for all blocks. Replace these edges simultaneously by their insertion paths. The resulting path or cycle has vertex set exactly V. A path retains its original endpoints; a cycle retains its cyclic order on the old vertices.

Distinct host edges need not be vertex-disjoint. Even adjacent edges may both be replaced: their common old vertex still has one incidence from each replacement path. The requirements are edge distinctness, disjoint omitted blocks, and omission of every inserted vertex from P. No further conflict graph is needed for these fixed local insertion paths. If several candidate gadgets overlap in their inserted vertices because the block partition is not fixed, that more general selection problem needs additional conflict constraints.

If a protected collection F of old backbone edges must remain literal edges, simply remove F from the right side before testing Hall's condition. Every other old edge is retained, except the selected hosts. This is the explicit prescribed-seed version of the lemma.

The inequalities 49L<=N or 196L<=N ensure the numerical budget for the displayed gadgets on omitted labels at most L. They do not ensure backbone membership or Hall expansion.

## 4. An exact interval-pair packing

Let s>=3 be odd, S=s², and h=(S−1)/2. For each index 1<=x<=h pair x with S−x. These pairs partition the interval I=[1,S−1]. The symmetric host pair is

    u_x=25S−x,       v_x=24S+x.

The host pairs partition J=[24S+1,25S−1], and u_x+v_x=49S. Thus all insertions for I have disjoint hosts in J.

Put q=14s. Whenever indices stay in range, the additional edges

    v_x–u_(x−q−1),       v_x–u_(x+q−1)

have square sums (7s+1)² and (7s−1)² respectively. The first type alone joins each residue class modulo q+1 into an alternating path, so J has a path cover by at most min(h,q+1) paths, every host pair edge included. Absorbing all pairs gives a path cover of I union J with the same number of components.

There is a sharper exact cycle on every window of 2q consecutive indices. For a window W=[b,b+2q−1], let

    F(x)=b+((x−b+q−1) mod 2q).

Then F(x)−x is either q−1 or −q−1, so v_x–u_F(x) is one of the preceding square edges. Since q is even,

    gcd(q−1,2q)=1,

and F is a single cycle on W. Reading u_x,v_x and then the next index F(x) gives a square-sum cycle of length 4q on the host vertices in W. Replacing every host edge gives a cycle of length 8q on these hosts and their omitted pairs. In the original variables this is 56s host vertices and 112s total vertices. A full window fits when s>=57. Any subset of its paired omissions may instead be inserted, retaining all the other host matching edges.

These cycles leave the middle interval [S,24S] untouched. They do not themselves prove a spanning cycle on a consecutive interval, nor supply a universal merger of the different windows.

## 5. A connector-energy obstruction to small splices

Every oriented connector v_x–u_y with square root 7s+k satisfies

    y−x=−qk−k²,       q=14s.

An oriented closed index cycle therefore satisfies

    q sum(k)+sum(k²)=0.

If all k are ±1, its length is a multiple of 2q: all steps are q−1 modulo 2q, and that residue is coprime to 2q. More generally, two equal-size connector sets with the same total head and tail indices satisfy

    q Delta(sum(k))+Delta(sum(k²))=0.

If the absolute energy change is less than q, it must vanish. If the old set has only |k|=1 and all new connectors have k!=0, vanishing energy forces every new |k|=1. Thus a bounded same-orientation swap with this energy bound cannot introduce higher square-root offsets. This is an obstruction for the stated local move class, not for arbitrary rearrangements.

A two-edge splice reversing the orientation of one component also has a precise condition. Delete v_x–u_a and v_z–u_b, having square roots 7s+k_1 and 7s+k_2, and reconnect v_x–v_z and u_a–u_b. If their new square roots are r,t, then

    r²+t²=(7s+k_1)²+(7s+k_2)²,
    r²=48S+x+z,       t²=50S−a−b.

So such a splice needs a compatible sum-of-two-squares representation, as well as the index equations. For k_1=k_2=−1, the equal-root representation r=t=7s−1 requires x+z=S−14s+1; since x,z<=h, both indices must lie within 14s−2 of h. This explains why that immediate reversal splice is localized near the top of the index range.

## 6. Reproduction and scope

From the repository root run `python 05-knowledge/experiments/verify_prescribed_edge_absorbers_20261019.py` and repeat with `python -O`. The standard-library-only verifier checks the Pell identities, all singleton absorbers in a substantial integer range, exact host interval partitions, rotation cycles and their absorbed forms, simultaneous partial absorption with protected matching edges, and connector-energy identities. It uses explicit exceptions rather than removable assertions. These controls supplement the proofs; they are not a finite-search substitute for Hall expansion or universal window merging.

The remaining useful question is whether one can arrange eligible hosts while retaining designated seed or matching edges, perhaps within Gerbicz's known interval recursion. Neither an arbitrary prescribed-edge extension theorem nor an equivalence with Collatz or graceful tree labeling is asserted.

## 7. What an affine-lift recursion must certify to preserve edges

Here is a conditional structural lemma about the stated affine-lift architecture; it does not require importing or reproducing Gerbicz's glue tables.

Suppose a construction joins intact blocks obtained from input paths by

    a_j -> 49a_j + (−1)^j c,       −24<=c<=24,

possibly reversing whole blocks, together with singleton connectors at most 24. Suppose all input-path endpoints are at most K. Then every endpoint of an inter-block or singleton glue edge is at most B=49K+24. Every remaining edge is internal to a lifted block, and its sum is 49 times the corresponding input-edge sum. Therefore:

* Any output edge with an endpoint greater than B has square-root divisible by 7.
* If such a fixed edge has root not divisible by 7, it cannot appear in any non-base output of this architecture.
* A fixed edge whose root has finite 7-adic valuation may be traced through the uniquely determined centered residue decomposition as long as it is internal. After at most that many descents it must become an edge in the bounded glue region, or belong to an actual terminal base output. Merely adding it to finitely many base certificates does not make it persist at arbitrary recursion depth.

For the endpoint normalization {1,3,8} described in Gerbicz's discussion, K=8 and B=416. The edge 482–3362 has root 62, not divisible by 7, and lies outside this glue region. Thus the unchanged bounded-endpoint lift cannot preserve this edge in its non-base outputs. Its symmetric absorber host is 95618–92738. The unique centered decomposition is

    95618=49·1951+19,       92738=49·1893−19.

The descended edge 1951–1893 again has root 62 and lies outside the glue region. Hence this fixed host cannot persist once the relevant input is itself a non-base output. This is a structural obstruction to this proposed preservation argument, not to Gerbicz's Hamiltonicity theorem or to a modified recursion with additional states.

A valid protected-edge induction can instead be certified by explicitly listing, for each recursion state, its protected input edges and protected output edges; showing that every protected output edge is either a retained internal lift of a protected edge in whichever input path is selected, or a required glue edge; and checking compatible protected base cases. A finite collection of host edges with no such state transition proof is insufficient. Enlarged endpoint sets, extra protected glue states, or scale-dependent host choices remain possible.

The exact receipt is [prescribed_edge_absorber_controls.json](prescribed_edge_absorber_controls.json).
