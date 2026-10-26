# Square-sum arrival anatomy and endpoint-aware Hamiltonian kernels

Session 32, 2026-10-23. Let G_n be the **simple undirected** graph on {1,...,n}, with x adjacent to y iff x != y and x+y is a square.

## Status and inheritance

**PROVED (elementary):** arrival-degree formula; global classification of degree-two births; component/cycle-rank update; connectivity threshold; endpoint-aware chain kernel; marked-ear/path-cover correspondence. These are elementary graph/arithmetic facts, not priority claims.

**INHERITED PROVED, rechecked here:** the small obstructions, unique path at15, and endpoint-aware obstruction at24. In particular, the degree-two endpoint error in earlier research remains retracted: degree2 alone does not force both edges of a Hamiltonian path. A vertex may be an endpoint.

**FINITE-EXACT:** independent direct vertex-order DFS and endpoint-marked chain-kernel edge-subset enumeration agree for every 1<=n<=25. The script outputs all paths and endpoints through25, the displayed cut certificates for18–22, and arrival data through36. Exact checks are validation, not proofs of an infinite range.

**CITED:** existence of Hamiltonian paths for every n>=25 and cycles for every n>=32 is Gerbicz's result. This session does not independently audit his large finite base/glue tables.

Inheritance sources in the repository:

- `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`: audited degree/leaf law and first path.
- `05-knowledge/results/collatz_mod6_20260922_w6_summand_closure_square_filter.md`, §§2–5: arrival degrees, component founders, connectivity and leaf windows.
- `05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md`, §§4,6,7: corrected cut certificates, explicit25-fold family, and Gerbicz citation scope.
- `05-knowledge/results/decoder_prime_square_20260925.md`, §§4–5: endpoint-aware obstruction at24 and compulsory ear at25.
- `05-knowledge/results/crossroads_family_20260926_squares.md`, §§1–2: complete unicyclic criterion, nearby endpoint census, and warning against endpoint-blind suppression.

Closest mechanism: square reflections x -> s²-x. Canonical hostile: the valid G23 path ending at22, although22 has graph-degree2. Corrected near miss: forcing both edges at every graph-degree2 vertex. Useful sidecar: the endpoint marks retained by the kernel and the marked non-square edge retained by ear suppression.

Concept board: arrival interval; component/cycle rank; endpoint obligations; marked chain kernels; parity-compatible dilation. The first four interact exactly below; the last is the extra ingredient in the cited all-order construction, not an inference from small graphs.

## 1. Exactly which edges arrive

The new neighbors of n are s²-n with

    1 <= s²-n <= n-1,

or equivalently n < s² <= 2n-1. Therefore

    d_new(n) = floor(sqrt(2n-1)) - floor(sqrt(n)).       (1)

The right endpoint 2n-1 excludes the forbidden loop n+n. The left endpoint is open: zero is not a vertex. Equivalently floor(sqrt(2n))-floor(sqrt(n)) minus the indicator that2n is square.

For an old vertex x<=n the complete degree is

    deg_Gn(x) = floor(sqrt(x+n)) - floor(sqrt(x))
                - indicator(2x is a square).          (2)

Arrival degree is not the later graph degree. In particular a vertex born with one edge can gain more later.

### All degree-two arrivals, not just the first36

Put r=floor(sqrt(n)). Equation(1) equals2 iff

    max(r², ceil(((r+2)²+1)/2)) <= n
       <= min((r+1)²-1, floor((r+3)²/2)).              (3)

For r>=8, 2r²-1 >= (r+3)², because r²-6r-10>=0; hence every such n already has arrival degree>=3. Evaluating r=1,...,7 in(3) gives exactly

    [13,15] union [19,32] union [36,40] union [49,50]. (4)

The nonempty r-bands are r=3:[13,15], r=4:[19,24], r=5:[25,32], r=6:[36,40], r=7:[49,50]. Thus the requested list through36 is13,14,15;19,...,32;36. At33,34,35 the arrival degree is3;36 drops back to2 because its lower square boundary advances from25 to36 before the next upper square arrives.

For comparison, degree0 births are exactly1,2,4, and degree1 births are exactly3,5,...,12,16,17,18. Both are inherited classifications. At25 the window replaces square25 by square49 without changing its size: n24 sees25,36, whereas n25 sees36,49. So25 is not an arrival-degree threshold.

## 2. Components versus cycles: why13,14,15 do different things

Let r_n be the number of old components touched by the new neighbors of n, c_n the component count, and beta_n=|E(G_n)|-n+c_n the cycle rank. Then

    c_n = c_(n-1)+1-r_n,
    beta_n = beta_(n-1)+d_new(n)-r_n.                 (5)

Proof: adjoin an isolated vertex, then connect it to the r_n distinct old components. Exactly r_n component mergers occur; all d_new edges count toward |E|. This is an exact update, not an asymptotic heuristic.

At12 the components really are the paths

    2–7–9;
    12–4–5–11;
    10–6–3–1–8.

At13 the new edges are13–3 and13–12, so the last two old components merge. The merged component is NOT a path: vertex3 now has degree3. At14 the edges14–2 and14–11 merge the final two components. G14 is a tree with three leaves8,9,10, not a Hamiltonian path. Thus the component counts at12,13,14 are3,2,1. The two-component stage occurs only at13.

At15 the edges15–1 and15–10 touch the same existing component. Cycle rank rises from0 to1. The resulting unique cycle is1–3–6–10–15–1, with pendant paths at adjacent cycle vertices1 and3. Deleting their common edge1–3 gives the unique Hamiltonian path. Its unused **square sum** is4; vertex4 itself is visited and both its incident edges are used.

Connectivity is permanent from14. For n>=6, the first square above n is at most(sqrt(n)+1)²<2n, since n>3+2sqrt(2); therefore it gives an older neighbor. At5 there is also the older neighbor4. No later isolated vertex is born, so induction from connected G14 proves connectedness for every n>=14. The finite earlier edge lists show G_n disconnected for2<=n<=13. G1 is the trivial connected exception.

## 3. What a degree-two arrival actually preserves

Let H be any graph and adjoin v with neighbors a,b. A Hamiltonian path of H+v is of exactly one of two types:

1. v is an endpoint. Deleting v leaves a Hamiltonian path of H whose endpoint a or b was extended.
2. v is interior. Deleting v leaves exactly two disjoint spanning paths of H, with a and b on different paths and exposed as their joining endpoints. Singleton paths are allowed.

Conversely both constructions yield Hamiltonian paths, and deletion recovers the data uniquely. This is an exact marked two-path-cover correspondence. The state lost by mere path counts is the set of usable endpoints, and then the joining ports of spanning two-path covers.

If H has no Hamiltonian path, type1 is impossible and both edges a–v–b are compulsory. Suppressing v is then a bijection with Hamiltonian paths of H plus a **marked edge**a–b that use the marked edge. The marked edge can be non-square; deleting its mark or treating the quotient as another square-sum graph is invalid. If a–b already exists, retain a distinct marked parallel edge as the record of the subdivision.

Applications:

- At15, the compulsory ear is1–15–10.
- At16 and17, each new vertex merely extends the unique available endpoint, giving one path at each size.
- At18, the new leaf attaches at7, which is interior in the unique old path. It cannot be appended; three leaves give a certificate of failure.
- At23, the compulsory ear is2–23–13. It joins two spanning path pieces that G22 cannot join directly.
- At24, the old G23 endpoints are18 and one of2,9,22. Neither new neighbor1,12 is ever an old endpoint, so even though G23 is traceable, endpoint extension is unavailable. Interior insertion fails for the independent reason in§5.
- At25, the compulsory ear is11–25–24. It supplies the missing connection to escape the proper cycle obstruction in G24. This uses the new square49 for the first time, since distinct positive summands<=24 sum to at most47.

This explains why equal arrival degrees do not imply equal Hamiltonian behavior. Connectivity depends on which components are touched; a spanning path depends on exact exposed ports and compatibility with all forced interiors.

## 4. Endpoint-aware chain kernel: an exact finite reduction

Fix distinct candidate endpoints s,t in a finite connected graph G with at least2 vertices. Let B be {s,t} together with all vertices whose G-degree differs from2. Replace each maximal chain whose internal vertices lie outside B by one distinguishable multigraph edge, retaining its full vertex list. Parallel edges and loops must not be identified or discarded.

Every internal vertex of such a chain is a certified nonendpoint of degree2, so both incident edges, and hence the entire chain, are forced. Therefore:

- a kernel edge representing a chain of length>=2 original edges is mandatory;
- an edge representing one original edge is optional;
- selected kernel degrees must be1 at s,t and2 at all other vertices of B;
- the selected expanded subgraph must be connected.

These conditions are necessary and sufficient. Necessity is forced-chain propagation. Conversely the expansion is a connected spanning graph with two degree1 vertices and all remaining degrees2, hence a spanning path. A mandatory kernel loop creates a cycle and cannot satisfy the connected spanning-path conditions. There is no loss of information because the chain records specify a unique expansion.

Importantly, endpoint candidates of original degree2 are included in B, splitting their chains. Suppressing them as unmarked interiors reproduces the old invalid argument. A pure cycle causes no exception: its chosen endpoints supply two boundary vertices. Disconnected graphs are rejected separately.

`verify_square_sum.py` uses this theorem as one enumerator and a separate direct vertex-order DFS as the other. Their complete path sets agree, not merely their total counts, for every n<=25.

## 5. Human-checkable failures and successes

The general cut obstruction is

    Hamiltonian path in G implies c(G-S)<=|S|+1.

Deleting |S| vertices from a path leaves at most|S|+1 pieces, while every remaining component must contain a piece. The certificates are:

- n13: disconnected.
- n14: the tree has three leaves8,9,10.
- n18,n19: deleting S={7} leaves3 components.
- n20,n21: deleting S={5,7} leaves4 components.
- n22: deleting S={5,7,14} leaves5 components.

The exact components are included in `certificates.json`; no exhaustive Hamiltonian search is needed for these failures.

### The endpoint-aware obstruction at24 (inherited proof)

Since18 is a leaf, one endpoint is18. Write e for the other.

If e!=9, then9 is interior with neighbors7,16. Its edge to7 and the leaf edge7–18 saturate7, excluding2–7. If e is neither11 nor22, both are interior and force14–11 and14–22, saturating14 and excluding2–14. Since2 has exactly the neighbors7,14,23, these implications prove

    e in {2,9,11,22}.

Vertex20 is therefore interior and forces5–20. If e!=11, interior11 also forces5–11. If e=11, the saturated7 forces interior2 to use2–14 and2–23; interior22 forces14–22, so14–11 is excluded and endpoint11 must use11–5. Thus in every remaining case5 uses20 and11, excluding4–5.

All vertices of the proper cycle

    1–8–17–19–6–10–15–21–4–12–24–1

are interior by the endpoint gate. Degree-two vertices8,17,19,10,21,24 force all its edges except4–12. Interior4 has neighbors5,12,21, and4–5 is already excluded, so4–12 is forced too. A path cannot contain that proper cycle. This is a structural proof, not a census shortcut.

At25, vertex25 has neighbors11,24, so the exact marked-ear bijection of§3 applies. An explicit path is

    18,7,9,16,20,5,11,25,24,12,4,21,15,10,6,19,17,8,1,3,13,23,2,14,22.

### Exact small counts

Count paths up to reversal; retain the trivial singleton at1.

- n2,...,14:0.
- n15,16,17:1 each, with endpoint pairs{8,9},{8,16},{16,17}.
- n18,...,22:0.
- n23:3, with endpoint pairs{18,2},{18,9},{18,22}.
- n24:0.
- n25:10. Each has18 as one endpoint; the other is exactly one of2,3,4,8,9,10,11,13,22,23, one path per choice.

The n23 path ending at22 is the crucial hostile to unmarked degree-two forcing:

    18,7,9,16,20,5,11,14,2,23,13,12,4,21,15,10,6,19,17,8,1,3,22.

## 6. All orders: what is cited, and what is not proved by the small cases

Robert Gerbicz's own January21,2018 contribution to [OEIS A090461](https://oeis.org/A090461) records the all-n result: paths exist for every n>=25, and cycles for every n>=32. The primary [Mersenneforum thread](https://www.mersenneforum.org/node/17331/page2) explains the extra invariant: pairs of consecutive-length chains in which shared values occupy positions of equal parity, followed by square-scaled, alternating-shifted blocks and finite gluing tables. This is stronger data than merely one path at each size.

The OEIS page was read directly this session. Direct opening of the forum returned an access error, but its indexed primary author text was available through targeted search, including the parity-preservation explanation and distinction between moduli25 for paths and49 for cycles. The existing repository note had separately read the archived proof description. This is a **CITED construction**, not a fresh verification of its full base/glue tables. Finite path computations here do not prove the infinite statement.

The ordinary interval-existence problem is therefore not an open conjecture. Open research must ask for an additional invariant, prescribed endpoints/edges, a simpler independent induction, or a different graph family. Nothing above implies a Collatz theorem or a repeated three-component phase at larger scales.

## Reproduction and checks

    cd square-sum  # from the published packet root
    python3 verify_square_sum.py > certificates.json
    python3 -O verify_square_sum.py > certificates.optimized.json
    cmp certificates.json certificates.optimized.json

The script uses explicit checks, not Python assertions. Universe: all n<=25 for complete path-set agreement, n<=36 for the arrival ledger, and n<=100000 only as a supplementary numerical control of the independently proved degree-two classification. No external dependencies, Git mutation, or reliance on floating point.

Useful next tests, without confusing them with the result: endpoint-conditioned transfer across the user's25-vertex blocks; all-n preservation of a chosen pair of ports; whether a requested boundary invariant survives Gerbicz's gluing. These are stronger problems than ordinary path existence.
