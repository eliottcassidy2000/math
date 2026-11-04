# Square-sum arrival and deletion response through component ports

**Session 190, 4 November 2026. Status: PROVED elementary transfer, response, and horizon statements; FINITE-EXACT bounded square-sum data. No literature-priority or all-order endpoint theorem is claimed.**

## 1. What is inherited, and what is added

`../session190-planar/PROOF.md`, Sections 7–8, points from the negative hypotraceability bridge to endpoint/deletion response. Its Section 3 preserves component connectivity rather than only degrees. `../session32-square-sum/README.md`, Sections 3–4, already proves the degree-two marked-ear / two-path-cover correspondence and endpoint-aware chain suppression. We do not claim either again as new. The repository's THM-4510 reflection description supplies the arrival neighbors but not their compatible path connections.

This packet adds a closed transfer on all component-port states, its exact deletion-card update, a sharp finite-horizon truncation, and an explicitly invertible three-query response for four distinct ports. The finite experiment computes complete endpoint-count matrices, not just existence or scalar path counts. No planarity or noncrossing restriction is imposed on square-sum states.

## 2. The closed state and its exact transfer

Let G be a finite simple graph. A spanning linear forest is a spanning edge subset with no cycle and maximum selected degree at most two. Isolated vertices are allowed. Represent each positive-length path component by its unordered pair of endpoint labels (a,b), and a singleton v by (v,v). The state P is the unordered collection of these pairs, with a canonical sorted encoding. Different components have disjoint endpoint labels. Labels absent from P have selected degree two. The pair (v,v) is a singleton marker, **not a graph loop or two distinguishable stubs**.

Let F_G(P) count edge subsets realizing P, with components and paths unordered and paths unoriented. All component numbers are retained. In particular, the one-component entries are the full Hamiltonian endpoint matrix H_G(a,b); H_G(v,v)=1 only when G is the singleton v.

Adjoin z, whose old neighbors are S. For an endpoint vertex a in a component pair, write mate_P(a) for its other endpoint; mate_P(a)=a for a singleton. The complete transfer T_(z,S) consists of:

1. **No new edge.** Adjoin (z,z) to P. Component number increases by one.
2. **One new edge za.** Choose an endpoint vertex a in S, remove its pair, and insert (mate_P(a),z). Component number is unchanged.
3. **Two new edges za,zb.** Choose an unordered pair of endpoint vertices a,b in S belonging to distinct components. Remove their two pairs and insert (mate_P(a),mate_P(b)). Component number decreases by one.

Each choice contributes its old coefficient once to its target state. Distinct choices that have the same target are added. A singleton supplies one possible endpoint vertex, not two choices.

**Theorem 1.** F_(G+z)=T_(z,S) F_G.

**Proof.** In a linear forest the new vertex has selected degree zero, one, or two. At degree zero it is a singleton. At degree one, its old neighbor must have old selected degree at most one, precisely an encoded endpoint vertex. At degree two, both old neighbors must be encoded endpoints. If they lie in the same old component the two new edges close a cycle, which is forbidden; distinct components join to one path. These are exactly the three operations. Removing z and its incident selected edges recovers the unique previous forest and chosen neighbor set. Thus there is no orientation, component-order, or singleton multiplicity factor. ∎

This is a finite exact state description, not a bounded-size or polynomial-time claim. Keeping only one- and two-component entries is sufficient for **one specified arrival**, but is not a closed update for an unlimited succession of arrivals: two-component entries can receive contributions from three-component entries, and so on.

### Explicit endpoint formulas

For old s and new z,

    H_(G+z)(s,z) = sum_(a in S) H_G(s,a).

For distinct old s,t, define J_G^(s,t)(a,b) as the sum of F_G(P) over two-component states P in which a,b are exposed on different components and the two remaining endpoints are s,t. Singletons are permitted, using mate_P(a)=a. Then

    H_(G+z)(s,t) = sum_({a,b} subset S) J_G^(s,t)(a,b).

These equations distinguish endpoint extension from an interior join. If e_G(a) counts Hamiltonian paths having endpoint a, and j_G(a,b) counts compatible two-component covers with joining vertices a,b, then

    h(G+z) = sum_(a in S) e_G(a) + sum_({a,b} subset S) j_G(a,b).

For a singleton G, e_G(a)=1, not 2. The old graph is nonempty in these formulas; the empty initial graph is handled by the singleton transition.

## 3. A complete endpoint matrix still loses necessary information

Take the two labeled path graphs

    A: 0–1–2–3–4,
    B: 0–1–3–2–4.

Both have exactly the same Hamiltonian endpoint matrix: its only nonzero entry is H(0,4)=1. Add vertex z=5 with neighbors {1,2}. Neither old path can be extended at an endpoint. In A, delete edge 1–2 and insert 1–5–2, giving the Hamiltonian path

    0–1–5–2–3–4.

In B, vertices 0 and 4 are forced endpoints. Vertex 1 must use its leaf edge to 0 and vertex 2 its leaf edge to 4. If z is interior, its two incident edges saturate both 1 and 2, leaving vertex 3 disconnected from the path. Hence B+z is nontraceable.

Thus **even the full old Hamiltonian endpoint matrix fails to determine whether the next graph is traceable**. This is stronger than scalar-count insufficiency. Both old graphs are connected and traceable; neither uses an endpoint-at-degree-two error.

## 4. Four distinct ports: three probes recover all three pairings

Fix four distinct vertices a,b,c,d of G. Let

    x1 = F_G({(a,b),(c,d)}),
    x2 = F_G({(a,c),(b,d)}),
    x3 = F_G({(a,d),(b,c)}).

Every other vertex, if present, is internal in one of the two components. These states have no singleton component. Define R_ab as the number of Hamiltonian paths with endpoints c,d after adjoining a new vertex adjacent **only** to a,b. Define R_ac and R_ad by the analogous complementary endpoint prescriptions.

**Theorem 2.**

    R_ab = x2+x3,
    R_ac = x1+x3,
    R_ad = x1+x2.

Consequently

    x1 = (R_ac+R_ad-R_ab)/2,
    x2 = (R_ab+R_ad-R_ac)/2,
    x3 = (R_ab+R_ac-R_ad)/2.

**Proof.** The new vertex cannot be an endpoint because the prescribed endpoints are two distinct old vertices. Its two neighbors must belong to different components after deleting it. Pairing ab|cd is excluded by the ab probe, while the other two pairings each have exactly one valid join. Apply the same argument to the other probes. ∎

The response matrix is the adjacency matrix of K_3, with determinant 2, rational rank 3, and rank 2 over F_2. Its modulo-two kernel is generated by (1,1,1). Every response triple has even sum and satisfies the triangle inequalities, since its inverse has nonnegative integer coordinates. Conversely every nonnegative integral response triple with even sum and these inequalities has a nonnegative integral inverse **as a formal coefficient vector**, not necessarily as one graph profile.

Each basis vector is physically realizable by a graph consisting of the corresponding two disjoint edges. Moreover G=K_4 realizes x=(1,1,1), producing the even nonzero response (2,2,2). This last physical example concerns the selected four-distinct-port sector; K_4 also has other forest states. No claim is made that its entire forest profile has only those three entries.

This is the unrestricted **two-path-cover / one-new-vertex** response. It is not the planar cycle-gluing matrix from the earlier packet. The latter first loses parity injectivity at six planar ports; this observer loses it at four unrestricted ports. Both statements retain named pairing information and concern the entire specified response operator, not the false claim that parity decides existence.

## 5. Exact deletion-card update

The endpoint-pair state does not expose the internal position of a deleted old vertex. Therefore a deletion card is maintained separately, not obtained by erasing its label from the undeleted state.

Maintain the deck consisting of F_G and F_(G-v) for each old v. Under the same arrival z with neighbors S, the exact update is

    undeleted:       F_(G+z)     = T_(z,S) F_G,
    old deletion v:  F_(G+z-v)   = T_(z,S minus {v}) F_(G-v),
    new deletion z:  F_(G+z-z)   = F_G.

**Theorem 3.** These equations update the complete one-vertex deletion deck without losing endpoint or component-pairing information.

**Proof.** Deleting an old vertex commutes with adjoining z after its neighbor set is restricted. Theorem 1 applies to the resulting graph. Deleting z recovers G exactly. ∎

For any fixed old deletion set D the same equation holds with S minus D and F_(G-D). This is a graph/deck identity, **not** a claimed linear erasure map on an aggregated undeleted port vector. Extracting the one-component entries gives every deletion-card endpoint matrix and hence every traceability decision. Empty/singleton card conventions must still be specified when relevant; all reported square-sum deletion cards have at least 23 vertices.

## 6. A sharp finite-horizon rule and sound completion pruning

**Theorem 4.** If exactly r vertices remain to be added, a current forest with more than r+1 components cannot become a Hamiltonian path. During a target computation, at a stage with t additions remaining, discarding states with more than t+1 components is exact for the final Hamiltonian endpoint matrix.

**Proof.** Each added vertex lowers component number by at most one, as Theorem 1 shows. Reaching one component in t steps requires at most t+1 current components. Equivalently, deleting t vertices from a final path leaves at most t+1 nonempty path pieces. ∎

**Strictness, with precisely r additions.** For every r>=1, take r+1 disjoint edges u_i–v_i, i=0,...,r. Add exactly r new vertices z_i, i=1,...,r, with neighbor sets {v_(i-1),u_i}. The final graph is the single alternating concatenation of the r+1 old edges and the r joining vertices, hence has a unique Hamiltonian path. Every old spanning forest has at least r+1 components. Thus discarding that component level would erase every witness before the first arrival. This proves worst-case necessity of level r+1 for an r-addition horizon, not a claim that every fixed square-sum instance needs it, nor a general minimal-state theorem for arbitrary nonlinear encodings.

For a known target graph and vertex order, call a current endpoint *retired* if it has no neighbor among vertices still to arrive. The implementation also discards:

- a retired singleton, when the target has more than one vertex;
- a state with more than two retired degree-one endpoints;
- a component both of whose endpoints are retired while any target vertex remains to arrive.

These are sound: retired singleton vertices can never acquire an edge; retired degree-one vertices must remain final endpoints; and a component with no available endpoint can never join the unprocessed vertices. These tests depend only on the recorded state and target graph, so all forests aggregated into a state are treated identically. They accelerate the exact recurrence but do not certify small state spaces in general.

For fixed G, all endpoint responses and all deletion-card responses are monotone under enlargement of S, coefficient by coefficient. That is immediate from their nonnegative extension/join sums. This monotonicity concerns **the same old graph and same new vertex with more allowed incident edges**. It is not monotonicity in the square-sum interval size.

## 7. Square-sum results beyond the inherited scalar census

For Q_n on 1,...,n, the arriving neighbor set is exactly

    S_n = {q²-n : n < q² <= 2n-1}.

The open lower bound excludes zero, and 2n-1 excludes the forbidden loop. We apply the unrestricted transfer with these actual neighbors. Full endpoint matrices agree with an independent vertex-order DFS for every 1<=n<=30. Complete single-vertex deletion matrices agree by both methods for Q24, Q25, Q26, Q29 and Q30, 134 deletion cards in total.

The new interval data are:

| n | paths up to reversal | endpoint pairs | endpoint extensions | interior joins |
|---|---:|---:|---:|---:|
| 25 | 10 | 10 | 0 | 10 |
| 26 | 12 | 9 | 2 | 10 |
| 27 | 35 | 16 | 3 | 32 |
| 28 | 52 | 19 | 5 | 47 |
| 29 | 19 | 15 | 0 | 19 |
| 30 | 20 | 15 | 2 | 18 |

Thus Q29 gives a concrete later compulsory-interior arrival even though Q28 is traceable: neither arriving neighbor 7 nor 20 is an endpoint of a Q28 Hamiltonian path. All 19 paths instead use the joining ear 7–29–20. The count drops from 52 to 19, so increasing interval size does not preserve path count or endpoint-pair support size.

Deletion responses are nonmonotone even when the undeleted path count rises: Q25 minus vertex 3 has one Hamiltonian path, while Q26 minus vertex 3 has none, although the full path counts rise from 10 to 12. Here the loss also has a short structural proof, not just a search certificate.

**Proposition 5 (structural deletion-support loss).** Q25−3 is traceable, but Q26−3 is not.

**Proof.** A positive path in Q25−3 is

    18,7,9,16,20,5,11,25,24,1,8,17,19,6,10,15,21,4,12,13,23,2,14,22.

In Q26−3, vertices 18 and 22 are leaves, so they must be the two endpoints of any Hamiltonian path. Vertex 9 is interior with neighbors {7,16}; hence 7–9 and the leaf edge 7–18 saturate 7, excluding 2–7. Interior vertices 13 and 26 have neighbor sets {12,23} and {10,23}, respectively. Their forced edges 13–23 and 26–23 saturate 23, excluding 2–23. But interior vertex 2 has precisely the neighbors {7,14,23}, leaving it only one usable edge. Contradiction. ∎

`verify_structural.py` checks all eight local adjacency lists used above and every edge of the positive path. This is exactly the endpoint-certified forcing inherited from the earlier chain-kernel method, applied to a new deletion transition. All remaining negative decisions in this section are FINITE-EXACT; none is an infinite claim. Full endpoint profiles, positive witnesses, and zero decisions are retained in RESULTS.json.

For n=23,...,30, the output additionally retains the complete predecessor port slice that can finish after the specified arrival. Its one- and two-component entries independently reproduce the extension/join split. The target-aware pruning never changes a retained state's coefficient: every discarded earlier state was already unable to complete under that fixed target.

## 8. Reproduction and limits

Run from this directory:

    python3 verify.py > RESULTS.json
    python3 -O verify.py > RESULTS.optimized.json
    cmp RESULTS.json RESULTS.optimized.json
    python3 verify_structural.py > STRUCTURAL.json

`verify.py` implements the component-pair transfer. `direct.py` is a separate vertex-order DFS that does not call that transfer. Additional edge-subset enumeration compares the **entire all-component profile** for every labeled simple graph of orders 1 through 5: 1,099 graphs. The same universe compares the pruned endpoint matrix with DFS. There are 4,306 deletion-update tests and 192 physical four-port response queries. Explicit exception checks, not Python assertions, guard correctness.

The results extend the prescribed-endpoint/deletion direction identified in the earlier packet. They do not audit Gerbicz's infinite base/glue construction, settle all-order prescribed endpoints or deletion robustness, transfer planar coloring reducibility, or derive anything about Collatz. The useful next arithmetic step is to test whether a specific square-sum block construction preserves these named port/deletion responses; ordinary path existence alone cannot supply that invariant.
