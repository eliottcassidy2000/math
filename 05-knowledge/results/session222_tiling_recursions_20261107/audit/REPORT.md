# Session222 independent mathematical audit

Status: PROVED audit deductions and FINITE-EXACT independent checks. No priority claims.

## Decision

Accept the shell, nested response doubling, and lattice traversal constructions with their stated type distinctions. Correct THM-452 before treating its weight and affine-coset wording as compatible with the standard backward-bit tiling gauge. The audit found two substantive additional corrections, not merely an indexing typo: the purported coset representative c_n actually lies in the blue subspace, and the singleton input is an exception to the nonblue assertion.

## 1. Shell geometry and Hamilton counts

For m=n−3, the staircase Q_m has binom(n−1,2) free cells. Translating it by (1,1) in Q_(m+2) leaves an L of 2n−1 cells; retaining it anchored leaves diagonals of lengths n−1,n. The producer's piecewise map is a bijection and not a rigid isometry or tournament-edge embedding. Its leg inclusion preserves old arcs; its diagonal inclusion changes (x,y) to (x+2,y), so it need not preserve H.

The zero L-boundary adds a source and sink: removal gives an exact Hamilton-path bijection for every old tournament. The sole order-three reversed arc yields H=3, zero-L order-five H=3, zero-diagonal order-five H=9. This minimal informative counterexample is decisive.

Independent endpoint subset-DP checks enumerate **all 1,098** fixed descending-path tilings for n=3,4,5,6 (2,8,64,1024). Every zero-L extension preserves H. Zero-diagonal extension changes H in respectively 1,7,63,1023 cases: all nonzero assignments in this finite universe. This is NOT an all-order monotonicity theorem. Geometric maps independently checked at every n=3..100.

Parity checked separately on **all 1,099 unfiltered labelled tournaments** orders1..5. H is positive odd in each case, with maxima1,1,3,5,15. The all-order statement is inherited classical tournament parity, not proved by enumeration. It prohibits H'=2H for nonempty tournaments; it does not prohibit even differences between odd counts, nor doubling H−1.

## 2. Physical selector doubling

The nested lift is valid. At a marked10 with i≥1 and zero right charge, insert0 immediately before position i−1. Old contributing left bits gain one place of separation and double their charge; the inserted0 contributes1; the old immediate neighbor stays immediate and contributes0. Hence L'=2L+1 and N'=2N. Deleting the new vertex restores both endpoint graphs and the same marked-arc transposition, so induced containment holds for the literal one-arc-flip pair, not just its word encodings.

The strongness condition starts0/ends1 is maintained. The special N=1 seed counts9/11→19/23 are correct. The universal binary seeds for N≥2 are correct: binary(2N−1)=binary(N−1) followed by1, so complemented digits acquire a terminal0.

Independent reverse-recursive path enumeration checks **354 eligible marked pairs** in all words of lengths3..8, including31 pairs strong at both endpoints; all satisfy gap2N→4N. This tests all eligible marked locations, not just the producer's canonical seed family.

## 3. THM-452 correction and safe patch

Read `01-canon/theorems/THM-452-doubling-tiling-coset-law.md` against the backward-bit convention in THM-280/554 and the skew block construction.

### Gauge and weights

In the canonical path P followed by reversed second-copy P, every twin arc is forward. Exactly one twin is a path edge and is not a free tile; the other n−1 twins are free, forward tiles. Thus §1's twin bits1 and §6's recurrence use **forward-bit1**, the complement of the standard backward-bit1 convention.

Each same-copy block has weight w in its own input gauge. Off-twin cross tiles form binom(n,2) complementary pairs and contribute exactly binom(n,2). Therefore

- backward gauge: w_b(2n)=2w_b(n)+binom(n,2);
- forward gauge: w_f(2n)=2w_f(n)+(n−1)+binom(n,2).

At orders1,2,4,8,16,32, backward weights are0,0,1,8,44,208 and forward weights0,0,2,13,61,257. The original five weights were attached to four orders, omitting order1.

A minimal safe patch explicitly declares that the existing §1 formulas use forward bits, corrects §6's order list, and gives the backward equivalent. If instead rewriting §1 in the global backward convention, keep the two identity/mirror copies but replace cross A by1−A and twins1 by0. Do not change only the weight recurrence while leaving the setup's gauge implicit.

### Affine fiber, not translation by the syndrome

Let B=ker(1+σ), c_n the off-twin cross indicator. Then σ(c_n)=c_n, so c_n∈B. Consequently B+c_n=B. The existing phrase that the nonblue D-image lies in `(blue subspace)+c_n` is false for n≥2.

The correct statement is

    t(D(x)) ∈ {t : (1+σ)t=c_n} = t0+B,

where t0 is any one D-image, for example the double of the transitive framed input in the declared gauge. The syndrome is c_n; a representative is t0, not c_n. The dimension and relative-density formulas remain valid: dimB=n(n−1), input dimension binom(n−1,2), density2^(−(n−1)(n+2)/2).

### Singleton boundary

For n=1, the outer tournament has two vertices and zero tiles. The syndrome is0 and the unique tiling is blue. Thus §2's “NO tournament” and §3's “never itself a D-image” require **n≥2**. For n≥2 the syndrome is nonzero. No claim here re-audits the unrelated exhaustive anti-automorphism statements in §4 or the Walsh formula in §5.

Independent direct **arc-set construction**, not imported signed-matrix code, verifies these claims for all **1,100 framed inputs of orders1..6**, and verifies the tower through order32. Explicit checks remain active under Python `-O`.

Recorded correction mechanism: an undeclared global bit complement made a correct forward-weight law appear to be the canonical backward-weight law; a linear-map output/syndrome was confused with an affine-fiber representative; an all-input nonblue quantifier overlooked the zero-dimensional singleton boundary. Preserve the original forward formulas as correctly scoped content rather than deleting them.

## 4. Odd-chain plane

Unique factorization n=2^k(2j+1) labels the quadrant grid bijectively. Doubling advances k; odd-core successor advances j and adds2^(k+1), not uniformly2. The producer's shell/rank and plane/rank formulas are inverse bijections. Independent shell enumeration by concatenated edge lists verifies100 complete shells and10,000 distinct vertices in each grid, with every consecutive visit at Manhattan distance1.

Transporting labels along the two traversals does not preserve all grid edges. Our independent witness is quadrant neighbors(1,1),(1,2), whose whole-plane images are(1,1),(0,−1), distance3. The producer gives a separate horizontal witness. The quadrant corner has degree2 while every plane grid vertex has degree4, forbidding a full graph isomorphism.

The literal undirected operations +2 and ×2 contain the five-cycle1–3–6–4–2–1, so all such edges cannot embed injectively as unit square-grid edges, since that grid is bipartite.

A path through every integer lattice center is not a continuous space-filling surjection onto R². Unit squares centered at those vertices tile the plane; their countable union of center-joining segments has area0 and does not fill it. Use “lattice Hamiltonian traversal” and “cell tiling” precisely. No missing Peano-style continuity theorem is implied.

## Reproduction

From the packet directory:

    python3 -O audit/audit.py
    python3 -O audit/doubling_audit.py
    python3 -O audit/selector_audit.py
    python3 -O audit/plane_audit.py

Outputs: RESULT.json, DOUBLING_RESULT.json, SELECTOR_RESULT.json, PLANE_RESULT.json. All counters and universal ranges above are explicit. No producer implementation is imported. The shell program uses endpoint subset DP, doubling uses directed-edge sets, selector uses memoized remaining-set recursion, and plane uses concatenated shell-edge lists. These finite checks supplement the symbolic proofs, not replace them.

## 5. Subsequent positive continuous extension audit

The odd-plane producer subsequently supplied REPORT.md §6, an explicit continuous filler, beyond the polygonal traversal. This extension is **accepted**. On the Cantor set, split odd/even ternary0/2 digits into binary coordinates. Restricted ternary expansions are unique; if the first difference is at l, the input separation is at least3^(−l), proving continuity with output coordinate oscillation at most2^(−floor(m/2)) for matching first m digits. Every square point has an interleaved preimage. Affine interpolation over deleted intervals stays in the square and is continuous at their endpoints and accumulation points: shrinking gaps have both endpoints tending to the target; gaps bounded below have a fixed-interval subsequence. This covers every possible discontinuity subsequence.

For each visited closed unit square, the three stages are the incoming edge-midpoint connector, translated square-filler, and outgoing edge-midpoint connector. The displayed coefficients match at1/3,2/3, and the next integer time. Connectors remain in the same convex square; the middle stage covers it. Locally finite gluing yields a continuous G:[0,∞)→R², with G([t,t+1]) equal to the whole visited square. No continuity at infinity is required. The square enumeration makes the union all of R².

The extension is substantive, noninjective, and distinct from joining lattice centers by straight segments. It preserves the cell schedule, not a bijection of curve points with integers or the global ×2 adjacency. The producer explicitly marks the Cantor construction as elementary, without a novelty claim. No blocker remains for this scoped theorem.
