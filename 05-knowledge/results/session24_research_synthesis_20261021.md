# Session 24: converse exceptions, path-count multiplicities, and metagraph transport

Research date: 2026-10-21. All results below are scoped to their stated objects. Proof packets include independent reviews and reproducible exact controls. No claim of literature priority is made. The classification uses cited classical group theory; the all-order low-count inventory uses Busch's published lower bound.

## Main result: the converse-symmetry exceptions are completely classified

A finite tournament of order at least two has no vertex subset with trivial setwise stabilizer under its automorphisms and anti-automorphisms **if and only if** it is

\[
X[C_3^{[k]}],\qquad X\in\{C_3,R_5,P_7\},\quad k\ge0.
\]

Here R5 is the regular cyclic five-vertex tournament, and the bracket operation substitutes a directed triangle tower at every vertex. This closes the repository's HYP-9353, not merely its primitive case.

The missing structural argument has five steps:
1. If ordinary automorphisms have several vertex orbits, complement a regular subset on one odd orbit to break every possible anti-symmetry. Hence exceptions are vertex-transitive and strong.
2. An arbitrary permutation-block system has an odd-cardinality majority tournament, even before its blocks are known to be modules.
3. Regular-subset palettes and the primitive classification force minimal nontrivial blocks to be triangles.
4. Let each triangle carry a singleton/doubleton marker. Failure of every marker choice forces the full independent triangle-rotation kernel into the actual automorphism group.
5. Independent rotations make those blocks genuine modules. Induction gives exactly the three seed families.

The proof does not assume the finite census establishes completeness. Historical premature completeness prose was corrected, and its later proof is identified as a separate argument.

Proof and review: [classification packet](converse_exception_classification_20261021/session24-structure/STRUCTURE.md).

## Exact generating functions and labelled counts

There is exactly one exceptional isomorphism class at each order in

\[
\{3^{k+1},5\cdot3^k,7\cdot3^k:k\ge0\}.
\]

Its ordinary generating function satisfies

\[
E(z)=z^3+z^5+z^7+E(z^3).
\]

For a seed of order d and automorphism order a, the member at order n=d3^k has

\[
|\mathrm{Aut}|=a3^{(n-d)/2},\qquad
\#\text{labelled copies}=\frac{n!}{a3^{(n-d)/2}},
\]

where (d,a)=(3,3),(5,5),(7,21). The labelled exponential generating function satisfies

\[
L(z)=z^3/3+z^5/5+z^7/21+L(z^3/3).
\]

Both generating functions are non-D-finite: their coefficient sequences have arbitrarily long zero runs followed by nonzero terms, incompatible with a fixed polynomial-coefficient recurrence beyond its finitely many leading-coefficient roots.

The exact twisted palette transfers also give a quantitative distinction: with any fixed number c≥3 of named colors, almost every ordinary-regular coloring is signed-regular along these towers, yet the probability that an arbitrary coloring is ordinary-regular decays exponentially in the number of vertices. All three seed families share the same negative exponential rate per vertex.

Details: [palette dynamics](exception_palette_dynamics_20261021/session24-palette-dynamics/PALETTE_DYNAMICS.md).

## How often does a Hamiltonian-path count occur?

Let N_n(H) count ordinary isomorphism classes of order n with exactly H directed Hamiltonian paths; S_n(H) counts the self-converse ones, and M_n(H)=(N_n(H)+S_n(H))/2 counts converse-merged nodes. Every attained H occurs infinitely often across unbounded orders by adding sources or sinks, so order-indexed multiplicity is the useful question.

For strong atoms define A_q(z)=sum z^|T| over strong classes with H(T)=q. Unique ordered SCC decomposition gives

\[
F_H(z)=\sum_{q_1\cdots q_r=H}\frac{A_{q_1}(z)\cdots A_{q_r}(z)}{(1-z)^{r+1}}.
\]

The factors satisfy q_i>1; the empty product supplies F_1=1/(1-z). A self-converse SCC word has two converse halves and an optional self-converse center, giving an exact square-divisor formula for S_n(H).

Consequences:
- Fixed-H ordinary multiplicities are eventually polynomials, not exponential sequences.
- Their exact degree is the maximum number of nontrivial strong SCC factors with product H.
- Busch's bound H(T)≥5^((|T|−1)/3) makes the polynomial valid already for n≥floor(3 log_5 H). Symmetry/merged counts are eventual period-dividing-two quasipolynomials with an equally explicit logarithmic onset.
- Thus the entire fixed-H problem reduces to strong cores of order O(log H). This is a constrained enumeration improvement, not a general speedup for unrestricted A000568.

A complete census through order8, independently checked against the labelled universe, gives a globally complete strong-atom inventory for every odd H≤73. Selected exact laws:

| H | N_n(H) | S_n(H) | Validity |
|---:|---:|---:|---|
| 3 | n−2 | 1 for odd n, otherwise0 | n≥3 |
| 5 | n−3 | 1 for even n, otherwise0 | n≥4 |
| 15 | (n−4)(n−3) | 2 for odd n, otherwise0 | n≥5 |
| 19 | 2(n−5) | 0 | n≥6 |
| 33 | n²−5n−1 | 1 | n≥6 |
| 35 | 2(n−6) | 0 | n≥7 |
| 39 | (n−6)(n−5) | 0 | n≥7 |
| 41 | 6n−35 | 1 | n≥6 |
| 49 | 8n−49 | 1 | n≥7 |
| 63 | 2(n−7) | 0 | n≥8 |

The 39 row must include the SCC factorization3·13; simply padding strong H39 atoms gives the wrong answer. At order8, H345 is the most frequent count:49 ordinary classes,25 merged classes, one self-converse class. This is finite-exact, not an assertion of persistent unimodality.

Proof and full tables: [actual-H multiplicities](actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).

## The {7,21} gap has an exact structural effect

For a≥2,

\[
N_n(7\cdot3^a)\sim\frac{2}{(a-2)!}n^{a-1},
\qquad N_n(3^a)\sim\frac{n^a}{a!}.
\]

The first degree is two below the naive prime-factor count a+1. The unique factor carrying7 must also carry at least two3s, because7 and21 are unavailable. Exactly two strong H63 atoms realize the top-length factorization, both automorphism-free and mutually converse.

This also changes the typical metagraph degree. For fixed H=3^a and growing order, almost all merged nodes are black-degree-one leaves. For fixed H=7·3^a, almost all have black degree63. These are the fundamental paired-black multigraph degrees, not distinct-neighbor counts.

This is a proved connection from missing path counts to SCC structure and metagraph statistics. It does **not** prove that every other odd H occurs, nor explain the numerical coincidence between exceptional vertex orders7,21 and missing path counts7,21.

## Exact Hamilton and converse counts of the exceptional family

Ordinary path-cover substitution was already proved in THM-3121. Its new reversal-fixed refinement retains the local cover involution, allowing exact ordinary and converse-fixed orbit counts without enumerating all vertex permutations.

If h=H/|Aut| and b is the number of converse-fixed ordinary path orbits, the self-converse node's fundamental paired-black degree is (h−b)/2.

| Tournament | Order | H | h | b | Black degree |
|---|---:|---:|---:|---:|---:|
| C3 | 3 | 3 | 1 | 1 | 0 |
| R5 | 5 | 15 | 3 | 3 | 0 |
| P7 | 7 | 189 | 9 | 3 | 3 |
| C3[C3] | 9 | 3159 | 39 | 11 | 14 |
| R5[C3] | 15 | 197030475 | 162165 | 873 | 80646 |
| P7[C3] | 21 | 123156134945235 | 2681562805 | 71083 | 1340745861 |
| C3^[3] | 27 | 417382500592116023859 | 261792936934433 | 40255495 | 130896448339469 |

The triangle transfer uses O(m⁴) integer operations and O(m³) storage for a full cover vector at order3m; its fixed-cover refinement is smaller. Exact tower Hamilton/fixed counts were computed through order243, with independent vertex-level checks at the feasible small orders. Fixed outer templates give polynomial arithmetic algorithms throughout all three exception families. Bit costs are documented separately in the packet.

Details: [path-cover transfers](exception_pathcover_transfers_20261021/session24-pathcovers/PATHCOVERS.md).

## Two genuine metagraph connections, and one obstruction

Current repaired canon and the new classification imply

\[
\{\text{signed exceptions}\}\cap\{\text{path-rigid}\}=\{C_3\},
\]
\[
\{\text{signed exceptions}\}\cap\{\text{pure-blue}\}=\{C_3,R_5\}.
\]

Triangle substitution injects all converse-merged node classes across orders n to3n. The block-contiguous Hamilton-path lift commutes with converse reversal, preserves self-converse/paired status, and never decreases fundamental black degree. Consequently no signed exception is a black leaf: C3 and R5 have degree0 in the black subgraph, and every other exception has degree at least3.

But this lift does not commute with fixed-path complementation. The smallest explicit test gives

\[
H(\kappa(C_3[C_3]))=2983,
\qquad H(\kappa(C_3)[C_3])=27.
\]

At each boundary, the fixed Hamilton path preserves one port arc while all other cross-block arcs reverse. The original module structure is destroyed. Actual entry/exit ports are therefore indispensable for a proposed wiring transport. This does not refute every possible unmarked adjacency correspondence or classify nonleaf bridges.

Proof and independent permutation check: [metagraph connection corollaries](metagraph_connection_corollaries_20261021/session24-connections/CONNECTIONS.md).

## The 341 connection preserves actual arithmetic structure

The older trunk map R(n)=4n+1 is conjugate modulo341 to multiplication by4 through L(n)=3n+1. The newer golden quotient sends its generator to81 on the same carrier. Their joint multiplier group is exactly QR11×QR31, of order75.

This yields exactly16 compatible Cayley connection sets and four tournament isomorphism classes: two lexicographic types and two prime axis-reversed types. An all-order theorem proves the axis-reversed construction prime for arbitrary skew Cayley factors of coprime odd cyclic orders m,n>1. Coprimality matters: a diagonal module survives at3×3.

For the prime341 types, an intrinsic common-neighbor relation recovers the CRT blocks; minority matchings force the full automorphism group, of order25575. The three-corner subset {0,155,187} is a minimum-size signed-regular subset. Every member's Hamiltonian-path count is divisible by25575, but its actual count has not been computed.

The bracket prime set{2,3,11} was also recovered in its original meaning. Its squareclasses provide a Fano/Hamming labeling, but any three independent primes give that same construction. No special bracket-to-Hamilton-gap mechanism is claimed.

Details: [arithmetic clocks and tournaments](arithmetic_clock_tournaments_20261021/session24-arithmetic/README.md).

## Constructive three-point Paley certificates

For every prime p≥11 with p≡3 mod4, the same subset {0,1,3} breaks every automorphism and anti-automorphism of P_p. Three is the minimum possible size: every one- or two-point set has an affine anti-reflection stabilizer.

The proof identifies the full signed group as AGL(1,p), using Müller's cited prime-degree theorem, then normalizes a triple to {0,1,t}. The only failures are midpoint configurations t∈{−1,2,1/2} or the quadratic t²−t+1=0. At t=3 the latter becomes7=0. This gives an elementary direct explanation of the prime Paley exceptions3 and7; by complementing subsets, the zero regular-triple count at7 rules out every subset size.

For p>3, the exact number of signed-regular triples in the full affine action is

\[
\frac{p(p-1)}6\bigl(p-5-2\,[p\equiv1\pmod3]\bigr).
\]

The fixed certificate must not be extrapolated to nonprime fields: Frobenius fixes prime-subfield points. A separate full-degree construction repairs this. Write q=p^e≡3 mod4, e>1 odd, and choose t of degree exactly e. If3 does not divide e, every such t works. If3 divides e, put r=p^(e/3); precisely the two equations t^r=1/(1−t) and t^r=(t−1)/t obstruct regularity. A root-count bound proves that a full-degree parameter avoiding both always exists.

Consequently the minimum signed-regular subset size is exactly three for EVERY Paley tournament except P3 and P7, where no regular subset exists. A Möbius-inversion formula counts the full semilinear triple orbits exactly: the numbers at q=27,243,343 are1,8,18. This uses the full affine-semilinear group from a separately verified classical theorem, not the invalid replacement of that group by its affine subgroup. At q27, full degree alone can fail through an order-three semilinear stabilizer, providing a decisive hostile example.

Details: [prime certificates](paley_prime_certificates_20261021/session24-paley-certificates/RESULT.md) and [prime-power certificates and exact counts](paley_prime_power_markers_20261021/session24-paley-powers/RESULT.md).

## Concurrent-source reconciliation

Incoming audited commit071409c4 and closing letterd2dd623e were preserved, including MISTAKE-650 and the distinction between triangle substitution and ordered stacking. Concurrent S50 checkpointae56f34f also develops the SCC multiplicity law and a pre-audit order-nine census. The independently checked census in this report remains through order8; that incoming order-nine work is not silently promoted by this report. Reserved theorem/hypothesis identifiers were preserved, and these packets use descriptive result paths rather than claiming those namespaces.

## Remaining boundaries

Still not established here:
- the assertion that{7,21} is the complete forbidden odd Hamilton-count set;
- complete metagraph wiring or all bridges, including nonleaf bridges;
- every large-H strong-atom inventory;
- a general faster enumeration algorithm for unrestricted A000568;
- any transfer proving Collatz, the lonely-runner conjecture, or Hadwiger–Nelson.

The next focused mathematical direction is port-aware transport: combine the exact cover vectors with the boundary defects, rather than attempt to infer adjacency from counts or symmetry labels alone.
