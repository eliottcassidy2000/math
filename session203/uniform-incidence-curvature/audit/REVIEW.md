# Independent audit of uniform-incidence highest curvature

6 November 2026. **PASS in the stated uniform-incidence scope.** The all-order degree theorem, exact 2-adic valuation, and even-multiplier obstruction are correct. No general nonuniform-incidence obstruction or odd-spectrum obstruction follows.

## Frozen producer inputs

- PROOF.md: 5b7fc292400fbf3193fb1bf520288206dff3733cfdc05239bc36f062c030a5bd.
- check.cpp: 8267c6877873d1347e1d6f9f73ecdbca100e05bb1d1d002975b5b6e5016e7fc0.
- NORMAL.txt: b2d4c0248197a362eedf8a5779f8101f72fb2ee2228efcbd78dd17d79477f3e5.

The reviewed proof is preserved as REVIEWED_PROOF.md. The producer code and output are also frozen here. No mathematical correction was required. A status-only promotion is approved, provided the context assumptions and scope paragraphs stay unchanged; the promoted hash should be recorded separately.

## All-order argument

The bit-dependent edge indicators in any Hamilton word have total number at most 2k, because only edges incident to the k sensitive vertices depend on bits. Boolean multilinear reduction cannot increase that number. Thus the Boolean degree is at most 2k.

For a term to involve all 2k distinct free bits, every sensitive vertex must be internal and adjacent on both sides to free spine vertices. All 2k incidences must use distinct free vertices. A repeated variable, sensitive-sensitive adjacency, endpoint sensitive vertex, or sensitive-to-fixed-core adjacency uses capacity without supplying a new required bit. Hence every top-degree term has a unique decomposition into k disjoint triples (free, sensitive, free), separated or joined among the remaining vertices W.

Each triple contributes top coefficient -1 regardless of the sensitive vertex's copy/complement polarity. This establishes a common sign (-1)^k; there is no hidden cancellation from the sensitive core edges, which are not used in these terms.

For an ordering of W with r occupied gaps, the k blocks form r nonempty runs. There are binomial(k-1,r-1) choices of their run lengths and k! assignments of sensitive labels. The remaining constraints on the 2k free labels are k-r forward spine inequalities between neighboring blocks. The pairs of slots in these inequalities are disjoint: a block's first and last free slots are different. Independent swaps therefore divide the (2k)! assignments by exactly 2^(k-r), with no additional restriction from W because the free positions are consecutive and W is uniform toward that block.

This is a bijective count. The W ordering, occupied gaps, run lengths, sensitive labels, and free labels can all be recovered from the Hamilton word. It yields exactly the producer's formula

    Delta = (-1)^k [k!(2k)! / 2^(k-1)] Q_k,
    Q_k = sum_(r=1..k) a_r binomial(k-1,r-1) 2^(r-1).

The identification a_1=H(W+x) is exact, including W empty. The standard tournament Hamilton-path oddness theorem gives a_1 positive and odd. All other terms in Q_k are nonnegative and even, so Q_k is positive and odd. The general dependency is precisely that standard oddness theorem; no strong-order bound, exhaustive tournament census, or conjectured closure is used.

Since

    k!(2k)! / 2^(k-1) = 2 (k!)^2 (2k-1)!!,

its exact 2-adic valuation is 1+2v2(k!). Nonvanishing supplies degree at least 2k and thus exact degree 2k whenever m>=2k. Fixing arbitrary exterior bits is harmless: the exterior spine vertices lie entirely before or after the consecutive block and enter W with uniform incidences. Edges from sensitive vertices to those fixed exterior vertices are irrelevant to the top-degree terms.

## Even-multiplier consequence and boundaries

For two contexts on the same m-bit spine, with at least one positive sensitive count and m>=2 max(k_P,k_C), unequal sensitive counts give unequal Boolean degrees. A context with k=0 is a nonzero constant function, so this also covers the case where only one count is zero.

If the positive counts agree, the same consecutive-block derivative has exact valuation 1+2v2(k!) in both contexts. Multiplication by a nonzero even integer changes that valuation, while addition of a bit-independent constant does not change the derivative. Thus H_C=a H_P+c is impossible on the full Boolean cube for such an even a. Negative even integers are included. The constant may depend on the chosen contexts or length, but not on the bits.

This remains valid for length-dependent contexts only at lengths satisfying the displayed size inequality. It says nothing about contexts with more than m/2 sensitive vertices, restricted bit languages, position-dependent incidence tables, nonlinear word-dependent incidences, or changed spine edges. Both-sensitive-counts-zero is explicitly excluded. Odd multipliers are not ruled out; indeed prefixing any context by an ordered directed triangle gives H_C=3H_P with the same sensitive count and otherwise allowed uniform core attachments.

## Independent controls

The independent count.cpp constructs each tournament literally and counts Hamilton paths by a starting-at-vertex subset recurrence. For every graph of order at most 8 it also performs a separate recursive path enumeration without memoization. It imports no producer code. The producer instead uses outgoing push transitions for paths ending at a vertex.

The independent verify.py computes the entire Boolean polynomial by a Mobius transform of each literal cube. It enumerates the W gap coefficients by classifying gaps as allowed and mandatory and summing binomial(allowed-mandatory,r-mandatory), rather than multiplying the producer's gap polynomials. It also checks a_1 directly by literal permutation enumeration of W+x.

The retained normal and optimized-Python/debug-C++ receipts agree byte-for-byte:

- 88 contexts with k=0,1,2,3,4 and varied core edges, polarities, uniform incidences, and exterior spine vertices;
- 4,256 literal graph counts, including 768 also checked by recursive path enumeration;
- 324 consecutive-block mixed-difference checks, including every assignment of the selected exterior bits;
- exact full Boolean degree 2k (or 0 when k=0) in every context;
- all predicted coefficients, positive odd Q_k values, and exact valuations.

Two hostile controls enforce the scope:

1. With one sensitive vertex and four spine positions, freeing nonconsecutive positions 0 and 3 while fixing the middle bits to zero gives coefficient -4, of valuation 2 rather than the consecutive-block value 1. The theorem does not assert its valuation formula for scattered free positions.
2. With two sensitive vertices and four spine positions, let both sensitive vertices copy all bits except that the second complements position 0 only. The literal cube has degree 3 and fourth coefficient zero, rather than degree 4. This violates uniform polarity across positions and is outside the theorem. All 16 counts are independently checked by literal permutations.

An additional four-assignment positive control verifies the allowed odd multiplier 3 by literal permutations before and after prefixing a directed triangle. It is not an attempted even-multiplier counterexample.

The producer's original 38-context, 1,520-graph output was replayed exactly. These computations are controls for the all-order bijection, not substitutes for it.

## Reproduction

Run from this audit directory:

    g++ -O2 -std=c++17 count.cpp -o /tmp/uniform_curvature_independent
    python3 verify.py /tmp/uniform_curvature_independent > replay.normal.json
    cmp replay.normal.json NORMAL.json
    g++ -O0 -std=c++17 count.cpp -o /tmp/uniform_curvature_independent_O0
    python3 -O verify.py /tmp/uniform_curvature_independent_O0 > replay.debug.json
    cmp replay.debug.json OPTIMIZED_PY_DEBUG_CPP.json
    g++ -O3 -std=c++17 producer_check.cpp -o /tmp/uniform_curvature_producer
    /tmp/uniform_curvature_producer > replay.producer.txt
    cmp replay.producer.txt producer_NORMAL.txt
    sha256sum -c SHA256SUMS

All validations use explicit exceptions or failure branches and survive removal of Python assertions. The folder is portable and needs no sibling scripts or banks.
