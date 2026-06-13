# The Fractal Tournament Codec

**Session**: kind-pasteur-2026-03-24-S20cq

## The Idea

A tournament on n vertices contains n sub-tournaments on n-1 vertices (delete each vertex). This recursive structure can be exploited for compression: encode T(n) as T(n-1) + residual, where the residual describes how the deleted vertex connects to the rest.

## The Key Discovery: Vertex Selection Is Everything

### Theorem (Zero-Savings for Fixed-Index Deletion)

If the deleted vertex v is chosen by FIXED INDEX (e.g., always vertex 0), the conditional entropy of the residual given the sub-tournament's isomorphism class equals the naive entropy:

    H(residual | class(T\v)) = n-1 bits = H(residual)

**Proof**: For any (n-1)-class C and any residual r (an (n-1)-bit string), the number of n-tournaments T where T\v is in class C and the residual is r equals |fiber(C)| -- the number of labeled tournaments in class C. This count is INDEPENDENT of r because any residual can be freely paired with any labeled sub-tournament. Therefore P(r|C) = 1/2^{n-1} for all r. The conditional distribution is uniform. QED.

**Verified computationally**: n=3 through n=7, all show exactly 2^{n-1} distinct residuals per sub-class, all uniformly distributed.

### Score-Based Deletion Breaks the Uniformity

When the deleted vertex is chosen based on the tournament's STRUCTURE (e.g., the vertex with minimum out-degree), the conditional entropy drops significantly:

| n | Naive (n-1 bits) | H(resid \| class, min-score) | Savings |
|---|---|---|---|
| 3 | 2 | 1.06 | 47% |
| 4 | 3 | 1.75 | 42% |
| 5 | 4 | 2.37 | 41% |
| 6 | 5 | 3.23 | 35% |

### The Mechanism: Score Dominates

Investigation 6 reveals that the savings decompose as:

| n | H(r) | H(r\|score) | H(r\|class) | H(r\|score,class) |
|---|---|---|---|---|
| 4 | 3.00 | 0.79 | 1.75 | 0.77 |
| 5 | 4.00 | 1.38 | 2.37 | 1.38 |
| 6 | 5.00 | 2.04 | 3.23 | 2.03 |

**The score constraint provides ~97% of the savings.** Knowing the score of the deleted vertex fixes the Hamming weight of the residual to C(n-1, score) possibilities instead of 2^{n-1}. The isomorphism class adds only ~0.01 extra bits of prediction.

This makes perfect sense: the score IS the marginal statistic of the residual. Knowing the score tells you how many 1-bits the residual has. The class only adds information about WHICH positions have 1-bits, and this is nearly uniform within a weight class.

## The Six-Level Compression Architecture

| Level | Method | Bits | Mechanism |
|---|---|---|---|
| 0 | Raw matrix | n^2 | Full adjacency matrix |
| 1 | Antisymmetry | C(n,2) | A[i][j] = 1 - A[j][i] |
| 2 | Isomorphism | log2(A000568(n)) | Burnside counting |
| 3 | Base path | C(n-1,2) | Fix a Hamiltonian path |
| 4 | Band-limited | C(n-1,2)/2 | Krawtchouk low-pass (S305-S307) |
| 5 | Fractal score | ~0.65*C(n,2) | Recursive score-conditioned coding |

Levels 4 and 5 are potentially COMBINABLE: at each recursive level, the residual can be further compressed using band-limitedness. Estimated combined: ~0.325 * C(n,2), or about 3:1 over naive.

## Why Fractal Recursion Still Matters (Despite Modest Ratio)

The raw compression ratio of the fractal codec (~1.5x with score-based deletion) is modest compared to the isomorphism limit (~2.4x at n=7). But the fractal structure provides:

1. **Progressive refinement**: decode levels 2, 3, 4, ... to get successively larger sub-tournaments. Like JPEG2000 progressive mode.

2. **Random access**: to get the k-vertex sub-structure, decode only levels 2..k.

3. **Natural integration with deletion-contraction**: THM-082 says H(T) = H(T\e) + H(T/e). The fractal codec's recursive structure mirrors the DC tree.

4. **Score sequence as side channel**: the score sequence is a byproduct of the fractal encoding (you learn one score at each level). The score sequence is a powerful invariant (~97% of H variance).

5. **Bridge to spectral methods**: the band-limitedness (Level 4) can be applied at each recursive level, multiplying the compression gains.

## Connection to Known Compression Paradigms

### Fractal Image Compression (Jacquin 1992)
- Images partitioned into range blocks, each matched to a scaled domain block
- Tournament analog: sub-tournament = domain block, residual = range-domain transform
- Key difference: tournament self-similarity is EXACT (same combinatorial type), not approximate

### Graph Compression (Choi-Szpankowski 2009)
- Vertex-by-vertex peeling with arithmetic coding
- Structural entropy: C(n,2)*h(p) - n*log(n) + O(n)
- The n*log(n) savings = the isomorphism factor
- Our result: score-based peeling gives better entropy than random peeling

### Spectral Graph Wavelets (Hammond et al. 2011)
- Wavelets defined via graph Laplacian eigenvectors
- Multiresolution on graphs: low-pass = smooth, high-pass = detail
- Tournament analog: Krawtchouk spectral decomposition IS the wavelet basis
- Band-limitedness = tournament structure is a low-pass signal

### Combinatorial Object Compression (Steinruecken 2016)
- Arithmetic coding with structural models for permutations, graphs, multisets
- Near-optimal for symmetric combinatorial objects
- Tournament analog: our conditional model P(residual | score, class) is exactly this

### Modular Decomposition
- Tournaments decompose into prime (indecomposable) components
- The modular decomposition tree IS a fractal structure
- Non-prime tournaments can be encoded much more compactly
- Most large tournaments are prime, limiting this approach

## The Deep Insight

The fractal codec reveals a fundamental truth about tournament information:

**The SCORE SEQUENCE carries almost all the compressible information.**

When you recursively build a tournament by adding vertices one at a time, the score of each new vertex is the dominant predictor of its connections. The specific wiring pattern (which of the C(n-1, score) possibilities) is nearly uniform.

This connects to the shadow principle (S150): the score sequence is a "shadow" of the tournament that captures ~97% of the structural information. The remaining ~3% is in the exact wiring -- the "cusp form" in the shadow metaphor.

In information-theoretic terms:
- **Score sequence**: n values, each in {0, ..., n-1}, with sum = C(n,2). About n*log(n) bits.
- **Wiring residual**: about C(n,2) - n*log(n) + O(n) bits. This is the SAME as the isomorphism gap!
- **The fractal codec separates these two contributions naturally.**

The score sequence is the SLOW VARIABLE (captures the class structure). The wiring is the FAST VARIABLE (scrambled by isomorphism). Compressing tournaments = compressing the slow variable efficiently + accepting the fast variable as incompressible noise.

## Scripts and Data

- `/04-computation/fractal_codec_fast_s20cq.py` -- fast level-by-level analysis
- `/04-computation/fractal_vertex_choice_s20cq.py` -- vertex selection investigation
- `/05-knowledge/results/fractal_codec_fast_s20cq.out` -- compression measurements
- `/05-knowledge/results/fractal_vertex_choice_s20cq.out` -- vertex choice results
- Related: `/04-computation/fractal_codec_s20fy.py` (prior session's prototype)
