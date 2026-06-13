# Even Graphs Through the Metagraph Lens

**Session**: opus-2026-03-23-S260

---

## The Three Meanings of "Even Graph"

There are three distinct notions of "even" for graphs, and confusing them leads to errors:

### 1. Degree-Even (= Eulerian = Cycle Space = Two-Graph)
A graph where **every vertex has even degree**. These form the **cycle space** of K_n over GF(2), with dimension C(n-1,2). The number of iso classes is **A002854**: 2, 3, 7, 16, 54, 243, 2038, ...
(Note: A003049 counts the CONNECTED version; A002854 counts all Euler graphs.)

### 2. Royle-Even (= Automorphism-Orientation Compatible)
A graph where every orientation and every non-identity automorphism reverses an **even** number of arcs. The number of Royle-even iso classes equals **A000568** (tournament count): 2, 4, 12, 56, 456, 6880, ...

### 3. Edge-Count Even
A graph with an even number of edges. Not useful for our purposes.

**The Royle equinumerosity** (arXiv:2204.01947, Devillers-Freedman-Glasby-Praeger-Royle 2022): #(Royle-even graphs on n vertices) = #(tournaments on n vertices) = A000568(n). This is **NOT** about degree-even graphs.

**Key caution**: The Royle equinumerosity is NOT simply "split Burnside by permutation order." The Burnside split gives:
- Σ_{|g| odd} 2^{c(g_E)} / n! = A000568 (tournaments)
- Σ_{|g| even} 2^{c(g_E)} / n! = A000088 - A000568 (NOT equal to A000568 for n ≥ 4)

The Royle-even partition is based on an **intrinsic graph property** (no automorphism reversing an odd number of oriented arcs), not on the Burnside split. That this partition yields A000568 is a deep combinatorial fact.

The conjecture was originally made by Pontus von Bromssen, who computed even graphs on ≤ 10 vertices and noticed the match with A000568.

---

## The GF(2) Parity Obstruction (NEW THEOREM)

**Theorem**: Over GF(2), the edge space of K_n decomposes as Cut ⊕ Cycle **if and only if n is odd**.

**Proof sketch**: The cut space C and cycle space Z of K_n satisfy C ∩ Z = {vectors that are both cuts and have all even degrees}. A cut δ(S) gives degree |V\S| to vertices in S and |S| to vertices in V\S. For K_n:

- If n is odd: the only cut with all even degrees is δ(∅) = 0. So C ∩ Z = {0}, and Edge = C ⊕ Z. ✓
- If n is even: δ(S) with |S| = n/2 gives degree n/2 (even when n ≡ 0 mod 4) to every vertex. More generally, ANY cut with even-sized S gives all-even degrees when n ≡ 0 mod 4. The intersection dim(C ∩ Z) = n - 2. ✗

**Consequence**: The tournament → even graph projection via cycle space is only well-defined for **odd n**.

**Proof at n=4**: The three 4-cycles on K_4 are simultaneously cuts AND cycles. They generate a 2-dimensional intersection. The cut space has dimension 3, the cycle space has dimension 3, and they overlap in 2 dimensions, spanning only 3+3-2 = 4 dimensions of the 6-dimensional edge space.

---

## The Cycle-Space Projection (Odd n Only)

For odd n, every tournament T has a unique decomposition:

    T = T_cut ⊕ T_cycle

where T_cut encodes the **score hierarchy** and T_cycle encodes the **cyclic structure** (a degree-even graph).

**The projection formula**: T_cycle = (I + L(K_n)) · T mod 2, where L(K_n) is the adjacency matrix of the line graph of K_n.

**This is a single matrix multiplication** — no iteration or search needed.

### Results at Odd n

| n | V_tournament | V_even (degree) | Ratio | Max fiber | H variance: cycle | H variance: score |
|---|---|---|---|---|---|---|
| 3 | 2 | 2 | 1.000 | 1 | 100.0% | 0.0% |
| 5 | 12 | 5 | 2.400 | 4 | 85.3% | 14.7% |
| 7 | 456 | 38 | 12.000 | 65 | 72.7% | 27.3% |

### The Information Inversion

At the **labeled** level: scores carry 97% of H information (cut space dominates).

At the **iso-class** level: cycle structure carries 73-85% of H variance (cycle space dominates).

**Why the inversion?** Scores are "macroscopic" — highly informative per bit but redundant across relabelings. Cycle structure is "microscopic" — individually weak but collectively dominant after quotienting by S_n. The act of taking isomorphism classes **strips away the score redundancy** and reveals the cycle structure as the true differentiator.

This is analogous to thermodynamics: temperature (macroscopic) determines most of the energy at the particle level, but the microstate (microscopic) determines the identity of the system at the statistical level.

---

## Fiber Structure

The map T → T_cycle sends multiple tournament iso classes to the same even graph iso class. The "fiber" over each even graph class consists of tournaments with the same cyclic structure but different score assignments.

### Fiber Properties
- **Not injective** for n ≥ 4 (even in the odd-n projection)
- **Fiber sizes grow rapidly**: max fiber is 1, 4, 65 at n=3,5,7
- **Fibers contain multiple score types**: up to 20 distinct score sequences in a single fiber at n=7
- **SC tournaments are NOT enriched** in large fibers (SC rate ≈ overall rate)
- **~15% of metagraph edges stay within the same fiber** (same even graph class)

### The Largest Fiber at n=7
65 tournament classes project to a single even graph with 10 edges, degree sequence (4,2,2,2,4,4,2). These 65 classes span 20 distinct score types and H values from 15 to 111. The even graph is a "universal backbone" carrying many hierarchical instantiations.

---

## Metagraph Edge Projection

When two tournament classes are connected by an arc flip in G_n, their even graph projections usually differ:

| n | Tournament edges | Same even class | Different even class |
|---|---|---|---|
| 3 | 1 | 0 (0%) | 1 (100%) |
| 5 | 30 | 7 (23%) | 23 (77%) |
| 7 | 4086 | 615 (15%) | 3471 (85%) |

**~85% of arc flips change the cyclic structure.** Only 15% preserve it (changing only the score assignment). This means the metagraph's edge structure is primarily about cyclic topology changes, not score reshuffling.

---

## Open Questions

1. **Royle-even characterization**: What is the explicit characterization of Royle-even graphs? The automorphism-orientation condition is computationally expensive (exponential in edge count). Is there a polynomial criterion?

2. **V_even at even n**: Without the cycle-space projection, how do tournament classes relate to degree-even graph classes at even n? The n-2 dimensional intersection creates ambiguity.

3. **Asymptotic ratio**: V_tourn/V_even(degree) grows as 1, 2.4, 12 at odd n. What's the growth rate? Does it approach C · n! / 2^n for some C?

4. **The 15% stabilization**: Does the fraction of edges preserving even class stabilize? At what value?

5. **Even graph metagraph**: The even graph flip graph (edge toggle) has FEWER edges than the tournament metagraph (arc flip). Why? Is this because arc flips are "richer" than edge flips for cycle-space modification?
