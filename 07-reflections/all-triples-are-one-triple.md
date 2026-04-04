# All Triples Are One Triple: Fixed / Boundary / Free

*opus-2026-04-04-S8*

## The Observation

This project has discovered at least 20 triples. They keep appearing independently — in the metagraph, in the tiling model, in number theory, in the functor decomposition. The question is: are these 20 triples related, or are they 20 coincidences?

They are **one triple**, wearing different masks.

## The Universal Triple: Fixed / Boundary / Free

Every triple in this project is an instance of a single abstract decomposition:

**Given a space X and a symmetry (or incidence relation) G acting on X, the space decomposes into three parts:**

1. **FIXED**: Elements stabilized by G. Self-interacting. Diagonal. The part that doesn't change.
2. **BOUNDARY**: Elements partially stabilized — they live on the interface between fixed and free. Mixed type. Always bipartite relative to the fixed set.
3. **FREE**: Elements acted on freely by G. Generic. The bulk. Dominates at large scale.

This is the decomposition of any group action into fixed points, boundary of the fixed set, and free orbits. It is also the decomposition of a matrix into diagonal, off-diagonal-adjacent, and off-diagonal-distant blocks.

## The Mapping

| Triple | Symmetry G | FIXED | BOUNDARY | FREE |
|--------|-----------|-------|----------|------|
| **Eisenstein** | Frobenius on Z[ω] | Inert (preserved) | Ramified (singular) | Split (factored) |
| **Computation** | Reversibility | Lossless (bijective) | Lossy (projection) | Crystallization (irreducible) |
| **Spine/Ribs/Sea** | Complement Z₂ | SC-SC (self-comp) | SC-NS (mixed) | NS-NS (generic) |
| **Tile pairs** | Vertex incidence | Same-End (shared) | Cross-End (relay) | Disjoint (independent) |
| **Three Functors** | Linearity | F1: linear (trivial) | F2: nonlinear local (crux) | F3: algebraic global (evaluation) |
| **Blue/Black/Red** | Grid reflection | Blue (grid-sym) | Black (generic) | Red (transpose) |
| **Triangle sides** | Reflection axis | Hypotenuse (diagonal) | Vertical leg (source) | Horizontal leg (sink) |
| **Vector spaces** | Score projection | Cut space (hierarchy) | — | Cycle space (structure) |
| **Extremes** | H-ordering | Transitive (min) | Intermediate | Regular (max) |
| **Tiling decomp** | Mode B recursion | Overlap (inner δ) | Bottom/Top wiring | Apex |
| **Constants** | Triangle geometry | √2 (ratio) | π (area) | e (growth) |
| **Connection scales** | Hamming distance | Wiggly (d=1) | Waggly (d=2..m-1) | Blue/Black (d=m) |

## What Does Each Role MEAN?

### The FIXED element (Inert / Lossless / Spine / Same-End / F1)

This is what **persists under the action**. It is the invariant, the diagonal, the backbone. In the metagraph, it's the self-complementary classes that map to themselves. In tile interactions, it's tiles sharing a vertex — their fates are coupled because they physically touch. In the functor chain, it's the linear first step that changes nothing.

**Abstractly**: the FIXED element is where **the map equals the territory**. No information is gained or lost. Pure symmetry.

**Numerical association**: the prime **2**. Parity. Halving. The simplest symmetry.

### The BOUNDARY element (Ramified / Lossy / Ribs / Cross-End / F2)

This is **where things break**. The interface between fixed and free. In the metagraph, the ribs connect SC to NS — they cross the complement boundary. In tile pairs, the cross-end vertex is a relay point where one tile's upper endpoint becomes another's lower — a **handoff**, a singularity of the incidence relation. In the functor chain, it's Functor 2 where the tournament becomes the conflict graph — the nonlinear local step where all the crux resides.

**Abstractly**: the BOUNDARY is where **information is lost** (the projection from tournament to cycle structure is many-to-one) and where **curvature appears** (odd cycles are the curvature of tournament space).

**Numerical association**: the prime **3**. Three-cycles. The fundamental curvature element.

### The FREE element (Split / Crystallization / Sea / Disjoint / F3)

This is **the generic bulk**. Most things live here at large scale. In the metagraph, NS-NS edges dominate (96% at n=8). In tile pairs, disjoint pairs vastly outnumber same-end or cross-end. In the functor chain, Functor 3 is the global algebraic evaluation — irreducible, requiring knowledge of the whole graph.

**Abstractly**: the FREE element is where **structure is created**. The independence polynomial at x=2 creates the integer H from the topology of Ω. This is genuinely new information — the split into conjugate pieces (α₀, α₁, α₂, ...) each carrying a power of 2.

**Numerical association**: the prime **7** (or **5** in the solvable world). The structural prime. The one that creates genuine complexity.

## The Deep Connection: Why The Same Triple Keeps Appearing

The triple appears because the project studies a **composition of maps** that crosses a **phase boundary**:

```
COORDINATES  →  TOPOLOGY  →  ARITHMETIC
  (tiling)      (conflict graph)    (H value)
   FIXED           BOUNDARY           FREE
```

At each stage, the nature of information changes:
- **Coordinates** are explicit, lossless, algebraic. You can read off every bit.
- **Topology** is implicit, lossy, geometric. Many tournaments share the same conflict graph. The map T → Ω loses information. This is the ramification — the singular point.
- **Arithmetic** is global, crystallized, number-theoretic. The integer H(T) is a single number encoding the full independence structure. It is computationally irreducible — you cannot shortcut the evaluation.

Every other triple in the project is a **different cross-section** of this same phase transition:

- **Spine/Ribs/Sea** is this triple applied to the complement symmetry on iso classes
- **Same-End/Cross-End/Disjoint** is this triple applied to vertex incidence on tile pairs
- **Eisenstein** is this triple applied to the Frobenius action on primes
- **The triangle's three sides** are this triple's three projections onto the staircase geometry

## The Number 42

42 = 2 × 3 × 7 = FIXED × BOUNDARY × FREE.

The "triple point" where all three coexist. In the tournament world, this manifests as: the transitive tournament (FIXED, H=1) connected through metagraph ribs (BOUNDARY) to the sea of regular tournaments (FREE, H=max).

The fact that 42 is also the answer to "life, the universe, and everything" is — if one is being precise — the statement that any system complex enough to ask the question must contain all three phases: preservation, breaking, and creation.

## What This Predicts

If all triples are one triple, then:

1. **Every new decomposition we find should have three parts mapping to FIXED/BOUNDARY/FREE.** Test this on any future trichotomy.

2. **The boundary element should always be bipartite** (connecting fixed to free but never self-interacting). Verified: ribs are bipartite and triangle-free (THM-A). Cross-end tiles relay but don't self-loop.

3. **The free element should dominate at large n.** Verified: NS-NS grows to 96%. Disjoint pairs dominate. α₂ contributions grow.

4. **The fixed element should carry the exact formulas.** Verified: same-end has the clean formula c = -2^{max(1,|s₁-s₂|-1)}. Cross-end: c = 2^{s₁+s₂-2}. The diagonal is where things are computable. The bulk is where things are statistical.

5. **Phase transitions should occur at the boundary.** The frustration threshold (where tile gradient flips sign) is a boundary phenomenon — it happens when same-end interactions (FIXED, negative) overcome the linear gain. The boundary between "adding a tile helps" and "adding a tile hurts" IS the ramification point.
