# Wiggly Lines and the Twisted Bundle

**Session:** kind-pasteur-2026-03-23-S20eo/S20ep
**Date:** 2026-03-23

## The Fiber Bundle Picture

The tiling space decomposes as a direct product:

```
{0,1}^{C(n,2)} = Boundary × Overlap = {0,1}^{2n-3} × {0,1}^{C(n-2,2)}
```

- **Boundary** (2n-3 bits): bottom wiring (vertex 0 to inner), top wiring (vertex n-1 to inner), apex (0 to n-1)
- **Overlap** (C(n-2,2) bits): sub-tournament on inner vertices {1,...,n-2}

**Wiggly lines** = single tile flips within the overlap fiber. Each tiling has C(n-2,2) wiggly neighbors. The fiber of 2^{C(n-2,2)} tilings forms a hypercube.

## The Three Key Discoveries

### 1. The Isotropy Theorem

Every metagraph edge of G_n/Z_2 is reachable by EVERY type of tile flip: overlap, bottom, top, or apex. The translucent (overlap) subgraph = the opaque (boundary) subgraph = the full metagraph, as GRAPHS.

The only difference is line multiplicity: each edge gets translucent/total ratio exactly C(n-2,2)/C(n,2) = (n-2)(n-3)/[n(n-1)]. This ratio is the SAME for every edge (perfect uniformity).

**Proof:** S_n symmetry. Any arc flip witnessing an edge can be carried to any other arc position by a permutation.

### 2. The Non-Factoring Theorem

The fiber map φ_B: overlap → G_n/Z_2 does NOT factor through the overlap iso class map. Two isomorphic overlap sub-tournaments, placed in the same boundary frame, can produce different full-tournament iso classes.

**Factoring rates:**
- n=4: 37.5% of fibers factor
- n=5: 6.2%
- n=6: 1.6%
- Trend: → 0% as n → ∞

This means: the overlap's LABELED structure matters, not just its iso class. The boundary "sees" the specific labeling of the inner vertices.

### 3. The 2×2 Mixing Table

|  | G_n edge | G_n self-loop |
|--|----------|--------------|
| **Sub edge** (overlap class changes) | 96.2% | 3.8% |
| **Sub self-loop** (overlap class preserved) | 78.6% | 21.4% |

(Values at n=6)

The astonishing entry: **sub SL → G_n edge = 78.6%**. Two overlaps in the SAME G_{n-2} iso class, with the SAME boundary, give DIFFERENT G_n classes 78.6% of the time. The boundary *amplifies* intra-class relabeling into class-level distinctions.

The sub-edge → G_n SL entry (3.8%) shows the reverse: the boundary can *absorb* an inter-class change. But this is rare and getting rarer.

## The Twisted Bundle Interpretation

The fiber bundle Tiling → G_n/Z_2 is **genuinely twisted**:

- The **structure group** acts on the fiber by overlap relabeling (the subgroup of S_n fixing the boundary)
- The **transition functions** depend on the boundary → this is the twist
- The **obstruction to triviality** is the non-factoring rate (→100%)

In a trivial bundle, knowing the boundary and the overlap ISO CLASS would suffice to determine the full class. In reality, you need the full labeled overlap.

## The Forgetting Map

Each G_n class C decomposes over G_{n-2} classes:

- The transitive class (cn=0) uses only 1 overlap class (minimum entropy)
- Generic classes (|Aut|=1) use ALL V(n-2) overlap classes (maximum entropy)
- The dominant overlap class has ~40-67% of tilings (not uniform but not concentrated)

This is the "forgetful functor" from G_n to G_{n-2}: forget the boundary, remember only the overlap. It's many-to-many: one G_n class can contain multiple G_{n-2} classes, and one G_{n-2} class contributes to multiple G_n classes.

## What This Means for the Mode B Recursion

**G_n cannot be recovered from G_{n-2} × boundary data alone.** The bundle is twisted, meaning the "transition functions" (how fiber maps change as you vary the boundary) encode genuinely new information that grows with n.

However:
1. All G_n EDGES are visible in every fiber (isotropy)
2. The mixing rates (sub edge → G_n edge) are HIGH and increasing (→100%)
3. The forgetting map has bounded entropy per class

This suggests: G_n is "almost" determined by G_{n-2} plus a small twist correction. The twist correction is the obstruction to the Mode B recursion being a simple product. Quantifying this twist may yield the missing formula for E(G_n).

## Connection to Earlier Framework

- **Spine (SC-SC)**: The transitive class uses 1 overlap class (rigid spine)
- **Sea (NS-NS)**: Generic classes use all overlap classes (fluid sea)
- **Ribs (SC-NS)**: Intermediate — partial overlap class usage
- **The principal line**: Traces the path where overlap entropy increases from 0 (transitive) to maximum (regular)

The wiggly line structure is PERPENDICULAR to the earlier blue/black decomposition:
- Blue/black distinguishes grid-symmetric vs non-grid-symmetric tilings
- Wiggly/opaque distinguishes overlap vs boundary tile flips
- These two decompositions are INDEPENDENT (by the isotropy theorem)
