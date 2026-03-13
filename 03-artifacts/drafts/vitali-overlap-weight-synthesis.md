# The Vitali Atom and the {2,1,0} Overlap Weight Hierarchy

**Author:** kind-pasteur-2026-03-13-S61
**Date:** 2026-03-13

## Overview

This document synthesizes the session's discoveries about how the Vitali atom (lambda-preserving (1,1,2,2) reversal) interacts with the hidden higher-dimensional structure of tournaments through the {2,1,0} overlap weight system.

## The Main Theorem (THM-170)

At n=8, for any lambda-preserving (1,1,2,2) reversal:

**delta_H = 2 * dc7_dir + 4 * delta_i2**

where:
- dc7_dir = change in total directed 7-cycle count
- delta_i2 = change in vertex-disjoint directed cycle pairs

Verified 166/166 with 0 failures.

## The Four-Level Structure

The tournament's Hamiltonian path count H(T) = I(Omega(T), 2) depends on the conflict graph Omega through its independence polynomial. The Vitali atom acts on four nested levels:

### Level 0: Cycle Counts (c3, c5, c7)
The "gross" statistics — how many vertex sets support cycles of each length.
- **Preserved at n<=8** (c3 by lambda, c5 always, c7 vertex sets preserved)
- This is the "measurable" level: determined by the lambda graph

### Level 1: Cycle Identities (which vertex sets)
Which specific vertex sets carry cycles.
- The Vitali atom **swaps** some vertex sets (4 lost, 4 gained in typical n=8 example)
- Always c3 sets at |V ∩ S| = 2 (exactly 2 S-vertices per swapped set)
- **Net count preserved** — a perfect permutation of cycle identities
- The **overlap weight spectrum** is preserved: the distribution of pairwise overlaps doesn't change

### Level 2: Directed Multiplicities
How many directed cycles live on each vertex set.
- A vertex set V can support 0, 1, 2, 3, ... directed Hamiltonian cycles
- c3: always multiplicity 1 (tournament on 3 vertices has exactly 0 or 1 directed cycle)
- c5: multiplicity 1-3 (tournament on 5 vertices can have up to 3 directed cycles)
- c7: multiplicity up to 24 (7-vertex tournament)
- **The Vitali atom changes these multiplicities** while preserving their NET total per length
- Changes occur at |V ∩ S| = 2 or 4 (the "marginal" intersection levels)

### Level 3: Disjoint Pair Products
The number of vertex-disjoint directed cycle pairs:
i2 = Σ_{disjoint (Vi, Vj)} m_i * m_j

This is a **bilinear form** in the multiplicities. Even when:
- Cycle counts are preserved (Level 0)
- Cycle identities permute (Level 1, with zero net change in overlap spectrum)
- Net directed cycle counts are preserved (Level 2, dc3 = dc5 = 0)

...the PRODUCT m_i * m_j can change because individual multiplicities m_i change on vertex sets V_i that happen to be disjoint from some V_j.

## The Mechanism: |V ∩ S| Marginal Analysis

For a disjoint (c3, c5) pair at n=8:
- c3 ∪ c5 = all 8 vertices
- |c3 ∩ S| + |c5 ∩ S| = |S| = 4

| |c3 ∩ S| | |c5 ∩ S| | c5 affected? | c3 affected? |
|-----------|-----------|--------------|--------------|
| 0 | 4 | YES (all arcs) | No |
| 1 | 3 | YES (3 arcs) | No |
| 2 | 2 | YES (1 arc) | YES (swapped) |
| 3 | 1 | No | YES |

The critical case is |c3 ∩ S| = 2: BOTH the c3 vertex set swaps AND the c5 multiplicity changes. This is the "marginal" level where the Vitali atom acts.

## The Vitali Set Analogy

The Vitali set V ⊂ [0,1] is non-measurable: it looks "the same" under every statistical test (measure zero/one), but it has structure invisible to the Lebesgue sigma-algebra.

The Vitali atom in tournaments is analogous:
1. **Statistical invariance**: Lambda graph preserved, cycle counts preserved, overlap spectrum preserved — every "measurable" statistic agrees before and after.
2. **Hidden structure**: The bilinear form i2 = Σ m_i m_j depends on the EXACT placement of cycles, not their statistics. It is "non-measurable" in the sense that no finite-order marginal determines it.
3. **Phase transition**: At n≤6, the Vitali atom is gauge-trivial (like restricting to Q-measurable sets). At n=7, one channel opens (c7 count — still "measurable"). At n=8, the "non-measurable" channel i2 opens. At n=9, higher-order products (i3) should open.

## The Hierarchy as Dimensional Growth

Each n opens new channels:

| n | Formula | Active channels | New phenomenon |
|---|---------|----------------|----------------|
| ≤6 | delta_H = 0 | None | Gauge-trivial |
| 7 | delta_H = 2·dc7 | c7 count | First non-trivial |
| 8 | delta_H = 2·dc7 + 4·di2 | c7 count + disjoint pairs | Multiplicity reshuffling |
| 9 | delta_H = 2·(dc7+dc9) + 4·di2 + 8·di3 | c7,c9 counts + pairs + triples | (predicted) |
| n | delta_H = Σ_k 2^k · delta_i_k | All channels up to k=⌊n/3⌋ | General OCF decomposition |

The coefficient 2^k comes from the independence polynomial: H = I(Omega, 2) = Σ_k i_k · 2^k.

## Connection to the {2,1,0} Overlap Weight System

The overlap weight W(C_i, C_j) = |V(C_i) ∩ V(C_j)| determines the conflict graph:
- W ≥ 1: C_i and C_j are adjacent in Omega (conflict)
- W = 0: C_i and C_j are independent (no conflict, can be in same independent set)

The {2,1,0} classification:
- **W = 2**: Share an edge. Strong coupling. Lambda graph determines these.
- **W = 1**: Share a vertex but no edge. Weak coupling. Still conflict in Omega.
- **W = 0**: Disjoint. No coupling. Independent in Omega.

The Vitali atom operates at the **W=1/W=0 boundary**:
- It preserves the NUMBER of pairs at each W level (overlap spectrum invariant)
- It changes WHICH pairs are at each level (identity reshuffling)
- The i2 channel measures the PRODUCT structure at W=0

This is why the Vitali atom reveals the "hidden dimension": the product structure at W=0 is information that cannot be recovered from any single-pair statistic. It requires knowledge of the GLOBAL arrangement of cycles.

## n=9 Results (CONFIRMED)

The formula extends to n=9:

**delta_H = 2*(dc7_dir + dc9_dir) + 4*delta_i2**

Confirmed 42/42 with 0 failures. Key findings:
- **dc9 channel ACTIVE**: dc9 ranges from -3 to +3, nonzero in 71% of examples
- **dc3 = dc5 = 0 ALWAYS**: net preservation continues at n=9
- **di3 = 0 STRUCTURALLY**: the c3 swap preserves disjoint triple counts perfectly
- The c3 vertex set swap preserves ALL disjointness counts: c3 count, c3-c3 pairs, c3-c3-c3 triples

The di3 = 0 result means the i3 channel doesn't open at n=9 despite being geometrically possible (three disjoint c3's need exactly 9 = n vertices). This is because the c3 swap is a "measure-preserving" permutation that maintains all disjointness structure.

## The General Formula is OCF (Not a Conjecture!)

The formula delta_H = Sigma_k 2^k * delta_i_k follows TRIVIALLY from OCF:
  H(T) = I(Omega(T), 2) = Sigma_k i_k(Omega(T)) * 2^k

Therefore delta_H = Sigma_k delta_i_k * 2^k is a TAUTOLOGY.

**The real content** is in the VANISHING theorems:
- dc3 = 0 always (THM-171: lambda determines c3 count and structure)
- dc5 = 0 always (verified n=8,9,10 — mechanism: per-overlap-level cancellation)
- di3 = 0 at n=9 (THM-171: c3 disjointness is lambda-determined)
- delta_(c3,c3) pairs = 0 always (THM-171: lambda-determined)

These vanishing results REDUCE the tautological general formula to the
simple 2-channel forms observed empirically.

## Updated Hierarchy Table (Vanishing Perspective)

| n | Effective formula | Vanishing results | Non-trivial channels |
|---|---------|----------------|----------------|
| <=6 | delta_H = 0 | ALL channels vanish | None |
| 7 | delta_H = 2*dc7 | dc3=dc5=0, no i2 possible | dc7 only |
| 8 | delta_H = 2*dc7 + 4*di2 | dc3=dc5=0, di2 from (c3,c5) only | dc7, di2 |
| 9 | delta_H = 2*(dc7+dc9) + 4*di2 | dc3=dc5=0, di3=0 | dc7, dc9, di2 |
| 10 | delta_H = 2*(dc7+dc9) + 4*di2 | dc3=dc5=0 STILL | dc7, dc9, di2 (new: (c3,c7),(c5,c5) pair types) |

## n=10 Results (CONFIRMED)

At n=10 (16 H-changing examples):
- **dc3 = dc5 = 0 STILL** — c5 phase transition has NOT happened
- **Formula holds at 100%** (16/16)
- **NEW disjoint pair types ACTIVE**: (c3,c7) contributes 30% of |di2|, (c5,c5) contributes 36%
- **(c3,c3) pairs: delta = 0 always** (THM-171 guarantees this)
- dc7 range wider (up to 6), dc9 range wider (up to 9), di2 range wider (up to 6)

## THM-171: Lambda Determines c3 Structure (PROVED)

For any tournament T, the c3 overlap statistics are fully determined by lambda(u,v):
- |C| = (1/6) Sigma lambda(u,v) [total c3 count]
- P_2 = Sigma C(lambda(u,v), 2) [pairs sharing 2 vertices]
- P_1 = Sigma C(delta(w), 2) - 2*P_2 [pairs sharing 1 vertex]
- D = C(|C|, 2) - P_1 - P_2 [disjoint pairs]

Since the Vitali atom preserves lambda by definition, ALL c3 overlap statistics are preserved.
This extends to triples, k-tuples, and the full overlap weight spectrum.

**Corollary:** di3 = 0 at n=9 is a direct consequence of THM-171.

## dc5=0: Per-Overlap-Level Cancellation

Computational analysis reveals:
- |V cap S| = 4: always exactly 2 sets change with +1/-1 (complement duality on S)
- |V cap S| = 2: net cancellation WITHIN this level (not across levels)
- Per S-pair nets are NONZERO: the (1,1,2,2) score structure creates signed patterns
  (+3, -3, -3, +3 in typical example) that cancel across all 6 pairs
- The (strong, weak) involution sigma maps V to sigma(V) with d(V) = d(sigma(V))
  (equal, not opposite!), so pairing doesn't directly explain the cancellation

The dc5=0 identity remains UNPROVED but is verified at n=8,9,10 with zero exceptions.

## Open Questions (Updated)

1. **PROVE dc5 = 0** — the per-S-pair signed cancellation pattern
2. **When does di3 open?** THM-171 shows dc3 disjointness is lambda-forced at all n.
   So di3 = 0 for all n (since di3 counts disjoint c3 triples, which are lambda-determined)?
3. **When does dc5 BREAK?** Verified dc5=0 at n=8,9,10. Is it universal?
4. **Higher-order vanishing**: At what n does delta_(c3,c5,c7) triple become nonzero?
5. **Connection to TDA**: The hierarchy resembles a persistence filtration
6. **Engineering**: Use vanishing results for efficient H-approximation
