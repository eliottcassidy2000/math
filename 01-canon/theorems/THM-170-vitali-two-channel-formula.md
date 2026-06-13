# THM-170: Vitali Atom Two-Channel Formula at n=8

**Status:** CONFIRMED (166/166 verified, 0 failures)
**Author:** kind-pasteur-2026-03-13-S61
**Date:** 2026-03-13

## Statement

For any tournament T on n=8 vertices, let S be a 4-subset with score sequence (1,1,2,2), and let T' be obtained by reversing all arcs within S. If the lambda graph is preserved (lambda(T) = lambda(T')), then:

**delta_H = 2 * delta_c7_dir + 4 * delta_i2**

where:
- delta_c7_dir = (# directed 7-cycles in T') - (# directed 7-cycles in T)
- delta_i2 = (# vertex-disjoint directed cycle pairs in T') - (same in T)

## Key Structural Facts

1. **Net c3/c5 directed changes are zero:** The reversal reshuffles WHICH 3-cycles and 5-cycles exist (swapping some vertex sets), but the NET directed cycle count at each length is preserved for c3 and c5.

2. **Two independent channels:**
   - **c7 channel:** delta_c7_dir in {-3, -2, -1, 0, 1, 2, 3}. The 7-cycle vertex set uses 7 of 8 vertices, so every other cycle conflicts with it. Adding/removing a c7 directed cycle only changes i_1 (alpha_1), contributing 2 per cycle.
   - **i2 channel:** delta_i2 in {-2, -1, 0, 1, 2}. The reshuffling of c3/c5 cycle identities changes which directed cycle pairs are vertex-disjoint. Each such pair contributes 4 to H via I(Omega, 2).

3. **Hierarchy from n=7 to n=8:**
   - At n=7: delta_H = 2 * delta_c7_dir (only c7 channel, delta_i2 = 0 always)
   - At n=8: Both channels active, independently varying
   - Prediction for n=9: delta_H = 2*(dc7+dc9) + 4*di2 + 8*di3

## Proof Sketch

**Step 1: Clique Perturbation.** All directed cycles on the SAME vertex set form a clique in Omega (they share all vertices). Changing multiplicity m_A -> m_B on a vertex set V is equivalent to resizing a clique K_{m_A} -> K_{m_B}.

**Step 2: c7 Universality.** A 7-cycle uses 7 of 8 vertices. Any other cycle (on 3, 5, or 7 vertices) shares at least 2 vertices with this set. Therefore the c7 clique is UNIVERSALLY adjacent in Omega. Resizing it changes only i_1, contributing 2 per added/removed vertex.

**Step 3: c3/c5 Net Cancellation.** The (1,1,2,2) reversal preserves lambda, which determines c3 vertex set counts (and c5 counts, verified invariant through n=8). The reversal swaps some 3-cycle vertex sets (4 lost, 4 gained in typical example) but preserves the total. The net directed cycle change at c3 and c5 is always zero.

**Step 4: i2 Channel.** The swapped vertex sets change which specific pairs of directed cycles are vertex-disjoint. Since c3 (3 vertices) + c5 (5 vertices) = 8 = n, a (3,5) disjoint pair uses ALL n vertices. Whether a specific 3-cycle is disjoint from a specific 5-cycle depends on their vertex sets. Swapping changes this count.

**Step 5: Higher-Order Terms Vanish.** At n=8, vertex-disjoint triples require at least 9 vertices, so i_3 = 0. This means the formula is exact with just the two terms.

## Evidence

| Property | Value |
|----------|-------|
| Tested examples | 166 |
| Failures | 0 |
| delta_c7_dir range | {-3,...,3} |
| delta_i2 range | {-2,...,2} |
| delta_H range | {-6,...,6} |

## Connection to OCF

This decomposition follows from H(T) = I(Omega(T), 2) where I is the independence polynomial evaluated at x=2. The coefficients are:
- i_0 = 1 (always)
- i_1 = total directed odd cycles
- i_2 = vertex-disjoint directed cycle pairs
- i_k = 0 for k >= 3 at n=8 (need 9+ vertices)

So H = 1 + 2*i_1 + 4*i_2, and delta_H = 2*delta_i_1 + 4*delta_i2.

Since delta_i_1 = delta_c3_dir + delta_c5_dir + delta_c7_dir = 0 + 0 + delta_c7_dir, we get the formula.

## Related

- THM-169: Vitali atom characterization (WHEN H changes)
- THM-002: OCF (H = I(Omega, 2))
- HYP-830: n=8 two-channel formula (this theorem)
- HYP-831: c3/c5 net zero reshuffling
- HYP-832: Overlap weight spectrum preservation
