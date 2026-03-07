# THM-081: Walsh Spectrum of Directed Cycle Counts

**Status:** PROVED (analytical derivation + exhaustive verification at n=5)
**Author:** opus-2026-03-07-S35c8
**Date:** 2026-03-07
**Dependencies:** OCF (THM-002)

## Main Theorem

For a directed k-cycle count t_k(T) on n-vertex tournaments:

**hat{t_k}[S] = (1/2^k) * sum_{directed k-cycles C : S ⊂ edges(C)} (-1)^{asc(S,C)}**

where asc(S,C) counts the number of edges in S that cycle C traverses from smaller to larger vertex.

### Proof

A single directed k-cycle C has indicator:

I_C(T) = prod_{(u,v) in C} A[u][v]

where A[u][v] = 1 if u->v is an arc. In Walsh variables:
- A[i][j] = 1/2 - chi_{ij}/2 for i < j (ascending edge)
- A[j][i] = 1/2 + chi_{ij}/2 for i < j (descending edge)

The Walsh coefficient of I_C at monomial S:

hat{I_C}[S] = E_T[I_C(T) * chi_S(T)]

**If S contains an edge not in C:** The edge variable is independent of I_C, and E[chi_e] = 0, so hat{I_C}[S] = 0.

**If S ⊂ edges(C):** Factor the expectation over edges:
- Edges NOT in S: E[A_{dir}] = 1/2, contributing (1/2)^{k-|S|}
- Edges in S, traversed ascending (i->j, i<j): E[A[i][j] * (-1)^{T_e}] = -1/2
- Edges in S, traversed descending (j->i, i<j): E[A[j][i] * (-1)^{T_e}] = +1/2

Product: hat{I_C}[S] = (1/2)^{k-|S|} * prod_{e in S} (-1)^{asc_e(C)} * (1/2)
= (1/2)^k * (-1)^{asc(S,C)}

Summing over all directed k-cycles: hat{t_k}[S] = sum_C hat{I_C}[S]. QED.

## Properties

1. **Support:** hat{t_k}[S] = 0 unless S is a subset of some k-cycle's edge set.
2. **Degree range:** Nonzero coefficients have 0 ≤ |S| ≤ k.
3. **Uniform amplitude per cycle:** Each cycle contributes exactly ±1/2^k.
4. **Sign rule:** Two directed cycles using the same edges but opposite orientation on S give opposite signs. Two cycles with same orientation on S give the same sign.
5. **Complement behavior:** t_k is complement-invariant (even Walsh support), since reversing all edges maps each directed cycle to the reverse cycle. Degree-parity: only even |S| survive (see below).

### Parity of Support

For a k-cycle, the number of ascending edges has parity depending on cycle direction. Swapping direction flips all k ascents, changing asc count by k (mod 2). For odd k, the two directions of a cycle differ in asc count by an odd number, so their contributions to even-|S| monomials add up (same sign) and to odd-|S| cancel. This gives t_k even Walsh support for odd k.

## Walsh-Domain OCF

Combining with OCF (H = 1 + 2*alpha_1 + 4*alpha_2 + ...):

For S ≠ ∅:
**hat{H}[S] = sum_{k odd} (2/2^k) * sum_{C: S⊂edges(C)} (-1)^{asc(S,C)} + 4*hat{alpha_2}[S] + ...**

At n=5 (where alpha_2 = 0):
hat{H}[S] = (1/4) * sum_{3-cycles C ⊃ S} (-1)^{asc(S,C)} + (1/16) * sum_{5-cycles C ⊃ S} (-1)^{asc(S,C)}

## Verification

**n=5, t3:** 1024 tournaments, all 2^10 monomials. Max error: 0.0 (EXACT).
**n=5, t5:** 1024 tournaments, all 2^10 monomials. Max error: 0.0 (EXACT).

## Degree-2 Specialization

For |S| = 2 (pair of edges {e1, e2} forming a P2):

hat{t3}[{e1,e2}] = (1/8) * sum_{3-cycles C containing e1,e2} (-1)^{asc({e1,e2},C)}

A P2 pair within a triangle is contained in exactly 2 directed 3-cycles. Both cycles traverse the P2 with the SAME number of ascents on {e1,e2} (since the third edge flips direction). So:

hat{t3}[P2 in triangle] = 2 * (±1)/8 = ±1/4

This matches THM-077 for H at degree 2.

## Connection to THM-077 (H Walsh Formula)

THM-077 gives: hat{H}[S] = (-1)^{asc(S)} * 2^r * (n-|S|)!/2^{n-1}

THM-081 gives: hat{H}[S] = sum_k (2/2^k) * sum_{C⊃S} (-1)^{asc(S,C)}

Equating these provides a **counting identity**:

sum_k (1/2^{k-1}) * sum_{C⊃S} (-1)^{asc(S,C)} = (-1)^{asc(S)} * 2^r * (n-|S|)!/2^{n-1}

This relates the number of directed cycles containing a given path pattern to factorials. Proving this identity algebraically would give a new proof of OCF.
