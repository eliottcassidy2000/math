# E[H(T)] = n! / 2^{n-1}

**Session:** opus-2026-04-05-S24
**Status:** THEOREM (proved)

## Statement

For a uniformly random labeled tournament T on n ≥ 1 vertices:

**E[H(T)] = n! / 2^{n-1}**

Equivalently: the total sum of H over all labeled tournaments is

**W(n) = Σ_T H(T) = n! × 2^{C(n-1, 2)}**

## Proof

By **linearity of expectation**. There are n! labeled Hamiltonian paths (directed paths visiting all n vertices). Each path P is present in T iff all n-1 arcs of P point in the right direction, which happens with probability (1/2)^{n-1} independently. Therefore:

E[H(T)] = E[#{paths P in T}] = Σ_P P(P ⊂ T) = n! × (1/2)^{n-1} = n!/2^{n-1}

And W(n) = E[H] × 2^{C(n,2)} = n!/2^{n-1} × 2^{n(n-1)/2} = n! × 2^{n(n-1)/2 - n + 1} = n! × 2^{(n-1)(n-2)/2} = n! × 2^{C(n-1,2)}.  □

## Verification

| n | n!/2^{n-1} | W(n) = n! × 2^{C(n-1,2)} | Verified |
|---|-----------|--------------------------|----------|
| 1 | 1         | 1                        | ✓        |
| 2 | 1         | 2                        | ✓        |
| 3 | 3/2       | 12                       | ✓        |
| 4 | 3         | 192                      | ✓        |
| 5 | 15/2      | 7,680                    | ✓        |
| 6 | 45/2      | 737,280                  | ✓        |
| 7 | 315/4     | 165,150,720              | ✓        |

## OCF Decomposition

Since H(T) = 1 + 2α₁ + 4α₂ + 8α₃ + ..., we have:

E[H] = 1 + 2·E[α₁] + 4·E[α₂] + ...

From the formula E[H] = n!/2^{n-1}, and E[α₁] computed exactly, the relation:

**n!/2^{n-1} = 1 + 2·E[α₁] + 4·E[α₂] + ...**

At n ≤ 6, the truncation to degree 2 is EXACT (Σα₃ = 0 for all labeled tournaments). This means the expected number of size-3 independent sets in Ω(T) is zero through n=6.

## E[α₁] Formula

E[α₁] = Σ_{L=3,5,7,...} C(n,L) × (L-1)! / 2^L

The L=3 term is C(n,3)/4, which dominates only for small n. For large n, longer cycles dominate: by n=8, the L=3 term is only 14% of E[α₁].

## Growth

E[H] grows like n!/2^{n-1} ≈ √(2πn) × (n/2e)^n × √2 → super-exponential.

The ratio E[H(n)]/E[H(n-1)] = n/2, an arithmetic progression.

## Connection to Previous Work

The coefficient formula a(n) = (n-2)!/2^{n-4} from S22 is the LINEAR regression coefficient in ΔH ≈ a(n)·Δc₃. The ratio a(n)/E[H(n)] = (n-2)!/2^{n-4} / (n!/2^{n-1}) = 2^3/(n(n-1)·2^4) = 1/(2n(n-1)), which goes to 0. The linear approximation captures a vanishing fraction of the total H at large n.

## Note on W(n) Definition

The project previously tracked "W(n) = 1,2,8,32,158,928,6350,49752" which is a DIFFERENT quantity (likely per-isomorphism-class or per-scoring-rule sum). The W(n) here is the raw labeled sum.
