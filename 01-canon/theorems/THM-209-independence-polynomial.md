# THM-209: H = Independence Polynomial of Odd-Cycle Disjointness Graph at x=2

**Status:** VERIFIED (exhaustive n=3..7, sampled n=8: 5K, n=9: 5K, n=10: 200), PROOF NEEDED
**Found by:** opus-2026-03-14-S89b
**Verified in:** `04-computation/crown_jewel_verify_89b.py`, `04-computation/ip_verify_n7n8.c`, `04-computation/ip_sample_n9_large_89b.py`, `04-computation/ip_sample_n10_fixed_89b.py`

## Statement

For a tournament T on n vertices, let OddCyc(T) denote the set of all directed odd cycles in T (as vertex sets). Define the **odd-cycle disjointness graph** G(T):
- Vertices = directed odd cycles in T
- Edges = cycles sharing at least one vertex

Then:

**H(T) = IP(G(T), 2)**

where IP(G, x) = Σ_{k≥0} i_k(G) · x^k is the **independence polynomial** of G, and i_k(G) is the number of independent sets of size k in G (i.e., collections of k pairwise vertex-disjoint odd cycles).

Equivalently:

**H(T) = Σ_{S ⊆ OddCyc(T), pairwise disjoint} 2^|S|**

## Explicit Formulas by n

The formula simplifies because vertex-disjoint cycle collections are limited by n:

| n | Maximum |S| | Formula | New terms vs n-1 |
|---|---------|---------|-----------|
| 3 | 1 | H = 1 + 2·t₃ | t₃ |
| 4 | 1 | H = 1 + 2·t₃ | (same, t₅ impossible) |
| 5 | 1 | H = 1 + 2·t₃ + 2·t₅ | t₅ |
| 6 | 2 | H = 1 + 2·t₃ + 2·t₅ + 4·d₃₃ | d₃₃ (disjoint 3-pairs, need 6 verts) |
| 7 | 2 | H = 1 + 2·(t₃+t₅+t₇) + 4·d₃₃ | t₇ |
| 8 | 2 | H = 1 + 2·(t₃+t₅+t₇) + 4·(d₃₃+d₃₅) | d₃₅ (need 8 verts) |
| 9 | 3 | H = 1 + 2·Σt_{2k+1} + 4·Σd_{ij} + 8·d₃₃₃ | t₉, d₃₃₃ (need 9 verts) |

where:
- t_k = number of directed k-cycles
- d_{ij} = number of vertex-disjoint (i-cycle, j-cycle) pairs
- d₃₃₃ = number of vertex-disjoint 3-cycle triples

## Coefficient Structure

The coefficients follow a **power-of-2 hierarchy**:

| Level k | Coefficient 2^k | Meaning |
|---------|-----------------|---------|
| 0 | 1 | Empty set (base Hamiltonian path) |
| 1 | 2 | Each individual odd cycle |
| 2 | 4 | Each disjoint pair of odd cycles |
| 3 | 8 | Each disjoint triple of odd cycles |

## Key Properties

1. **H is always odd** (sum of powers of 2, with the 2^0=1 term always present).
2. **H ≥ 1** always (Rédei's theorem: every tournament has at least one Hamiltonian path).
3. **H = 1 ⟺ T is transitive** (no odd cycles ⟺ only the empty set contributes).
4. **The formula is a polynomial in cycle counts** with integer coefficients that are powers of 2.

## Connection to Rédei's Theorem

This formula provides a **constructive explanation** of Rédei's theorem: H(T) ≡ 1 (mod 2) because every term 2^k with k ≥ 1 is even, so H = 1 + (even) ≡ 1 (mod 2).

More precisely, H(T) ≡ 1 (mod 2) follows immediately from the independence polynomial structure, since all terms except the empty set contribute even values.

## Connection to OCF (Odd Cycle Collection Formula)

The Odd-Cycle Collection Formula expresses H in terms of independent sets of odd cycles in the cycle-clique hypergraph. THM-209 gives the EXACT weight: each independent set of k pairwise vertex-disjoint odd cycles contributes exactly 2^k.

## Verification Data

| n | Tournaments | All match? | Method | New terms tested |
|---|-------------|------------|--------|------------------|
| 3 | 8 | YES | Exhaustive | t₃ |
| 4 | 64 | YES | Exhaustive | - |
| 5 | 1,024 | YES | Exhaustive | t₅ |
| 6 | 32,768 | YES | Exhaustive | d₃₃ (level 2) |
| 7 | 2,097,152 | YES | Exhaustive (C) | t₇ |
| 8 | 268M | 5000/5000 | Random sample | d₃₅ |
| 9 | 69B | 5000/5000 | Random sample | t₉, **d₃₃₃ (level 3!)** |
| 10 | 35T | 200/200 | Random sample | d₃₇, d₅₅ |

## Open Questions

1. **Proof?** The formula is exact for all tested cases (n=3 through n=10). What is the proof mechanism? The connection to the OCF and the multilinear expansion of H on the tournament hypercube may provide the key.
2. **Connection to chromatic symmetric function?** Via Mitrovic-Stojadinovic, the chromatic function of the incomparability graph equals ω applied to the Rédei-Berge symmetric function. Does IP(G,2) arise from an e-positivity/h-positivity statement?
3. **Asymptotic behavior?** As n→∞, how does the independence polynomial grow?
4. **Connection to statistical mechanics?** H = Z(G, 2) is the hard-core lattice gas partition function at fugacity 2. What thermodynamic quantities arise from this?
