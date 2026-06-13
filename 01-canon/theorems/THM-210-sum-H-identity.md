# THM-210: Sum of H over All Tournaments

**Status:** PROVED
**Found by:** opus-2026-03-14-S89b
**Verified in:** `04-computation/sum_h_proof_89b.py`

## Statement

For all n ≥ 1:

**Σ_{T ∈ T_n} H(T) = n! · 2^{m-n+1}**

where T_n is the set of all 2^m tournaments on n vertices, m = C(n,2).

Equivalently: **E[H] = n!/2^{n-1}**.

## Proof

The sum counts pairs (T, π) where T is a tournament and π is a Hamiltonian path in T:

Σ_T H(T) = Σ_T Σ_π 1_{π is a path in T}

Swap summation:

= Σ_π Σ_T 1_{π is a path in T}

For a fixed permutation π = (π(1), π(2), ..., π(n)), the path π is Hamiltonian in T iff T contains arcs π(1)→π(2), π(2)→π(3), ..., π(n-1)→π(n). This constrains n-1 of the m arcs; the remaining m-(n-1) arcs are free. So:

Σ_T 1_{π is a path in T} = 2^{m-n+1}

There are n! permutations, so:

Σ_T H(T) = n! · 2^{m-n+1}. □

## Consequence: Master Equation

Combined with THM-209 (H = IP(G(T), 2)):

n!/2^{n-1} = E[H] = 1 + 2·Σ_k E[t_{2k+1}] + 4·E[Σ d_{ij}] + 8·E[Σ d_{ijk}] + ...

where:
- E[t_k] = C(n,k) · (k-1)!/2^k = C(n,k) · (k-1)!/(2^k)
- E[d₃₃] = C(n,3)·C(n-3,3)/2 · (1/4)²

## Verification

| n | m | n!·2^{m-n+1} | Verified |
|---|---|---------------|----------|
| 3 | 3 | 12 | YES |
| 4 | 6 | 192 | YES |
| 5 | 10 | 7,680 | YES |
| 6 | 15 | 737,280 | YES |
