# THM-011: General Block-Counting Formula

**Type:** Theorem (derived, needs full verification for n >= 6)
**Certainty:** 3 -- VERIFIED for small n, formula derived
**Status:** CONJECTURE (pending full verification)
**Added by:** opus-2026-03-05-S1
**Source:** file.txt (inbox contribution)
**Tags:** #block-counting #3-cycle #contiguous-block #formula

---

## Statement

For a tournament T with a directed 3-cycle C = (u -> w -> x -> u) and complement vertex set W = V \ {u, w, x}, the number of Hamiltonian paths of T where {u, w, x} form a contiguous block is:

    H_C^+(T) = sum_{Q in Ham(T[W])} f_C(Q)

where for each Ham path Q = (Q_1, ..., Q_m) of T[W]:

    f_C(Q) = 3 + sum_{k=2}^{m} |p(Q_k)| - sum_{k=1}^{m-1} p(Q_k) * sigma(p(Q_{k+1}))

Here:
- p(z) = (T(z,u), T(z,w), T(z,x)) is the "pattern" of z relative to C
- |p(z)| = number of vertices in C that z beats
- sigma is the cyclic permutation matching the 3-cycle orientation

---

## Key Properties

1. f_C(Q) >= 1 for all Q (conjectured, verified for small n)
2. For n = 4 (m = 1, single vertex z): f_C(Q) = 3 (recovers THM-010)
3. The "discordance" disc(Q,C) = sum p(Q_k) * sigma(p(Q_{k+1})) measures how consecutive pairs in Q interact with cycle C

---

## Relationship to OCF

H_C^+(T) counts only paths where the 3-cycle vertices are contiguous. The full Claim A requires accounting for ALL Hamiltonian paths, including those where C's vertices are interleaved with other vertices. The gap between H_C^+(T) and the Claim A contribution of C is the main difficulty.
