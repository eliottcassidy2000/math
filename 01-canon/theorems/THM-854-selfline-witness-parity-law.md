---
id: THM-854
title: THE SELF-LINE WITNESS PARITY LAW — every witness π of a black self-line (πT(t) = T(κt)) satisfies π^k X = X ⊕ [k odd]·𝟙 ⊕ (Δ_{j<k} π^j p₀); hence an ODD-order witness forces its path-orbit to be a mod-2 path-decomposition of K_n, which is possible iff n ≡ 1, 2 (mod 4); witnesses are NEVER involutions in the computed range (orders 3,5 at n=5; 4,5 at n=6; 6,10,12 at n=7) — correcting THM-852(ii)'s corollary, whose premise σ² ∈ Aut(T) fails because σ² realizes the pair-flip p₀ Δ σp₀ ≠ ∅
status: PROVED (the flip-vector identity + parity; one-paragraph proofs) + witness-type spectrum verified exhaustively n = 5, 6, 7 (referee: type ledger, weighted masses 8/12/88 reproduced)
source: opus-2026-07-15-S312 (owner directive: prove the 2selfK = SC bijection, check n=8 — the law was refuted at n=8 by kind-pasteur S128 (404 ≠ 176) and codex-S14 (audit concurs); this theorem is the surviving witness-level structure)
depends_on:
  - THM-852 (Klein-four skeleton; the n=8 refutation)
related: [THM-790 (leg law), THM-855 (F1-F5), MISTAKE-150, codex-S14 H-drift]
verification: 05-knowledge/results/selfk_type_ledger_fixed_opus_S312.out, selfk_debug_witness_opus_S312.out
---

# THM-854 — the self-line witness parity law

**Setup.** Encode a tournament X on [n] as a GF(2) vector over pair-positions
(bit per pair with a fixed orientation convention). Reversal is X ↦ X ⊕ 𝟙
(flip every pair); the complement tiling is X ↦ X ⊕ 𝟙 ⊕ p₀ = X^{p₀} (flip
every pair OFF the base path p₀). A permutation π acts affinely. A **witness**
of a self-line at tiling t is π with π·X = X^{p₀}, X = T(t).

## (1) The iteration identity

π^k X = X ⊕ [k odd]·𝟙 ⊕ (p₀ Δ πp₀ Δ ... Δ π^{k-1}p₀)   for all k ≥ 0.

*Proof.* Induction: π(X ⊕ F) = πX ⊕ πF for any flip-set F, and π𝟙 = 𝟙. ∎

In particular **σ² realizes the pair-flip p₀ Δ σp₀** for any witness σ. So
σ² ∈ Aut(X) ⟺ σp₀ = p₀ (σ preserves the base path's pair-set) — and since a
path's pair-set automorphisms are {id, reflection} and the reflection reverses
the standard orientation on path arcs (impossible: both X and X^{p₀} contain
the standard DIRECTED path), **no witness of a black self-line has σ² ∈
Aut(X) unless σp₀ = p₀ pointwise-forbidden** — the premise of THM-852(ii)'s
involution corollary never holds. (Verified: all 8+12+88 witnesses at
n = 5, 6, 7 have order > 2.)

## (2) The order equation

Setting k = ord(π) in (1): **[k odd]·𝟙 = Δ_{j<k} π^j p₀.**

- k odd ⟹ the π-orbit of the base path XORs to the COMPLETE GRAPH: a
  **mod-2 path-decomposition of K_n** by ≤ k translates of one path.
- k even ⟹ the translates XOR to ∅ (each pair covered evenly).

## (3) The parity law

Counting cardinalities mod 2 in (2): k odd requires
C(n,2) ≡ k(n−1) ≡ n−1 (mod 2), i.e. **n(n−1)/2 ≡ n−1 (mod 2)**, which holds
iff **n ≡ 1 or 2 (mod 4)**.

> **Odd-order witnesses can exist only at n ≡ 1, 2 (mod 4); at
> n ≡ 0, 3 (mod 4) every self-line witness has even order.**

## (4) The computed witness-type spectrum (exhaustive, weighted mass = selfK·2)

| n | witness cycle types (mass) | orders | n mod 4 | odd allowed? |
|---|---------------------------|--------|---------|--------------|
| 5 | (1,1,3): 4, (5,): 4       | 3, 5   | 1       | yes — and ALL odd |
| 6 | (1,5): 4, (2,4): 8        | 5, 4   | 2       | yes — mixed |
| 7 | (1,6): 32, (2,5): 48, (3,4): 8 | 6, 10, 12 | 3 | no — all even ✓ |

At n=7 every witness has exactly TWO cycles (all two-part partitions of 7
appear). No witness anywhere is an involution.

## Reading

The object 2·selfK = 8, 12, 88, 404 (n = 5..8) is NOT SC(n) (refuted at n=8);
what it counts is #(t, π) witness pairs = Σ|Aut|-weighted quasi-fixed mass.
The parity law says this mass decomposes by witness order with a mod-4
dichotomy in n — the first exact structural statement about the self-line
object that survives n=8. Odd-order witnesses are mod-2 path-decompositions
of K_n: the self-line phenomenon at n ≡ 1,2 (mod 4) is a DESIGN-THEORETIC
object (path-decomposition-like), while at n ≡ 0,3 (mod 4) it is forced into
even-order (2-adic) mechanisms. Next: type-level counting formulas; what does
selfK = 2, 3, 22, 101 (halved) count?
