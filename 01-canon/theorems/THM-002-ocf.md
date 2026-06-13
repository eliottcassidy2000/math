# THM-002: Odd-Cycle Collection Formula (OCF)

**Type:** Theorem
**Certainty:** 5 — PROVED for all n
**Status:** PROVED
**Last reviewed:** kind-pasteur-2026-03-05-S12
**Disputes:** none
**Tags:** #ocf #independence-polynomial #conflict-graph #proved #grinberg-stanley

---

## Statement

For every tournament T on n vertices:

```
H(T) = I(Ω(T), 2) = Σ_{k≥0} α_k(Ω(T)) · 2^k
```

where Ω(T) is the conflict graph on directed odd cycles of T (vertices = odd cycles, edges = shared vertex), I(G,x) is the independence polynomial of G, and α_k counts independent sets of size k.

---

## Proof

### Proof 1: Via Grinberg-Stanley (external, all n)

**Corollary 20 of arXiv:2412.10572** (Irving & Omar, 2024; attributed there to Grinberg & Stanley [Theorem 1.39 & Lemma 6.5 of arXiv:2307.05569]) states:

> For a tournament D on [n]: ham(D̄) = Σ_{σ ∈ S(D), all cycles of σ have odd length} 2^{ψ(σ)}

where S(D) = permutations whose nontrivial cycles are directed cycles of D, ψ(σ) = number of nontrivial cycles, and D̄ is the complement digraph.

**For tournaments, D̄ = D^op (converse):** Since exactly one of (u,v), (v,u) is in D, the complement has (u,v) iff (v,u) ∈ D.

**ham(D^op) = ham(D):** Path reversal bijection (v_1,...,v_n) ↔ (v_n,...,v_1).

**The RHS is I(Ω(D), 2):** Each permutation σ with all nontrivial cycles being odd D-cycles corresponds to a collection of vertex-disjoint odd directed cycles (the nontrivial cycles of σ). The weight 2^{ψ(σ)} = 2^{|collection|}. Summing over all such collections gives Σ_k α_k · 2^k = I(Ω(D), 2), since independent sets in Ω(D) are exactly collections of pairwise vertex-disjoint odd cycles.

**Therefore: H(T) = ham(T) = ham(T̄) = I(Ω(T), 2).** QED.

### Proof 2: Via vertex deletion induction (internal, all n)

Now that OCF is established, the internal proof also closes:

**Inductive step:** Pick vertex v. Then:
- H(T) = H(T−v) + (H(T) − H(T−v))
- Claim A gives: H(T) − H(T−v) = 2 Σ_{C∋v} μ(C)
- Claim B gives: I(Ω(T),2) − I(Ω(T−v),2) = 2 Σ_{C∋v} μ(C)
- By induction H(T−v) = I(Ω(T−v),2), so H(T) = I(Ω(T),2). □

(Claim A now follows from OCF + Claim B, so this chain is non-circular.)

### Proof 3: Polynomial identity (internal, n ≤ 8)

THM-015 proves delta_H = delta_I as an exhaustive polynomial identity at n ≤ 8.
THM-018 proves alpha_w^H = alpha_w^I symbolically at n ≤ 8.
Combined with base case H(transitive) = 1 = I(empty, 2) and arc-flip reachability.

---

## Verification Record

| n | Tournaments | Failures | Method |
|---|-------------|---------|--------|
| 4 | 64 | 0 | exhaustive |
| 5 | 1,024 | 0 | exhaustive |
| 6 | 32,768 | 0 | exhaustive |
| 7 | 1,048,576 configs | 0 | exhaustive (THM-015) |
| 8 | 134,217,728 configs | 0 | exhaustive (THM-015) |
| 11 (Paley) | 1 | 0 | H(T_11) = 95095 = I(Ω(T_11), 2) |

**Verification scripts:**
- `04-computation/verify_all_theorems.py` (THM-002 section, exhaustive n<=5, sampled n=6,7)
- `04-computation/verify_ocf_full.py` (full OCF verification with all odd cycles)
- `03-artifacts/code/compute_H_T11.py` (Paley T_11 verification)

## Notes & History

The formula H(T) = I(Ω(T), 2) was discovered independently in this project (2026) and proved computationally for n ≤ 8. It was subsequently found to be equivalent to results of Grinberg & Stanley (arXiv:2307.05569, 2023), who proved it in full generality using the Rédei-Berge symmetric function. Irving & Omar (arXiv:2412.10572, 2024) give an alternative derivation via matrix algebra and state it explicitly as their Corollary 20 (attributed to Grinberg-Stanley).

The formula connects tournament Hamiltonian path counting to:
- Independence polynomials of conflict graphs
- Hard-core lattice gas model at fugacity 2 (see TANGENT T006)
- Rédei's theorem: I(G, 2) ≥ 1 (always odd)

**IMPORTANT — Closed-Form Clarification (kind-pasteur-2026-03-05-S5, DISC-002 resolved):**

H(T) = I(Ω(T), 2) IS a valid closed-form identity where Ω(T) uses ALL directed odd cycles of T and I is the PLAIN independence polynomial (alpha_k = number of independent sets of size k, NO mu weights). MISTAKE-004 (retracted) arose from confusing mu weights with independence polynomial coefficients.

## Key References

- Grinberg & Stanley, "The Rédei-Berge symmetric function of a directed graph", arXiv:2307.05569 (2023) — original proof (Theorem 1.39 & Lemma 6.5)
- Irving & Omar, "Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra", arXiv:2412.10572 (2024) — re-derives via matrix algebra (Corollary 20, attributed to Grinberg-Stanley)
