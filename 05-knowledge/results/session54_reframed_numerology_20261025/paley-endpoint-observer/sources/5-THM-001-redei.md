# THM-001: Rédei's Theorem

**Type:** Theorem
**Certainty:** 5 — PROVED
**Status:** PROVED
**Last reviewed:** SYSTEM-2026-03-05-S1
**Disputes:** none
**Tags:** #redei #hamiltonian-paths #parity #classical

---

## Statement

For every tournament T on n vertices, H(T) ≡ 1 (mod 2).

Equivalently: every tournament has an **odd** number of directed Hamiltonian paths.

---

## Proof / Proof Sketch

The paper gives four independent proof routes:

**Route A (Q-Lemma):** Fixed-point-free involution on two-block path decompositions. Conditional on the per-3-cycle parity conjecture (verified exhaustively for n≤4). The involution itself (toggle lemma) is unconditional and is the k=2 case of Forcade's F₂-invariance.

**Route B (β-involution):** Reduction to a smaller tournament via a β-involution. Unconditional.

**Route C (automorphism action):** Free automorphism-group action. Unconditional.

**Route D (OCF):** H(T) = I(Ω(T), 2). Now unconditional — Claim A is PROVED (see CONJ-001, THM-002). Since I(Ω(T), 2) ≥ 1 (the empty independent set) and satisfies the deletion-contraction recurrence, parity follows.

The classical proof (Rédei, 1934) proceeds by induction on n. The base case n=1 is trivial (H=1). The inductive step uses vertex deletion and a careful parity argument.

## Verification Record

Verified computationally for all n ≤ 6 (exhaustive) and n = 7 (random sampling).

**Verification script:** `04-computation/verify_all_theorems.py` (THM-001 section)

## Notes & History

Rédei's 1934 result. This is the foundational theorem everything else is built on. Routes A–D in the paper give new structural explanations for why it holds.
