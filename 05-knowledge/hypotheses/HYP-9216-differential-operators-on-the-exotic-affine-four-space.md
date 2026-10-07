---
id: HYP-9216
title: "Differential operators on openai/math #047's exotic affine 4-space: is D(X) isomorphic to the Weyl algebra A_4? (X x A^1 = A^5 but X is not A^4, per #047; D(X) tensor A_1 = A_5 and gr D(X) = C^[8] hold, and K-theory and Hochschild homology agree with A_4)"
status: >
  OPEN question (no conjectured answer). CONDITIONAL on openai/math #047 (unrefereed): X = Spec A, A[w] = C^[5], A not C^[4].
  PROVED from that: projective O(X)-modules are free (Quillen-Suslin on A[w]), so T*X = A^8 as a variety,
  D(X) tensor A_1 = A_5, gr D(X) = C^[8]; K-theory and Hochschild homology of D(X) agree with A_4 (H_dR(X) = C).
source: mac-mini-2026-10-06-oaimath2, 05-knowledge/results/oai2_openai_math_second_reading_20261006.md (section on openai/math #047)
related:
  - THM-1300 (the explicit Dixmier counterexample on A_3); THM-4562 (slice lemma)
  - PROBLEM-LEDGER section A (Jacobian / Dixmier ledger lines A2, A3)
---

# HYP-9216 — does the exotic 4-space have the Weyl algebra as its ring of differential operators?

**Setting (CONDITIONAL on #047).**
* `A = C[p, s, u, F, J]/(H)`, with `x = s² + u³ + p²F` and `H = x²F − (1 + 2sx)J − p²J² − pu`.
* `X = Spec A` satisfies `X × A¹ ≅ A⁵` and `X ≇ A⁴`.
* The non-isomorphism comes from a Derksen-type invariant, `DK(A) ⊆ F_0 A`, proved with Mason–Stothers on a (3,2)-cusp.

**Question.** Is `D(X) ≅ A_4`?

**What is known.**
* `D(X) ⊗ A_1 ≅ D(X × A¹) = A_5`.
* `gr D(X) ≅ O(T*X) ≅ C^[8]`, because projective modules over `A` are free.
* Algebraic K-theory and Hochschild homology cannot tell `D(X)` from `A_4`.

**Either answer is interesting.**
* **No:** cancelling a Weyl tensor factor would fail. That is a noncommutative cancellation counterexample built from a commutative one.
* **Yes:** two non-isomorphic smooth affine varieties would have isomorphic rings of differential operators (and `D(X)` would be a Weyl algebra on a non-affine-space variety).

**Heuristic (session reader).** Through the mod-`p` centres (Tsuchimoto; Belov-Kanel–Kontsevich), "yes" forces `(T*X, ω)` to be symplectomorphic to the standard `A⁸` mod `p` for almost all `p`. The repo's Dixmier work (THM-1300's explicit transfer) is the natural toolbox.
