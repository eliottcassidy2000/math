---
id: THM-4559
title: "Conditional on the tensor classification of finite Euclidean Ramsey sets (openai/math #172), every finite spherical set with algebraic coordinates is Ramsey: the separability idempotent of a number field turns the sphere equation into the required tensor certificate (Graham's spherical conjecture holds for algebraic configurations; Palvolgyi's non-Ramsey heptagon needs a transcendental radius)"
status: "CONDITIONAL on openai/math #172 (unrefereed preprint, 2026-09; its release ships a Lean development of the classification, whose comparator statement ComparatorChallenges/EuclideanRamsey.lean was read here and matches the paper; the development was not built here). The deduction from #172 is PROVED (three lines). Exact certificates checked by the session's reader for two configurations. Not found stated in #172 or in Palvolgyi's abstract; literature search not exhaustive. Independent audit: see the results note."
session: mac-mini-2026-10-06-oaimath2 (found by the session's reader of openai/math #172; re-derived here)
source: 05-knowledge/results/oai2_openai_math_second_reading_20261006.md
related:
  - openai/math #172, A classification of finite Euclidean Ramsey configurations (2026-09-23), Theorem 1.1
  - D. Palvolgyi, A cyclic non-Ramsey heptagon, arXiv:2609.23327 (2026): non-Ramsey for every transcendental r > 2
  - Erdos-Graham-Montgomery-Rothschild-Spencer-Straus (1973): Ramsey implies spherical
  - Kriz (1991), Frankl-Rodl (1990), Leader-Russell-Walters (2012): classical positive results and the subtransitive conjecture
  - THM-440, THM-431, THM-4558 (the unit-distance configurations this settles)
---

# THM-4559 — algebraic spherical sets are Ramsey (conditional on #172)

## The criterion (#172, Theorem 1.1; CONDITIONAL)

Let `A = {a_1, …, a_s} ⊂ R^d` (`s ≥ 2`) span `R^d` affinely. Put `p_i = (1, a_i)` and `F = Q(coordinates)`, and let `m : F ⊗_Q F → F` be multiplication.

`A` is Ramsey **iff** some `P ∈ Mat_(d+1)(F ⊗_Q F)` satisfies:
* `(p_i ⊗ 1)^T P (1 ⊗ p_i) = 0` for every `i`;
* `m(P_(αβ)) = δ_(αβ)` on the spatial block `1 ≤ α, β ≤ d`.

## Theorem (PROVED from the criterion)

If `A` is spherical and its coordinates are algebraic, then `A` is Ramsey.

*Proof.*
1. Since `F` is a number field, `F/Q` is finite separable. So `F ⊗_Q F ≅ F × (other factors)`, with `m` the projection onto the first factor.
2. Let `e` be the corresponding idempotent. Then `m(e) = 1` and `e·(x ⊗ 1 − 1 ⊗ x) = 0` for every `x`, hence `e·(x ⊗ y) = e·(xy ⊗ 1)`.
3. Sphericity gives `H ∈ Mat_(d+1)(F)` with spatial block `I_d` and `p_i^T H p_i = 0`. This is the sphere equation, solvable over `F` because the center is determined by linear equations over `F`.
4. Put `P = e·(H ⊗ 1)`. Then `(p_i ⊗ 1)^T P (1 ⊗ p_i) = e·Σ (p_(iα) H_(αβ) ⊗ p_(iβ)) = e·(p_i^T H p_i ⊗ 1) = 0`, and `m(P_(αβ)) = H_(αβ) = δ_(αβ)`. ∎

**Algebraic squared distances suffice.** A configuration with algebraic squared distances has a congruent copy with real-algebraic coordinates (Cholesky over the real algebraic numbers). The criterion does not depend on the representative.

## Consequences (CONDITIONAL on #172)

* **Graham's spherical conjecture holds for algebraic configurations.** Every non-Ramsey spherical set is transcendental.
  * This matches Pálvölgyi's heptagon (non-Ramsey for every transcendental radius `r > 2`) and #172's own negative examples (nine algebraically independent parameters; a Liouville angle).
  * Both use derivations or algebraic independence, and derivations vanish on number fields.
* **Our configurations.**
  * Every spherical configuration in the repo's unit-distance work (Eisenstein, Moser field, Heegner fields, `P_7` heptagon subsets) is Ramsey.
  * The THM-431/THM-440 extremal configurations are not spherical, so they are not Ramsey (unconditional, EGMRSS 1973).
  * This has no bearing on `u(22)`.
* **Check (FINITE-EXACT, session reader).** The certificate `P = e·(H ⊗ 1)` was verified exactly in `Q[x, y]/(f(x), f(y))` for two sets:
  * the 18 points `ζ_6^j {1, ω_3, ω̄_3}`;
  * the twelve points "three squares rotated by `±arccos(5/6)`", the algebraic twin of #172's Liouville example.
  * In both, the evaluations do not vanish without `e`.

## Not claimed

* Nothing about transcendental spherical sets beyond #172 itself.
* Nothing unconditional: if #172's sufficiency direction fails, this theorem fails with it.
* Whether 6 generic concyclic points are Ramsey is decidable by the criterion but not settled here.
  * At most 5 concyclic points are Ramsey (#172).
  * Generic 7-point sets are not, per the appendix of Pálvölgyi's paper (a claim he marks as unchecked).
  * Any superset of a non-Ramsey set is non-Ramsey.
