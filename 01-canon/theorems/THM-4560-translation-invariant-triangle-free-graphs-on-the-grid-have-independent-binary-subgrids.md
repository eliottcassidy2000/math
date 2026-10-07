---
id: THM-4560
title: "Every translation-invariant (gap-determined) triangle-free graph on N^n has an independent binary subgrid, for every n: a chain of minimal idempotent ultrafilters q_1 <= ... <= q_n on the level semigroups chooses the gaps of any prescribed tree shape so that every contiguous gap sum avoids the sum-free gap set; hence t_dead(Finv) < infinity for every n, and the strong-Specker barrier of THM-521 D holds unconditionally (Erdos 592)"
status: "PROVED (idempotent-ultrafilter argument; Hindman-Strauss Ch. 1-4 facts cited). Found by the session's reader of openai/math #164 (Hindman's FS+FP conjecture), whose own proof uses no idempotents; the method here is the classical Galvin-Glazer one. FINITE-EXACT companions: fully gap-determined Q_inv(2,3) SAT, Q_inv(2,4) UNSAT (independent CaDiCaL run here); Q_inv(3,t) SAT for t = 4, 5, 6 (reader's two runs, witnesses brute-verified). Independent audit: see the results note."
session: mac-mini-2026-10-06-oaimath2
source: 05-knowledge/results/oai2_openai_math_second_reading_20261006.md
scripts:
  - 04-computation/experiments/oai2_20261006_fields_ramsey_checks.py (+ .out; part 3)
related:
  - THM-453 (the tree-grid reframe; binary subgrids, Q(n,t), strong witnesses)
  - THM-470 (feature algebras; König equivalence A2; coarsening collapse A3; MISTAKE-577 for the Finv label)
  - THM-469 (sum-free gradings), THM-521 D (the strong-Specker barrier), HYP-2396 (R(n,2) = 2n+1, still open), HYP-2558
  - N. Hindman, D. Strauss, Algebra in the Stone-Cech Compactification (2nd ed.), Thm 1.60 (minimal idempotents below a given idempotent), Cor 4.18 (closures of ideals), Thm 5.8 (Galvin-Glazer)
---

# THM-4560 — invariant triangle-free grids always leave a binary subgrid independent

## Setting (THM-453 C/D, THM-470)

* `N^n` is ordered lexicographically. `P_n` is the set of lex-positive vectors of `Z^n`, and `lev(v)` is the index of the first nonzero coordinate.
* A graph on `N^n` is **gap-determined** if `x ~ y` (`x <_lex y`) iff `y − x ∈ E`, for a fixed `E ⊆ P_n`.
* It is triangle-free iff `E` is **sum-free**: there are no `d₁, d₂ ∈ E` (equal allowed) with `d₁ + d₂ ∈ E`.
* A **binary subgrid** is the 2^n-leaf set of a height-`n` tree that picks two children at every node (THM-453 C).

## Theorem (PROVED)

Let `E ⊆ P_n` be sum-free and `w ∈ [n]^k` any word. Then there are `δ_1, …, δ_k ∈ P_n` with `lev(δ_i) = w_i` such that **every contiguous sum `δ_i + … + δ_j` avoids `E`.**

*Proof.*
1. **Level semigroups.** Let `S_ℓ` be the vectors of level `≥ ℓ` (a subsemigroup of `P_n`) and `P^(ℓ)` those of level exactly `ℓ`, a two-sided ideal of `S_ℓ`. Its closure is an ideal of `βS_ℓ`, so it contains the smallest ideal `K(βS_ℓ)`.
2. **A chain of idempotents.**
   * Take a minimal idempotent `q_n` of `βS_n`.
   * For `ℓ = n−1, …, 1`, `q_(ℓ+1)` is an idempotent of `βS_ℓ ⊇ βS_(ℓ+1)`. Take a minimal idempotent `q_ℓ ≤ q_(ℓ+1)` of `βS_ℓ` (H–S Thm 1.60).
   * Then `P^(ℓ) ∈ q_ℓ`, and `q_a + q_b = q_b + q_a = q_(min(a,b))` because the chain is `≤`-ordered.
3. **The complement is large.** A sum-free set lies in no idempotent: by Galvin–Glazer, every member of an idempotent contains some `x, y, x + y`. So `F := P_n ∖ E ∈ q_ℓ` for every `ℓ`.
4. **The induction.** Choose the `δ`'s left to right, keeping `C_0 = F` and `C_(j+1) = F ∩ (−δ_(j+1) + C_j)`. So `C_j = {y ∈ F : δ_i + … + δ_j + y ∈ F for all i ≤ j}`.
   * *Invariant:* `C_j ∈ q_b` for every `b = min(w_(j+1), …, w_m)`, `m > j`. These are the levels of the sums that will start at `j + 1`.
   * *Step.* Let `a = w_(j+1)` and let `b′` range over the levels of sums starting at `j + 2`.
     * If `b′ ≥ a`, then `q_a + q_(b′) = q_a`, and `C_j ∈ q_a` by the invariant.
     * If `b′ < a`, then `q_a + q_(b′) = q_(b′)`, and `C_j ∈ q_(b′)` because `b′ = min(a, b′)` is itself an invariant level.
   * Either way `{x : −x + C_j ∈ q_(b′)} ∈ q_a`. Hence `C_j ∩ P^(a) ∩ ⋂_(b′) {x : −x + C_j ∈ q_(b′)} ∈ q_a`, and any `δ_(j+1)` in it keeps the invariant. ∎

**Corollary 1 (binary subgrids).**
* Take the ruler word `w_i = n − v_2(i)`, `1 ≤ i < 2^n`, and the prefix sums `s_0 = (M, …, M)`, `s_i = s_0 + δ_1 + … + δ_i`, with `M` large.
* Two prefix sums `s_i, s_(i′)` differ by a contiguous sum whose level is the first bit in which `i` and `i′` differ, counted from the top. So the `s_i` are the leaves of a binary subgrid.
* All their pairwise gaps are contiguous sums, so the subgrid is **independent**.
* Base-`b` ruler words give `b`-ary subgrids; any finite tree shape works.

**Corollary 2 (Erdős 592).**
* No gap-determined strong witness exists on `N^n`, for any `n`.
* By König (THM-470 A2), **`t_dead(Finv) < ∞` for every `n`**. By the coarsening collapse (THM-470 A3), every gap-determined feature algebra (`(sign, v₂)`, jets, cross-gap, leading digit, …) dies at a finite `t`.
* **THM-521 D (the strong-Specker barrier) holds without HYP-2396 and without "invariant witnesses = valuation gradings".** A strong witness, if one exists, must use value-dependent features (Larson partial sums; HYP-2558 stays open).

**Corollary 3 (the role of the seam).** Sum-freeness is exactly what keeps `E` out of every idempotent ultrafilter. So the 2-adic seam (THM-469) can delay the death of an invariant witness but never prevent it.

## Finite data (FINITE-EXACT)

* **`n = 2`, fully gap-determined:** `Q_inv(2,3)` SAT, `Q_inv(2,4)` UNSAT, so the cutoff is 4 (CaDiCaL here; the reader agrees).
  * The *row-invariant* game of THM-453 F/G has cutoff 5 (SAT at 4, UNSAT at 5, recomputed here).
  * These are different games (MISTAKE-577).
* **`n = 3`, fully gap-determined:** `Q_inv(3,t)` is SAT for `t = 4, 5, 6`, with `|E| = 34, 70, 109`.
  * Two independent runs produced different witnesses.
  * Every witness was brute-checked, including all `1.7·10^8` binary subgrids at `t = 6`.
  * `(3,7)` is undecided. By the theorem, `t_dead(Finv)` at `n = 3` lies in `[7, ∞)` and is finite.

## Not claimed

* No explicit bound on `t_dead(Finv)`. The argument is non-constructive (ultrafilters); a Milliken–Taylor-style finite version would give one.
* Nothing about value-dependent witnesses, the free game `Q(3,7)`, or `R(n,2)` itself (HYP-2396).
