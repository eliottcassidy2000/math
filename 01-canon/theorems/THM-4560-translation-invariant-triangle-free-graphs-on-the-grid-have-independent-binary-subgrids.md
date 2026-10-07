---
id: THM-4560
title: "Every translation-invariant (gap-determined) triangle-free graph on N^n has an independent binary subgrid, for every n: a chain of minimal idempotent ultrafilters q_1 <= ... <= q_n on the level semigroups chooses the gaps of any prescribed tree shape so that every contiguous gap sum avoids the sum-free gap set; hence t_dead(Finv) < infinity for every n, and the strong-Specker barrier of THM-521 D holds unconditionally for fully gap-determined witnesses (Erdos 592; row-invariant witnesses are not covered)"
status: "PROVED (idempotent-ultrafilter argument; Hindman-Strauss Ch. 1-5 facts cited). Found by the session's reader of openai/math #164 (Hindman's FS+FP conjecture), whose own proof uses no idempotents; the method here is the classical Galvin-Glazer one, and the chain-of-minimal-idempotents device is standard in Milliken-Taylor / variable-word proofs (e.g. Bergelson-Blass-Hindman 1994). FINITE-EXACT companions: Q_gap(2,3) SAT, Q_gap(2,4) UNSAT (this session, the reader, and the audit: two solvers plus exhaustive enumeration of all 2^24 gap sets); Q_gap(3,t) SAT for t = 4, 5, 6 (witnesses brute-verified twice, independently). INDEPENDENTLY AUDITED 2026-10-06 (audit B: every step of the proof checked; PASS WITH CORRECTIONS, applied: the scope of the THM-521 D corollary, notation, wording)."
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
* **`Q_gap(n,t)`** is the finite game for gap-determined rules on `[t]^n`: the rule must be triangle-free and hit every binary subgrid. This is THM-470's Finv instance.
  * **Warning.** It is *not* THM-453 G's `invQ(n,t)`, which is the *row-invariant* game (`R_a = R`, `B_(a,a') = B_(a'−a)`, arbitrary column relations). See MISTAKE-577.
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
* **THM-521 D (the strong-Specker barrier) holds unconditionally for fully gap-determined witnesses**, including every valuation grading of the gap vector, without HYP-2396.
  * For *row-invariant* witnesses (THM-453 F/G, including the dyadic `B_(v₂ g)` family, whose `n = 2` cutoff is `5 = 2n+1`) the barrier is not covered here and remains conditional/open.
  * A strong witness, if one exists, is not fully gap-determined. Value-dependent features such as Larson's partial sums are one possibility; HYP-2558 (the strong-Specker barrier entry) stays open.

**Corollary 3 (the role of the seam).** Sum-freeness keeps `E` out of every idempotent ultrafilter (the exact criterion is that `E` contains no IP set). So the 2-adic seam (THM-469), one source of sum-free gradings, can delay the death of a gap-determined witness but never prevent it.

## Finite data (FINITE-EXACT)

* **`n = 2`, fully gap-determined:** `Q_gap(2,3)` SAT, `Q_gap(2,4)` UNSAT, so the cutoff is `4 = 2n`.
  * Confirmed by CaDiCaL here, by the reader, and by the audit (CaDiCaL and Glucose, plus exhaustive enumeration: 75 winning `E` of `2^12` at `t = 3`, none of `2^24` at `t = 4`).
  * The gap-determined `(sign, v₂)` algebra also has cutoff 4 (audit, exhaustive).
  * The *row-invariant* game of THM-453 F/G has cutoff 5 (SAT at 4, UNSAT at 5, recomputed here and by the audit). These are different games (MISTAKE-577).
* **`n = 3`, fully gap-determined:** `Q_gap(3,t)` is SAT for `t = 4, 5, 6`, with `|E| = 34, 70, 109`.
  * This was already implied by THM-470 B (F2J SAT at `(3,4..6)`, A3), so these runs are confirmations.
  * Two independent runs produced different witnesses.
  * Every witness was brute-checked twice, independently, including all `1.7·10^8` binary subgrids at `t = 6`.
  * **`Q_gap(3,7)` is UNSAT (2026-10-07), so `t_dead(Finv) = 7` at `n = 3`.**
    * How it was certified UNSAT: a CEGAR loop that also adds the 7 coordinate-reflected copies of each subgrid clause returned UNSAT after 467 iterations (CaDiCaL 1.9.5, 2,864 s, deterministic rerun); the dumped 669,086-clause set was re-solved UNSAT from scratch by MapleChrono (1,668 s, reader) and by CaDiCaL 1.5.3 (1,379 s, this session, `04-computation/experiments/oai2_20261006_finv37_resolve.py`); every clause was validated as a genuine constraint (340,299 realizable triangles, 328,787 binary-subgrid clauses) by the reader's audit script and independently by `04-computation/experiments/oai2_20261006_finv37_validate.py`; no DRAT proof.
    * The certificate is in `04-computation/experiments/oai2_20261006_readers/hindman_ramsey/finv_unsat_n3_t7.{cnf,leaves}.gz`.
    * Earlier runs without the reflection images, and the repo's 2-hour run, had timed out.
* **Summary.** The fully gap-determined game dies at `t = 3, 4, 7` for `n = 1, 2, 3`. This matches `2n+1` at `n = 1, 3` but not at `n = 2` (NUMEROLOGY for now). It says nothing about the free game `Q(3,7)` (HYP-2396).

## Not claimed

* No explicit bound on `t_dead(Finv)`. The argument is non-constructive (ultrafilters); a Milliken–Taylor-style finite version would give one.
* Nothing about value-dependent witnesses, the free game `Q(3,7)`, or `R(n,2)` itself (HYP-2396).
