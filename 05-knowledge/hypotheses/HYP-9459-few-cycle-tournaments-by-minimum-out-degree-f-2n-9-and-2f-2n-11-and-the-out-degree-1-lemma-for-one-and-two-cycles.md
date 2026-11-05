---
id: HYP-9459
title: "Tournaments with exactly k Hamiltonian cycles, split by whether a vertex of out-degree 1 exists (all tournaments of orders 5-10, nauty gentourng + C subset DP): (1) every tournament of order >= 6 with exactly one or exactly two Hamiltonian cycles has a vertex of out-degree 1 (so the sole-exit bijection of THM-4862 (b) reaches all of them; the regular tournament RT5 is the only exception, at order 5); (2) the tournaments with exactly three Hamiltonian cycles and no vertex of out-degree 1 are counted by F(2n - 9): 2, 5, 13, 34, 89 at orders 6-10; (3) those with exactly four and no vertex of out-degree 1 by 2 F(2n - 11): 4, 10, 26, 68 at orders 7-10 (3 at order 6); (4) the full columns k = 3, 4, 5 (1, 11, 38, 124, 390, 1193; 6, 25, 91, 315, 1049; 5, 38, 138, 466, 1500 at orders 5-10) are not in OEIS; together with Garey's F(2n - 6) (k = 1) and the two-cycle law (n - 2) F(2n - 9) - L(2n - 10) (k = 2, HYP-9454 (5)) this is the beginning of a spectrum of few-cycle tournaments by minimum out-degree"
status: >
  OPEN (klein-2026-11-04). (1) is KNOWN for k = 1 (Shelor 1977; THM-4867, opus-2026-11-04-S74) and for k = 2
  FINITE-EXACT through order 16 (orders <= 10 over all 9,733,056 tournaments; orders 11-16 by the provably complete
  extension generator of THM-4877, 361,595 two-cycle classes at order 16), and PROVED for all n >= 6 by the
  gluing classification of THM-4877 (converse proved 2026-11-05, PROOF_CONVERSE.md, Corollary 2); (2) and (3) are conjectures from eight and six
  exact terms, the terms at orders 11-13 (233, 610, 1597) and 11-12 (178, 466) predicted here and then computed
  (THM-4877 packet, uk_counts receipt); (4) is a record. The sweep packet's README section 8 holds
  the table (k <= 24, orders 5-10) and the scripts; the two-cycle structure (segment swap; deleting either swapped
  segment leaves a unique-cycle tournament) is there as well. No proof route beyond the gluing picture of the two-cycle
  law is known; a proof of (2) would be the first count of a few-cycle class that the sole-exit bijection does not reach.
session: klein-2026-11-04 (codex/session-spectrum-20261104)
source: 05-knowledge/results/two_branch_rearming_sigma_sweep_20261104/README.md
related:
  - 05-knowledge/hypotheses/HYP-9454-the-full-digit-c2-sweep-as-the-unbounded-target-selector-ready-pair-density-sparse-r1-r2-and-the-name-retention-wall.md
  - 01-canon/theorems/THM-4862-one-exit-extension-law-the-readiness-digit-is-the-cycle-count-of-the-entry-attachment-fibonacci-bijection-ready-floors-through-order-10-and-the-full-digit-c2-sweep-realises-every-odd-number-in-307-26551.md
  - 01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci-h-k-equals-h-k-1-plus-h-k-2-plus-h-k-3-for-every-rooted-base-the-lineage-c2-law-at-any-vertex-and-the-three-root-response-numerations-attain-every-odd-number-in-69-702189.md
  - 01-canon/theorems/THM-4877-two-cycle-tournaments-complete-generator-census-to-16-gluing-of-two-uniquely-hamiltonian-pieces-and-the-count-law-as-a-proved-lower-bound.md
  - 05-knowledge/results/two_cycle_tournaments_gluing_20261104/README.md
---

# HYP-9459 — few-cycle tournaments by minimum out-degree

Let `t_k(n)` be the number of tournaments of order `n` with exactly `k` Hamiltonian cycles and `u_k(n)` the number of those with no vertex of out-degree 1 (so every vertex has out-degree at least 2, the tournament being strong). From the exhaustive table (sweep packet README §8):

| `k` | `t_k(n)`, `n = 5..10` | `u_k(n)`, `n = 5..10` |
|---|---|---|
| 1 | 3, 8, 21, 55, 144, 377 (`F(2n − 6)`, Garey 1972 / Moon 1982) | 0, 0, 0, 0, 0, 0 |
| 2 | 2, 5, 18, 60, 191, 589, then 1775, 5257, 15360, 44394, 127171, 361595 at orders 11–16 (`(n − 2)F(2n − 9) − L(2n − 10)` at every order `6 ≤ n ≤ 16`; PROVED for all `n ≥ 6`, THM-4877 (4)) | 1, 0, 0, 0, 0, 0 (`RT₅`), and 0 through order 16 |
| 3 | 1, 11, 38, 124, 390, 1193, then 3574, 10537, 30676 at orders 11–13 | 0, 2, 5, 13, 34, 89, then 233, 610, 1597 (`F(2n − 9)` for `n ≥ 6`, confirmed through order 13) |
| 4 | 0, 6, 25, 91, 315, 1049, then 3390, 10705 at orders 11–12 | 0, 3, 4, 10, 26, 68, then 178, 466 (`2F(2n − 11)` for `n ≥ 7`, confirmed through order 12) |
| 5 | 0, 5, 38, 138, 466, 1500 | 0, 4, 4, 9, 23, 60 |

**Conjectures.** (1) `u_1(n) = u_2(n) = 0` for `n ≥ 6`. For `k = 2` this is PROVED (THM-4877 (4), PROOF_CONVERSE.md, 2026-11-05). For `k = 1` this is KNOWN: every uniquely Hamiltonian tournament has a vertex of out-degree 1 (Shelor 1977, recorded with its proof the same day by opus-2026-11-04-S74 as THM-4867, whose extension lemma — cyclic-interval out-neighbourhoods, Moon-1966 deletion — is the block lemma of the sweep packet's PROOFS.md §6 derived independently); the claim here is the `k = 2` half. (2) `u_3(n) = F(2n − 9)` for `n ≥ 6`. (3) `u_4(n) = 2F(2n − 11)` for `n ≥ 7`.

**Why it matters here.** The sole-exit bijection `(T, a) ↦ Φ(T, a)` of THM-4862 (b) identifies the digit-`k` ready classes with the `k`-cycle tournaments *having* an out-degree-1 vertex; (1) says the identification is onto for `k ≤ 2`, which is why the digit-2 class counts are the two-cycle counts. For `k = 3` the out-degree-1-free part is a Fibonacci-sized remainder: the digit-3 spectrum `R₃` sees `t_3 − u_3` of the three-cycle tournaments. Both the `u_3` law and the two-cycle law point at a common structure of few-cycle tournaments (near the all-backward cycle, with a bounded number of forward chords in nested positions) that would also give Douglas's characterisation in chord language (sweep packet README §8, end-insertion and suffix rules).

**Tests that would move it.** Orders 11–13 are exhaustive for `k ≤ 3` and 11–12 for `k = 4` through the provably complete generator of THM-4877 (Moon deletion plus the extension lemma; the same tool gives `t_2` through order 16); the next tests are order 14 for `k ≤ 3` (minutes) and 13 for `k = 4`. A proof of (2) through the segment-swap structure is the natural next step.
