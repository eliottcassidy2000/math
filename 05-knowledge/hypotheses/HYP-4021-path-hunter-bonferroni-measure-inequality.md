---
id: HYP-4021
title: THE PATH-HUNTER (BONFERRONI) MEASURE INEQUALITY -- the combinatorial heart of the LRC(14) 7-wall crossing, formalized sorry-free. The union bound mu(union D_i) <= sum mu(D_i) DIES at j=7 (each danger arc has measure 2*(1/14)=1/7, so seven tile the circle). Hunter's inequality (2nd-order path-Bonferroni) recovers a pair credit: mu(union_{i<n} A_i) + sum_{i=1}^{n-1} mu(A_i cap A_{i-1}) <= sum_{i<n} mu(A_i). PROVED (LRCHunterLedger.lean, general Mathlib measure, additive form = NO ENNReal subtraction) by pure disjointification induction: measure_union_add_inter (mu(S cup T)+mu(S cap T)=mu S+mu T) + monotonicity (A_{i-1} subset of union_{j<i}, so the path term A_i cap A_{i-1} is dominated by the running-intersection term). Plus the LEDGER COEFFICIENT: 1 - c/7 + (c-1)/49 = (48-6c)/49, strictly positive for every block size c <= 7 (48-6*7=6>0, boundary at c=8) -- the exact arithmetic of the 7-wall crossing designed in kps-S22 (LRCFatBlockChain, HYP-3979). This is the input kps's path-Bonferroni ledger needs; combined with the pair-floor mu(I cap D_i cap D_{i+1}) >= |I|/49 - err (mac-mini's JointRateCore per-cell obligation) it gives good >= |I|*(48-6c)/49 - fees > 0 for c<=7, crossing past the union-bound wall. Absent from the repo before (LRCBonferroniMeasure had only the PAIRWISE Bonferroni mu(A cap B)>=mu A+mu B-1)
status: VERIFIED (Lean, sorry-free, foundational-axioms-only; #print axioms path_hunter_add_le = [propext, Classical.choice, Quot.sound]; registered in the root module, corpus green). path_hunter_add_le proved for a general MeasureTheory.Measure by disjointification induction (no subtraction, only measurable A_i); ledger_coeff + ledger_coeff_pos by ring/linarith. HONEST: this is the COMBINATORIAL / measure-theoretic heart of the 7-wall crossing -- the general Hunter inequality + the ledger arithmetic. It does NOT close the crux by itself: the ANALYTIC pair-floor (mu(I cap D_i cap D_{i+1}) >= |I|/49 - err for near-equal runners, the joint-rate per-cell obligation) is mac-mini's active work, and wiring {Hunter + singles-upper-bound + pair-floor} into the concrete danger-region (goodRegion/teeth) framework is the remaining assembly. Provides the reusable Hunter lemma the ledger consumes.
source: klein-2026-07-02-S116
depends_on:
  - HYP-4020   # S115: the far-count-7 dispatch (the >=7-far leg hge7 that this ledger targets)
  - HYP-3979   # kps-S22: the 7-wall crossing DESIGN (path-Bonferroni ledger + pair-floor interface)
related:
  - HYP-3874   # mac-mini JointRateCore: the pair-floor / per-cell obligation
  - HYP-2832   # LRCBonferroniMeasure: the pairwise Bonferroni (this is the multi-set generalization)
  - OPEN-Q-108 # the j=7 union-bound wall
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCHunterLedger.lean
---

# HYP-4021 — the path-Hunter (Bonferroni) measure inequality

## The inequality (LRCHunterLedger.lean, sorry-free)
For a general measure `μ` and measurable `A : ℕ → Set α`:
> `μ (⋃ i ∈ range n, A i) + Σ_{i ∈ Ico 1 n} μ (A i ∩ A (i-1))  ≤  Σ_{i ∈ range n} μ (A i)`.

This is **Hunter's inequality** — the second-order Bonferroni bound along a spanning **path** — in
**additive form** (no ENNReal subtraction). Proof: disjointification induction. The step uses
`measure_union_add_inter` (`μ(S∪T) + μ(S∩T) = μ S + μ T`) and monotonicity: since `A_{i-1} ⊆ ⋃_{j<i} A_j`,
the *path* term `μ(A_i ∩ A_{i-1})` is dominated by the *running-intersection* term `μ((⋃_{j<i})∩A_i)`, which
`measure_union_add_inter` converts exactly. No measurability beyond each `A_i`.

## The ledger coefficient (the 7-wall arithmetic)
`ledger_coeff : 1 − c/7 + (c−1)/49 = (48 − 6c)/49`, and `ledger_coeff_pos : c ≤ 7 → 0 < (48−6c)/49`. At
`c = 7` the credit is `6/49 > 0`; boundary at `c = 8` (`= 0`). This is exactly the arithmetic of kps-S22's
designed crossing.

## Why this crosses the 7-wall
The union bound gives `good ≥ |I| − Σ μ(I∩D_i) ≥ |I|(1 − c/7) − fees`, which is `≤ 0` for `c ≥ 7` (seven
danger arcs of measure `1/7` tile). Hunter adds the pair credit: `good ≥ |I| − [Σ μ(I∩D_i) − Σ_path μ(I∩D_i∩
D_{i+1})]`, and with the singles bound (`μ(I∩D_i) ≤ |I|/7 + fee`) and the **pair-floor** (`μ(I∩D_i∩D_{i+1}) ≥
|I|/49 − err`, mac-mini's per-cell obligation), this is `≥ |I|(1 − c/7 + (c−1)/49) − fees = |I|(48−6c)/49 −
fees > 0` for `c ≤ 7`. So near-equal blocks of up to 7 runners are lonely — past the wall.

## Honest scope
This formalizes the **combinatorial / measure-theoretic heart** — the general Hunter inequality (new to the
repo; `LRCBonferroniMeasure` had only the pairwise case) plus the ledger arithmetic. It does **not** close the
crux alone: the **analytic pair-floor** (near-equal runners' danger overlap `≥ |I|/49 − err`) is mac-mini's
active JointRateCore per-cell work, and wiring `{Hunter + singles + pair-floor}` into the concrete
`goodRegion`/teeth framework is the remaining assembly. This provides the reusable Hunter lemma that ledger
consumes.
