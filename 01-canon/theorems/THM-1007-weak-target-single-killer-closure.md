---
id: THM-1007
title: THE WEAK-TARGET SINGLE-KILLER CLOSURE — THM-724's near-tight residual is EMPTY at M > 1/14, closed by the balance lemma ALONE (no shallow witness, no census). STATEMENT: every primitive 13-set S = C ∪ {v_f} with |C| = 12 and v_f > 13·max(C) (single-killer) satisfies M(S) > 1/14, UNCONDITIONALLY. PROOF (three lines): let t₀ be a core optimum and μ := M(C) ≥ 1/13 by LRC(13) (settled). (i) If ‖v_f t₀‖ ≥ 1/13 then M(S) ≥ min(μ, 1/13) = 1/13 > 1/14 (THM-724 Case 0). (ii) Otherwise the killer is (near-)resonant and THM-724's Lemma 1 (balance, unconditional) gives M(S) ≥ μ·v_f/(v_f + s), where s ≤ max(C) is the fastest descending binding runner. Since x ↦ v_f/(v_f + x) is decreasing and max(C) < v_f/13 (single-killer), M(S) ≥ μ·v_f/(v_f + max C) > (1/13)·v_f/(v_f + v_f/13) = (1/13)·(13/14) = 1/14. ∎ WHY THIS CLOSES THE RESIDUAL: THM-724's residual was exactly the configs where the balance value μ·v_f/(v_f+s) falls SHORT of the sharp target 14/183 — requiring s > v_f(183μ−14)/14, a "fast binding runner", for which a shallow-type witness had to be exhibited and was only closed by census (2336+3234 configs, no counterexample). At the WEAK target the condition becomes s < v_f(14μ−1), and since μ ≥ 1/13 gives 14μ−1 ≥ 1/13 while single-killer gives s ≤ max C < v_f/13, the inequality holds AUTOMATICALLY for every configuration — the residual has no members. The 7% gap between 14/183 = 0.07650 and 1/14 = 0.07143 is exactly what converts an empirical census into a two-line proof. QUANTITATIVE: the worst case (μ = 1/13, s = max C = m, v_f = 13m+1) gives balance = (1/13)(13m+1)/(14m+1) and margin exactly 1/(182(14m+1)) over 1/14 — thin but strictly positive (m = 12 ⟹ 1/30758); and μ = 1/13 forces C to be a dilated AP (prime-13 tightness, HYP-4382), which is THM-724's Case 1/2 with better bounds, so the true residual margin is strictly larger. NOTE: the proof does NOT use the covering hypothesis — it holds for every primitive single-killer 13-set, covering or not
status: PROVED (three lines from THM-724's Lemma 1 + Case 0 + LRC(13); no new machinery) + VERIFIED exactly (120 single-killer configs: μ ≥ 1/13, max C < v_f/13, balance > 1/14, ZERO violations; deep-well bound 0.072165 > 1/14 with actual M = 14/183; worst-case margin formula confirmed exactly at m = 1,5,12,25,100). Inherits THM-724's handling of the near-resonant sub-case (their stated "O(‖v_f t₀‖) shift only helps")
source: kind-pasteur-2026-07-18-S128 (cont.52; owner: run THM-724's near-tight residual at the weak target)
depends_on:
  - THM-724   # Lemma 1 (balance, unconditional), Case 0 (killer-safe), the residual this empties
  - THM-523   # the q-witness/covering reduction framing
  - LRC(≤13)  # settled by owner directive: μ = M(C) ≥ 1/13 for |C| = 12
related:
  - THM-995 (XI)(XII) — the weak-target localization that predicted this; THM-1007 executes it
  - THM-726 — the multi-killer half, whose \|P\| ≤ 8 strata are empirically 1.91× loose at the weak target (the analogous closure is the named next)
  - HYP-4382 (prime-13 tightness: μ = 1/13 ⟺ dilated AP), HYP-2566 (the covering-set hard core)
script: 04-computation/thm724_weak_target_kps_S128c52.py -> 05-knowledge/results/thm724_weak_target_kps_S128c52.out
---

# THM-1007 — the weak-target single-killer closure

## Statement

Let S = C ∪ {v_f} be a primitive 13-set, |C| = 12, with v_f > 13·max(C). Then **M(S) > 1/14**.

## Proof

Let t₀ attain μ := M(C). By LRC(13) (settled), μ ≥ 1/13.

- **Killer safe.** If ‖v_f t₀‖ ≥ 1/13, evaluating at t₀ gives M(S) ≥ min(μ, ‖v_f t₀‖) ≥ 1/13 > 1/14.
- **Killer (near-)resonant.** THM-724's Lemma 1 gives M(S) ≥ μ·v_f/(v_f + s) with
  s ≤ max(C) the fastest descending binding runner. The map x ↦ v_f/(v_f+x) decreases, and
  single-killer gives max(C) < v_f/13, so

  M(S) ≥ μ·v_f/(v_f + max C) > (1/13)·v_f/(v_f + v_f/13) = (1/13)·(13/14) = 1/14. ∎

## Why the residual vanishes

THM-724's residual = {balance value < 14/183} = {s > v_f(183μ−14)/14}. At the weak target the
requirement flips to s < v_f(14μ−1), and

> μ ≥ 1/13 ⟹ 14μ − 1 ≥ 1/13,  single-killer ⟹ s ≤ max C < v_f/13,

so s < v_f/13 ≤ v_f(14μ−1) always. **No configuration can sit in the residual.** The census
(2336 + 3234 configs) is not needed at this target.

## Quantitative margin

Worst case μ = 1/13, s = max C = m, v_f = 13m+1: balance = (1/13)(13m+1)/(14m+1), margin
exactly **1/(182(14m+1))** — verified at m = 1, 5, 12, 25, 100. (And μ = 1/13 forces a
dilated-AP core — THM-724 Case 1/2 — so the genuine residual margin is larger.)

## Named next
- **The multi-killer analogue:** run THM-726's \|P\| ≤ 8 open strata at the weak target by
  iterating the balance lemma over ≥2 killers. THM-995(XII) measured those strata at 1.91×
  the threshold empirically, so the analogous closure should be comfortable — that would
  close "covering ⟹ M > 1/14" outright and with it the equality horn.
- Lean rendering: the chain is one inequality over ℚ given Lemma 1 as an interface.
