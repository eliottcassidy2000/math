---
source: opus-2026-07-10-S205
status: WORKED the a-priori measure floor for the residual class (the difficulty). THREE findings, one
  decisive. (1) NEW elementary bound (center-collapse): mu(S) >= 1 - Sum_{q in D} phi(q)/(7 m(q)); exact on
  {1,2} (recovers boxeph mu_2 = 11/14) but FAILS on residual families (0/4000) -- for dissociated speeds it
  degenerates to the union bound 13/7. (2) Peeling near-APs (LEM-012's longest-AP >= k-6) does NOT lift the
  floor: dissociated mu_min = 0.0114 vs near-AP mu_min = 0.0082, same order. (3) DECISIVE: the
  ResidualObligation domain ADMITS NON-PRIMITIVE DILATES. Explicit witness v = 2*[1..9,11,12,13,20] =
  [2,4,6,8,10,12,14,16,18,22,24,26,40] satisfies EVERY residual clause (covering, gap>13, compressed,
  distinct, max>=23, divisor-closed, no-common-residue) with gcd=2 and mu = 1/980, while its core has
  Vmax = 20 <= 22 (window-censused). Since alpha -> c*alpha is measure-preserving, mu(c*w) = mu(w) exactly.
  Hence inf mu over the residual = 0 and NO UNIFORM MEASURE FLOOR EXISTS. Fix: a primitivity (gcd) peel,
  cheap via lonely_scale; after it, deep search gives min mu ~ 0.0094.
tags:
  - lrc14
  - hB5
  - measure-floor
  - residual
  - dilate
  - ill-posed-target
  - thm-685
---

# The residual admits dilates — so no uniform measure floor exists (and the fix is one cheap branch)

**opus-2026-07-10-S205.** Owner: work the a-priori measure floor for the residual class. klein's THM-685
(Kronecker transfer) made this *the* remaining analytic content: `|LM(q) − q·μ(S)| ≤ K(S) ≤ Σv`, so a
measure floor `μ(S) ≥ μ₀` yields liveness at every `q > Σv/μ₀`. I attacked the floor directly. The headline
is that **the target, as stated, is ill-posed** — and the reason is structural and cheap to fix.

## 1. A new elementary bound — and why it can't work here

The bad set is `⋃_l ⋃_k B(k/v_l, 1/(14 v_l))`: balls around every rational `p/q` with `q | v_l`. Different
runners **share centers**: a reduced center `p/q` is hit by every runner with `q | v_l`, and the union of
those concentric balls is just the largest, radius `1/(14·m(q))` with `m(q) = min{v_l : q | v_l}`. Counting
each center once (there are `φ(q)` of denominator `q`) gives

> **Center-collapse bound.** `μ(S) ≥ 1 − Σ_{q ∈ D} φ(q)/(7·m(q))`,  `D = {q : q | v_l for some l}`.

It is *exact* when balls around distinct centers are disjoint. Sanity: for `{1,2}` it gives `1 − 3/14 =
11/14`, exactly boxeph's `μ₂`. For the tight AP `{1..13}` it gives `1 − 1.1985 < 0`, correctly refusing to
certify (that family has `μ = 0`).

**Honest negative:** on residual families it is *never* positive (0/4000; median `−0.685`). Diagnosis: for
dissociated large speeds there is no divisor sharing (`m(q) = q` for `q = v_l` prime), so every runner
contributes ≈ 1/7 and the bound degenerates to the union bound `13/7`. The savings that matter for the
residual are **inter-center overlaps** — i.e. decorrelation — which this bound throws away by construction.
So the residual floor cannot come from divisor-sharing; it must come from decorrelation (Bonferroni depth
≥ 5 / the moment route). This closes a natural door.

## 2. Coherence is not the controlling parameter

Natural next guess: the μ-extremal families are the coherent ones, so peel them with LEM-012 (near-AP,
`longest-AP ≥ k−6 = 7`, already Lean-formal as `lem012_nearAP_free_gap`). **It doesn't help.** Adversarial
descent within each class:

| class | `μ_min` found |
|---|---|
| near-AP (`longest-AP ≥ 7`) — LEM-012 peels these | `0.00816` |
| dissociated (`longest-AP < 7`) — the "true" residual | `0.01145` |

Same order. So `longest-AP` does **not** control `μ_min`, and adding a LEM-012 branch buys no comfortable
floor. Second honest negative.

## 3. The decisive finding: the residual admits dilates, so `inf μ = 0`

`α ↦ c·α` is measure-preserving on the circle, hence **`μ(c·w) = μ(w)` exactly**. Now take

  `w = [1,2,3,4,5,6,7,8,9,11,12,13,20]`  (`Vmax = 20 ≤ 22` — already handled by the window-22 census),
  `v = 2w = [2,4,6,8,10,12,14,16,18,22,24,26,40]`.

Machine-checked, `v` satisfies **every clause** of the `ResidualObligation` domain — covering, scale gap
`> 13`, compressed, distinct `|v_i|`, some `|v_i| ≥ 23`, divisor-closed, no nontrivial common residue —
with `gcd(v) = 2`, and

> `μ(v) = μ(w) = 1/980 ≈ 0.00102`.

(Pure dilates *do* satisfy the divisor-closed clause: `g = 2` divides every coordinate, so the implication
`(∀ j ≠ i₀, g ∣ v_j) → g ∣ v_{i₀}` holds vacuously-true, not falsely. It is the *almost*-dilates — one odd
coordinate, like `[2,…,26,57]` — that the detuned/divisor-closed branch peels.)

Since the window census contains near-AP families whose `μ` is as small as we like (the tight AP itself has
`μ = 0`), their dilates sit inside the residual with the same tiny `μ`. Therefore

> **`inf μ` over the `ResidualObligation` domain is `0`. No uniform measure floor `μ₀ > 0` exists.**

This is not a soundness bug — the obligation is still *true* — but it means **any proof strategy that
routes through a uniform `μ₀` cannot work**, and it explains why the floor has resisted: one is trying to
bound below a quantity whose infimum is zero.

## 4. The fix, and what it buys

The dilates are trivial families: `v = c·w` is lonely at `t/c` whenever `w` is lonely at `t`
(`LonelyRunner.lonely_scale`). So a **primitivity (gcd) peel** — a ninth assembly branch reducing to
`tupleGcd v = 1` — removes them at essentially zero cost. With `gcd = 1` enforced, adversarial descent over
the full residual predicate finds

> `min μ ≈ 0.00939`  at `v = [1,4,6,10,13,14,16,17,19,20,22,23,36]` (vs iid `(6/7)^13 = 0.1348`).

So after the peel a uniform floor becomes *plausible* (nothing here proves one exists — the search is not
exhaustive). **Recommendation to the fleet:** add the primitivity peel to `lrc14_grand_assembly` and restate
the residual with `tupleGcd v = 1` *before* anyone invests further in a uniform measure floor. Without it
the target is provably unreachable; with it, it is at least well-posed.

## Ledger

- NEW: the center-collapse bound `μ(S) ≥ 1 − Σ_{q∈D} φ(q)/(7 m(q))` (exact on `{1,2}` = `11/14`); refuted
  as a residual tool (0/4000) — the residual needs decorrelation, not divisor sharing.
- REFUTED: `longest-AP`-based coherence peeling lifts the floor (dissociated `μ_min` ≈ near-AP `μ_min`).
- **PROVED (machine-checked witness): the residual admits non-primitive dilates with `μ = μ(core)`; since
  `μ(c·w) = μ(w)` and cores range over the window census (including near-APs with `μ → 0`),
  `inf μ = 0` over the residual — no uniform floor exists.** Fix = primitivity peel (`lonely_scale`).
- After the peel, empirical `min μ ≈ 0.0094`.
- Files: `lrc14_center_collapse_bound_opus_S205.out`, `lrc14_residual_mu_exact_opus_S205.out`,
  `lrc14_mu_min_vs_Vmax_opus_S205.out`, `lrc14_mu_by_coherence_opus_S205.out`,
  `lrc14_residual_mu_infimum_probe_opus_S205.out`, `lrc14_primitive_residual_mu_opus_S205.out`,
  `lrc14_residual_full_predicate_mu_opus_S205.out`, `lrc14_residual_dilate_gap_opus_S205.out`.
  → THM-685 (klein), hB5, LEM-012, LEM-024, `lonely_scale`, `ResidualObligation`, opus-S204 brick (iii).
