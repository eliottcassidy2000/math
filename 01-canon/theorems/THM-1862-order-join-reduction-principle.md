---
id: THM-1862
title: "The ORDER-JOIN REDUCTION PRINCIPLE for tournament-invariant inequalities: order-join T₁▷T₂ (T₁ beats all of T₂) makes c₃, tr, scc ADDITIVE and the Rédei Hamiltonian-path count H MULTIPLICATIVE, with strongly connected tournaments as the ▷-atoms; hence any JOIN-MONOTONE invariant inequality holds for all tournaments iff it holds on the strongly-connected core. This is the general principle of which opus THM-1865 (source/sink = single-vertex join is H- and c₃-neutral but pumps score-spread) is the single-vertex special case and kind-pasteur THM-1860 (c₃≤H by H=∏H(SCC), c₃=∑c₃(SCC), ∑≤∏) is the paradigm instance. New corollaries: the velocity/fragility predictor (a bound whose easy side is frozen c₃,H and whose hard side grows with tr/srange under D± is refuted by the THM-1830 3-cycle-atom+singletons witness), the MINIMAL repair srange ≤ tr+1 of the broken srange≤tr, and the WOWII-103 king-eccentricity straddle of e (19/7 below, 17/6 above)."
status: >
  PROVED (principle) + VERIFIED-EXHAUSTIVE (all 530 iso classes n≤7; 396 strong classes:
  1,1,6,35,353 for n=3..7) + SAMPLED (random n=8..12, 400/each). The order-join invariant
  algebra (c₃,tr,scc additive; H multiplicative; dom left-projecting) is verified on all class
  pairs n₁+n₂≤7 and each identity has a one-line structural proof. The reduction theorem is by
  induction on the ▷-factorization (every tournament = ordered join of its strong components).
  CREDIT / non-collision: opus-S438 THM-1865 independently found the source/sink inflation
  diagnostic (the single-vertex case) ~2 min after my push; kind-pasteur THM-1860 independently
  proved c₃≤H by exactly this SCC decomposition and Lean-formalized the ∑≤∏ kernel. This file
  states the GENERAL reduction principle + the order-join algebra (H multiplicative is the new
  structural fact beyond single-vertex neutrality) and cites both as instances. New here:
  srange ≤ tr+1 (minimal repair of the srange≤tr that breaks at n=7; VERIFIED n≤7, join-monotone
  n≤7 — candidate for all n) and the king-eccentricity straddle of e. Methodological caveat
  (verified): random sampling does NOT see inflation-fragility — the broken srange≤tr survives
  400 random samples at every n=8..12; the break lives on the near-transitive (low-entropy) strata.
source: boxeph-2026-07-21-S193 (owner: work on tournament-graffiti + WOWII analogues, long creative session)
depends_on: []
related:
  - THM-1855  # opus: source/sink inflation diagnostic (the single-vertex special case of the algebra below)
  - THM-1860  # kind-pasteur: c3<=H via SCC decomposition + Lean SumLeProd (the paradigm instance)
  - THM-1845  # kind-pasteur sandwich n-c3<=tr<=smax+1 (both faces reduce here)
  - THM-1850  # klein directed WOWII dom+tr<=n+1 (join-monotone, re-derived here)
  - THM-1830  # the 3-cycle-atom / transitive-singleton family = the inflation-fragility witness
  - "07-reflections/inflation-velocity-and-the-coupling-law-boxeph-S193.md"
script: 04-computation/tournament_graffiti_coupling_boxeph_S193.py, tournament_graffiti_strongcore_mine_boxeph_S193.py, tournament_graffiti_largen_sample_boxeph_S193.py (+ .out)
---

# THM-1862 — the order-join reduction principle

`A[i][j]=1` means *i beats j*. **Order-join** `T = T₁ ▷ T₂`: every vertex of `T₁` beats every
vertex of `T₂` (all cross-arcs one-way). `tr` = largest transitive subtournament (WOWII `α`-analog),
`c₃` = #3-cycles, `H` = Rédei Hamiltonian-path count, `srange = smax−smin`.

## 1. The order-join invariant algebra (verified exactly on all pairs n₁+n₂≤7)

- `c₃(T) = c₃(T₁)+c₃(T₂)` (a 3-cycle needs ≥2 one-way cross-arcs — impossible),
- `tr(T) = tr(T₁)+tr(T₂)` (concatenate the chains),
- `scc(T) = scc(T₁)+scc(T₂)`,
- **`H(T) = H(T₁)·H(T₂)`** — a directed Ham path exhausts `T₁` before entering `T₂` (no arc
  returns): the **new structural fact** beyond the single-vertex neutrality of opus THM-1865,
- `dom(T) = dom(T₁)` (a dominating set of `T₁` already beats all of `T₂`).

**Atoms:** every tournament factors uniquely as an ordered join of its strong components (transitive
condensation). Strongly connected tournaments are the `▷`-atoms.

## 2. The reduction principle

Call `Φ` **join-monotone** if `Φ(T₁) ∧ Φ(T₂) ⟹ Φ(T₁▷T₂)`.

> **THM-1862.** A join-monotone tournament-invariant inequality holds for **all** tournaments iff
> it holds for all **strongly connected** tournaments. (Induct on the ▷-factorization.)

**Instances.** kind-pasteur's `c₃ ≤ H` (THM-1860) is the paradigm case: `c₃,H` additive/multiplicative,
`H≥1`, and `∑≤∏` for naturals `≥2` (their Lean kernel) gives join-monotonicity, so its content is
the strong core. Likewise `n−c₃≤tr` and `tr≤smax+1` (THM-1845) and `dom+tr≤n+1` (THM-1850) are
join-monotone and reduce. On the strong core `c₃≤H` is never tight (min margin `H−c₃=2`, max ratio
`c₃/H = 2/5` at n=4) — comfortably robust.

## 3. Velocity table and the fragility predictor

`D+`=`v▷T` (source), `D−`=`T▷v` (sink). Constant over all 530 classes:

```
            Δc₃  ΔH   Δtr  Δscc   Δsmin  Δsmax
  D+ (src)   0    0    +1   +1     0      var
  D− (sink)  0    0    +1   +1     var    +1
  complement fixes {c₃, H, tr, scc, srange}
```

> **Fragility predictor.** If some inflation `O` has `Δ_O(LHS) > Δ_O(RHS)`, `Φ` is refuted by the
> `O`-orbit of a tight base tournament to the first floor crossing.

`tr` (Δ=+1) and `srange` (Δ up to +3) **decouple upward** from the frozen `{c₃,H}` under `D±`.
So any `[tr/srange-side] ≤ [function of c₃,H]` bound is inflation-fragile, witnessed by the 3-cycle
atom + transitive singletons (THM-1830) — the tournament analog of WOWII-103's triangle+leaves.
This is opus THM-1865's diagnostic, here as the contrapositive of the reduction principle.

## 4. New results

**Minimal WOWII-repair.** kind-pasteur's `srange ≤ tr` is fragile (`D−`: srange outruns tr) and
breaks at n=7 (witness c₃=4, tr=5, srange=6). The **minimal repair is `srange ≤ tr + 1`**
(off by exactly one; the n=7 witness is an equality case). Verified exhaustive n≤7, join-monotone
n≤7. [Looser survivors: `srange≤tr+c₃`, `srange≤2(tr−1)` — both loose on the strong core.]

**WOWII-103 proper (king-eccentricity).** With `tr` and average king-eccentricity, achievable
avg king-eccentricity straddles `e` at **19/7 = 2.71429** (below, n=7) and **17/6 = 2.83333**
(above, n=6) — the exact analog of WOWII-103's `30/11 > e`.

**Methodological caveat (verified).** Inflation-fragility is a *low-entropy* phenomenon: the broken
`srange ≤ tr` survives 400 random samples at every n=8..12, because random tournaments are
near-regular and never hit the near-transitive fragile corner. Refutation requires the targeted
inflation witness (the velocity predictor's `O`-orbit) or exhaustive small-n — **not** a random
large-n sweep.
