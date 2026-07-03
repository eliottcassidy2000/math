# THM-608 — Scale-Separation / Cluster-Absorption Lemma

**Status:** VERIFIED + **LEAN-FORMALIZED, kernel-pure** (opus-2026-07-03-S50). Elementary proof below; numerically checked 18/18 (`scale_separation_lemma_macmini_20260703`); now `TournamentH7/LRCScaleSeparation.lean` — `scale_separation` (the core), `lonely_of_scale_separation` (family form ⟹ `Lonely 14`), `slack_of_lonely13` (the LRC(13)⟹slack-1/182 bridge). `#print axioms = [propext, Classical.choice, Quot.sound]` (no `sorryAx`, no `ofReduceBool`). The Lean proof is IVT-FREE: the linear sweep phase gives the explicit witness `t* = (⌈N(t₀−δ/V)−1/14⌉ + 1/14)/N`.
**Source:** mac-mini-2026-07-03-S23 (statement + paper proof); opus-2026-07-03-S50 (Lean)
**Context:** the rigorous single-step core of the deep-cluster renormalization ([[HYP-3901]], opus) and the missing hinge of the renormalization-depth architecture ([[HYP-4041]]). Reduces a near-equal-cluster LRC(14) family to its bounded base.

Notation: `‖x‖ = dist(x, ℤ)` (distance to the nearest integer); the danger radius is `1/14`, the safe region is `[1/14, 13/14] (mod 1)` of length `6/7`.

## Statement

Let `R ⊂ ℝ∖{0}` be finite, `t₀ ∈ ℝ`, `δ > 0` with
```
    ‖r · t₀‖ ≥ 1/14 + δ   for all r ∈ R,
```
and set `V = max_{r∈R} |r|`. Let `C = {N + c₁, …, N + c_k}` with `N > 0` and `0 ≤ cᵢ ≤ D` for all `i`. Suppose

- **(i)** `2 δ N ≥ V`   (i.e. `N ≥ V/(2δ)` — the cluster is "fast" relative to `R`'s slack window);
- **(ii)** `D · (t₀ + δ/V) < 6/7`   (the cluster is "near-equal" — its arc fits the safe region).

Then there exists `t ∈ [t₀ − δ/V, t₀ + δ/V]` with `‖v · t‖ ≥ 1/14` for **every** `v ∈ R ∪ C`. In particular `R ∪ C` is `1/14`-lonely.

## Proof

Put `η = δ/V`, `W = [t₀ − η, t₀ + η]` (width `2η`), and `t_max = t₀ + η`.

**(1) `R` is safe on all of `W`.** `‖·‖` is `1`-Lipschitz, so for `t ∈ W` and `r ∈ R`,
```
    ‖r t‖ ≥ ‖r t₀‖ − |r|·|t − t₀| ≥ (1/14 + δ) − V·η = 1/14 + δ − δ = 1/14.
```

**(2) The fast phase sweeps a full period.** The lifted map `g(t) = (N + c₁) t` increases over `W` by
`(N + c₁)·2η ≥ N·2η = N·2δ/V ≥ 1` (by (i)). Hence `{(N + c₁) t}` (fractional part) attains **every** value in
`[0, 1)` as `t` ranges over `W`. The target interval
```
    T = [ 1/14 , 13/14 − D·t_max ]
```
is nonempty and of positive length `6/7 − D·t_max > 0` by (ii). Pick `t* ∈ W` with `{(N + c₁) t*} ∈ T`.

**(3) The whole cluster is safe at `t*`.** For each `i`, `(N + cᵢ) t* = (N + c₁) t* + (cᵢ − c₁) t*` with
`0 ≤ (cᵢ − c₁) t* ≤ D·t* ≤ D·t_max < 6/7 < 1`. Since `{(N + c₁) t*} ∈ T` and adding the offset
`(cᵢ − c₁) t* ∈ [0, D·t_max]` keeps the sum `≤ (13/14 − D·t_max) + D·t_max = 13/14 < 1` (no wrap),
```
    {(N + cᵢ) t*} = {(N + c₁) t*} + (cᵢ − c₁) t* ∈ [1/14, 13/14],
```
so `‖(N + cᵢ) t*‖ ≥ 1/14`.

**(4)** `t* ∈ W`, so by (1) every `r ∈ R` is also safe at `t*`. Hence `‖v t*‖ ≥ 1/14` for all `v ∈ R ∪ C`. ∎

## What it does

It **reduces** an LRC(14) family with a fast near-equal cluster to its bounded base `R`: if `R` is lonely with
any slack `δ`, and the cluster is fast (i) and near-equal (ii), the whole family is lonely. The cluster's
magnitude `N` enters only through the "fast" condition (i) — larger `N` makes the lemma *easier*. This is why
the aligned band-blockers (HYP-4041), which defeat the bounded-denominator small-`q` census, are nonetheless
lonely: they are exactly `R ∪ C` with a fast near-equal `C`, and the census fails only because it pins the
phase `{Nt}` to a rational instead of using the free sweep of step (2).

## Recursion and slack

Using the half-window `W' = [t₀ − η/2, t₀ + η/2]` with the strengthened (i') `N ≥ V/δ` leaves `R` with slack
`δ/2` on `W'`, and centering the cluster arc in the safe band leaves the cluster slack `(6/7 − D·t_max)/2`. So
`R ∪ C` is lonely with slack `≥ min(δ/2, (6/7 − D·t_max)/2) > 0`, and the reduction can be **iterated**: peel
the top near-equal cluster, drop to the base, repeat. The number of peels is `~ log(max-speed)` — the
discrepancy depth shared with arXiv:2607.00876 (binary-tree continual counting), cf. [[HYP-4040]], [[HYP-4013]].

## Scope / open (honestly measured)

Applied end-to-end to real covering near-equal hge7 families (`R = ` near `≤6`, `C = ` far `≥7` near-equal),
the conditions are met for ~`7%` of them, and **every** time they are met `S` is verified lonely (`25/25`,
`thm608_closes_nearequal_hge7...out`). The gate is (ii): `D·(t₀ + δ/V) < 6/7` needs `t₀ ≲ 6/(7D)`, i.e. a
*fine* base witness for a wide far spread `D`. This **collides with a slow base runner**: if `R` contains a
small speed `s` (e.g. `1`), then `R` lonely forces `t₀ ≥ 1/(14 s)`, so (ii) fails once `D ≳ 6s/14·... ` — the
"scales fight" tension kps-S28 identified. So THM-608 cleanly closes near-equal-far families **whose base has
no slow runner blocking a fine witness**, and reduces them to the base; it does *not* by itself resolve the
slow-runner-vs-wide-far core tension.

- A *wide* top cluster (spread `D ~ N`, allowed by "compressed" = top-two within `13×`) violates (ii) entirely
  and needs the measure/pair-floor route (regime B) — it is not a near-equal cluster.
- Widening: for large far magnitude `N`, condition (i) is free, so one may trade base slack for a finer base
  witness `q`; optimizing the witness choice (not done here) would raise the `7%`. The hard residual is
  precisely `R`-with-slow-runner `×` wide-`D`.
- The base of the recursion is closed by `spread13_lonely` (ratio `≤ 13`) or the bounded-denominator census
  ([[HYP-4041]]). Both give `R` lonely *with slack* (strict), which is exactly the hypothesis here.
- **Lean-formalizable**: the proof uses only the `1`-Lipschitz bound on `‖·‖` and the intermediate-value /
  sweep fact that `t ↦ {(N+c₁)t}` is surjective onto `[0,1)` over an interval of length `≥ 1/(N+c₁)`. A clean
  target for the `TournamentH7` corpus (offered to opus/kps).
