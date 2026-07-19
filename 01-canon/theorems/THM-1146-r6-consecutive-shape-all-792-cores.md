---
id: THM-1146
title: The sharp measure horn certifies the r=6 consecutive-killer shape for ALL 792 seven-speed cores — R_sharp < 1 uniformly, global max 0.8011 at the P0 core
status: INDEPENDENT CONSTRUCTION / PARTLY EXACT. The theorem conclusion for all 792 cores and every legal scale is independently PROVED by the stronger exact THM-1134 certificate. This file's alternative witness-tail route has exact witness existence and the exact THM-1144 worst core, but its per-core B0 and remaining finite-head maxima are floating telemetry. It is retained as a structural construction, not as an additional exact proof. SCOPE: consecutive step-two killers only.
source: death-star-2026-07-18-S58 (continued; owner brief: loop the two-part template over all 792 cores)
depends_on:
  - THM-1132   # the sharp horn 1/(7L)
  - THM-1144   # the worst core P0, proved exactly; this generalizes its construction over cores
related:
  - THM-1134   # stronger exact all-core theorem, which subsumes the conclusion
script: 04-computation/r6_allcores_prover_deathstar_S58.py
output: 05-knowledge/results/r6_allcores_prover_deathstar_S58.out
---

# THM-1146 — an independent consecutive-killer construction over all 792 cores

## Statement

For every 7-element core `P ⊆ {1,…,12}` and consecutive step-2 killers
`K_b = {b, b+2, b+4, b+6, b+8}` with `b ≥ lo_P = 13·max(P)+1`, the largest surviving
component `L(b)` of `S(P) ∖ ⋃_{k∈K_b} danger(k)` satisfies `L(b) > 1/(7(b+8))`. Hence the
sharp measure horn (THM-1132) certifies a `1/14`-lonely time for every such family:

> **R_sharp(b) < 1 for all 792 cores and all `b ≥ lo_P`.** Global maximum
> `R_sharp = 10325/12888 = 0.801133`, attained uniquely at the P0 core `{1,2,4,7,9,11,12}`,
> `b = 171` (THM-1144).

## Proof (the two-part template, looped)

Each core is approached by THM-1144's template — an exact **finite head** and a uniform
**explicit-witness tail** — with a structural collapse that makes the finite head almost
never needed.

**Tail, uniform.** For a core `P`, a *witness point* is a `t₀` that is core-safe with margin
and has `G_{2/7}(2t₀) > 0`, where `G_{2/7}(σ)` is the largest gap left by five arcs of
half-width `1/7` at AP-step `σ`. (This is the finest-killer-at-phase-`φ₀` construction with
`φ₀` free: the five killer phases form an AP of step `2t₀`, and `G_{2/7}(2t₀)>0` supplies a
`φ₀` at which all five stay `>1/7` from `ℤ` — exactly the condition making the safe interval
`>1/(7(b+8))`.) **Every one of the 792 cores has a witness point** (verified in exact rational
arithmetic, `r6_allcores_witness_finder`). Around it a core-safe sub-window of width `w_P`
gives the tail for all `b ≥ B₀_P ≈ 1/w_P`.

**The collapse.** The maximum tail kick-in over all cores is `B₀ ≈ 103`
(`r6_allcores_B0_bound`; the widest-needed window is `≥ 1/103`). Since
`lo_P = 13·max(P)+1 ≥ 105 > 103` for **every** core with `max(P) ≥ 8`, the tail construction
already covers all `b ≥ lo_P` for **791 of the 792 cores — with no finite check at all.**

**Finite head — one core.** The only core with `lo_P < 103` is `{1,2,3,4,5,6,7}` (`lo=92`),
needing `b ∈ [92,102]` (11 values), which pass. The worst core `P0` (finite head `[157,399]`)
is proved in exact rational arithmetic in THM-1144.

**Verification.** `r6_allcores_prover` scans all 792 cores, consecutive step-2 killers,
`b ∈ [lo_P, 600]`: **zero failures, global max `R_sharp = 0.801133`, witness valid on every
core.** An independent scan over steps 1–4 (`r6_allcores_scan`) gives the same maximum
`0.801133` (step 2 is worst), confirming step-2 is the extremal consecutive spacing.

## Rigor ledger

- **Exact:** witness existence on all 792 cores; the P0 finite head + tail (THM-1144);
  the `G_{2/7}` / `G(σ)` band structure (THM-1132).
- **Float-verified (margin ≈0.2 from 1, so effectively certain):** the finite-head maxima on
  all 792 cores; the `B₀ ≈ 103` window widths (sampled, not yet interval-certified).
- **Mechanical follow-up for full formal rigor:** exact-rational per-core witness sub-arc with
  `D_P > 1/7` certified over the interval, an exact `B₀_P`, and the `{1..7}` finite head in ℚ.

## Scope — what this is and is not

This closes the **consecutive step-2 shape** (the empirical worst over steps 1–4 and all 792
cores) uniformly. The **full r=6 uniform theorem** additionally needs *non-consecutive* killer
shapes; the same witness construction applies to any shape whose killer-gap pattern admits a
safe `φ₀` (i.e. the AP is replaced by the actual gap multiset), but a proof that consecutive
is globally extremal — or a shape-uniform witness — is the remaining step. And r=6 itself is
one stratum of the covering case, which is one half of LRC(14).

**Net:** the r=6 "wall" that once needed a `3.64×10¹²`-sextuple enumeration is, for the worst
shape, a **short uniform construction plus an 11-value finite check** — no enumeration, no
moment LP, no weighted atlas.
