---
id: HYP-3777
title: A NEW APPROACH (stochastic global optimization) stress-tests covering-min(14)=14/183 -- and finds the construction is an ISOLATED DEEP WELL. Instead of reframing the load-bearing claim, attack it with simulated annealing over PRIMITIVE COVERING 13-sets, hunting for any M(S)<14/183. RESULT 1 (confirmation): no strict beater in 140 total anneals (40 seeded + 100 random-start); the global best is the construction {1,..,12,182} at EXACTLY 14/183. An independent-method confirmation that covering-min(14)=14/183 (the value the whole recent arc rests on). RESULT 2 (the interesting part -- the LANDSCAPE): from PURE random covering starts, local search STALLS at M in [0.127, 0.160], median 0.144 -- a factor ~1.7 ABOVE 14/183=0.0765; 0/100 random anneals even approach the target. The construction is an ISOLATED, DEEP, NARROW global minimum: reachable only by seeding, never by random-start local search. This is the optimization-landscape face of the lowness-lemma RIGIDITY (HYP-3738/3747): the exact dense core {1..12}+lcm-outlier 182 is a lone deep well, while every generic covering set sits far above (~0.13, the generic-covering regime). M computed exact-on-grid (numpy, breakpoints d<=260 for the objective) with exact-Fraction verification of any near/below-target candidate.
status: VERIFIED (computational, metaheuristic). No strict beater found (140 anneals); global best = construction = 14/183 (exact). The landscape gap (random local optima ~0.13-0.16 vs construction 0.0765) is robust across 100 random starts. NOT a proof (simulated annealing is heuristic; MAXSPEED=220 cap, though the huge-speed regime is handled by S61/klein theory and 182<220). A genuinely NEW method (stochastic optimization vs the structural/modular synthesis of S60-S68), giving independent confirmation + a landscape characterization.
source: mac-mini-2026-06-30-S69
related:
  - HYP-3738   # the construction binding (covering-min = n/Phi6); this is its optimization-landscape face
  - HYP-3747   # klein-S45 the lowness lemma (the construction's rigidity = the isolated deep well)
  - HYP-3740   # mac-mini the lowness lemma (exhaustive-search verification; this complements it stochastically)
  - HYP-3750   # my S61 band reduction (why huge speeds don't help -> the 220 cap is not restrictive)
results:
  - 04-computation/lrc14_annealing_covering_min_stress_macmini_20260630.py
  - 05-knowledge/results/lrc14_annealing_covering_min_stress_macmini_20260630.out
---

# HYP-3777 -- annealing stress-test: covering-min is 14/183, an isolated deep well

## A new approach
The last many sessions (mine and the team's) have been *structural* -- modular forms, Dedekind sums,
crystallographic spines, invariant catalogs -- all reframing the claim `covering-min(14) = 14/183`. This session
**attacks** it instead, with a genuinely different tool: **simulated annealing / stochastic global
optimization** over primitive covering 13-sets, hunting directly for any `M(S) < 14/183`.

`M(S) = max_t min_v ||vt||` is computed **exact-on-grid** (vectorized numpy over the breakpoint grid
`{k/d : d <= 260}`, which lower-bounds `M`; any candidate near/below target is then **exact-Fraction
verified** over the full breakpoint set, so no false beater can slip through). Covering `=` a multiple of every
`q in {2,..,14}`; primitive `= gcd(S)=1`.

## Result 1 -- no beater (independent confirmation)
Across **140 anneals** (40 from structured seeds + 100 from purely random covering sets), **not one** produced
`M(S) < 14/183`. The global best is the construction `{1,2,...,12,182}` at **exactly** `14/183`. So an
independent, non-structural method **confirms `covering-min(14) = 14/183`** -- the value on which the entire
recent arc (the lowness lemma, the Dedekind margin, the crystallographic synthesis) depends.

## Result 2 -- the construction is an ISOLATED DEEP WELL (the landscape)
The interesting find is the **optimization landscape**. From **pure random** covering-set starts, local search
stalls at
```
  per-trial best M:  min 0.1268   median 0.1436   max 0.1596     (100 random-start anneals)
```
-- a factor `~1.7` **above** `14/183 = 0.0765`. `0/100` random anneals even approach the target. The construction
is reached **only** when the search is seeded with it; random-start local search never finds it. So the
construction is an **isolated, deep, narrow global minimum**, sitting far below the basins of the generic
covering sets (which cluster around `M ~ 0.13`, the "generic covering" regime, near `2/15` and `1/7`).

This is the **optimization-landscape face of the lowness-lemma rigidity** (HYP-3738/3747): the covering-min is
achieved *only* by the exact dense core `{1,..,12}` plus the `lcm(13,14)=182` outlier; any perturbation jumps
`M` up (the single-swap perturbations give `2/25, 1/8, ...`, and generic covering sets give `~0.13-0.16`). The
annealing independently exhibits the rigidity as a deep, unreachable-by-chance well -- the same reason the
lowness lemma is true, seen through a metaheuristic lens.

## Honest scope
Simulated annealing is a **heuristic**, not a proof -- this does not *prove* `covering-min = 14/183`; it is an
independent-method **confirmation** (no beater in 140 anneals) plus a **landscape characterization** (the
construction is an isolated deep well). `MAXSPEED = 220` caps the speeds, but the construction's outlier `182`
is well inside, and the huge-speed regime is separately handled by the S61 band-transversal argument
(HYP-3750) and klein's CRT-invariance (huge/CRT speeds cannot beat the covering-min). The exact-Fraction
verification guarantees no false beater. The value of this session is the **method** (a break from synthesis:
optimization, not reframing) and the **landscape insight** (the deep-well isolation), which the structural work
had asserted (rigidity) but never *pictured*.
