---
id: HYP-3099
title: Applying the tournament-as-proof catalog to three targets — cap-optimality is a NON-transitive improvement tournament (bounded finite check); baby-Hodge holes are c5-spectral NOT forbidden-H (α₂) instances; apex-7 ↔ forbidden-H=7 is a numerical coincidence, the genuine lever being the order-2 antipodal obstruction
status: VERIFIED (three applications, each with an independent verification script; two correct seductive over-claims). Diagnostic, not new bounds.
source: mac-mini-2026-06-27-S65
reflections:
  - tournaments-as-proof-engines-a-generative-catalog
applies:
  - THM-029   # forbidden H=7 (completeness-forcing)
  - THM-079   # forbidden H=21 (i2-jump)
  - THM-200   # 3 triangles force directed C5 (the pentagon)
  - THM-499   # two-layer baby-Hodge dichotomy H = 1 + 2(c3+c5) + 4 alpha2
  - THM-576   # cap = pairwise avoidance (minimizers)
  - THM-577   # symbolic overlap (apex-14 threshold)
related:
  - HYP-2907  # H=7 realizability guardrail
  - HYP-2908  # forbidden-H=7 state-lift (the OPEN first arrow)
  - HYP-2605  # winding tournament
  - OPEN-Q-099 # baby-Hodge holes
  - OPEN-Q-108
---

# HYP-3099 — Three applications of the tournament-as-proof catalog

Generated the catalog (`tournaments-as-proof-engines-a-generative-catalog`, 12 techniques, 4 novel) and
applied three of its techniques to three targets. Each application carries an independent verification
script; two of the three **correct a seductive over-claim with evidence** — the honest, valuable mode.

## App A — Technique #6 (exchange / improvement tournament) → cap-minimizer optimality
`04-computation/lrc_cap_optimality_exchange_macmini_S65.py`. The cap optimality (THM-576/577 value symbolic,
optimality open) attacked via the improvement tournament (orient config `P→P'` if a single speed-swap lowers
`meas(lonely)`).
- **VERDICT (diagnostic):** the minimizer is **bounded** (identical over `{1..13},{1..16},{1..20}`) — so
  optimality is a **finite check**; but the single-swap improvement digraph is **NON-transitive** — 4
  spurious local minima at j=3 (`{1,12,13}`=global, plus `{2,9,11},{3,10,11},{5,8,9}`), and greedy descent
  gets **stuck** at j=5 (k=8) on `{1,10,11,12,13}` vs the true `{1,5,7,8,9}`.
- **What the technique bought:** "why is cap-optimality hard" is answered precisely — *the config improvement
  tournament has 3-cycles (non-transitive), so no clean exchange/greedy proof exists; the k=8 'break' IS that
  non-transitivity.* Optimality reduces to a finite, bounded check (rule out the finitely many local minima),
  not an unbounded search.

## App B — Technique #2 (realizability-hole = forbidden-H) → baby-Hodge cycle-count holes
Background scout + `04-computation/baby_hodge_forbidden_h_crosscheck_opus.py` (independent 2^15 sweep).
- **VERDICT (corrects an over-claim):** baby-Hodge holes are **NOT** forbidden-H instances. The flagship hole
  `(c3,c5)=(8,10)` at n=6 is PROVED unrealizable and moment-feasible (`50 = ⅓·40 + ⅔·55`, skew-Hankel PSD),
  but its neighbors carry `H=41,43` — **both realizable, neither forbidden**. By THM-499's dichotomy
  `H = 1 + 2(c3+c5) + 4α₂`, the `(8,10)` hole is a **c5/spectral exclusion** (Landau score-stratification),
  *orthogonal* to the α₂/conflict-graph H-gap that THM-029 completeness-forcing rules out.
- **What the technique bought:** it cleanly **separates the two hole-layers** — forbidden-H (THM-029/079)
  closes the α₂ sub-family `{7,21}`; the c5-family needs per-`(n, score-class)` certificates and the holes
  *migrate* (`c5=10` forbidden at n=6 is realized at n=7), so no single uniform forcing lemma. The deep
  single-object bridge (THM-200's directed-C5 = the E₇ metagraph pentagon) is the real open form, at
  `H∈{7,21}` not at the c5-hole. (Caught a subtlety: `(c3,c5)` does NOT determine `H` in general — 6 fibers
  split at n=6.)

## App C — Technique #1/#5 (forbidden-H transfer / winding tournament) → apex-7 ↔ H=7
Background scout + `04-computation/apex7_vs_forbiddenH7_bridge_audit.py`.
- **VERDICT (corrects a slogan):** `apex-7 ↔ forbidden-H=7` is a **numerical coincidence**, not a bridge.
  (i) The winding tournament `T(t)` (HYP-2605) avoids `H=7` **vacuously** — identical spectrum `[1,9,11,15]`
  for tight / GW / non-tight configs (H=7 is forbidden for *every* tournament, THM-200), so it carries **zero**
  discriminating info about tightness. (ii) The real apex-7 is `n/2 = 7` antipodal diameters at `t=1/14` —
  a **triangle-free perfect matching** (verified count `m/2` for all even m), structurally the **opposite** of
  K₃. "Apex tied count 7" and "I(K₃,2)=7" are two different sevens colliding at denom 14. (iii) THM-577's
  overlap-triangles (`p+q>14`) are *undirected speed* triangles and run the *good* (loneliness-helping)
  direction — no functor to `Ω=K₃`. The project's own guardrails (HYP-2907/2908, THM-572) exist to block this
  over-claim; the load-bearing "first arrow" `apex over-cover → Ω=K₃` is unproved and the route is retired.
- **The genuine lever** (the honest residue): the **order-2 antipodal obstruction** (S48) — a tournament has
  **odd |Aut|** (no order-2 automorphism, since a pair-swap reverses an arc), so the symmetric extremal at
  `t=1/14` (the 7 tied diameters) must resolve, each of the `2^7` resolutions giving `M≥1/14`. Correct and
  tournament-native, but reduces to the **consec-maximizes crux** that all LRC routes hit.

## Meta-finding
Across the three, the tournament techniques mostly **diagnose** rather than resolve, and twice **refute** a
seductive coincidence with evidence (the value of the catalog is as much in killing false bridges as in
building true ones). The reusable LRC tournament levers that survived: the **order-2 antipodal / odd-|Aut|**
obstruction (App C) and the **non-transitivity** characterization of cap-optimality (App A) — both point at
the same consec-maximizes crux from different sides. The forbidden-H technique is real but its LRC target is
the α₂-layer (`{7,21}` / the THM-200 pentagon), NOT the c5/spectral or apex/measure phenomena.

## Next
1. App A's bounded-finite-check: enumerate the finitely many single-swap local minima per j and prove the
   top-cluster/break beats each (a finite certificate for cap-optimality).
2. App C's order-2 obstruction × App A's non-transitivity: both reduce to consec-maximizes — is the odd-|Aut|
   resolution argument exactly the certificate that rules out App A's spurious local minima?
3. THM-200 pentagon single-object identification (the genuine forbidden-H ↔ E₇ bridge), OPEN-Q-099 c5-family.
