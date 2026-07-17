# THM-951 — The decidable census pipeline (death-star-2026-07-17-S43)

**Status:** PROVED (Lean — `TournamentH7/LRCDeepCertificate.lean`; pipeline theorems
kernel-pure standard trio; the demo uses kernel-executed `decide`; verify the build
report in the session log). Source: HYP-7204. The integration of the census funnel.

## Statement

1. `instDecidableCoverageCapped`: the coverage cap DECIDES (finite multiplier
   range, computable bandCount).
2. `lonely_of_census` (**the pipeline**): for nonzero speeds,
   `CoverageCapped v q 6` + `#{bandCount = 6} < liveCount` ⟹ `∃ t, Lonely 14 v t`
   — THM-945's exact identity `B5 = live − deepSix`, `mreach_ge_of_B5_pos`, and
   `lonely_of_Mreach_ge` composed. `lonely_of_census_capped5`: at cap 5 one live
   multiplier suffices (`B5 = liveCount`).
3. `census_demo`: the family `v_i = i + 2` is certified lonely END-TO-END at
   `q = 31` with cap and census both by `decide` — **the first loneliness proof in
   the corpus produced by the B5 funnel** rather than a direct witness or a pack
   row. Concrete families now acquire machine-checkable loneliness certificates
   with no witness search and no analysis.

## Companion

`00-navigation/LRC14-FORMALIZATION-PICTURE-2026-07-17.md` — the cohesive map:
trunk (cite + residual), the shrink table (THM-937→942B + codex/opus wires), the
B5 funnel column (THM-940/942/943/944/945/950 → this pipeline), the killer/witness
arc (THM-883 → 947 → 949 → 950 → codex's Plücker layer), the parallel analytic
arcs, and the honest three-item open list.
