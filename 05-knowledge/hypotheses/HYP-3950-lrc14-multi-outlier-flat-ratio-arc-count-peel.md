# HYP-3950: LRC(14) r=2 residual, multi-outlier 11-cores — the flat-ratio ledger + uniform arc-count peel

**Status:** CLAIMED (stub) — kind-pasteur-2026-07-01-S28, in progress this session.
**Numbering note:** adopting klein-S68's per-machine HYP-block proposal (kps = 3950+) to stop the collision cascade.

## Context
The r=2 residual of the multi-far moment relaxation reduces LRC(14)'s open covering branch to:
**inf meas(L_C) >= 1/36 over all 11-speed cores C** (band 1/14). kps-S27 closed the *bounded* part
(15472-core census, min = pentagon 313/9702 = 0.032261, ALL >= 1/36, margin 1.161x). The single-outlier
part (10 bounded + 1 far) has min {1,2,3,5,7,8,9,11,12,13} + w=50 at 0.03231 = the (6/7) x (10-core min)
equidistribution limit. **The residual is the multi-outlier case (>= 2 far speeds).**

## Claim (to be tested/proved this session)
1. **FLAT-RATIO LAW:** min over bounded m-cores of meas(L_B) ~ (7/6)^(11-m) x 0.0323 (each bounded slot
   costs one equidistribution factor 6/7, same as a far slot). Ledger min_m for m=7..11 computable exactly.
2. **UNCONDITIONAL UNION FLOOR:** meas(L_{B u W}) >= meas(L_B) x (1 - r/7) for r = |W| outliers
   (elementary union bound; no spread condition). With the ledger this closes r <= 3 outliers if
   min_m x (1-(11-m)/7) >= 1/36 holds for m >= 8.
3. **UNIFORM ARC-COUNT PEEL** (the sharpening for r >= 4): L_S is a union of <= 1 + sum(v in S) arcs;
   peeling the largest outlier w gives meas(L ∩ safe(w)) >= (6/7) meas(L) - 2 N(L)/(7w), N(L) = arc count;
   N grows by <= w per peel. Spread outliers peel at ~(6/7)^r; clustered outliers handled by pair-overlap
   (integer combs cannot be disjoint — Weyl fills the torus; pair overlap >= 1/49 - D_{p,q}-type error).
4. Multi-outlier min stays >= 1/36 (target: verified census + assembled peel skeleton).

## Evidence / artifacts (to be filled)
- Script: 04-computation/lrc14_multi_outlier_flat_ratio_peel_kps.py
- Results: 05-knowledge/results/lrc14_multi_outlier_flat_ratio_peel_kps.out

## Depends on / relates to
THM-522, THM-523, HYP-3787 (O(1/w) far-comb), kps-S27 census, THM-563/564 (pair Dedekind machinery),
mac-mini HYP-3824 (sub-critical slope clarification), OPEN-Q-108.
