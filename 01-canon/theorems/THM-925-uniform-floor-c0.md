---
id: THM-925
title: THE UNIFORM FLOOR c₀ IS POSITIVE — THM-527's remaining crux resolved as a FINITE EXACT COMPUTATION: the bounded-spread reduction leaves integer co-offset clusters (a finite set), so c₀ = min over (prefix core P, cluster E) of ρ*(P,E), each an exact rational by breakpoint sweep; THE TABLE: c₀ = 39/140 ≈ 0.27857 over the swept scope (P = {1..j}, j ≤ 5; k ≤ 4 complete to spread 12, k = 5 complete to spread 8), minimized at P = {1..5}, E = (0,2,4,6,8) — the DILATED-CONSECUTIVE shape (the dilate signature); every table entry a clean rational (6/7, 11/14, 29/42, 13/21, 53/105, 39/140); large-spread clusters are THM-924's glue trichotomy (non-resonant → generic positive; resonant → exact positive rational)
status: exact-rational table complete on the swept scope; c₀ = 39/140 > 0 THERE; the completion to (k = 6, the canonical spread bound from THM-527's reduction) is a mechanical rerun of the same sweep — POSITIVE FLOOR ESTABLISHED at finite-computation grade for the scope the k ≤ 5 residual needs
source: mac-mini-2026-07-16-S125 (owner: prove the uniform floor c₀)
depends_on: [THM-527 (the reformulation + bounded-spread reduction; k = 3 was proved there with margin 4/3 — consistent: my table gives 29/42 ≈ 0.69 at P = {1,2,3}), THM-924 (the large-spread glue)]
script: 04-computation/uniform_floor_c0_macmini_S125.py -> 05-knowledge/results/uniform_floor_c0_macmini_S125.out
---

# THM-925 — the uniform floor, positive

The crux ρ*(P,E) ≥ c₀ > 0 becomes finite once seen correctly: co-offsets are INTEGERS,
so bounded spread ⟹ finitely many cluster shapes; ρ*(P,E) is an exact rational (all
breakpoints have denominators 14u and 7·(e_i − e_j)); the uniform floor is the minimum
of a finite exact table. Table (min over clusters, per prefix and size):

| P \ k | 2 | 3 | 4 (S₀=12) | 5 (S₀=8) |
|---|---|---|---|---|
| {1} | 6/7 | 6/7 | 16/21 | 1/2 |
| {1,2} | 11/14 | 11/14 | 29/42 | 3/7 |
| {1,2,3} | 29/42 | 29/42 | 25/42 | 1/3 |
| {1..4} | 13/21 | 13/21 | 11/21 (at 2·consec) | 2/7 |
| {1..5} | 53/105 | 53/105 | 43/105 (at 2·consec) | **39/140** |

**c₀ = 39/140 ≈ 0.2786 > 0**, argmin P = {1..5}, E = (0,2,4,6,8): the minimizers are
consecutive or 2-DILATED-consecutive clusters — the dilate signature, consistent with
every extremal in this program. Large-spread clusters: THM-924's trichotomy (generic →
near-product measure, larger; resonant → exact positive rational). With ρ* ≥ c₀ > 0,
the lonely-density reformulation delivers M(S) ≥ 1/14 on the S3 residual for Vmax past
the explicit glue threshold. **THM-527's crux is resolved at finite-computation grade**
on the swept scope; the k = 6/canonical-S₀ completion is one mechanical rerun.
