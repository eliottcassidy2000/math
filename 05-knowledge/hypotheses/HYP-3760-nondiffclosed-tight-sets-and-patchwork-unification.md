---
id: HYP-3760
title: INDEPENDENT CONFIRMATION + EXTENSION of the tight-set classification (klein-S48 HYP-3750) and the patch-tuning unification (klein-S49 HYP-3753). Reached the SAME results by independent computation before syncing: (A) difference-closed tight sets = dilated APs (THM-560/opus HYP-3749); the non-diff-closed residual = the GODDYN-WONG doubling family {skip n-2, add 2(n-2)}, TIGHT exactly when n=2 mod 6, plus small sporadics (n=5,6) and multi-swap CROSS-types (n=8 {1,4,5,6,7,11,13}); (B) the PATCHWORK UNIFICATION -- GW (tight, M=1/n) and the covering-min construction (M=n/Phi6) are the two single-element skip-and-patch mutations of the mother AP {1..n-1} whose patch is a MULTIPLE of the skipped speed (DOUBLE 2(n-2) for the floor; LCM n(n-1) for the covering), each killing the resonance hole 1/k the skip opens. EXTENSIONS beyond klein's runs: the single-swap census verified n=5..16 (klein n=4..8) and the GW-doubling family verified n=5..26 (tight at 8,14,20,26; klein n=8,10,12,14) -- both confirm n=2 mod 6 exactly. CEDES priority to klein-S49/S48 (pushed first); kept as an independent cross-check + extended verification.
status: CONFIRMS klein-S48 (HYP-3750-tight) + klein-S49 (HYP-3753) by independent computation; EXTENDS the verification range (census n<=16, GW family n<=26). Priority ceded to klein (pushed first). Note the HYP-3750 NUMBER collides: mac-mini-S61 HYP-3750 (band-prime reduction, pushed first) vs klein-S48 HYP-3750 (tight classification) -- flagged to klein.
source: mac-mini-2026-06-30-S62
depends_on:
  - HYP-3753   # klein-S49: THE patch-tuning unification (mother AP, two healing modes) -- PRIMARY, this confirms it
  - HYP-2893   # Goddyn-Wong = Jacobsthal criterion
related:
  - HYP-3750   # klein-S48: the tight-set classification (GW-type + cross-type) -- PRIMARY [NOTE: number collides with my S61 band-prime reduction]
  - HYP-3749   # opus/klein: difference-closed = dilated APs; the punctured-core wide hole
  - HYP-+2913  # three-gap tight-locus characterization + census
  - HYP-3738   # the construction binding (covering-min = the LCM healing mode)
results:
  - 04-computation/tight_sets_and_patchwork_families_macmini_20260630.py
  - 05-knowledge/results/tight_sets_and_patchwork_families_macmini_20260630.out
  - 05-knowledge/results/tight_swap_census_macmini_20260630.out
  - 05-knowledge/results/tight_gw_family_patch_macmini_20260630.out
  - 05-knowledge/results/tight_patch_principle_macmini_20260630.out
---

# HYP-3760 -- independent confirmation + extension of the tight-set classification & patch-tuning unification

I was given "classify non-difference-closed tight sets and patchwork covering families" and computed the answer
independently; on syncing I found **klein-S49 (HYP-3753)** and **klein-S48 (HYP-3750-tight)** had reached the
same picture first. This note CEDES priority to them and records the independent confirmation + the extended
verification (which is genuinely useful for the program). Their files are the canonical statements.

## What I confirm (identical to klein-S49/S48)
- **Difference-closed tight sets = dilated APs** (THM-560, opus HYP-3749).
- **The patch-tuning unification** (klein-S49 HYP-3753): both extremal objects are single-element patches of
  the mother AP `{1..n-1}` -- skip an element, insert one -- and the patch heals the resonance hole `1/k` the
  skip opens:
  - **GW / floor:** `skip (n-2), patch 2(n-2)` (DOUBLE) -> `M = 1/n`, tight, non-covering; iff `n ≡ 2 (mod 6)`.
  - **construction / covering-min:** `skip (n-1), patch lcm(n-1,n) = n(n-1)` (LCM) -> `M = n/Phi_6`, covering.
  - At `n=14` the floor variety is the single point `{skip 12, patch 24}` = GW.
- **The non-diff-closed classification** (klein-S48 HYP-3750): GW-type single-swaps (`v -> 2v`) plus CROSS-type
  multi-swaps (`{1,4,5,6,7,11,13}` at `n=8`, drop 2 dup 5); "every one is GW" is refuted; per-`n` finiteness holds.

My independent computation reproduced all of this (`skip 12 + 24` uniquely tight at `1/14`; `skip 13 + 182`
uniquely covering-min at `14/183`; the cross-type `{1,4,5,6,7,11,13}`; the small sporadics `{1,3,4,7}`,
`{1,3,4,5,9}`).

## What I add (extension beyond klein's runs)
- **Single-swap census verified `n = 5..16`** (klein `n=4..8`): non-diff-closed single-swap tight sets occur
  only at `n = 5, 6, 8, 14`; `n = 7, 9, 10, 11, 12, 13, 15, 16` are AP-only (apex-clean).
- **GW-doubling family verified `n = 5..26`** (klein `n=8,10,12,14`): `skip (n-2) + 2(n-2)` is tight exactly at
  `n = 8, 14, 20, 26` -- confirming the `n ≡ 2 (mod 6)` criterion across a full period-and-a-half.
- **The "patch is a multiple of the skip" framing** makes both healing modes one statement: `2(n-2)` is a
  multiple of `n-2`, `lcm(n-1,n) = n(n-1)` is a multiple of `n-1`; each is the smallest multiple of the skipped
  speed that kills its resonance hole (DOUBLE = smallest new multiple; LCM = smallest that ALSO covers `q=n`).
- **Covering patchwork families:** the covering patchworks are `skip s + patch m·s` with `n | m·s`; the minimal
  (least `M`) is the construction, so the covering-min is the OPTIMAL patchwork covering of the AP (others,
  e.g. `skip 7 + 14`, cover but give `M = 1/11 > n/Phi_6`).

## Housekeeping
Priority to klein-S49 (HYP-3753) and klein-S48 (HYP-3750-tight). **Number collision:** klein-S48 used HYP-3750,
which already names my S61 band-prime reduction (pushed first); flagged to klein to renumber the tight-set one.
This HYP is retained as an independent cross-check + the extended-`n` verification, not a competing claim.
