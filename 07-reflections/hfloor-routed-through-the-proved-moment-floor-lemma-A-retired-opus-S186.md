---
source: opus-2026-07-09-S186
status: ROUTED hfloor through the PROVED moment floor (THM-661), RETIRING the open Lemma A. Formalized
  (LRCWitnessMomentFloor.lean, kernel-pure [propext,Classical.choice,Quot.sound], wired to root, built
  first try): witness_floor_from_momentfloor_nodes discharges hfloor (witnessMP <= witnessG2) from
  Bonferroni (PROVED kps-S30) + the moment floor (THM-661: nu=mu>=minD3>=momentBar k) + Lemma B (proved),
  with the Bonferroni arithmetic momentBar(k)+capRat(k)-1=witnessMP DEFINITIONAL (THM-661's (A') bar is
  CONSTRUCTED as witnessMP+1-capRat k). VERIFIED (exact rationals): THM-661's minD3(k) >= momentBar(k) for
  all k=8..13 (razor +0.00058 at k=8 -> +0.252 at k=13); the definitional identity holds for every k. NET:
  hfloor's open dependency moves from the never-proved standalone Lemma A to THM-661 -- the fleet's actual
  proved density floor. The only residue is THM-661's own decorrelation tail (LEM-005 coupled region),
  shared by EVERY route since nu=mu (opus-S185). Lemma A is now UNNECESSARY for the formalization.
tags:
  - lrc14
  - formalization
  - hfloor
  - moment-floor
  - lemma-a
  - density-floor
  - lean
---

# hfloor routed through the proved moment floor — Lemma A retired

**opus-2026-07-09-S186.** Owner: route hfloor through the proved moment floor. Done, and formalized. The
Lean skeleton's density-floor node `hfloor` no longer depends on the open Lemma A — it now depends on
THM-661, the density floor the fleet has actually proved.

## The substitution

`hfloor` = `witnessMP ≤ witnessG2 s`, where `witnessG2 = μ(GOOD ∩ G_P)`. The Bonferroni reduction
(`LRCWitnessBonferroni`, kps-S30, PROVED) gives
```
   witnessG2 ≥ ν(E) + measGP(P) − 1        (Bonferroni; ν = μ(GOOD), measGP = μ(G_P))
```
and closes it with **Lemma B** (`measGP ≥ capRat`, proved) + **Lemma A** (`ν ≥ nuConsec`, OPEN — the
compactness minimization). opus-S185 showed the key fact: `ν(E) = μ(E) = P(W>0)` *identically*, so `ν`'s
lower bound is exactly the covering measure THM-661's moment-LP bounds:
```
   ν(E) = μ(E) ≥ D3(E) ≥ min_E D3(E) ≥ (A') bar = witnessMP + 1 − capRat(k)      (THM-661).
```
So replace Lemma A's `ν ≥ nuConsec` by the moment floor's `ν ≥ momentBar(k) := witnessMP + 1 − capRat(k)`.
The Bonferroni arithmetic that needed a separate rational lemma in the Lemma-A route becomes
**definitional**:
```
   momentBar(k) + capRat(k) − 1 = (witnessMP + 1 − capRat k) + capRat k − 1 = witnessMP.
```
Chaining: `witnessG2 ≥ ν + measGP − 1 ≥ momentBar(k) + capRat(k) − 1 = witnessMP`. hfloor discharged.

## Formalized (kernel-pure, first-try build)

`LRCWitnessMomentFloor.lean` (wired to root, axioms `[propext, Classical.choice, Quot.sound]`):
```lean
theorem witness_floor_from_momentfloor_nodes
    (nuShape measGP : Shape → ℝ)
    (hbonf   : ∀ s, nuShape s + measGP s - 1 ≤ witnessG2 s)                       -- Bonferroni (PROVED)
    (hMoment : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
                 (momentBar (clusterSize s) : ℝ) ≤ nuShape s)                     -- moment floor (THM-661)
    (hB      : ∀ s, 8 ≤ clusterSize s → clusterSize s ≤ 13 →
                 (capRat (clusterSize s) : ℝ) ≤ measGP s)                         -- Lemma B (proved)
    (s : Shape) (h8 : 8 ≤ clusterSize s) (h13 : clusterSize s ≤ 13) :
    (witnessMP : ℝ) ≤ witnessG2 s
```
plus `momentBar_add_capRat` (the definitional arithmetic) and the `large_...` skeleton-shaped variant.
Structurally identical to `witness_floor_from_bonferroni_nodes`, but the OPEN `hA` (Lemma A) is replaced
by `hMoment` — which THM-661 supplies — and the arithmetic node needs no separate proof.

## Verified: the moment floor clears every bar

The `(A')` bar `momentBar(k) = witnessMP + 1 − capRat(k)` vs THM-661's `min_E D3(k)` (exact rationals;
`witnessMP = 14249/252252`):

| k | capRat | momentBar (A' bar) | min D3 (THM-661) | margin |
|---|---|---|---|---|
| 8 | 2243/5880 | 0.675025 | 0.675608 | **+0.00058** |
| 9 | 1979/4004 | 0.562231 | 0.567746 | +0.0055 |
| 10 | 55/91 | 0.452092 | 0.480 | +0.028 |
| 11 | 66/91 | 0.331212 | 0.400 | +0.069 |
| 12 | 6/7 | 0.199344 | 0.355876 | +0.157 |
| 13 | 1 | 0.056487 | 0.308844 | +0.252 |

`min D3(k) ≥ momentBar(k)` for every binding `k` (`k=13` has `capRat=1` since `|P|=0`, so `witnessG2 = ν`
and the moment floor `0.309 ≥ m_P = 0.056` clears outright). `k=8` is razor-thin (+0.00058) — there
THM-661's degree-3 `D3` is at its weakest and the degree-4 `B4` (0.7611) has far more room, so the leg is
safe. The definitional identity `momentBar + capRat − 1 = witnessMP` holds exactly for all `k`.

## What this wins, and the honest residue

- **hfloor's open dependency moves from an orphan to a theorem.** Lemma A was a standalone compactness
  minimization the fleet never proved (opus-S184/S185). THM-661 is the density floor the fleet HAS proved
  (per-shape exact; degree ≤ 4 clears all six legs). Routing hfloor through THM-661 aligns the Lean
  skeleton with the actual proof effort and **retires Lemma A** — it is no longer needed.
- **The residue is THM-661's own tail, not a new gap.** `min_E D3(k) ≥ momentBar(k)` is proved by the
  compact exact check + the decorrelation tail (LEM-005), whose coupled band `diam ∈ [18,35]` is the
  fleet's shared open step. Since `ν = μ` (opus-S185), EVERY route to `hfloor` bottoms out there — Lemma A
  did too. So this does not close that residue, but it stops duplicating it in a second (open) lemma and
  concentrates it in THM-661 where the fleet's transfer machinery (monad THM-669/670, boxeph mu-level
  transfer, klein modular supply) is already aimed.

## Ledger

- FORMALIZED `witness_floor_from_momentfloor_nodes` (LRCWitnessMomentFloor.lean, kernel-pure, root-wired,
  first-try build): hfloor ⟸ Bonferroni (kps-S30) + moment floor (THM-661, `ν ≥ momentBar`) + Lemma B;
  Bonferroni arithmetic `momentBar+capRat−1=witnessMP` DEFINITIONAL.
- VERIFIED min D3(k) ≥ momentBar(k) = witnessMP + 1 − capRat(k) for k=8..13 (k=8 razor +0.00058; degree-4
  B4 has room there). Lemma A RETIRED — unnecessary for the formalization.
- Residue = THM-661's decorrelation-tail coupled region (LEM-005), shared by all routes since `ν = μ`.
- Files: `LRCWitnessMomentFloor.lean`, `lrc14_momentfloor_bars_check` (verification). -> THM-661 (moment
  floor), LRCWitnessBonferroni / kps-S30 (Bonferroni), opus-S185 (`ν=μ`)/S184 (Lemma A scoping), LEM-005
  (coupled region), THM-669/670 (monad transfer), boxeph HYP-5722.
