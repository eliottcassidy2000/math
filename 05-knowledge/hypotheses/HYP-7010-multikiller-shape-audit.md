# HYP-7010 — Multi-killer covering-min: the SHAPE GAP audit (non-interval small parts)

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S17)
**Verify-first:** anyone building on this should re-run the audit script before citing; numbers
below are stub-stage until the results file lands.

## What is claimed

THM-726 (multi-killer covering-min rigidity, `≥2` outliers `≥13` ⟹ `M ≥ 1/13`) rests on a
finite check (`lrc14_multikiller_landscape_macmini_S70.py`) that enumerated **interval cores
`{1..k}`, `k = 9,10,11` only**, with every outlier a carrier multiple (`≤ 220`). Two shape
classes are NOT in that basis:

1. **Non-interval small parts** `P ⊊ {1..12}`, `|P| ≤ 11`, e.g. `P = {1..12}\{6}` — these admit
   **free outliers** (outliers carrying no missing modulus), impossible for interval cores:
   e.g. `S = {1..12}\{6} ∪ {w, 182m}` is primitive covering multi-killer for EVERY `w ≥ 13`.
2. **Low-`|P|` shapes** where outliers carry compound moduli (`|P| ≤ 8`, `j = 13−|P| ≥ 5`).

The S111 low-M rigidity assembly ALSO leans on THM-726 with a definitional seam (outlier `≥13`
in THM-726 vs `>14` in the assembly wording) — the named "outlier-threshold bookkeeping".

## This session's work

(a) Exhaustive exact audit of ALL small-part shapes `P ⊆ {1..12}` with `j = 2,3,4` outliers
(adaptive caps), early-exit `1/13`-witness search + exact Farey fallback; boundary cases
`{1..14}\{a}` and the free-outlier family computed exactly. (b) Either the shape gap closes
(THM-726 Step 2 basis upgraded to all shapes) or a counterexample is filed (court case).
(c) The outlier-threshold bookkeeping written out precisely: covering + `≤1` element `≥13`
⟹ `S = {1..12, 182m}` (deep-well ladder, THM-724 Case 1) — so the covering-min assembly
needs multi-killer ONLY in the `≥2`-outlier-`≥13` sense; the `>14` wording is superfluous.

## Known at claim time

- `{1..13}\{6} ∪ {182}`: `M = 2/23 ≈ 0.087 ≥ 1/13` ✓ (HYP-6660, known).
- `{1..11,13,84}`: `M = 7/89` (THM-726 min, interval core k=11).
- `{1..14}\{6}` (near-AP, outliers `{13,14}`): M unknown at claim time — computed this session.

-> THM-726, THM-724, mac-mini-S111 assembly, MISTAKE-126 (fixed-window monotonicity refutations).
