---
source: kind-pasteur-2026-07-24-S145 (Opus 4.8)
status: INDEPENDENT AUDIT of the THM-2941/THM-735 foundation (the h/7 wall and the finite-tail seal),
  reproduced from scratch in exact rational arithmetic; PLUS an honest boundary on what I can and cannot
  census. Confirms the base of the current live obligation; does NOT claim census progress on k=2,3.
tags: [lrc, lonely-runner, audit, THM-2941, THM-735, h7-wall, six-body-carrier, census, honest]
related: [THM-2941, THM-2923, THM-735, THM-2928, kps-S144]
---

# Independent audit of the h/7 wall, and where my reach ends on the k=2,3 census

Task: dent the live obligation — the `k=2,3` six-body census (THM-2941). I re-derived the foundation from
scratch and verified it; I also hit, and am recording, the exact point past which I cannot reliably contribute.

## 1. Setup (as I now understand it, from THM-2941 §1)
Six-body carrier `E ∈ C({1,…,14}, 6)`; `D_w = {t : ‖wt‖ < 1/14}`; `C_E = T \ ⋃_{e∈E} D_e`; `h = μ(C_E) > 0`.
External labels `w ≥ 15`; `c(w) = μ(C_E ∩ D_w)`. A 13-speed cover is `E ∪ {w_1,…,w_7}`; LRC(14) here asks that
the seven externals never cover `C_E`.

## 2. A misread, corrected (recorded per the thread's method-discipline)
I first read "the wall `h/7`" as an **upper bound** `c(w) ≤ h/7` and my computation "refuted" it (`c(w)/(h/7)`
≈ 1.4–1.8 for small `w`). That reading was **wrong**. `h/7` is the **equidistribution limit**:
`μ(D_w) = 2·(1/14) = 1/7`, so as `w → ∞`, `c(w) → μ(C_E)·μ(D_w) = h/7`. Small/resonant `w` cover *more* than
`h/7` — which is precisely why the census must examine *small* externals (`z_1 ≤ 1742/312`), not large ones.
(Fifth "assume-hard-behaves-like-easy" slip of this thread, caught by computing.)

## 3. The audit (independently confirmed, exact arithmetic)
For carriers `E = {1,…,6}` (`h = 16/35`, `r = 12` components) and `E = {9,…,14}` (`h ≈ 0.35491`, `r = 38`):
- **Wall = equidistribution limit.** `c(w) → h/7` monotonically-in-scale: e.g. `E={1..6}` gives
  `c = 0.06619, 0.065143, 0.065405, 0.065343` at `w = 100, 500, 2000, 10000` → `h/7 = 0.065306`. ✓
- **Finite-tail seal (THM-735).** `c(w) < h/7 + γ/w` with `γ = 99r/490` holds in **every** sampled case
  (`w = 23,100,500,2000,10000`, both carriers). ✓
So the base THM-2941/THM-735 rests on — `c(w) → h/7` with an explicit `γ/w` finite-`w` excess — is
**independently reproduced.** This validates the wall and the tail that make the seven-slot balance critical
(`7 · h/7 = h`): coverage is only possible if resonant small externals supply the excess that aligned large
ones cannot.

## 4. Honest census boundary (what I cannot reliably do)
Censusing `k=2,3` needs the drift decomposition: which externals sit "at the wall" (aligned) vs "drift", the
pair functional `B_2 = max_{u<v} μ(C_E ∩ (D_u ∪ D_v))`, and the Gram/spanning-tree Hunter scalar `(h/49)(7I−J)`
that yields the `z_1 ≤ 1742/312` caps. The aligned externals are **not** bounded by the cap (only the first
drift is), so the bank is not a simple bounded box, and I do not have an independent, error-free handle on its
enumeration. Given this thread's error rate in unfamiliar apparatus, I will **not** assert a census result I
cannot verify.

## 5. What I CAN offer the census owners (reliable, direct)
- **A direct config-checker.** For any specific 13-set `E ∪ {w_i}` the fleet flags as a candidate, I can compute
  its exact gap (`max_{a/q} min_v ‖v·a/q‖`) and certify `≥ 1/14` (loose) or `< 1/14` (cover) — no drift
  machinery needed. This is an independent second-opinion on any zero-margin terminal or suspicious bank row.
- **The audited foundation (§3)** as a reusable exact-arithmetic module (`/tmp/wall2.py`): `C_E` construction,
  `c(w)`, and the `γ/w` tail — handy for spot-checking the seal on carriers beyond the two here.
- The **max-deletion gate** (kps-S144) still applies: any candidate whose top-speed-deleted 12-set has
  `δ·max > 2h` is loose without further work.

## 6. Honest status of task "dent the live obligation"
Delivered: an independent audit of the foundation (§3) and a corrected understanding (§2). NOT delivered:
census progress on the `k=2,3` banks — that requires the drift/Hunter apparatus I flagged in §4. The most
useful next step for me specifically is to serve as an independent gap-checker for the banks' critical rows
rather than to re-derive the drift reduction.

Files: `/tmp/{carrier,wall2}.py`.
