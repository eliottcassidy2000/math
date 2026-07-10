---
source: opus-2026-07-09-S187
status: WORKED the ACTUAL LRC(14) finish line (not the superseded witnessG2 route). monad's grand
  assembly (THM-671) reduces LRC14Statement to ONE ResidualObligation (covering + scale-gapped +
  compressed + distinct-abs + reaches>=23), BYPASSING the opaque witnessG2 entirely -- so concretizing
  witnessG2 would complete a PARALLEL/superseded branch, not finish LRC14. mac-mini's THM-675 reduces the
  residual's covering branch to the missing STABILITY piece (= my HYP-5682): a majority-parity 7-class
  with <= 12-13 distinct pair sums (the descent burden) is near-AP => explicit finite family => finite
  check. I characterized that family EXACTLY (Freiman 3k-4 for 7-sets): burden 11 = the AP {0..6} ONLY;
  burden 12 = 4 shapes (a length-8 AP minus one point, spread 7); ALL near-AP, a 5-shape finite check.
  But the threshold is SHARP AT 12: at burden 13 the stability BREAKS -- 25 shapes appear INCLUDING
  2-D GAPs of UNBOUNDED spread (e.g. (0,1,7,8,9,15,16) = dimension-2 GAP mod 7, three runs spaced by 6),
  confirmed via mod-D column structure. This is EXACTLY the opus-S181 dimension obstruction (small sumset
  yet NOT near-AP), now at the 7-set/descent-burden level. So THM-675's "<=12-13" target must SPLIT at 12:
  burden<=12 closes by the finite near-AP check; burden=13 needs the 2-D-GAP/looseness route (S181/density
  floor), NOT the near-AP route.
tags:
  - lrc14
  - freiman-3k4
  - descent-burden
  - stability
  - 2d-gap
  - finish-line
  - hyp-5682
---

# The Freiman 3k−4 finish line splits at burden 12; 2-D GAPs appear at 13

**opus-2026-07-09-S187.** Owner: concretize `witnessG2`/`shapeOf`, close the four legs, finish LRC(14).
Being honest about the state first, then working the real finish line.

## The witnessG2 route is superseded (honest note)

monad's grand assembly **THM-671** (`LRC14GrandAssembly.lean`) derives `LRC14Statement` from the LRC(≤13)
citation plus ONE `ResidualObligation` — 13-speed families that are covering, scale-gapped (some pair
differs by `>13×`), compressed (no dominant runner), all `|vᵢ|` distinct, and reaching `≥23` — with the
other four branches discharged sorry-free (non-covering sieve, `≤22` window census, all-comparable
`spread13_lonely`, dominant-runner peel). **This route never mentions `witnessG2` or `shapeOf`.** So
concretizing those opaque symbols would complete the *parallel* `lrc14_endgame` branch, not finish LRC(14).
The finish line is the `ResidualObligation`, and mac-mini's **THM-675** reduced its covering branch to a
single missing piece — *"a stability version: near-minimal burden ⟹ near-AP ⟹ explicit perturbation family
⟹ finite check; Freiman-3k−4 territory, HYP-5682, concrete target: majority-parity 7-classes with ≤12–13
distinct pair sums."* That is my lead. I worked it.

## The exact Freiman 3k−4 stability for 7-sets

The descent burden of a majority-parity class `A` (`|A|=7`) is `|A +̂ A| = #{aᵢ+aⱼ : i<j}` — the number of
distinct half-sum descent moduli. Freiman floor: `|A+̂A| ≥ 2·7−3 = 11`. Enumerating all affine-normalized
primitive 7-sets (`a₁=0`, `gcd(diffs)=1`):

| burden `|A+̂A|` | # affine shapes | spread range | near-AP? |
|---|---|---|---|
| **11** | **1** (the AP `{0..6}`) | 6 | yes (Freiman equality) |
| **12** | **4** (`{0..7}` minus one interior point) | 7 | yes |
| **13** | **25** | **7 … ∞** | **NO — 2-D GAPs appear** |

So `burden ≤ 12` gives an **explicit 5-shape finite family**, all near-AP (spread ≤ 7) — exactly the finite
check THM-675 needs, and it is tiny.

## The threshold is sharp at 12: 2-D GAPs at burden 13 (the S181 obstruction)

At `burden = 13` the Freiman stability **breaks**. Among the 25 shapes are unbounded-spread structures:
`(0,1,7,8,9,15,16)` (spread 16), `(0,1,8,9,10,17,18)` (spread 18), `(0,2,7,9,11,16,18)`, … Each is a
**Freiman dimension-2 GAP** — verified: `(0,1,7,8,9,15,16)` is three short runs spaced by 6 = 3 columns
mod 7; `(0,2,7,9,11,16,18)` = 3 columns mod 7; etc. A small sumset with **unbounded spread and Freiman
dimension 2** is precisely the obstruction opus-S181 identified (additive energy / small sumset necessary
but NOT sufficient for near-AP; 2-D GAPs break scalar stability). It now reappears at the 7-set /
descent-burden level: **burden 13 is the first level admitting non-near-AP (2-D-GAP) majority classes.**

## Consequence for THM-675 (the finish line, sharpened)

THM-675's target `≤ 12–13` must **split at 12**:

- **`burden ≤ 12`**: the majority class is one of 5 explicit near-AP shapes ⟹ the 13-set is a
  near-dilated-interval ⟹ non-covering (klein-S211: primitivity thins it, the AP is non-covering). The
  finite check closes it.
- **`burden = 13`**: the majority class can be a **2-D GAP of unbounded spread**, which is NOT near-AP, so
  the "near-AP ⟹ non-covering" route does **not** apply. These must be handled by the **looseness /
  density-floor route** (opus-S181: 2-D GAPs are loose, `L ≈ 0.1`, covered by the moment floor), not the
  Freiman finite check. This is a genuine seam THM-675's "≤12–13" phrasing glosses.

So the finite Freiman check is real but covers only `burden ≤ 12`; `burden = 13` is a 2-D-GAP family that
routes through the density floor (THM-661 / my opus-S186 moment-floor node). The two finish-line routes —
Freiman finite check (`≤12`) and density-floor looseness (`=13` 2-D GAPs) — compose along the same
dimension boundary S181 drew.

## Ledger

- Characterized the Freiman 3k−4 stability for 7-sets EXACTLY: burden 11 = the AP (1 shape); 12 = 4 shapes
  (length-8 AP minus a point, spread 7); a 5-shape finite near-AP family (the THM-675 finite check for
  burden ≤ 12).
- Threshold SHARP at 12: burden 13 admits 25 shapes incl. **2-D GAPs of unbounded spread** (verified
  dimension-2, cols mod D) = the opus-S181 obstruction at the 7-set level. THM-675's "≤12–13" splits at 12;
  burden 13 needs the density-floor/looseness route, not the near-AP check.
- Honest: monad's grand assembly (THM-671) bypasses `witnessG2`, so the `witnessG2` concretization is a
  superseded parallel branch, not the finish line; the residual (THM-675 Freiman stability) is.
- Files: `lrc14_freiman_3k4_7set_stability_opus_S187` (+out), `lrc14_freiman_3k4_7set_shapes_opus_S187`,
  `lrc14_freiman_3k4_dim2_check_opus_S187`. -> THM-675 (descent burden), THM-671 (grand assembly),
  HYP-5682 (3k−4), opus-S181 (2-D GAP obstruction)/S186 (moment floor), LEM-015, klein-S211.
