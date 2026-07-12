---
source: opus-2026-07-11-S239
status: A REFRAMING (correcting S237/S238) + a UNIFICATION. The three-gap coverage advantage (mac-mini
  cont.44) shows the spread bulk is the FAVORABLE case: divisor-complete ⟹ spread ⟹ bad coverer ⟹ many
  clearing multipliers (min 8 over the window, verified). The genuinely hard object is the AP (good coverer,
  0 window-clearing), dispatched by t=1/14. Both LRC(14) routes reduce to the SAME crux — the AP inverse
  theorem ("the good coverer is the {kα} orbit").
tags:
  - lrc14
  - three-gap
  - coverage
  - spread
  - unification
  - inverse-theorem
  - reframing
---

# The spread bulk is favorable — both routes reduce to the AP inverse theorem

**opus-2026-07-11-S239.** Owner: work the three-gap disjunction. Reading mac-mini cont.44 (three-gap
regularity = why the AP is the best coverer) reframes my spread-bulk residual and unifies the two LRC(14)
proof routes.

## The reframing (corrects S237/S238)

Clearing at a modulus `q` (a lonely multiplier `p`, `bandCount = 0`) means the family **fails to cover** the
danger arc at `p`. mac-mini's three-gap coverage advantage: the AP (the Steinhaus `{kα}` orbit) is the
**unique best coverer** — it covers 7–24× more than iid because three-gap regularity is anti-clumping; every
non-AP (spread) family covers like iid, i.e. **badly**. A bad coverer therefore has **many** lonely/clearing
multipliers. And divisor-complete **forces** spread (S237: 99% longest-AP ≤ 7; structural — it needs multiples
of 8, 9, 11, 13, 14, incompatible with the tight AP). So:

> **divisor-complete ⟹ spread ⟹ bad coverer ⟹ many clearing multipliers ⟹ clears easily.**

**Verified** (window = non-14 moduli in [15,31]):

| family | window clearing-count | clears at |
|---|---|---|
| **AP `{1..13}`** (good coverer, bucket A) | **0** | no non-14 modulus (only mult. of 14, `t=1/14`) |
| **spread divisor-complete** (bad coverers) | **min 8, mean 31, max 56** | mean ~7 of 16 window moduli |

`0/800` spread families have window-clearing ≤ 3 — they clear **robustly**. So the S237/S238 framing ("the
spread bulk is the wall") was **backwards**: the spread bulk is the *favorable* case. The genuinely hard
object is the **AP** — the good coverer with **zero** window-clearing — and it is **not** divisor-complete
(bucket A), dispatched by the elementary `t=1/14` witness. (kps's "window-hard = spread" means spread families
clear at *varying* `q` — a wide window, no small shortcut, S238 — **not** that they are hard to clear.)

## The unification

Both LRC(14) proof routes now reduce to the **same** crux:

- **Density route** (mac-mini/klein): the moment-ladder base needs *consec/AP is the extremal* — the AP
  **maximizes** coverage (the good coverer). Hard part: proving the extremal **has** the `{kα}` structure.
- **Clean-ruler route** (kps/opus): the residual needs *divisor-complete (spread) clears* — spread = **bad**
  coverer ⟹ clears. Hard part: proving spread ⟹ **not** the AP good coverer.

These are the **same inverse theorem** — "the good coverer is the `{kα}` orbit" — approached from opposite
sides. The density route needs *the AP is the unique maximal coverer*; the clean-ruler route needs *everything
non-AP covers less* (the same statement, contrapositive). So the two routes the fleet has run in parallel for
months have **one shared wall**: the AP coverage-extremality / inverse theorem. This is mac-mini's own
verdict ("the open difficulty everywhere is the same: proving the extremal HAS the `{kα}` structure") now seen
to govern the clean-ruler route too.

## Honest scope

This is a **reframing + unification**, not a proof. The residual (spread ⟹ clears) is still the inverse
theorem — but now correctly oriented (the **favorable** direction: bad coverers clear) and **quantitatively
robust** (min 8 clearing multipliers over the window, not marginal). The value:

1. **Corrects** S237/S238 — the spread bulk is easy, not the wall.
2. **Identifies** the true hard object as the AP (good coverer), which is dispatched by `t=1/14` — so the
   clean-ruler route's difficulty is *not* in its own residual but in the shared AP inverse theorem, already
   handled on the clean-ruler side by the elementary bucket-A witness.
3. **Unifies** the density and clean-ruler routes onto one crux — so progress on mac-mini's three-gap
   coverage-extremality (proving the AP is the unique good coverer) closes *both* routes at once.

The concrete next target is therefore mac-mini's lane: **prove the AP is the unique maximal coverer** (the
`{kα}` inverse theorem). The clean-ruler route contributes the observation that its own residual is the
*favorable* direction of that theorem, robustly satisfied — so if the inverse theorem is proved for the
density route, the clean-ruler route's spread residual follows a fortiori.

→ mac-mini cont.44 (three-gap coverage advantage — the mechanism), THM-661 (AP best coverer), opus-S237/S238
(the spread bulk — reframed here as favorable), opus-S233/S235 (bucket A / band-edge — the AP dispatch),
LEM-010 (three-gap good-period), kps cont.36 (decoupling). Files:
`lrc14_threegap_reframing_opus_S239.py` (+`.out`).
