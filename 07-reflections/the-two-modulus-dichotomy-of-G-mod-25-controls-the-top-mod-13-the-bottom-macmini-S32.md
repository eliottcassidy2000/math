# The two-modulus dichotomy of (G): mod-25 controls the top, mod-13 the bottom — and the residual is exactly the full transversals

*mac-mini-2026-07-06-S32 (HYP-4622). Owner: creatively whittle away the crux.
Sharpens (and corrects) kps-S41's mod-25 covering reduction (HYP-4567) by pinning the
exact clearability characterization, and connects it to my own S7 "gap member = full
transversal mod 25." Verified:
`lrc_two_modulus_transversal_macmini_S32.out`,
`lrc_transversal_adversarial_macmini_S32.out`.*

## The setup

(C), the crux, is the per-order achievability gauntlet (kps HYP-4557, opus-S122): no
12-integer-family has `M ∈ (1/13, 2/25)`. The two endpoints are `1/(N+1) = 1/13` and
`2/(2N+1) = 2/25` at `N=12`. kps-S41 reduced the *lower-bound* half to a "mod-25
covering fact": a family clears to `M ≥ 2/25` at `t = c/25` when some unit `c` rotates
every speed off the forbidden band `{0, ±1} mod 25`. kps conjectured **every** family
without a multiple of 25 is clearable this way, and flagged it as a finite covering
fact to prove.

## The exact characterization (corrects kps's conjecture)

A clearing rotation exists ⟺ some unit `c` has `c·v_i ∉ {±1} (mod 25)` for all `i`.
The **forbidden** `c` for runner `i` is `{±v_i^{-1}}` — the `±`-class of `v_i^{-1}`.
So a clearing rotation exists ⟺ the forbidden set `⋃_i {±v_i^{-1}}` does **not** cover
`(ℤ/25)*`. Inversion permutes the ten `±`-classes `{[1],[2],[3],[4],[6],[7],[8],[9],
[11],[12]}` bijectively, so:

> **Characterization.** A no-clearing-rotation family is exactly one whose unit-speeds
> **hit all ten `±`-classes mod 25** — a **full transversal mod 25**.

Therefore kps's "every no-mult-of-25 family clears mod-25" is **false precisely for the
full transversals** (e.g. the AP `{1..12}` itself: its units `{1,2,3,4,6,7,8,9,11,12}`
hit all ten classes, so no rotation clears it — correctly, since `M(AP)=1/13`). The
genuine residual of kps's reduction is the full-transversal families — which is exactly
my S7 result "a gap member must be a full transversal mod 25," seen from the dual side.

## The proved half (explicit closed-form witness)

If the unit-speeds **miss** a `±`-class `[a]` (family not a full transversal, no
multiple of 25), then `t = a^{-1}/25` is a witness: every `c·v_i` with `c = a^{-1}`
lands in `[2,23] mod 25`, so **`M ≥ 2/25`**. (Multiples of 5 but not 25 are always safe:
`5m·c ∈ {5,10,15,20}`.) Verified: family missing `[7]` clears at `t = 18/25`
(`18 = 7^{-1}`), all distances `≥ 2/25`. This is a clean rational-point witness —
exactly kps's `LRCMod25Floor` / `rational_point_margin` atom at `s=25, μ=2`, now with
the hypothesis pinned to a **decidable residue condition** ("not a full transversal").
So **the non-transversal half of (G) is proved**, and Lean-ready.

## The residual, resolved computationally

The residual is: **full transversal mod 25 with `M < 2/25` ⟹ dilated AP.** Evidence:

- **10,685 structured `≥3`-defect families**: 0 in-gap; among all 306 full transversals,
  the dilated AP is the **unique** one with `M < 2/25` (at the boundary `1/13`) — every
  other full transversal has `M ≥ 2/25`.
- **4,000 adversarial full transversals** (random reps per class, heights to ~120): 0
  in-gap; all `M ≥ 2/25`.

So being a full transversal mod 25 **and** near-tight (`M < 2/25`) forces the AP.

## The frame: (G) is a two-modulus incompatibility

The two endpoint denominators split the work by modulus:

- **mod `2N+1 = 25` controls the TOP.** Miss a `±`-class ⟹ `M ≥ 2/25` (explicit
  witness). To *avoid* the top a family must be a **full transversal mod 25**.
- **mod `N+1 = 13` controls the BOTTOM.** The AP is pinned at `M = 1/13` (kps's `r=1`
  argument / my S7 THM-593A mod-13 pinning). To be *above* the bottom a family must be
  **non-AP**.

A gap member must be **both** a full transversal mod 25 (else cleared to the top) **and**
non-AP (else sitting at the bottom). The residual says these are **incompatible**: full
transversal mod 25 + `M < 2/25` ⟹ AP. So the two moduli pinch the gap shut. Half of
this pinch (the top, via mod 25) is now a proved explicit witness; the other half (full
transversal ⟹ AP, a mod-13 rigidity) is the sharp remaining core, with 0
counterexamples in ~15k structured + adversarial families.

## Net

- **Characterization**: mod-25 clearing rotation exists ⟺ family is **not** a full
  transversal mod 25. Corrects kps-S41's "all no-mult-of-25 clear."
- **Proved half of (G)**: non-transversal (no mult 25) ⟹ `M ≥ 2/25` via the explicit
  witness `t = a^{-1}/25` (missed class `[a]`). Lean-ready (kps's atom + a decidable
  residue-covering fact).
- **Residual pinned**: the full-transversal families, where the AP is the unique
  `M < 2/25` member (boundary) — a mod-13 rigidity, 0 counterexamples in ~15k families.
- **Frame**: (G) = two-modulus incompatibility — `2N+1=25` closes the top, `N+1=13`
  pins the bottom; a gap member must thread both and cannot.

## Pointers

- Scripts: `lrc_two_modulus_transversal_macmini_S32.py`,
  `lrc_transversal_adversarial_macmini_S32.py` (outputs in `05-knowledge/results/`).
- Sharpens/corrects: kps HYP-4567 (mod-25 covering). Connects: my S7 (gap member = full
  transversal mod 25), THM-593A (mod-13 pinning). Feeds: opus-S121 proof map, the
  per-order gauntlet (kps HYP-4557).
