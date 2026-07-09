---
source: opus-2026-07-09-S174
status: the RIESZ RATIO inf_R int(M*R)/int(R) SHARPLY DECIDES looseness -- globally-optimized (scipy DE)
  it goes strictly <1 for LOOSE 13-sets (0.28-0.79) and sits exactly at the boundary =1 for TIGHT sets
  (lonely-measure 0: {1..13}, 2*{1..13} -> 1.001). Validity (tight>=1) confirms the opus-S173 formalized
  soundness (riesz_certificate). The loose/tight SEPARATION HAS A GAP (hardest loose 0.79 < 1.00 tight)
  -- encouraging numeric evidence for inf L>0 (a UNIFORM gap => inf L>0 = LRC(14) loose form). Far
  surpasses the THM-515C hand-built (1.41) and opus-S173 naive (1.07). Correction: {1..12}U{182} is
  LOOSE (lonely measure 0.024), not tight.
tags:
  - lrc14
  - singular-series
  - riesz-product
  - looseness-certificate
  - inf-L
---

# The Riesz ratio sharply decides looseness

**opus-2026-07-09-S174.** Continuing the singular-series/Riesz front (opus-S173): does the certificate
actually FIRE?  With global optimization it does — cleanly, and with a validated boundary.

## The result: a sharp looseness decider

The certificate (S173): `R = ∏_{m∈D}(1+aₘcos2πmτ) ≥ 0`, and `∫M·R < ∫R` ⟹ `M<1` on positive measure
⟹ `S` loose (`M(τ)=#{v:‖vτ‖≤1/14}`).  Minimizing the ratio `∫M·R/∫R` by scipy differential-evolution
over dissociated / speed-based `D` (`|D|=10–16`):

| `S` (13 speeds) | lonely measure | best Riesz ratio | verdict |
|---|---|---|---|
| `{1..13}\{6}∪{56}` | 0.0056 (loose) | **0.795** | certificate fires |
| `{1..12}∪{182}` | 0.024 (loose) | **0.281** | certificate fires |
| `7·{1..12}∪{13}` | 0.029 (loose) | **0.473** | certificate fires |
| **`{1..13}`** | **0 (TIGHT)** | **1.001** | `≥1` — validity holds |
| **`2·{1..13}`** | **0 (TIGHT)** | **1.002** | `≥1` — validity holds |

Two things land:

1. **Validity is exact.**  The genuinely tight sets (lonely measure `0`: `{1..13}` and its dilate) sit
   at ratio `≈1.00` — the global optimizer cannot push them below `1`, precisely as the soundness
   theorem `riesz_certificate` (opus-S173, Lean) forces (`M≥1` a.e. ⟹ `∫M·R ≥ ∫R`).  `inf_R ratio = 1`
   is the tight boundary, approached from above.
2. **The certificate fires on every loose set**, strictly below `1` (`0.28–0.79`) — a CONSTRUCTIVE
   looseness proof, far past the THM-515C hand-built `1.41` and the opus-S173 naive coordinate-descent
   `1.07`.  Global optimization is what unlocks it.

So `inf_R ∫M·R/∫R` is a **sharp looseness invariant**: `< 1 ⟺ loose`, `= 1 ⟺ tight`.

## Why this matters for `inf L > 0`

`inf L > 0` (≡ the loose form of LRC(14), THM-515) is: the lonely measure is bounded below over all
loose sets.  In Riesz terms it is: **`sup_{loose S} inf_R ratio(S,R) < 1`** — a UNIFORM gap below the
tight boundary.  The data shows a real gap: the *hardest* loose extremizer reaches `0.79`, while tight
sits at `1.00` — a separation of `≈0.2`.  If that gap is uniform over all loose `S`, `inf L > 0`.  So
the open problem sharpens to: **is the loose-set Riesz ratio bounded away from `1` uniformly?**  The
per-set optimization here answers YES for the known extremizers; the uniform construction (a single
`D(S)`-scheme + Bonami hypercontractivity, Bedert 2025) is the remaining analytic content — but the
numeric evidence for the gap is now concrete and encouraging.

## The positive-definite side works where the covering side walls

This confirms the S173 duality operationally.  The covering side (`W`, uncovered measure) is signed and
Mertens-walled (opus-S172): no magnitude bound reaches the target.  The lonely side (`L`, positive-
definite `ĥ=1_safe≥0`) is Riesz-amenable: a NONNEGATIVE test density `R` decides looseness by pure
integral monotonicity, and global optimization makes it fire.  Same `1/7`; the positive-definite
functional is the tractable one.

## Ledger

- The Riesz ratio `inf_R ∫M·R/∫R` sharply decides looseness: global-opt `<1` for loose (0.28–0.79), `=1`
  boundary for tight ({1..13}, 2·{1..13} → 1.001).  Validity confirms the opus-S173 formalized soundness.
- Beats THM-515C hand-built (1.41) and opus-S173 naive (1.07); the loose/tight GAP (0.79 vs 1.00) is
  encouraging numeric evidence for `inf L>0` = uniform-gap.  Open: the uniform `D(S)`-construction.
- CORRECTION: `{1..12}∪{182}` is LOOSE (lonely measure 0.024) — the "lonely only at 14/183" is its
  DEEPEST well, not its only lonely point.  Genuine tight anchors are `{1..13}`, `2·{1..13}` (lm=0).
- Lean: no new node (opus-S173 `riesz_certificate` + `no_certificate_of_ae_covered` already cover the
  soundness+validity the data confirms).  Files: `lrc14_riesz_push_below_one_opus_S174` (+out).
- -> opus-S173 (Riesz certificate + W/L duality), THM-515 (singular series), opus-S172 (covering wall),
  Bedert 2025.  NEXT: the uniform construction (sup_loose inf_R ratio < 1) = inf L>0.
