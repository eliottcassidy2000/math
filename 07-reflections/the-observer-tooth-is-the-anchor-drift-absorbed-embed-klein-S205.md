---
source: klein-2026-07-09-S205
status: hembed's realization ADVANCED and FORMALIZED on the large-ruler regime. LRCDriftEmbed.lean
  (sorry-free, kernel-pure) proves the DRIFT-ABSORBED ruler embedding: a tooth-free gap `(a,a+g) ⊆ [0,1]`
  at ruler index `j` with `g > 1/7 + 2·spread·φ/Vmax` (`φ = a+g/2`) yields a lonely time
  `τ = (j+φ)/Vmax` with `minReach ≥ 1/14` — WITH the `e_i = Vmax − v_i` binding kps-S105 correctly
  flagged as missing from klein-S203/S204's abstract `hembed`. CREATIVE CORE: the drift is
  **φ-PROPORTIONAL** (`d_i = e_i·φ/Vmax`, not merely `≤ spread/Vmax`), and the observer's own tooth
  `c_0 = 0` (since `e_0 = Vmax − Vmax = 0`) is an **ANCHOR that never drifts** — so the drift-optimal gap
  is the one abutting it. Verified: the anchor gap wins in every tested cluster (`a = 0`, `g = 1/2`,
  `φ = 1/4`), lowering the firing threshold from `Vmax > 3.5·spread` to `Vmax > 1.41·spread` (2.4×).
  HONEST SCOPE: `1.41 > 7/6`, so the embed does NOT reach the hard regime; the residual is the BOUNDED
  window `Vmax ∈ (spread, 1.41·spread]` — a finite check, not open analysis.
tags: [lrc14, thm-527-part-a, hembed, ruler-embedding, drift, anchor, lean, slow-fast]
---

# The observer's tooth is the anchor — and the drift is proportional to the fast phase

**klein-2026-07-09-S205.** Owner: work toward the finish target creatively; mine past threads. The finish
target (klein-S204) is criterion-C's REALIZATION: exhibit a real `τ` whose fast phase lands in the
good-period gap. Triple convergence (klein-S204 / kps-S105 / opus-S176) had already reduced it to the
**tooth wobble**. This session turns that into a proved, formalized embedding — and finds the structure
that makes the wobble small.

## The setup, and the exact drift

Put the witness time on the `Vmax`-ruler with a FREE fast phase, `τ = (j + φ)/Vmax` (`j : ℤ`, `φ ∈ [0,1)`).
Then `frac(Vmax·τ) = φ` — the fast phase is ours to choose. With the co-offset binding `v_i = Vmax − e_i`,

> `v_i·τ = (φ − t_i − d_i) + j`,  `t_i := e_i·j/Vmax` (the tooth), **`d_i := e_i·φ/Vmax` (the DRIFT)**.

So the fast phase must clear each tooth *after* it has drifted. Since `nearInt` is `ℤ`-periodic and
characterised by `∀ n : ℤ, r ≤ |y − n|` (kps-S31 `GapReach`), placing `φ` at the midpoint of a tooth-free
gap of width `g` gives margin `g/2`, and the drift eats at most `|d_i|`. Hence, with **no Lipschitz
estimate**, just a triangle inequality:

> **`g > 1/7 + 2·max_i|d_i|` ⟹ every runner is `> 1/14` from the origin at `τ = (j + a + g/2)/Vmax`.**

Formalized: `LRCDriftEmbed.{minReach_ge_of_driftGap, Mreach_ge_of_driftGap}` (sorry-free, kernel-pure).
**It carries the `e_i = Vmax − v_i` binding** — the omission kps-S105 correctly called out in my S203/S204
abstract `hembed` (which, unbound, was false in isolation: pick `E` unrelated to `v`).

## The creative half: the drift is φ-proportional, and the observer is an anchor

The fleet's bound was `|d_i| ≤ spread/Vmax` (kps/opus's "tooth wobble"), obtained from `φ ≤ 1`. But the
drift is **exactly** `e_i·φ/Vmax`. Two consequences the coarse bound hides:

1. **Low fast phases drift less.** `|d_i| ≤ spread·φ/Vmax`. Placing `φ` low shrinks the drift *linearly*.
2. **The observer's tooth is an ANCHOR.** The observer is the runner `v = Vmax`, so its co-offset is
   `e_0 = Vmax − Vmax = 0`; its tooth sits at `c_0 = frac(0·j/Vmax) = 0` for EVERY `j`, and **it never
   drifts** (`d_0 = 0`). So `0` is a permanent tooth, and the gap immediately above it starts at `a = 0` —
   the lowest gap available, hence the least-drifting placement.

The observer's own safety forces `φ ≥ 1/14` (it must clear its own tooth at `0`), so the least achievable
drift is `spread/(14·Vmax)` — a **14× floor** below the naive `spread/Vmax`.

**Verified** (`lrc14_drift_optimal_gap_klein_S205`): the anchor gap wins in *every* tested cluster —
`a = 0.0000`, `g = 1/2`, `φ = 1/4` — and the firing threshold drops from

> `Vmax > 3.5·spread` (coarse) to **`Vmax > 1.41·spread` (sharp)** — a `2.4×` reduction,

with the constructed `τ` giving `minReach ≥ 1/14` in 100% of cases (`lrc14_drift_embed_verify_klein_S205`;
soundness of the Lean conclusion). The realized gain is `2.4×`, not `14×`, because the winning anchor gap
has `g = 1/2` (so `φ = 1/4`, not `1/14`) — the `14×` is the floor, attained only by a razor-thin anchor gap.

## Honest scope, and what remains

`1.41 > 7/6 = 1.167`, so **the drift embed does not reach the hard regime** (`Vmax ≤ 7·spread/6`, where
`j=1` fails). What it does deliver is the complementary, and for THM-527 the *relevant*, half: under the
**bounded-spread compact reduction** (THM-527's own title), `Vmax → ∞` with `spread` bounded, so
`Vmax > 1.41·spread` eventually holds and the embedding fires. The residual `hembed` corner is therefore
the **bounded window `Vmax ∈ (spread, 1.41·spread]`** — a finite check, corroborating kps-S105/opus-S176's
"formalization gap, not open analysis," and sharpening their cited `V* ≤ 234 / 1106 / 3^12` to an explicit
`⌈1.41·spread⌉` once the compact reduction bounds the spread.

## Mined threads (the tangential concepts that paid)

- **`ScaleSeparation.scale_separation_phase`** (the existing sorry-free embedding kps pointed to) places the
  fast phase at the band midpoint `1/2` and absorbs drift via `Δφ + Dd·(δ/V) ≤ 3/7`. Its condition fails
  exactly at `spread = 6Vmax/7` — mac-mini's knife-edge. My embed is the *gap-centred* generalisation
  (arbitrary `φ`, not `1/2`), which is what lets the anchor's low `φ` be exploited.
- **mac-mini's "0-neighbourhood" arc and the `j=1` wraparound.** The arc mac-mini kept finding at `0`
  (`maxgap = 1 − spread·x`, the `6/(7·spread)` arc) is *the anchor gap*. Its geometric cause is now named:
  `e_0 = 0` pins a tooth at the origin forever, so the origin-adjacent gap is a structural invariant of the
  co-offset picture, not an artefact.
- **The `13`-comb Eisenstein resonance** (`t* = 14/183`, `scale_separation_phase`'s docstring): phase spread
  tiny despite huge speed spread. Same phenomenon, dual form — here the *drift* is tiny despite huge spread,
  because `φ` (not `spread`) is the small factor.
- **THM-527's bounded-spread compact reduction** is exactly what makes `Vmax > 1.41·spread` generic.

Files: `LRCDriftEmbed.lean` (built, sorry-free, kernel-pure); `lrc14_drift_embed_verify_klein_S205`,
`lrc14_drift_optimal_gap_klein_S205` (+outs). Builds on `LRCCriterionC` (klein-S204), `LRCGapReach`
(kps-S31), `LRCMreachConcrete`. Fixes the binding kps-S105 flagged. Remaining: the bounded corner
`Vmax ∈ (spread, 1.41·spread]` (finite check) — and, for the good-period leg's hard regime, that corner IS
the whole hard regime, so the finite check is where the covering case's realization now lives.
