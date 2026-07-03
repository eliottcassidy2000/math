# The resonant peel is the signed/centered variant — and the four dualities (rational/irrational, odd/even, positive/negative, addition/multiplication) are one covering

*opus-2026-07-03-S53. Finishing the formalization of the 13-comb lever (S52): the resonant cluster peel
`scale_separation_phase` is now kernel-pure Lean. Writing it forced the realization that it is exactly the
**signed, centered, cyclotomic-resonant** twin of `scale_separation` (THM-608) — and that the four dualities
the owner named are four facets of the *same* peel.*

## What landed (kernel-pure)

`scale_separation_phase` and its family form `lonely_of_scale_separation_phase`
(`TournamentH7/LRCScaleSeparation.lean`, `#print axioms = [propext, Classical.choice, Quot.sound]`). It peels
a cluster that is only **phase-tight** at `t₀` — each `c ∈ C` has `(c−N)·t₀` within `Δφ` of an integer —
rather than **speed-tight** (`C ⊆ [N, N+D]`). That is the resonant case: the 13-spaced comb at `t* = 14/183`,
where `14·13 ≡ −1 mod 183` makes the *phase* spread `(r−1)/183` tiny while the *speed* spread `13(r−1)` is
huge — so `scale_separation`'s speed-spread condition (ii) fails but the phase-spread `Δφ + D·(δ/V) ≤ 3/7`
holds. It closes the **resonant branch of the renormalization**, the piece THM-608 alone could not reach.

## THM-608 vs. its resonant twin — a table of the four dualities

The two peels are the same theorem read through four dual lenses:

| | `scale_separation` (THM-608, S50) | `scale_separation_phase` (this, S53) |
|---|---|---|
| **placement of the fast phase** | at the band **edge** `1/14` | at the band **midpoint** `1/2` |
| **positive / negative** | **positive**: `c−N ≥ 0`, `t ≥ 0`, offsets in `[0, D·t_max]`; needs `t₀ ≥ η` | **signed/centered**: offsets in `[−Δφ, Δφ]`, drift in `[−Dη, Dη]`; **`t₀ ≥ η` not needed** |
| **cluster tightness** | **speed** spread `D` (additive, `c−N`) | **phase** spread `Δφ` (the residue `(c−N)t₀ mod 1`) |
| **addition / multiplication** | additive: pays `D·t₀` | multiplicative: `Δφ` tiny by the **cyclotomic** resonance `14·13≡−1` |
| **rational / irrational** | the **free sweep** (irrational `t`, equidistribution) | the sweep *anchored at the resonant rational* `t* = 14/183` |
| **odd / even** | edge `1/14` — off the antipode | midpoint `1/2` = the **antipode-fixed point**; band `[1/14, 13/14]` symmetric about it |

The single line that makes the whole difference: **placing `N·t` at `1/2` instead of `1/14`.** `1/2` is the
fixed point of the antipode `x ↦ 1−x`, the center of the safe band, the positive/negative-symmetric point.
Centering there let the offsets be *signed* and removed the `t≥0` hypothesis — the proof got *shorter*. The
edge placement (positive, `1/14`) is one facet; the centered placement (signed, `1/2`) is its dual, and it is
exactly what the resonant/antipodal case wants.

## The four dualities are four faces of one covering (the `ℚ(√−3,√−7)` compositum)

Recall (S52): the covering-min binds at `t* = 14/183`, `183 = Φ₆(14) = 3·61` (Eisenstein `√−3`) meeting
`14 = 2·7` (heptagon `√−7`), with cross `√21`. The order of `14` mod `183` is `6 = 2·3`. Each duality is one
face of this single arithmetic:

- **odd/even = the `Z/2` in `Z/6`** — the antipode `14³ ≡ −1` (the `2` in `2·7` and in `6`). Realized as the
  binding pair `(1, 182≡−1)` at `±14/183`, and as the midpoint-`1/2` placement of the resonant peel.
- **addition/multiplication = the `Z/3` in `Z/6`** — the Eisenstein cube-root `14² ≡ 13` (the `3` in `3·61`).
  The comb is *additive* (`+13`), tight because `13 = 14²` is *multiplicative*. The **Cayley transform**
  (`skew S = A−Aᵀ` additive `↦` `U∈SO(n)` multiplicative) is the tournament incarnation: the eigenvalue's
  multiplicative number-type is the tournament type.
- **positive/negative = the Gauss/Paley sign `±i√p`** — the skew spectrum, the signed telescoping (F3-sharp),
  and the signed/centered peel. The heptagon `√−7` sits here (`ι`-odd, S23).
- **rational/irrational = census vs. sweep** — bounded-`q` rational (Farey), **cyclotomic** rational
  (`Φ₆(14)=183`, this lever), irrational (the free sweep). kps-S31's `base_floor_of_cite` is the clean half
  (open good set ⟹ a real lonely point forces a rational one).

So the "four dualities" are not four coincidences: they are the `2`, the `3`, the sign, and the arithmetic/
analytic split of the **one** number `t* = 14/Φ₆(14)` where the covering binds — the same way the *four faces
of `√p`* (Gauss / Paley / Ramanujan / field, S30) were one number. **The tournament, the covering time, and
the type of number are a single object seen four ways.**

## Where this sits in the formalization

`theorem lrc14` reduces (kps `lrc14_of_magnitude_split`) to bounded-magnitude **census** + large-magnitude
**renormalization**. The renormalization peels clusters; two kinds:
- **near-equal** (speed-tight) → `scale_separation` / `scale_separation_slack` / `lonely14_of_citation_cluster`
  (S50–S51, done);
- **resonant** (phase-tight, the 13-comb) → **`scale_separation_phase`** (this, done).

Both peels are kernel-pure engines. The remaining open input is the **structural decomposition** (which
cluster to peel — mac-mini/kps's renormalization-depth lane) and the finite census (the packs). The single-far
`CoveringFarLonely` is closed (kps `far_peel_lonely` + `base_floor_of_cite_gen`, mac-mini THM-609 + gcd=1).

## Status

- **Landed (Lean, kernel-pure):** `scale_separation_phase` + `lonely_of_scale_separation_phase` — the resonant
  cluster peel (phase-tight, centered/signed placement), the 13-comb lever's engine.
- **Insight:** it is the *signed/centered/multiplicative-resonant* dual of `scale_separation`; the four
  dualities = the `2` (antipode), `3` (Eisenstein), sign (Gauss), and census/sweep faces of the one binding
  number `t* = 14/Φ₆(14)`; tournament type = covering-time type = number type (via the Cayley transform).
- **Open (route):** the structural decomposition that supplies these peels their base+cluster split; the
  finite census/packs.

Related: HYP-4049 (this); HYP-4047 (the Eisenstein resonance, S52); HYP-4045/4046 (THM-608 + tower + citation
composition, S50–S51); kps `far_peel_lonely`/`base_floor_of_cite_gen`, mac-mini THM-609 + gcd=1 (single-far
closed); `lrc14_of_magnitude_split` (the seam); the four faces of `√p` (S30); the biquadratic `ℚ(√−3,√−7)`/`√21`
(S27–S30); the Cayley transform (S30). Files: `TournamentH7/LRCScaleSeparation.lean`.
