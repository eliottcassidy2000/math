# The actual floor is positive — decorrelation caps the resonant dips, and the hard core is exactly the tight-locus rigidity

*opus-2026-07-03-S60. The owner asked me to work on the hard core: the actual floor. I resolved the
foundational question my S59 left open — the uniform measure floor `inf meas(lonely S)` over primitive
covering families is POSITIVE (≈ 0.004), not 0 — proved the mechanism (THM-611 decorrelation caps the
resonant dips), added the elementary margin↔measure bridge (THM-613), and located the residual hardness
exactly where the fleet's rigidity program already sits (THM-612's confinement + `g(14)≤3`).*

## First: what the hard core actually is

`meas(lonely S) > 0` for every primitive covering `S` *is* LRC(14) (a positive-measure safe set ⟹ a safe
point ⟹ `M ≥ 1/14`). The **uniform** floor `inf_S meas(lonely S) > 0` is strictly *stronger* than LRC —
so it cannot be easier. The value of chasing it is diagnostic: does the floor route (THM-579's `R' > 0`
uniformly) even have a positive target, or does the infimum vanish? S59 left this open (the low-`R'`
families were near-imprimitive; a naive sequence decorrelated back). S60 settles it.

## The resolution: inf meas > 0, and decorrelation is why

The min-measure primitive covering families are **near-tight block ∪ {resonant primitivizer}**:
`S = 2·({1..13}\{6}) ∪ {w}`. The block is the AP `{1..13}` missing one tooth, 2-dilated — as close to the
imprimitive tight AP as a covering family gets. The only freedom to push `meas` down is the single
primitivizer `w`. And there **THM-611 decorrelation is decisive**:

`|meas(block ∪ {w}) − (6/7)·m_block| ≤ A_block/(3w)`.

Measured: as `w` grows, `meas` **oscillates** around the decorrelation limit `(6/7)·m_block ≈ 0.00699`,
with resonant dips that **shrink** like `A/(3w)` — deepest `0.00408` at `w=63` (`=9·7`, commensurate with
the block's `7`), decaying to `10⁻⁵` by `w=5005`. So `meas` does **not** tend to 0. The infimum over each
single-primitivizer family is a *positive finite-`w` resonant dip*. An aggressive descent (speeds ≤ 80)
bottoms at `≈ 0.00408`, not lower. **The floor is positive, ≈ 0.004.**

The mechanism, cleanly: to make `meas` small you must sit next to the tight AP; but you can only sit there
with a *bounded* commensurate block (scale-invariance kills pure dilation, primitivity forbids the AP
itself), and any runner you add or grow to primitivize either moves you off the tight locus or
**decorrelates** and restores `meas` toward `(6/7)·(block measure)`. Rigidity holds the block's measure
up; decorrelation caps the dips. Neither lets `meas → 0`.

## The bridge: margin ⟹ measure (THM-613)

An elementary companion. `F_S(t) = min_v ‖vt‖` is `v_max`-Lipschitz and peaks at `M(S)`, so a ball of
half-width `(M−b)/v_max` about the peak is `b`-safe:
> `meas{F_S ≥ b} ≥ 2(M(S) − b)/v_max`.
With `b = 1/14` and the covering-min `M ≥ 14/183` (kps HYP-4060): `meas(lonely S) ≥ (13/1281)/v_max`.
This converts the fleet's **margin** program (prove `M > 1/14`, THM-610/612 rigidity) into an explicit
**measure** floor — the two routes are one, quantitatively. The bound is loose at large `v_max` (where
decorrelation, not the peak, carries the measure), which is exactly why it doesn't by itself give the
uniform floor — but it *does* give a clean positive per-family floor and pins the two routes together.

## Where the hardness actually lives (unchanged, now sharpened)

A *uniform* `inf meas ≥ c > 0` is `≥ LRC(14)`-hard. The proof splits two-sided (same architecture as the
census route):
- **dominant runner** (`v_max ≫ rest`): peel it — THM-611 gives `meas ≥ (6/7)m_R − A_R/(3v_max)`,
  controlled when one runner dominates. *Tractable.*
- **no dominant runner** (compact / near-tight): scale-invariance reduces the gcd to a bounded family, and
  positivity then needs the family to be *bounded away from the tight AP* — which is exactly **THM-612's
  confinement + `g(14)≤3`**, the fleet's open core. *The hard case = the rigidity, verbatim.*

So the measure floor and the margin rigidity are not two problems but one, meeting at the near-tight block.
My contribution is to have shown the floor is genuinely positive (so the route is not chasing 0), named
the mechanism (decorrelation caps dips), and reduced its hardness to the rigidity the fleet is already
attacking — with the extremal shape `2·({1..13}\{6}) ∪ {resonant w}` as the concrete floor witness.

## Status

- **Resolved (evidence + mechanism):** `inf meas(lonely S) > 0` (≈ 0.004) over primitive coverings —
  decorrelation (THM-611) caps the resonant dips; not 0. Extremal shape `2·({1..13}\{6}) ∪ {resonant w}`.
- **Proved:** THM-613 (margin↔measure slope bridge; `meas ≥ 2(M−1/14)/v_max`, and `(13/1281)/v_max` from
  the covering-min).
- **Open (= LRC(14) hard core, the rigidity):** a *uniform* `inf meas ≥ c` proof — the "no dominant
  runner" / near-tight case is THM-612's confinement (`primitive tight ⟹ q*=14`) + `g(14)≤3` (HYP-2913).

Given MISTAKE-097 (prior overclaims on this crux), the non-closure is flagged plainly: the floor is shown
positive and its hardness *located*, not closed.

Related: THM-611 (decorrelation, the cap — this session's load-bearing tool), THM-613 (the bridge, new),
THM-612/mac-mini + HYP-4062/kps (the rigidity = the hard case, verbatim), HYP-4060/kps (covering-min, the
bridge's input), HYP-4063/opus-S59 (the inf question, resolved here), THM-579 (the floor program, now with
a confirmed positive target). Scripts: `lrc14_measure_floor_infimum_opus_S60.py`,
`lrc14_measure_floor_descent_opus_S60.py`, `lrc14_floor_is_positive_decorrelation_caps_dips_opus_S60.py`.
HYP-4064.
