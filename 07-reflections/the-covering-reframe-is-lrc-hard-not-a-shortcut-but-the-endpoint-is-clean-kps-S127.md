# The covering reframe is LRC-hard, not a shortcut; but the endpoint structure is clean

*kind-pasteur-2026-07-11-S127 cont.64. Owner: "work the covering reframe: 12 smooth collars can't tile
[1/14,13/14]." I tested it for a genuine obstruction — measure budget, structure, a universal gap point — and it
has none: the reframe is equivalent to the `|core|=1` residual, as hard as LRC itself. The honest value is
ruling it out as an easier route, plus one clean local structure (the endpoint) and the confirmation that the
single-killer slice is already proved in covering language.*

---

## The reframe (cont.63)

`{1}∪B` is lonely (LRC(14)) ⟺ the twelve non-core collars `Cᵢ = {t : ‖bᵢ t‖ < 1/14}` leave a **gap** in runner-1's
good region `[1/14, 13/14]`. The hope: the smooth structure of the `bᵢ` (each divisible by a prime `≤ 13`)
obstructs the covering in a way a plain loneliness argument would miss.

## No obstruction — three ways it stays hard

**(1) Measure permits covering.** Each collar has measure `2/14 = 1/7`; twelve of them total `12/7 = 1.714`, vastly
more than the interval `|[1/14,13/14]| = 12/14 = 0.857`. Computed the actual covered/gap split inside the interval:
covered `≈ 0.83`, gap `≈ 0.005–0.03`. So the collars have a `2×` measure surplus — **no budget obstruction**; the
gap is a thin residue, not a mass deficit. (Same wall as the measure route, cont.63.)

**(2) The gap location varies — no universal gap point.** Single-killer `{1..12,182}`: gap at `14/183 ≈ 0.0765`
(just past `1/14`). Multi-killer `{1..11,13,84}`: gap at `37/89 ≈ 0.416` (mid-interval). So there is no fixed `t*`
that is uncovered for every covering family — the gap is genuinely family-specific, exactly as loneliness is.

**(3) The interval can't be shrunk uniformly.** Because the gap ranges over `[14/183, ~1/2]`, no proper
sub-interval of `[1/14,13/14]` contains all gaps; the covering must be defeated on the whole interval.

Together: the covering reframe is **equivalent to the `|core|=1` residual and as hard as LRC** — not an easier
formulation. The crux stays the fine loneliness / discrepancy that opus's Fourier route (S259–S262, mollified
Erdős–Turán, klein's Weyl) targets.

## The one clean structure — the endpoint

Just above the endpoint `1/14`, coverage has a rigid form: a runner `b` has an arc centred near `1/14` iff
`b·(1/14)` is near an integer, i.e. `b ≡ 0 (mod 14)`. So **near `1/14`, only the multiple-of-14 runner covers** —
verified: for the deep-well body only `182` is bad at `1/14` (reaching to `1/14 + (1/14)/182 ≈ 0.0718`); for the
`S₂` body only `364`. This is why the *single-killer* gap sits just past `1/14`: the lone mult-of-14 runner covers
a shrinking sliver `[1/14, 1/14 + 1/(196m)]` (`m` = the multiple), and immediately past it the gap opens (`14/183`
for `m=13`). It is a genuine, clean, *local* handle — but only local: the mult-of-14 arc is tiny, and for
multi-killer families the gap has moved to mid-interval where this structure says nothing.

## What is actually settled here — the single-killer slice, in covering language

The reframe does re-express one proved fact cleanly: for the single-killer ladder `{2,…,12, 182c}`, the collars
leave the gap at `14c/(182c+1)`, and this is exactly the machine-checked `reach ≥ 14c/(182c+1)` (Lean cont.60/61,
`LRCSingleKillerLadder`). So *"the single-killer body's collars cannot tile `[1/14,13/14]`"* is theorem, not
conjecture — the covering reframe's structured extremal slice is closed. The open part is the same as everywhere:
the general smooth body (multi-killer, arbitrary `13`-smooth non-core), where the gap is mid-interval and only the
analytic route bites.

## Net

The covering reframe is a faithful restatement of the `|core|=1` residual, not a shortcut: measure permits
covering (`2×` surplus), the gap point is family-specific (no universal obstruction), and only the tiny endpoint
neighbourhood has rigid structure. Honest outcome — it is **ruled out as an easier route** (saving the fleet a
plausible detour), the endpoint fact is documented as the clean local handle (explaining the near-`1/14`
single-killer gaps), and the single-killer slice is confirmed proved in covering language. The live route remains
opus's mollified-Fourier discrepancy over the smooth body.

*Files: lrc14_covering_reframe_kps_S127.py (+.out). Follows cont.63 (the reframe), cont.60/61 (single-killer Lean);
complements opus-S259–S262 / klein-S278 (the Fourier/Weyl route that this shows the reframe does not bypass).
HYP-6234.*
