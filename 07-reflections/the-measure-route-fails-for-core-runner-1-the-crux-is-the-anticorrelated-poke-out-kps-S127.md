# The measure route fails for core-runner-1: the smooth good set is anti-correlated and small; the crux is the poke-out

*kind-pasteur-2026-07-11-S127 cont.63. Owner: "attack the general core-runner-1 case with arbitrary smooth body."
I tried the cleanest possible sufficient condition — a measure bound that would sidestep opus's discrepancy
machinery entirely — and it fails, for an informative reason. The finding rules out the measure route, quantifies
why (the smooth body's good set is strongly anti-correlated and small), and reframes the crux as a covering
problem. The positive content: the difficulty is irreducibly a fine *location* fact, exactly opus's target.*

---

## The setup and the clean idea

A core-runner-1 covering family is `{1} ∪ B`, `B` = twelve non-core runners (each divisible by a prime `≤ 13`).
Runner 1 (speed 1) is bad on a **single arc** `D₁ = (−1/14, 1/14)`, measure `|D₁| = 1/7`. The body's good set is
`G' = {t : ‖bt‖ ≥ 1/14 for all b ∈ B}`. Since

> if `|G'| > |D₁| = 1/7`, then `G'` cannot fit inside `D₁`, so `∃ t ∈ G'∖D₁` where every runner (incl. runner 1)
> is lonely ⟹ `M({1}∪B) ≥ 1/14`,

a **measure bound `|G'| > 1/7`** would close the whole `|core| = 1` residual for LRC(14) — no Erdős–Turán, no
mollification. It is the first thing to try.

## It fails — and the reason is structural

Computed `|G'|` at level `1/14` for the smooth bodies (grid `6·10⁵`):

| body `B` | `|G'|` | `> 1/7 = 0.143`? |
|---|---|---|
| deep-well body `{2..12, 182}` | 0.0851 | **no** |
| `{2..12, 364}` | 0.0890 | no |
| `{2,4,…,14, 3,9,5,11,13}` | 0.0839 | no |
| `{4,6,8,9,10,12,14,15,21,22,26,33}` | 0.1077 | no |

**Every body has `|G'| ≈ 0.08–0.11 ≪ 1/7 = 0.143`** — the good set is far *smaller* than runner-1's arc. The
independent model would give `(6/7)¹² = 0.157 > 1/7`, but the real bodies sit at `≈ 0.54×` that: the twelve
`13`-smooth runners are **strongly anti-correlated** — their bad sets (centred at `k/bᵢ` for small primes)
interleave and cover *efficiently*, shrinking the good set well below independence. So the measure bound has the
sign **backwards**: `|G'| < |D₁|`, and `G'` fits inside the arc with room to spare, measure-wise. The measure
route is dead.

## What the failure reveals — the poke-out is a location fact

The good set is small and lies mostly *inside* runner-1's arc; loneliness of `{1}∪B` comes from the thin part of
`G'` that **pokes past** `‖t‖ = 1/14`. Quantifying `coreCover = |G' ∩ D₁| / |G'|` (want `< 1`) on covering
bodies: it runs `≈ 0.72–0.74` on the ones tested here (margin `1−coreCover ≈ 0.26`), up to `≈ 0.92` on the
tightest (mac-mini-S74). In *every* case `coreCover < 1` (LRC holds), but the surplus is a fraction of a good set
that is itself only `~0.085` wide — an absolute poke-out of order `10⁻²`. There is no measure cushion; the
positive margin is entirely about **where** the small good set sits relative to `1/14`.

This is precisely why the analytic route needs discrepancy/mollification (opus-S259/S260/S261) and not a measure
bound: proving `coreCover < 1` = proving the anti-correlated good set genuinely crosses `‖t‖ = 1/14` for **all**
smooth bodies, a fine location statement, not a mass statement.

## The reframe — a covering problem

Contrapositive of `coreCover = 1`: whenever `‖t‖ ≥ 1/14`, some `bᵢ` is bad. That is,

> **`coreCover = 1` ⟺ the twelve smooth bad sets `⋃ᵢ {‖bᵢ t‖ < 1/14}` COVER runner-1's good region
> `[1/14, 13/14]`; LRC(14) for `{1}∪B` ⟺ they leave a gap.**

So the `|core| = 1` crux is exactly: *twelve `13`-smooth-divisible runners cannot cover the interval
`[1/14, 13/14]` with their `1/14`-collars.* Total collar measure is `12·(1/7) = 12/7 = 2×` the interval, so
measure permits covering — the obstruction is arithmetic (the collars are centred at `k/bᵢ`, forced to pile up
near `0`, so their coverage of the *far* interval is inefficient). This is a covering-system-flavoured statement
(cf. the project's `q`-covering language), a possible alternative to the Fourier route, and it makes the smooth
structure the load-bearing hypothesis: the `13`-smooth centres cannot tile `[1/14,13/14]`.

## Net

The measure sufficient condition `|G'| > 1/7` is **refuted** (`|G'| ≈ 0.085 < 0.143` for every smooth body,
because the smooth good set is anti-correlated to `~0.54×` independent) — a clean dead end that saves the fleet
from a natural but doomed approach. The residual is irreducibly the fine location fact "the small anti-correlated
good set pokes past `1/14`," i.e. opus's discrepancy target, equivalently "twelve `13`-smooth collars cannot
cover `[1/14, 13/14]`." The structural half (the extremal deep well, the single-killer ladder) stays pinned and
Lean-formalized; this half is analytic and now has its natural measure shortcut ruled out.

*Files: lrc14_core1_good_measure_kps_S127.py (+.out), lrc14_core1_margin_kps_S127.out. Complements opus-S259/260/261
(discrepancy/mollification — now with the measure route eliminated + the covering reframe), mac-mini-S74/75
(runner-1), kps cont.62 (core-runner-1 residual). HYP-6232.*
