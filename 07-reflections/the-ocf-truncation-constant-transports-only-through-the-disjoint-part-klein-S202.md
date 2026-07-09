---
source: klein-2026-07-09-S202
status: Worked opus-S171's NEXT critical-path item — transport THM-076's clean-truncation constant to the
  covering-side residual R_abs, to turn kps-S96's E_grid good-period route a-priori. RESULT: the constant
  transports ONLY PARTIALLY. The tournament Walsh-OCF (THM-076) is a DISJOINT cycle-covering — it telescopes
  cleanly (constant-term identity) to 2^r(n-2k)!/2^{n-1}. The LRC covering is an OVERLAPPING covering
  (k·(1/7)=13/7=1.86>1, "barely-covers"). The constant transports through the DISJOINT/smoothness part —
  Φ(y)=uncovered-measure is continuous piecewise-linear ⇒ Φ̂(M)=O(1/M²), giving the RIGOROUS
  R_grid_abs ≤ TV(Φ')/(12Vmax²), TV(Φ')≈12·spread² — but this closes only the large-ruler regime Vmax≥3·spread
  (redundant with j=1, LEM-010(i)). It is ~30-50× too loose at the HARD boundary (spread≈Vmax): there the
  OVERLAP/barely-covers cancellation is essential — the multivariate absolute grid-residual sum OVERSHOOTS
  E[W] by 1.4-1.7× (tightAP 11-14×). Genuine contribution: the SOUND split R = R₀(signed=density floor) +
  R_grid(absolute), and the corrected averaging threshold E_grid[W] > (6/7)/Vmax (not >0). Corrects kps-S98
  total-absolute + opus-S171 R_abs<(6/7)^k over-claim; respects the tight AP (klein-S201).
tags: [lrc14, good-period, egrid, walsh-ocf, transport, barely-covers, disjoint-vs-overlapping, density-floor, cross-thread]
---

# The OCF truncation constant transports only through the disjoint part

**klein-2026-07-09-S202.** Owner: work the transport of the OCF/decay truncation constant, then formalize.
opus-S171's flagged NEXT: *"transport THM-076's truncation constant to R_abs ⟹ a-priori theorem"* — turn
kps-S96's `E_grid[W] = (6/7)^k + R > 0 ⟹ good period` from a verified bound into a proved one, by porting the
tournament Walsh-OCF factorization's clean truncation to the covering-side residual. I worked it to the end.
It transports **exactly as far as the covering is disjoint** — and no further.

## The two sides, and the one place they differ

opus-S171's cross-thread connection: the covering master law `𝒲̂(m) = Σ_{balanced n} ∏ĝ_{n_i}` (LEM-011/007)
IS the tournament Walsh–OCF factorization `|Î[S]| = 2^r·(n−2k)!/2^{n−1}` (THM-076) — both a Fourier
coefficient as a signed sum over coverings of a product. On the **tournament** side it truncates CLEANLY: the
constant-term identity `C(m,j)·j!·(m−j)! = m!` makes every covering term contribute equally, the sum
telescopes to a single product, higher Walsh degrees follow from degree-0. The reason is structural: tournament
cycle-coverings are **vertex-disjoint** — a genuine partition. The covering "count" has no overlap.

The **LRC** side covers `[0,1)` with `k = 13` arcs each of length `1/7`. Total mass `k/7 = 13/7 = 1.86 > 1`:
the arcs **cannot be disjoint** — this is the *barely-covers* regime. mac-mini-S57 proved (exhaustively) the
covering master law does NOT telescope: support-truncation GROWS (`Σ_m|Ŵ_r|² = 0.077, 0.226, 0.932,…`), the
sum converges only by extreme signed cancellation. opus-S171 named the gap exactly — *"the one place they
differ (tournament had no k/7>1 barely-covers obstruction) is where to look."* That is where the constant dies.

## What DOES transport: the smoothness (disjoint) half → a rigorous O(1/Vmax²)

The clean part of the OCF structure is the **piecewise-linear/telescoping smoothness**, and it does port. The
single-variable covering profile `Φ(y) = ` (uncovered measure of `{frac(e_i y)}`) `= Σ(gap−1/7)_+` is
**continuous piecewise-linear** in `y` (gaps vary continuously; slopes ±1). So its Fourier coefficients obey
the standard two-integration-by-parts bound `|Φ̂(M)| ≤ TV(Φ')/(2πM)²`, and the grid residual

> **`R_grid_abs := Σ_{m≥1} 2|Φ̂(mVmax)| ≤ TV(Φ')/(12·Vmax²)`**  (rigorous, `O(1/Vmax²)`),

with `TV(Φ') = ` total variation of `Φ'` = the explicit covering **truncation constant**, measured `≈ 12·spread²`
(`TV/spread² = 11.4, 11.8, 12.2, 13.9` for spread `12→129`; provably `O(k²·spread²)` from `O(k²·spread)`
breakpoints × `O(spread)` slope jumps). This turns opus-S171's *verified* `R_abs` into a *rigorous* decay bound
— the OCF constant, transported.

## Why it only closes the redundant regime

`TV(Φ') ≈ 12·spread²` makes the bound `R_grid_abs ≤ (spread/Vmax)²`. Setting it below the covering density
`E[W] ≈ 0.13` needs `Vmax ≳ spread/√0.13 ≈ 2.8·spread`. But **hard** clusters (where a good period is
non-trivial: `j=1` fails) have `spread ≥ 6Vmax/7`, i.e. `Vmax ≤ 1.17·spread`. So the transported bound closes
only `Vmax ≥ 3·spread` — which is `spread ≤ Vmax/3 < 6Vmax/7`, i.e. **`j=1` already works** (LEM-010(i)).
The rigorous transport is redundant with the trivial period. At the hard boundary (`spread ≈ Vmax`) the bound
is `≈ 1`, while the true `R_grid_abs ≈ 0.01–0.02` — **30-50× loose**. The looseness is precisely the
cancellation the TV bound cannot see: the multivariate ABSOLUTE grid-residual `Σ_{Vmax|n·e} 2|𝒲̂(n)|`
(exact LEM-011 magnitudes, no cancellation) **OVERSHOOTS `E[W]` by 1.4-1.7×** (7-struct `1.78`, dissoc `1.23`,
tightAP `11-14×`). The overlap's essential cancellation is exactly what a clean constant would have removed —
and on the LRC side it cannot be removed a-priori. **The transport is as partial as the covering is
non-disjoint.**

## The genuine contribution: keep R₀ signed (the density floor), take only R_grid absolute

kps-S98 tested the TOTAL absolute `Σ_{Vmax|n·e}|𝒲̂(n)|` and found it exceeds `(6/7)^k` at small/mid spread
(`1.55@s50`), reverting to LEM-013. That took EVERYTHING absolute — including the `n·e=0` exact relations,
whose absolute sum is the AP's own decorrelation and is genuinely large. **Split instead:**
`R = R₀ + R_grid`, `R₀ = Σ_{n·e=0}𝒲̂(n) = E[W]−(6/7)^k` (the exact relations = the **density floor**; keep
SIGNED, bounded because `E[W] ≥ 0`), `R_grid = Σ_{n·e=mVmax,m≠0}𝒲̂(n)` (pure grid; bound absolutely). Then,
accounting for the trivial shift `W(0) = 6/7` (`Σ_{j=1}^{V−1}W(j/V) = V·E_grid[W] − 6/7`), a good period
exists once

> **`R_grid_abs < E[W] − (6/7)/Vmax`**  (SOUND: verified 0 over-claims over ~1400 adversarial hard clusters).

This is the corrected averaging threshold — `E_grid[W] > (6/7)/Vmax`, **not** `E_grid[W] > 0` (which is
trivially true from the `j=0` term and was the latent flaw in the "E_grid>0 ⟹ good period" framing). It fires
iff a good period exists, and correctly does NOT fire for the tight AP `{0..12}@V=13` (`E_grid[W] = 6/7/13`
exactly, no good period — klein-S201). Keeping `R₀` signed sidesteps kps-S98's wall for the part that admits it.

## Triple convergence (independent, same day) — and this session's unique piece

Two other agents reached compatible conclusions the same day, independently:
- **opus-S172** (the transport's author) worked the SAME port and reached the SAME honest negative: the OCF
  clean truncation does NOT transfer to an a-priori `|R| < (6/7)^k`. opus measured `TV(W') ≈ 13·spread²`
  (exponent 2.03) — matching my `≈ 12·spread²` — so `TV/(12Vmax²) ≈ 1.1 ≈ 8·main`, useless; the actual
  `R_abs` is `40-200×` below by cancellation; `k/7 = 13/7 > 1` over-covering (`~26·spread` order-swap
  breakpoints, jumps `~spread/2`) IS the obstruction, and the tournament OCF truncates cleanly precisely
  because it does NOT over-cover. Identical diagnosis, independently derived (strong cross-check). opus also
  ran the L² energy bound `E[W'²] ≈ spread²` — also `200-400×` short. HYP-5610, Lean `resonant_tail_le`.
- **mac-mini-S64** derived the SAME `R₀`-signed / `R_grid`-absolute split below.

**This session's UNIQUE piece** (not in opus-S172 or mac-mini-S64): the corrected `j=0` **Lean lemma**
`exists_good_of_grid_residual` — the threshold `E_grid[W] > (6/7)/Vmax` (not `> 0`), concluding a NONTRIVIAL
`j ≠ 0`. It fixes the latent bug in kps-S96's `LRCEgridExistence.exists_good_of_residual_small`, whose
conclusion `∃ j ∈ Finset.range N, 0 < W j` is a TAUTOLOGY (`j = 0` always works, `W 0 = 6/7`) — the S201/S202
`j=0` blind spot, now formalized out. Plus the multivariate-overshoot quantification (`1.4-1.7×`) as a direct
measure of the over-covering.

## Convergence: mac-mini-S64 derived the same split (independently, same day)

The owner steered mac-mini-S64 to "think R₀-signed/R_grid-absolute split" the same day; we converged on the
identical decomposition and the identical `j=0` correction (`strict good period ⟺ V·(E_x + R_grid) > 6/7`).
mac-mini supplied the complementary MAIN-TERM fact I had left as "E[W] ≥ floor": **`V·E_x[W] ≥ 5.65·(6/7)` for
every valid `V`** (0 failures, grows linearly — the density floor in units of the collapse term `6/7`), so the
surplus over the trivial `W(0)=6/7` is `≥ 4.65` units, and the residual `R_grid` is a bounded `~0.5`-unit
correction — comfortable in the WIDE regime (`spread > 6V/7`), with the knife-edge `spread = 6V/7` being j=1's
non-strict job (mac-mini-S64). mac-mini also placed `R₀` on the "winning side": `R₀` correlates with the
additive energy / kissing number (AP-max), so the AP-extremal resonances sit inside the GROWING `E_x`, and
"AP maximizes kissing uniformly in Vmax" (kps's worry) evaporates for the RESIDUAL (wraparound-only, no additive
triples). **This session's distinct additions to that shared split:** (i) the RIGOROUS residual bound
`R_grid_abs ≤ TV(Φ')/(12Vmax²)` — mac-mini observed decay; this is the explicit `O(1/Vmax²)` transport of the
OCF/smoothness constant; (ii) the barely-covers QUANTIFICATION (multivariate absolute overshoots `E[W]` by
`1.4-1.7×`) diagnosing why that bound is `~30-50×` loose at the hard boundary; (iii) the corrected Lean
`exists_good_of_grid_residual` (the `j=0` threshold, sorry-free), fixing the vacuous `∃ j ∈ range N` conclusion
in kps-S96's `LRCEgridExistence`. The two sessions together are the complete honest state of the residual route.

## Where this leaves the good-period leg (honest)

The residual/`E_grid` route is an **abundance** result (it counts good periods when they exist), not a
universal existence proof — many hard clusters (the dense/AP ones) have NO good period and are the density
floor's / exact-check's job (the klein-S201 2×2: {near-AP, dissociated} × {`V≥Q+1`, `V≤Q`}). The genuine
closure of the good-period leg is LEM-012 (near-AP, `V≥Q+1`) + LEM-013 (dissociated) via the **maxgap margin**,
plus density-floor/exact-check for the rest. The OCF-transport does not change that; it makes opus-S171's
verified `R_abs` bound *rigorous in the large-ruler regime* and *correctly scoped* (SOUND split) elsewhere, and
it **quantifies the barely-covers obstruction** (`1.4-1.7×` overshoot) as the reason the constant cannot fully
transport. The cross-domain lesson: **a covering's clean-truncation constant transports through its disjoint
(smoothness) part and dies on its overlap; `k/7 > 1` is exactly how much of the LRC covering is overlap.**

Files: `lrc14_grid_residual_split_klein_S202`, `lrc14_TV_transport_klein_S202`, `lrc14_grid_sweep_klein_S202`,
`lrc14_multivariate_gridres_klein_S202` (+outs). Connects: opus-S171 (transport target), kps-S96 (E_grid route),
kps-S98 (total-absolute not uniform), mac-mini-S57 (LEM-007 barely-covers wall), THM-076 (Walsh-OCF clean
truncation), LEM-011 (exact 𝒲̂), klein-S201 (tight-AP / j=0 / density-floor 2×2).
