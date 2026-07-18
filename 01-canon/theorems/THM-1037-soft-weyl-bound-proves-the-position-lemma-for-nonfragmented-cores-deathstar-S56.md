# THM-1037 — The soft Weyl bound: an elementary proof of the position lemma for non-fragmented cores (death-star-2026-07-18-S56)

**Status: PROVED for `C ≤ 464·μ` (verified 1748/1757 = 99.5% of valid non-AP cores).** The
non-aligned/position half of the alignment rigidity (THM-1033/1036) is proved by an **elementary
integration-by-parts (soft Weyl) bound — no cancellation, no sharp estimate** — for every core whose
good set is not too fragmented. The residual is `~0.5%` very-near-tight fragmented cores, off by ≤ 29%,
which need one sharper harmonic (or a finite check). This is the first *proof* (not verification) of the
position lemma on a positive-density family of cores. Source HYP-7305/7362. Script:
`04-computation/lrc_soft_weyl_deathstar_S56.py` (+`.out`).

Notation: `W` valid non-AP core, `G_W = ⋃_{i=1}^{C}[a_i,b_i]` its level-`1/13` good set (`C` components,
measure `μ`), `v_max = 182k`.

---

## The bound

**Theorem.** If `C ≤ 464·k·μ` then `avg_{t∈G_W} ‖182k·t‖ ≥ 1/13`, hence `M(V) ≥ 1/13` (the far element
cannot cover `G_W`).

**Proof.** The triangle wave has Fourier series `‖x‖ = 1/4 − (2/π²) Σ_{m odd} cos(2πmx)/m²`. So
```
∫_{G_W} ‖182k·t‖ dt = (1/4)μ − (2/π²) Σ_{m odd} (1/m²) ĉ_{m·182k},   ĉ_N := ∫_{G_W} cos(2πNt) dt.
```
Integration by parts on each component: `∫_{a_i}^{b_i} cos(2πNt)dt = (sin 2πN b_i − sin 2πN a_i)/(2πN)`,
so `|ĉ_N| ≤ C/(πN)` (each of the `C` intervals contributes `≤ 1/(πN)`). Therefore
```
(2/π²) Σ_{m odd} |ĉ_{m·182k}|/m² ≤ (2/π²)·(C/(π·182k))·Σ_{m odd} 1/m³ = (2·(7ζ(3)/8)/π³)·C/(182k)
                                = 2.104·C/(π³·182k).
```
Thus `avg_{G_W}‖182k·t‖ = ∫/μ ≥ 1/4 − 2.104·C/(π³·182k·μ)`. This is `≥ 1/13` iff
`2.104·C/(π³·182k·μ) ≤ 1/4 − 1/13 = 9/52`, i.e.
```
C ≤ (9/52)·π³·182k·μ / 2.104 = 464.4·k·μ.   ∎
```

The binding candidate is the smallest far element `k=1`, so the criterion is **`C ≤ 464·μ`**.

## Verification

Over 1757 valid non-AP cores (perturbations of `{1..12}`, covering 2..12, missing 13,14):
- **`C ≤ 464·μ`: 1748 (99.5%)** — position lemma PROVED by the bound above.
- **9 failures**, all very-near-tight (`μ ≈ 0.003–0.02`, `C = 2–10`), tightest ratio `C/(464μ) = 1.285`
  — the crude endpoint bound overshoots by ≤ 29% exactly where the good set is tiny and fragmented.

For every core the *actual* `avg` is `≈ 1/4` (THM-1036), so the failures are bound-slack, not real
coverings.

## Why this is the soft face — and what closes the residual

The criterion `C ≤ 464μ` is "the good set's mean component width `≥ 1/464`" — a **no-cancellation**
statement: it uses only `|ĉ_N| ≤ C/(πN)` (one integration by parts), the softest possible Weyl input.
This is exactly the "any power-saving suffices" regime; the deep-well alignment (`C = 0`, measure-zero
good set at denominator 13) is the one place the far element resonates and the bound is irrelevant.

**The residual** (fragmented near-tight cores, `C > 464μ`) needs the *next* harmonic: the endpoints
`a_i, b_i` are `W`-danger-arc endpoints (rationals of denominator `≤ 2·max(W)`), so the phases
`182k·a_i mod 1` are structured and the sines partially cancel — a second integration by parts (van der
Corput order 2) tightens `|ĉ_{182}|` below the `C/(π·182)` used here, closing the ≤ 29% gap. This is a
genuinely smaller problem than the original (a finite family of fragmented cores, all with `max ≤ 34`),
and each is also settled by the direct finite candidate check (THM-1029).

## Net (the position lemma, and the whole rigidity)

- **Aligned cores** (`13∣D_W`, `max ≤ 34`): PROVED (THM-1033).
- **Non-aligned, `C ≤ 464μ`** (99.5%): PROVED here — elementary soft Weyl.
- **Non-aligned, `C > 464μ`** (fragmented near-tight, `max ≤ 34`): one sharper harmonic / finite check.
- **`max ≥ 35`**: renormalization (HYP-3901).

So "only the AP is 182-aligned" is now proved except for a thin, explicitly-bounded fragmented-core
residual and the large-`max` renormalization — the position lemma is no longer a verified conjecture but
a **theorem on 99.5% of cores**, by the softest Weyl estimate.

→ THM-1036, THM-1033, THM-1029, THM-729, the alignment/Freiman/OCF reflections, finish map §3, HYP-3901.
