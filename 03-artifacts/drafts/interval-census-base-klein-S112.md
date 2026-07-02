# The interval-census base: the covering band has a uniform spectral gap

klein-2026-07-02-S112 (HYP-4017). Companion to `04-computation/lrc14_band_margin_analysis_klein_S112.py`
(output in `05-knowledge/results/lrc14_band_margin_analysis_klein_S112.out`).

## The finding

For every covering 13-subset S of [1..22] (the censused band, 31471 families):

    max_t min_i ||s_i t||  >=  1/12   —   uniform margin 1/84 above the band radius 1/14.

**Zero tight rows.** The minimum optimum over the entire band is exactly 1/12, attained at
(1,2,3,4,10,11,12,13,14,15,16,17,18) (same extremal row at W = 20 and W = 22). Margin
distribution: min 1/84 ≈ 1.19e-2, median ≈ 5.4e-2, max ≈ 0.24.

## Why the classical tight family does not obstruct

The tight 13-runner family (1,2,...,13) has optimum exactly 1/14 — but it MISSES q = 14,
so it is NON-COVERING and never enters the census. It is handled by the denominator sieve
(`sieve_one_div` at t = 1/14), where its tightness is exactly the missing-modulus mechanism:
||s/14|| = 1/14 iff s ≡ ±1 (mod 14). **The extremals live on the easy leg.** Inside the
covering band — the leg that carries all the remaining open content — the loneliness is
uniformly strict.

## Consequence: the ε-schedule base is unobstructed

kps-S17's design note: the analytic peel (far-element transport) cannot consume bare
point-loneliness; the induction invariant must carry an interval ("length-positivity, the
epsilon schedule"). The worry was that tight base rows certify only isolated points.
There are none. Concretely, for every band row with optimum M ≥ 1/12 at witness t₀:

    ∀ t, |t − t₀| ≤ δ  ⟹  min_i ||s_i t|| ≥ M − W·δ ≥ 1/14   for δ = (M − 1/14)/W ≥ 1/1848,

i.e. every covering band family has an INTERVAL of 1/14-lonely times of length ≥ 1/924.
The repeat leg (LRC(≤13) citation) has its own margin 1/13 − 1/14 = 1/182, giving intervals
of length ≥ 2·(1/182)/22 = 1/2002 for bounded families with a repeated entry.

**So the full bounded window at W = 22 admits an interval-loneliness upgrade with the
uniform floor ε = 1/2002** — the exact base the length-invariant peel wrapper consumes.

## The Lean upgrade path (next session-scale piece)

1. `speedOKr` — the kernel gate at parametric radius: `(den) ≤ ρ·min((s·num) % den, den − (s·num) % den)`
   with ρ = 12 for the band rows (re-emit witnesses from the margin analysis: the
   first-passing witnesses in `LRCWindowData22.lean` were found at ρ = 14; the OPTIMAL
   witnesses achieving ≥ 1/12 are computed and stored by the margin script).
2. The interval transport lemma (pure abs algebra over ℝ):
   `Lonely-at-radius-1/12 at t₀ → ∀ t ∈ [t₀ − δ, t₀ + δ], Lonely 14 v t` for `δ ≤ (1/12 − 1/14)/max|v|`.
3. The interval form of `lonely14_of_repeat` (same margin argument through the citation).
4. `hwindow22_interval (cite) : ∀ v ≠ 0, |v| ≤ 22 → ∃ a b, b − a ≥ 1/2002 ∧ ∀ t ∈ [a,b], Lonely 14 v t`.
5. The length-invariant census+peel wrapper (kps's design note) consuming exactly this.

## What this does NOT solve

The analytic peel itself (the middle band 22 < N ≲ N*): the interval base makes the
far-runner argument *statable*, but the union bound still fails for f ≥ 7 far runners at
radius 1/14 (each far runner can be unsafe on 1/7 of any interval; 13·(1/7) > 1), and one
far bad-arc of width 1/(7N) can swallow a short interval when N is just above the cut.
The rate-lemma/F3-sharp machinery (joint rate_core) remains the intended engine there.
