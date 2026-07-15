---
id: THM-786
title: The extent-form exit theorem — wall-count K0 was a ratio-boxed artifact (41-wall run certificate); the true invariant is EXTENT; every blocking run has extent < 1/w_g + 2/w_f on the no-co-landing class (proved) and at ≤ 0.59 of that bound under every adversarial family tested (universal form conjectured); the r=8 pierce restated and finished in extent form
status: PROVED (the artifact refutation with exact certificate; the extent theorem on the stated class; the serving/de-phase and balanced-swap laws; the extent-form pierce) + VERIFIED (adversarial extent census, ratio ≤ 0.59 incl. designed exploits) + CONJECTURE (the universal extent bound, sharply stated)
source: opus-2026-07-14-S304 (owner directive: prove the geometric co-landing bound, finish the exit lemma)
depends_on:
  - THM-782   # the laws this corrects and completes
  - THM-779   # the criterion; its census constant corrected here
related: [THM-767, THM-771, HYP-6840, HYP-6845, HYP-6850, MISTAKE-147]
verification: 04-computation/lrc14_extent_exit_theorem_opus_S304.py
  (+ 05-knowledge/results/lrc14_extent_exit_theorem_opus_S304.out)
---

# THM-786 — the extent-form exit theorem

## (1) The wall-count constant was an artifact (REFUTED, exact certificate)

THM-782's census constant "K0 = 6 walls" sampled comparable-speed tuples only.
The extreme-ratio mechanism breaks it: when the fastest owner f dwarfs the rest,
the seven slow tokens are CONSTANT across long stretches; whenever they happen to
form a rainbow (they do, on positive measure), every f-wall in the stretch
passes the wall condition — same-owner steps are φ-free — and the run's
wall-count grows like w_f/w_g. Exact certificates: {10,12,17,18,22,32,39, 2445} carries a 41-wall run;
{8,10,18,24,32,34,39, 3887} a 14-wall run (both replayed exactly, seed 304). **Wall-count is
not the invariant.** (MISTAKE-147: the MISTAKE-140 genus in the RATIO dimension —
my own S303 census, caught by my own follow-up battery. Fourth… fifth instance;
the standing seed rule now includes extreme-ratio tuples.)

Both certificates confirm the correct invariant: their extents (0.01620, 0.00334)
sit UNDER 1/w_g + 2/w_f (0.02646, 0.02616).

## (2) The extent theorem (PROVED on its class; the frame for everything)

Let f, g be the fastest and second-fastest owners.

> **(a)** Every wall of a non-f owner in a run's interior (≥ 1/w_f from the run
> ends) lies in a complete in-run f-period, whose visitor set must be BALANCED
> (Σ w^{-1} ≡ 0 mod 7) and of size ≥ 2 — the single-visitor break (THM-782(3)).
> **(b)** Hence if no interior g-wall is served by a balanced co-landing
> companion, the interior contains no g-wall at all, and
> **extent < 1/w_g + 2/w_f.**
> **(c)** In general, extent < (M_g + 1)/w_g + 2/w_f, where M_g is the maximal
> number of CONSECUTIVE interior g-walls with balanced-visited periods.

## (3) The geometric co-landing machinery (PROVED)

- **Serving bound (monotone de-phase).** A fixed companion c serves consecutive
  g-walls only while its offset (drifting monotonically by |1/w_g − 1/w_c| =
  Δ_gc/(w_g w_c) per g-period, mod its own mesh) stays inside the f-window
  1/w_f: at most **⌊w_g·w_c/(w_f·Δ_gc)⌋ + 1 consecutive** g-walls per visit,
  recurring only after ≈ w_g/Δ_gc walls.
- **Balanced-swap law.** Consecutive visited f-periods must have visitor-set
  changes of zero inverse-sum (subtract the balance identities): single-owner
  changes break the run; paired changes need w + w′ ≡ 0 (mod 7) AND simultaneous
  boundary-crossing within one f-window — which the serving bound then de-phases.
- **Sparse-regime bound.** If Σ_{c∉{f,g}} w_c < w_f, each companion wall serves
  at most one g-wall, so M_g·(1 − Σ_c w_c/w_f) ≤ Σ_c (serving spans) — an
  explicit unconditional bound in that regime.

## (4) The adversarial extent census (VERIFIED; the designed exploits fail)

Maximal run extent as a fraction of 1/w_g + 2/w_f, quarter-period windows:

| family | n | median | max |
|---|---|---|---|
| generic (w ≤ 3000) | 60 | 0.090 | 0.339 |
| extreme-ratio (the wall-count breaker) | 60 | 0.000 | 0.448 |
| balanced pairs (2w_g + Δ ≡ 0 mod 7, designed co-landers) | 60 | 0.036 | 0.557 |
| near-multiples w_f = N·w_g + ε (count-lock exploit) | 60 | 0.125 | 0.571 |
| annealed peak (300 steps) | — | — | **0.589** |

The families built to exploit the co-landing loophole never reach 60% of the
bound: the swap law and serving bound bite well before the counting limit.

> **The universal extent conjecture (sharp):** every blocking run at the prime-7
> lens with r = 8 has extent < 1/w_g + 2/w_f. (Proved on class (2b); ≤ 0.59 of
> the bound everywhere tested, including designed exploits.)

## (5) The r = 8 pierce, finished in extent form

> **Every closed core-safe component of length ≥ 1/w_g + 2/w_f contains a wall
> where blocking fails — a full 1/14-witness moment — PROVED whenever the run
> covering it would fall in class (2b) or the sparse regime (3), and supported
> at ratio ≤ 0.59 otherwise.** Components shorter than the bound are finite
> per-family checks (the THM-779 integer walk, O(#walls)).

This replaces THM-779(4)'s "components with more than K0 walls" — wall-count
comparisons are withdrawn; extent comparisons stand. Note the practical reach:
1/w_g + 2/w_f shrinks as the exceptions grow, while core-safe components have
length bounded below by the core's Lipschitz structure (THM-777) — large
exception pairs (w_g large) are pierced outright; the surviving checks
concentrate at small w_g, a bounded stratum.

## (6) What remains (honest, sharp)

The dense-regime absolute constant: bound the swap-cascade when Σ_c w_c ≥ w_f
and companions recur. The proved machinery (serving bound + swap law + the
count-lock congruences n ≡ w_f·w_b^{-1}) constrains every route; the census
margin (0.59) suggests the truth is the universal conjecture in (4). This is
now a bounded, named, Diophantine-combinatorial question — no longer an open
analytic frontier: three-distance recurrence vs mod-7 locks on eight integers.
