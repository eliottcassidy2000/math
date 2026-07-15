---
id: THM-786
title: The extent-form exit study — wall-count K0 is a ratio-boxed artifact; the no-companion extent bound is proved and the universal bound is conjectural; the original serving/de-phase completion is withdrawn
status: PROVED (the wall-count artifact refutation with exact certificates; the no-companion extent theorem and its conditional pierce consequence) + REPORTED/NOT REPRODUCED BY THE STORED SCRIPT (finite adversarial extent census, stated ratio ≤ 0.59) + CORRECTED (the serving/de-phase and sparse-regime claims are REFUTED/WITHDRAWN by MISTAKE-148) + CONJECTURE (the universal extent bound) + OPEN (active-period/co-landing bound and core incidence)
source: opus-2026-07-14-S304 (owner directive: prove the geometric co-landing bound, finish the exit lemma)
depends_on:
  - THM-783   # corrected period-sum/single-visitor and conditional extent laws
  - THM-779   # the criterion; its census constant corrected here
related: [THM-767, THM-771, THM-784, THM-788, HYP-6840, HYP-6845, HYP-6850, MISTAKE-147, MISTAKE-148]
verification: 04-computation/lrc14_extent_exit_theorem_opus_S304.py
  (+ 05-knowledge/results/lrc14_extent_exit_theorem_opus_S304.out)
---

# THM-786 — the extent-form exit theorem

> **Correction (codex-S10 referee audit).** The no-companion extent theorem in
> §2 is sound. The serving/de-phase bound and sparse completion originally
> claimed in §3 are false as stated and are withdrawn; see MISTAKE-148. The
> census in §4 is finite evidence only. Consequently §5 does not finish the
> general r=8 pierce. THM-788 gives a sound replacement reduction through the
> number of active fastest periods.

## (1) The wall-count constant was an artifact (REFUTED, exact certificate)

THM-783's census constant "K0 = 6 walls" sampled comparable-speed tuples only.
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
> (Σ w^{-1} ≡ 0 mod 7) and of size ≥ 2 — the single-visitor break (THM-783(3)).
> **(b)** Hence if no interior g-wall is served by a balanced co-landing
> companion, the interior contains no g-wall at all, and
> **extent < 1/w_g + 2/w_f.**
> **(c)** In general, extent < (M_g + 1)/w_g + 2/w_f, where M_g is the maximal
> number of CONSECUTIVE interior g-walls with balanced-visited periods.

## (3) CORRECTION — the advertised geometric co-landing completion is withdrawn

- **Serving/de-phase is false as stated.** For `(f,g,c)=(9,8,6)`, four
  consecutive g-walls are served by c, while the displayed formula gives only
  three. Paired indices can cross, so their signed offset ranges across a
  window of width `2/f`, not `1/f`. General co-visits need not even advance
  both indices by one.
- **What survives algebraically.** Subtracting two zero-sum visitor identities
  shows that the multiset change has inverse-residue sum zero. A literal
  one-out/one-in change therefore pairs opposite residues. This does not prove
  simultaneous boundary crossing, fixed order, de-phasing, or a mandatory
  higher-cluster handover.
- **Sparse-regime claim withdrawn.** Its serving-span premise used the false
  bound and its displayed inequality was not derived. A valid sparse or dense
  result must instead control active f-periods or prove additional pairing/order
  hypotheses.

## (4) The adversarial extent census (REPORTED; stored script does not reproduce it)

Maximal run extent as a fraction of 1/w_g + 2/w_f, quarter-period windows:

| family | n | median | max |
|---|---|---|---|
| generic (w ≤ 3000) | 60 | 0.090 | 0.339 |
| extreme-ratio (the wall-count breaker) | 60 | 0.000 | 0.448 |
| balanced pairs (2w_g + Δ ≡ 0 mod 7, designed co-landers) | 60 | 0.036 | 0.557 |
| near-multiples w_f = N·w_g + ε (count-lock exploit) | 60 | 0.125 | 0.571 |
| annealed peak (300 steps) | — | — | **0.589** |

The session reported that the tested families built to exploit the co-landing
loophole never reach 60% of the proposed bound. The stored script prints this
summary but implements only the two exact wall-count certificates; it does not
reproduce the table. Treat the table as unreproduced evidence for the
conjecture, not evidence for the withdrawn swap/serving proof mechanism.

> **The universal extent conjecture (sharp):** every blocking run at the prime-7
> lens with r = 8 has extent < 1/w_g + 2/w_f. (Proved on class (2b); ≤ 0.59 of
> the bound everywhere tested, including designed exploits.)

## (5) The r = 8 pierce in the proved no-companion class

> **Every closed core-safe component of length ≥ 1/w_g + 2/w_f contains a wall
> where blocking fails — a full 1/14-witness moment — PROVED whenever the run
> covering it would fall in class (2b), and supported by the stated finite
> census otherwise.** Components shorter than the bound are finite
> per-family checks (the THM-779 integer walk, O(#walls)).

This replaces THM-779(4)'s "components with more than K0 walls" — wall-count
comparisons are withdrawn; the conditional extent comparison stands. Without
the no-companion hypothesis, neither the census nor the fact that the proposed
bound shrinks proves that every core-safe component is pierced.

## (6) What remains (honest, sharp)

Bound the number `A` of active fastest periods after empty fastest-owner blocks
are contracted. THM-788 proves that `A<=C` would give
`extent < (C+1)/w_g+(C+2)/w_f` and the corresponding ratio-sensitive wall
bound. The finite census margin (0.59) suggests the universal conjecture in
(4), but no valid serving/cascade proof is currently known. The remaining
problem is both Diophantine-combinatorial (varying-index co-landings versus
mod-7 balance) and geometric (incidence with the core-safe components).
