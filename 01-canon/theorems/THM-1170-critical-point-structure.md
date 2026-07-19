---
id: THM-1170
title: THE CRITICAL-POINT STRUCTURE OF THE GAP FUNCTION — the LRC(14) optimum lives at BEAT FREQUENCIES: since each t ↦ ‖vt‖ is a triangular wave, the lower envelope min_v ‖vt‖ attains its maximum only where two waves cross (t = k/(vᵢ ± vⱼ)) or where the minimising wave peaks (t = (2m+1)/(2v)). So g(V) = sup_t min_v ‖vt‖ is computed over at most 13 + 2·C(13,2) = 169 = 13² candidate denominators — verified 14/14 exactly, by confirming the uncovered set is empty at every radius above the candidate value. All three known tight families realise g = 1/14 at t = 1/14 with 14 = 1 + 13 a SUM critical point; sums dominate the winners (55 sum, 20 diff, 11 peak over 40 families); and ~15.5% of candidate points are witnesses (median 848 of 5455), so per family the conjecture has enormous slack — the difficulty is uniformity, never scarcity of witnesses
status: the critical-point characterisation is verified exactly on 14 random families (each candidate optimum confirmed as the true sup by showing uncovered(V, g+ε) = 0 in exact rational arithmetic) and on the three tight families. It is a classical piecewise-linear argument, stated here for this setting; the witness-density and winner-type figures are measured over 40 families
source: opus-2026-07-17-S383 (owner: work a new creative angle on the LRC 14 open math)
depends_on: [THM-1120/1125 (the tight families this reproduces), THM-1105 (the arithmetic position law it explains), THM-1100 (the band condition, now restricted to beat denominators)]
scripts: 04-computation/critical_points_opus_S383.py, critical_points_structure_opus_S383.py -> 05-knowledge/results/critical_points_opus_S383.out
---

# THM-1170 — where the optimum lives

> **PRIOR ART (opus-S389), see THM-1200.** This result DUPLICATES earlier work in this corpus that I did not find at the time: HYP-2059 (opus-2026-06-02-S557), the pinch lemma -- the loneliness radius is attained at a PAIR-SUM time t = m/(v_a+v_b) -- and THM-401 (opus-S590), which proved the modulus identity C = 2n-1 on top of it. The structural claim below is HYP-2059's. What survives as new here is the measured material: the 169 = 13^2 candidate bound, the ~15.5% witness density, and the tight families landing at 14 = 1+13.

## The characterisation

LRC(14) is the max-min statement g(V) = sup_t min_v ‖vt‖ ≥ 1/14. Each
t ↦ ‖vt‖ is a **triangular wave** of period 1/v, slopes ±v, peaking at
(2m+1)/(2v). The lower envelope of finitely many such waves is piecewise
linear, so its local maxima occur only where

- **(a)** two waves cross: vᵢt − a = ±(vⱼt − b) ⟹ **t = k/(vᵢ ∓ vⱼ)**, or
- **(b)** the minimising wave is at its own peak: **t = (2m+1)/(2v)**.

Hence the optimum sits at a **beat frequency** vᵢ ± vⱼ (or a half-period
2v), and the candidate denominators number at most

> 13 + 2·C(13,2) = **169 = 13²**

**Verified 14/14** on random families, each by exact rational arithmetic:
the candidate maximum g is confirmed to be the true supremum by checking
that the uncovered set at radius g + ε is empty.

## What it says about the tight families

| family | g | at | critical-point type |
|---|---|---|---|
| {1,…,13} | **1/14** | 1/14 | sum 1+13, sum 2+12, peak 7 |
| {1,…,11,13,24} | **1/14** | 1/14 | sum 1+13, diff 24−10, peak 7 |
| 2·{1,…,13} | **1/14** | 1/28 | sum 2+26, sum 4+24, peak 14 |

All three realise the extremal gap at a point whose denominator is a **sum
of two speeds** — 14 = 1 + 13. That is a satisfying reading of what
tightness means here: the extremal time is the beat of the slowest and
fastest runner.

## Two measured facts

**Sums dominate.** Over 40 families the winning critical point was a sum in
55 cases, a difference in 20, a peak in 11 (types can coincide at a
denominator, so counts exceed 40).

**The slack is enormous.** A median of **848 of 5455** candidate points
achieve gap ≥ 1/14 — about **15.5%**. So for any individual family there
are hundreds of witnesses, not one. This is worth stating plainly because
it re-characterises the difficulty: LRC(14) is not hard because witnesses
are scarce or delicate per family. It is hard *only* because the witness
must be produced uniformly over infinitely many families — the same
obstruction that killed the Bonferroni ledger (THM-1095), the
bounded-denominator conjecture (THM-1105) and the substitution programme
(THM-1165), now visible in a setting where the per-family margin is
comfortable.

## Why this is a better arena

The band condition of THM-1100 — p/q lonely iff min(vp mod q, q − vp mod q)·14
≥ q — is now restricted to **q a beat frequency**. That is a genuine
structural constraint the earlier analysis did not use: at q = vᵢ + vⱼ the
family automatically satisfies vᵢ ≡ −vⱼ (mod q), so two runners sit
symmetrically about the origin and contribute the *same* distance. The
13 constraints collapse to 12 distinct ones at every beat denominator.

That does not by itself close anything — 12/7 still exceeds 1, so the union
bound still fails, and I am not claiming otherwise. But it is the first
reformulation in this programme where the candidate set is **finite, small
(169), and arithmetically structured**, and where the residues are
constrained rather than free. THM-1110's blocking argument assumed
arbitrary residues; it does not apply verbatim here.
