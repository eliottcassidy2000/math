---
id: THM-1205
title: LRC(14) AS A PURE ACTIVE-PAIR RATIO — building on boxeph-S120's located-maximizer theorem, at the maximizer two runners straddle (vᵢt* = aᵢ + g, vⱼt* = aⱼ − g), so g·(vᵢ+vⱼ) = vᵢaⱼ − vⱼaᵢ = D, a positive integer. Hence g(V) = D/(vᵢ+vⱼ) and **LRC(14) ⟺ vᵢ + vⱼ ≤ 14·D at the active pair** — a purely arithmetic condition on ONE pair, with no measure anywhere, so it sits on the surviving side of the THM-1185 triage. Verified 18/18. The ratio (vᵢ+vⱼ)/D runs 4.75–6.50 (median 5.73) on random families, comfortably inside the threshold 14; and BOTH tight families give ratio **exactly 14**, at active pair (1,13) with D = 1 — the boundary is real and is attained precisely there
status: the formula is verified exactly 18/18 (rational arithmetic, maximizer located via the HYP-2059 critical-point set); the ratio statistics are measured over 12 random families; the tight-family values are exact. The reformulation is an EQUIVALENCE, not a proof — it restates LRC(14), it does not establish it
source: opus-2026-07-19-S390 (owner: pull from other agents and work a new improved angle)
depends_on: [boxeph-S120 / HYP-7745 (the located-maximizer theorem this is built on), HYP-2059 (pinch lemma, the critical-point set), THM-1185 (the measure-vs-pointwise triage that motivated seeking an arithmetic form), THM-1120 (the tight families)]
scripts: 04-computation/active_pair_ratio_opus_S390.py -> 05-knowledge/results/active_pair_ratio_opus_S390.out
---

# THM-1205 — the conjecture as one pair's arithmetic

## Derivation

boxeph-S120 located the maximizer at a straddling active pair. Writing that
out: at t* two runners satisfy

> vᵢt* = aᵢ + g  and  vⱼt* = aⱼ − g

Eliminating t* gives g(vᵢ + vⱼ) = vᵢaⱼ − vⱼaᵢ =: **D**, a positive integer.
So g(V) = D/(vᵢ+vⱼ), and the conjecture becomes

> **LRC(14) for V  ⟺  vᵢ + vⱼ ≤ 14·D  at the active pair.**

Verified exactly **18/18**.

## What the ratio looks like

| family | active pair | D | vᵢ+vⱼ | ratio |
|---|---|---|---|---|
| random | (34,14) | 8 | 48 | 6.00 |
| random | (25,13) | 8 | 38 | 4.75 |
| random | (12,27) | 6 | 39 | 6.50 |
| … | | | | median **5.73** |
| **{1,…,13}** | **(1,13)** | **1** | **14** | **14 exactly** |
| **{1,…,11,13,24}** | **(1,13)** | **1** | **14** | **14 exactly** |

Random families sit at ratio ~5–6.5, a factor of about 2.4 inside the
threshold. Both tight families sit **exactly on it**, with the same active
pair (1,13) and D = 1.

## Why this is the right kind of statement

It is **purely arithmetic and pointwise** — no measures, no integrals, no
densities. That matters because THM-1185 established that every
measure-based method (Bonferroni, the density bounds, the Delsarte LP) is
structurally blind to the tight families, since those have nothing to see at
the level of measure. This formulation sees them perfectly: they are exactly
the equality case.

It also explains the tight families in one line. They are tight because
their active pair is (1,13), whose sum is 14 and whose determinant is 1 —
the unique way to sit exactly at the boundary with D = 1.

## The concrete target it opens

A counterexample needs an active pair with vᵢ + vⱼ > 14·D. With D = 1 that
means a coprime active pair summing to **more than 14** — e.g. (1,14) would
give g = 1/15 < 1/14. So:

> **Is there a 13-family whose active pair has D = 1 and vᵢ + vⱼ ≥ 15?**

That is a sharp, finite-flavoured question about one pair, and it is the
whole of LRC(14) at D = 1. Higher D requires sum > 14D, which the measured
ratios (≤ 6.5) suggest is far off, but that is measurement, not proof.

## Status

An equivalence, not a proof. It restates LRC(14) in arithmetic terms and
locates the extremal families exactly, which the measure-theoretic
formulations could not do.
