---
id: THM-689
title: The dead-zone lemma — (A) k ≤ 7 RIGIDITY (PROVED): seven arcs of length 1/7 cover the circle only as a PERFECT NET (pairwise-disjoint open arcs with zero slack ⟹ centers equally spaced ⟹ finitely many α), so meas(D) = 0 and μ∞(P,E) > 0 unconditionally for k ≤ 7; (B) the MEASURE CRITERION + EXHAUSTIVE P-CENSUS: μ∞ = 0 needs meas(D(E)) ≥ m_P, and min m_P over ALL P ⊆ [1,13] (exhaustive per size) beats the consecutive dead zone at every k = 8..12 with 2–6× margins — the whole two-scale dead-zone question reduces to CONSECUTIVE-MAXIMALITY of meas(D) (one extremal lemma; hill-climb + battery + dilation-invariance evidence); (C) the FIRST-WINDOW THEOREM: W_P = [1/(14·pmin), 13/(14·pmax)] ⊆ G_P always, degenerating EXACTLY at the ratio-13 tight locus; (D) the hunt: zero classes with μ∞ = 0 (min found 0.0447, at the consecutive apex shape again)
status: PROVED ((A) complete proof below; (C) two-line inclusion; (B)'s criterion proved, the P-side census EXHAUSTIVE per size, the E-side maximality of the consecutive shape CONJECTURED with strong evidence — battery + hill-climbs at k = 8, 11, 12 all return consecutive; dilation invariance meas(D(cE)) = meas(D(E)) proved, explaining consec = AP_d exactly). With (A)+(B): the two-scale dead-zone lemma holds unconditionally at k ≤ 7 and modulo consecutive-maximality at k = 8..12 (|P| = 0, k = 13 is non-residual — ratio ≤ 13). Per-class positivity remains exactly decidable regardless (the evaluator), so THE WALL's per-class pipeline is unaffected by the remaining lemma.
source: klein-2026-07-10-S238 (HYP-5915; the named structural gap of THM-687/688; hypohamiltonian-criticality frame per owner directive — see the companion reflection)
depends_on:
  - THM-687 / THM-688   # the limit measures this closes the positivity of
related:
  - kps LRCStrictRuler (the wall chain), THM-667 Lemma A (m_P floors)
  - 07-reflections/hypohamiltonian-criticality-and-the-lonely-runner.md (the frame)
---

# THM-689 — the dead-zone lemma: rigidity, criterion, and consecutive-maximality

## Setting

D(E) := {α ∈ [0,1) : the arcs (eα − 1/14, eα + 1/14), e ∈ E, cover T}
= {α : every gap of the center multiset {eα mod 1} is ≤ 1/7}. μ∞(P,E) = 0
forces G_P ⊆ D(E) up to measure zero (m_E > 0 off D). 0 ∈ E throughout
(the cluster contains V_max).

## (A) The k ≤ 7 rigidity theorem (PROVED)

For |E| = k ≤ 6: k arcs of total mass k/7 < 1 never cover — m_E ≥ 1 − k/7
pointwise (THM-687(C)). For k = 7: if m_E(α) = 0 then meas(∪ arcs) = 1 =
Σ meas(arc), so all pairwise intersections are null; open intervals with
null intersection are disjoint; seven disjoint open arcs of total measure 1
leave exactly seven boundary points, so consecutive centers are exactly 1/7
apart: the centers form a perfect translate of (1/7)ℤ/ℤ. Then for any two
distinct e, e′ ∈ E, (e − e′)α ∈ (1/7)ℤ (mod 1), confining α to the finite
set (1/(7|e−e′|))ℤ. Hence {α : m_E(α) = 0} is finite, meas(D) = 0, and

> **μ∞(P,E) = ∫_{G_P} m_E > 0 whenever m_P > 0, for every E with |E| ≤ 7.**

(Confirmed numerically: meas(D) = 0 exactly on consecutive/AP₂/AP₃/random
k = 7 shapes.) ∎

## (B) The measure criterion and the census (k = 8..12)

**Criterion (proved):** μ∞(P,E) = 0 ⟹ meas(D(E)) ≥ m_P.

**P-side (EXHAUSTIVE):** min m_P over ALL P ⊆ [1,13] of size 13 − k:
|P| = 5: 0.381463 at {1,5,7,8,9}; |P| = 4: 0.494256 at {1,11,12,13};
|P| = 3: 0.604396 at {1,12,13}; |P| = 2: 0.725275 at {1,13};
|P| = 1: 6/7 at {1}.

**E-side (battery + hill-climb):** meas(D) for consecutive E = {0..k−1}:
k = 8: 0.059864, k = 9: 0.159864, k = 10: 0.224490, k = 11: 0.373696,
k = 12: 0.430097 (k = 13: 0.557514, non-residual). Dilation invariance
meas(D(cE)) = meas(D(E)) is proved (α ↦ cα is measure-preserving and maps
the gap structure of {ceα} to that of {eα}), so all AP shapes tie
consecutive exactly (observed). Random/detuned/designed shapes all fall
below; hill-climbs at k = 8, 11, 12 return consecutive.

**Consequence:** at every k = 8..12, min m_P EXCEEDS the consecutive dead
zone with margin 6.4× / 3.1× / 2.7× / 1.9× / 2.0× — so

> **the entire two-scale dead-zone lemma reduces to CONSECUTIVE-MAXIMALITY:
> meas(D(E)) ≤ meas(D({0,…,k−1})) for all |E| = k** (conjectured; strong
> evidence). Given it, μ∞(P,E) > 0 for every residual two-scale class.

## (C) The first-window theorem (PROVED)

For α ∈ W_P := [1/(14·pmin), 13/(14·pmax)]: every p ∈ P has pα ∈
[p/(14 pmin), 13p/(14 pmax)] ⊆ [1/14, 13/14] with no wrap, so **W_P ⊆ G_P
always**; W_P ≠ ∅ ⟺ pmax ≤ 13·pmin, degenerating to the single point 1/14
exactly at ratio 13 — the tight locus is the vanishing of the first safe
window (kps's knife-edge, seen from the window side). Containment tests
against W_P are the cheapest per-class disproof of μ∞ = 0: the designed
adversaries put ZERO dead mass in W_P for |P| ≥ 4.

## (D) The hunt

Over the hard grid (min-m_P and ratio-13 P-shapes × consecutive/AP/detuned/
random/window-designed E-shapes, all k): **zero classes with μ∞ = 0**;
minimum found μ∞ = 0.044710 at P = {1,2,3,4}, E = {0,…,8} — the consecutive
apex shape, again the extremal, again positive.

## Status of the wall

k ≤ 7 two-scale classes: closed UNCONDITIONALLY (A). k = 8..12: closed
modulo consecutive-maximality (B) — one clean extremal statement about
circle coverings by dilate-centered arcs. Everything remains per-class
decidable exactly (the evaluator), so the wall pipeline (THM-685 transfer +
kps strict chain) is unconditional per class today and wholesale modulo (B).

## Verification & files

`04-computation/lrc14_dead_zone_lemma_klein_S238.py` (+ `.out`): rigidity
checks, the census, the exhaustive min-m_P table, the hunt.
`05-knowledge/results/lrc14_dead_zone_hillclimb_klein_S238.out`: the E-shape
hill-climbs (consecutive stands at k = 8, 11, 12).

## Addendum — (A) FORMALIZED = the ≤7-arcs pigeonhole DE-CITED (klein-2026-07-11-S251, HYP-5985)

`LRCSevenGapRigidity.lean` (kernel-pure) proves the grand assembly's citation
`SmallClusterFull` outright — for every nodup cluster `E` with `3 ≤ |E| ≤ 7`,
`slowμ (goodSet E) = 1` — which is exactly (A)'s rigidity made pointwise:

- **`cgap_sum`**: the k cyclic gaps of a sorted phase enumeration telescope to 1
  (the wrap capped at `1 + f 0`; `Finset.sum_range_sub` does the telescope).
- **`cgap_witness` (the bridge)**: ANY gap `> 1/7` makes its left endpoint a
  goodSet witness. The key simplification found here: an interior gap at `s_i`
  forces `s_i < 6/7`, so every wrap term `1 + s_j − s_i` clears `1/7` for free —
  NO maximum-gap selection is needed, any big gap works; the wrap gap is
  witnessed at the top phase.
- **`gap_dichotomy`**: `k ≤ 6` ⟹ some gap `> 1/7` (mac-mini's pigeonhole,
  `LRCMaxGapPigeonhole`); `k = 7` with no big gap ⟹ ALL gaps exactly `1/7`
  (sum-erase bound) ⟹ the first sorted pair sits at `fract = 1/7` EXACTLY —
  the PERFECT NET.
- **`pair_fiber_countable`**: `fract(d·x) = 1/7`, `d ≠ 0` pins `x = (n + 1/7)/d`
  — a countable fiber; the ≤ 49 tooth pairs give a countable, hence null, bad set.
- **`smallClusterFull_proved : SmallClusterFull`** assembles: off the null set
  every time is good; `slowμ(goodSet) = 1`.
- **`lrc14_from_moment_and_supply`**: LRC(14) = [THM-661 moment floor] +
  [certificate supply]. **The grand assembly is down from THREE citations to TWO.**

Phases with collisions dedup through `Finset.image` (coincident phases only help:
fewer distinct points, bigger max gap), which is where nodup-of-values vs
distinct-teeth is handled — the theorem needs no arithmetic on the teeth at all.
