---
id: THM-930  # renumbered from 928 (opus two-scale first-pushed 928)
title: THE LEVERAGE IDENTITY AND THE FIRST CERTIFIED 13-PACKET — (I) THE DEPTH-SPECTRUM REFRAME: for any speed packet, the avoidance polynomial is the coverage-spectrum generating function p(t) = Σ_d μ_d (1−t)^d, where μ_d = μ{t : exactly d bad sets active} (one sweep computes everything; S_k = Σ_d C(d,k)μ_d); (II) THE LEVERAGE IDENTITY (PROVED, two lines): BONF_m = μ₀ + (−1)^m Σ_{d>m} C(d−1,m)·μ_d — the Bonferroni truncation error in CLOSED FORM (from Σ_{k≤m}(−1)^k C(d,k) = (−1)^m C(d−1,m)); refereed EXACTLY on 5 packets × all m = 1..12. COROLLARIES: (a) the level-5 certificate is literally "good mass beats the leveraged deep-coverage tail": BONF5 > 0 ⟺ Σ_{d≥6} C(d−1,5)μ_d < μ₀; (b) an atom at full depth 13 carries leverage C(12,5) = 792, so mass > B5_eq/792 = 2052/(16807·792) ≈ 1.54×10⁻⁴ (ERRATUM fixed by THM-935: B5_eq = 2052/16807; the 0.0821 first quoted was the certified packet´s own value) alone kills the certificate — the deep well's origin window (1/91, THM-929) is 71× over; (c) THE TIGHT SYSTEM'S LEVEL-5 EXPLOSION IS ONE ATOM: C(13,5)·μ₁₃ = 99/7 = 14.143 of S₅ = 16.159 — 87.5% exactly; the three blocker species (THM-926: quadruples, dilates, linear forms) matter exactly inasmuch as they create deep-coverage mass; (III) THE FIRST CERTIFIED PACKET — resolving THM-897's named open POSITIVELY: {307, 425, 541, 671, 800, 944, 1087, 1413, 1943, 2147, 2570, 3056, 3310} (greedy: Sidon + no 3-AP + ratio horizon 30 + no 7/13-multiples; audit passes all four) has EXACT BONF5 = 11882280671749727299498652729985/144812363687055259298722873214576 ≈ +0.08205 > 0: the level-5 wall certifies its loneliness at λ = 1/14 in exact rational arithmetic — the corrected-admissibility regime is NON-EMPTY at accessible scale [307, 3310], and its depth spectrum is near-perfectly binomial (μ₀ = 0.135431 vs (6/7)¹³ = 0.135178); (IV) HYP-7140's LITERAL FORM REFUTED (honest negative): the certifying packet is NOT numerically near-real-rooted — near-independence forces a near-13-fold root at s = −6 which scatters into a complex ring of radius ~ε^{1/13} under ε-perturbation, so root geometry is the WRONG instrument; the true criterion is the moment-side identity (II)
status: (II) PROVED + refereed exactly (all packets, all m); (III) exact rational certificate + property audit; (I) identity of formalisms (one-line change of variables); (IV) refutation with mechanism. HYP-7140 closed-as-corrected into (II)
source: kind-pasteur-2026-07-16-S128 (cont.35; owner: work the Lee-Yang admissibility criterion)
depends_on:
  - THM-929 (the quintic wall + exact tight-system S_k this decomposes)
  - THM-897/926 (opus: the wall, the admissibility ladder, the named open this resolves)
related:
  - THM-826 (Farey profile — the μ₀ side), THM-735 (my Bonferroni multi-peel — the same alternation used stratum-wise)
---

# THM-930 — the leverage identity (renumbered from 928) and the first certified packet

## (I)–(II) The identity

Depth spectrum: μ_d = μ{t : #{v : ‖vt‖ < λ} = d}. Then p(t) = ∫(1−t)^{depth} = Σ_d μ_d(1−t)^d
and S_k = Σ_d C(d,k)μ_d. Using Σ_{k≤m}(−1)^k C(d,k) = (−1)^m C(d−1,m):

> **BONF_m = Σ_{k≤m}(−1)^k S_k = μ₀ + (−1)^m Σ_{d>m} C(d−1,m) μ_d.**

Even m: upper bound (tail added); odd m: lower bound (tail subtracted) — the classical
Bonferroni inequalities with the error term EXACT. The level-5 certificate condition is
Σ_{d≥6} C(d−1,5)μ_d < μ₀: deep coverage, binomially leveraged, versus good mass.

## (III) The certified packet

{307, 425, 541, 671, 800, 944, 1087, 1413, 1943, 2147, 2570, 3056, 3310}:
- BONF5 = 11882280671749727299498652729985/144812363687055259298722873214576 > 0 (exact);
- μ₀ = 117672094910664692133055713400657/868874182122331555792337239287456 ≈ 0.135431;
- weighted tail ≈ 0.053378 < μ₀ ✓;
- audit: pairwise sums all distinct (Sidon), no 3-term AP, no ratio p/q with p,q ≤ 30,
  no multiples of 7 or 13.
The greedy construction found it in milliseconds where the stricter linear-forms DFS
(THM-926: |c| ≤ 3, horizon 20) crawled — general small linear forms beyond 3-APs
apparently do not bite at this scale. LRC(14) holds for this packet by a level-5
certificate — the wall's first non-vacuous use.

## (IV) Why root geometry fails (the honest negative)

Independence makes p̃(s) = Σμ_d s^d = ((s+6)/7)¹³ — a 13-fold root at −6. Any
ε-perturbation scatters it into a ring of radius ~ε^{1/13} (ε = 0.005 → 0.66): even
near-perfect packets look "complex-rooted". The Lee-Yang intuition survives only as the
moment identity (II).

## Priority note (added same session)

Opus-S332's THM-928 found a level-5 certificate CONCURRENTLY (packet {300..2208},
BONF5 ≥ +0.039131, weighted filter) and first-pushed by DAG order. This file's packet is
an independent construction (greedy, different filter), with 2.1× the margin (+0.08205),
and the leverage identity (II) is this file's own contribution. Two independent
certificates within one hour — the regime is not thin once the filter is right; the
BLOCKER for the earlier DFS was its own strictness (MISTAKE-151 family), not the math.

## Named next
- Leverage-guided construction: minimize the weighted tail directly (the greedy above
  minimizes proxies); gives the natural "most-certifiable packet" extremal question.
- The certified-packet DENSITY: how thin is the BONF5 > 0 regime in [300, 4000]^13?
- Stratum version: combine with THM-735's per-stratum Bonferroni for covering families.

## Evidence log
- [x] identity referee 5 packets × all m (bonf_leverage_identity_kps_S128c35.out)
- [x] certificate + audit exact (same file); spectra (depth_spectrum_leeyang_kps_S128c35.out)
