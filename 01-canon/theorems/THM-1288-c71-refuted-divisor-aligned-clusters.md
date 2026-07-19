---
id: THM-1288
title: SUNGKAWICHAI–TRAKULTHONGCHAI CONJECTURE 7.1 IS FALSE FOR EVERY k ≥ 2 — for any d with (k+1) ∤ d, the translated cluster V_d = {d+1,…,d+k} is coprime and non-tight (strict band witness at t* = 1/(2d+k+1), margin d(k−1)/((k+1)(2d+k+1))) yet has NO witness time in (1/d)ℤ, because at the d-grid its distances EQUAL the AP's (‖(d+i)j/d‖ = ‖i·j/d‖) and no j/d with (k+1) ∤ d is an AP-witness (pigeonhole forces ‖m·j/d‖ = 1/(k+1) exactly for some m ≤ k, whence (k+1) | d). Hence no uniform D exists: their conjectured universal witness denominator fails at every k ≥ 2, at arbitrarily large d
status: PROVED — self-contained three-step elementary argument (band witness; translation vanishes on the d-grid; pigeonhole + quantization), independent of the S-T paper's content beyond the conjecture's verbatim statement (arXiv:2604.23906 §7, extracted 2026-07-19: "There exists a constant D such that, for any integer d ≥ D, every non-tight speed tuple v ∈ ℤ_{>0}^k with gcd(v)=1 has a witness time in (1/d)ℤ"). Exact verification k=12,13 at d ∈ {25, 100, 997, 10007, 100003} plus (k+1)|d controls (which correctly HAVE witnesses). The only extraction-dependence is the conjecture's own wording (uniform D); if the intended statement lets D depend on the family, see the repaired-conjecture discussion in the S401 reflection — the mechanisms bounding per-family D(v) are cataloged there
source: opus-2026-07-19-S401 (owner: keep working creative angles like HYP-7920)
depends_on: [the pigeonhole three-gap seed (k+1 points, gaps sum to 1 — elementary, in-proof), repo translation lore (S399 synthesis §3 axis 1; THM-1225 genus), S-T Remark 3.2 (the AP's witness set {s/(k+1)} — rederived in-proof via the same pigeonhole, so not load-bearing as citation)]
scripts: 04-computation/lrc_c71_refutation_divisor_alignment_opus_S401.py -> 05-knowledge/results/lrc_c71_refutation_divisor_alignment_opus_S401.out
---

# THM-1288 — C7.1 refuted by divisor-aligned clusters

## Statement

For every integer k ≥ 2 there is NO constant D = D(k) such that every non-tight coprime
speed tuple v ∈ ℤ_{>0}^k has a witness time in (1/d)ℤ for all integers d ≥ D.
Counterexamples: for any d with (k+1) ∤ d, the cluster V_d = {d+1, …, d+k}.

## Proof

**(1) V_d is coprime and non-tight.** Consecutive integers share no common factor. At
t* = 1/(2d+k+1) the positions (d+i)t*, i = 1..k, lie in the band
[(d+1)/(2d+k+1), (d+k)/(2d+k+1)], symmetric about 1/2 since (d+1)+(d+k) = 2d+k+1; the
minimum circle distance is (d+1)/(2d+k+1) = 1/(k+1) + d(k−1)/((k+1)(2d+k+1)) > 1/(k+1)
strictly (k ≥ 2, d ≥ 1). Strict witnesses exist, so V_d is non-tight.

**(2) The d-grid cannot see the translation.** For t = j/d and each i:
(d+i)·j/d = j + i·j/d, so ‖(d+i)·j/d‖ = ‖i·j/d‖. On (1/d)ℤ the cluster's distance profile
is identical to the AP {1,…,k}'s.

**(3) No grid point is an AP-witness when (k+1) ∤ d.** Suppose min_{1≤i≤k} ‖i·j/d‖ ≥
1/(k+1). The k+1 points {0, t, 2t, …, kt} (t = j/d) partition the circle into k+1 arcs
summing to 1, so two of them lie within 1/(k+1); their difference gives 1 ≤ m ≤ k with
‖m·t‖ ≤ 1/(k+1), hence ‖m·j/d‖ = 1/(k+1) exactly. But ‖m·j/d‖ is a multiple of 1/d′ with
d′ = d/gcd(j,d) | d, so (k+1) | d′ | d — contradiction.

By (2)+(3), no j/d is a witness for V_d; by (1), V_d is a legitimate non-tight coprime
family. Since every D admits some d ≥ D with (k+1) ∤ d, the conjecture fails. ∎

## Remarks

- **Mechanism (repo lore):** the (1/d)ℤ-restriction is blind to translations by multiples
  of d. V_d is the AP translated by d: its true M jumps to ≈ 1/2 − k/(4d), but the grid
  problem is exactly the AP's — whose witness set is the measure-zero {s/(k+1)}, invisible
  to the grid unless (k+1) | d. Translation-invariance blindness (S399 triage axis 1),
  fifth appearance, first strike outside the repo.
- **What survives:** per-family D(v) < ∞ always. Cataloged blowups: translation-alignment
  (this proof, D(v) ≳ v_min-scaled divisor structure) and near-floor width collapse
  (S400 probe, D(v) ~ v_max/(M − 1/(k+1))) — the latter constrained by HYP-7930's
  accumulation finiteness. A repaired conjecture must normalize by the family (e.g.
  D(v) ≤ poly(v_max)) — that version is open and is where the repo's gap program engages.
- Lean-ready: steps (2),(3) are decidable-per-instance arithmetic; the general statement
  is a five-line kernel target (pigeonhole via Finset).
