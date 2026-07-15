---
id: THM-874
title: THE MÖBIUS-LOG² GRAMMAR OF THE FAREY LADDER — the depth-layer constants of the interval-core measure profile assemble into Σ_d (μ(d)/d²)·log²(1/(1−x^d)), i.e. [x^s] = (2/s)·H*(s) for every depth s; prime depths are the pure-log² terms (THM-819's harmonic law), composite depths are Möbius-corrected; COROLLARY: the first composite corridor law m({1..13}; λ) = (H*(14)/7)·(1 − 14λ) on [1/15, 1/14] — the LRC(15) deep-well corridor, constant 1666/6435, coprime filter {1,3,5,9,11,13} live
status: PROVED (two-line GF proof below) + machine-exact (identity to s = 40 in ℚ; Farey corridor constants k = 3..13; LRC(15) corridor at 3 rational points)
source: mac-mini-2026-07-15-S111 (owner: LRC(14) through the new lenses, think roots of unity)
depends_on:
  - THM-826 (the Farey profile ladder m = Σ_gaps (1 − λ(i+j))₊/(ij) — this file is its grammar)
  - THM-819 (the primitive harmonic law = the prime-depth case)
related:
  - THM-853(II) (the k=12/prime-13 corridor, subsumed as the depth-13 coefficient)
  - THM-873 kps cont.25 (Ramanujan-Fourier expansion of the good SET — the Fourier face of the same coprime structure; this file is the measure-PROFILE face, exact in ℚ)
  - THM-868 (contrast: the figurate ladders are geometric in ONE substituted variable; the Farey ladder's grammar needs the full Möbius scale tower — see Reading)
script: 04-computation/lrc14_lenses_moebius_rigidity_macmini_S111.py -> 05-knowledge/results/lrc14_lenses_moebius_rigidity_macmini_S111.out
---

# THM-874 — the Möbius-log² grammar

**Theorem.** Let H*(s) = Σ_{i<s, gcd(i,s)=1} 1/i (the primitive harmonic number). Then

> Σ_{d≥1} (μ(d)/d²) · log²(1/(1−x^d))  =  Σ_{s≥2} (2·H*(s)/s) · x^s,

and (2/s)H*(s) is exactly the depth-s layer constant of THM-826's profile ladder: for the
interval core {1..k}, the Farey gaps with i+j = s contribute total length
A_s-layer = Σ_{i+j=s} 1/(ij) = (2/s)H*(s) whenever s ≤ k+1 (all coprime pairs realized),
each dying linearly at λ = 1/s. In particular the FIRST SEGMENT of every profile is
m({1..k}; λ) = (2H*(k+1)/(k+1))(1 − (k+1)λ) on [1/(k+2), 1/(k+1)] — the general corridor
law, prime case = THM-853(II).

*Proof.* log²(1/(1−y)) = (Σ_{i≥1} y^i/i)² has coefficients [y^s] = Σ_{i+j=s} 1/(ij) =
(2/s)H_{s−1} (split 1/(ij) = (1/s)(1/i + 1/j)). Filtering pairs by gcd(i,j) = 1 via
1_{gcd=1} = Σ_{d | gcd} μ(d) rescales (i,j) = (di′, dj′): Σ_d μ(d)/d² · [coefficient at
s/d], which is the displayed GF; on the coefficient side the same filter turns H_{s−1}
into H*(s) (i coprime to s ⟺ coprime to j = s−i). The Farey identification is THM-826's:
gaps of F_k ↔ ordered coprime pairs (i,j), i+j > k, realized once each. ∎

**Corollary (the LRC(15) corridor — first composite instance).** On [1/15, 1/14]:

> m({1..13}; λ) = (H*(14)/7)·(1 − 14λ),  H*(14) = 1 + 1/3 + 1/5 + 1/9 + 1/11 + 1/13 = 11662/6435,

machine-exact at λ = 1/15, 15/211, 1/14. The Möbius filter (depth 14 = 2·7 composite)
removes the even and 7-divisible chairs — the first corridor where THM-819's coprimality
is visible (all previous deep-well corridors had prime depth 13). Möbius form:
H*(m) = Σ_{d|m} (μ(d)/d)·H_{m/d−1} (verified m < 40).

## Reading (the grammar taxonomy, THM-868 contrast)

THM-868's figurate ladders are geometric series in ONE substitution variable (u or v) —
"solvable" towers with product-form tails. The Farey ladder REFUSES a single substitution:
its grammar is a Möbius-weighted sum of log² over ALL scales d — the scale tower is
irreducible. Same lens, opposite outcome, and the outcome is the right statement of why
interval-core arithmetic is hard: the profile is not one geometric ladder but a ζ/μ-dressed
family of them (log² = the depth-2 harmonic ladder; μ(d)/d² = the scale dressing). kps's
THM-873 (Ramanujan-Fourier of the good set) is the same taxonomy statement on the Fourier
side; the two faces should be merged in the eventual write-up.
