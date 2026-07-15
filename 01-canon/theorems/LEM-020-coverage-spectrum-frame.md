---
id: LEM-020
title: THE COVERAGE-SPECTRUM FRAME for the LRC(14) covering case — six exact identities on the multiplicity spectrum μ_k(x) of the 1/7-arcs: sum rules, S_d = Σ C(k,d)μ_k, the truncation-deficit identity (B_D deficit = Σ C(k−1,D)μ_k ⟹ B_D EXACT ⟺ max multiplicity ≤ D), the covering-slack identity (covering x ⟹ overlap excess EXACTLY 6/7), the RIGIDITY of the minimal covering spectrum ((μ₁,μ₂) = (1/7, 6/7) is forced — the LRC twin of the K(A₅) multiplicity-≤2 Kakeya witness), and the MID-BAND LAW (maxC(AP, 1/q) = min(13, ⌊q/7⌋+[7∤q]): B5's exactness zone is the mid-band; the far tail is the 13-deep single-cluster stack)
status: PROVED (each identity 1–3 lines; all verified exactly in rational arithmetic on 4 instances × 7+ rational x, plus q-sweeps 7..200; 8/8 checks)
source: klein-2026-07-15-S313 (cont.4); the (r,g) truncation grammar (HYP-6946) and the Kakeya tightness mechanism (THM-870) transferred to the covering case
depends_on: [THM-604 (the truncation identity in use), THM-663 (covering-case frame), THM-867/868/869 (the corona grammar)]
related: [THM-863 (ρ pair floor = d=2 summand floors), THM-667 (mid-band realization), LEM-A/THM-710 (far-tail transfers), THM-856 (Kakeya-rank first-moment frame)]
verification: 04-computation/coverage_spectrum_bonferroni_corona_klein_S313.py -> 05-knowledge/results/coverage_spectrum_bonferroni_corona_klein_S313.out
---

# LEM-020 — the coverage-spectrum frame

Cluster E = {e_1,…,e_j} (j ≤ 13), slow variable x; arcs A_i = [frac(e_i x), frac(e_i x)+1/7);
C(t) = #{i : t ∈ A_i}; **spectrum** μ_k(x) = meas{t : C(t) = k}. Everything below is exact.

**(R0)** Σ_k μ_k = 1.  **(R1)** Σ_k k·μ_k = j/7 (arc mass).
**(R2)** μ₀ = (7−j)/7 + Σ_{k≥2}(k−1)μ_k — *uncovered = overlap excess − (j−7)/7*.
 ⟹ **covering-slack identity**: at covering x (μ₀ = 0, j = 13) the overlap excess is
 EXACTLY 6/7. The covering adversary has a hard overlap budget; loneliness ⟺ making
 overlaps exceed budget somewhere.

**(S)** S_d := Σ_{|T|=d} meas(∩_{i∈T} A_i) = Σ_k C(k,d)·μ_k — subset sums are spectrum moments
(no subset enumeration ever needed; cross-checked literally at d = 2, 3).

**(T)** Σ_{d≤D}(−1)^d S_d = μ₀ + (−1)^D Σ_{k≥1} C(k−1, D) μ_k. For odd D this is THM-604's
Bonferroni with the deficit now EXACT and nonneg:
> **deficit(B_D) = Σ_{k≥D+1} C(k−1,D) μ_k, so B_D is EXACT iff max multiplicity ≤ D.**
The whole depth ladder B₃/B₅/B₇ reads binomial-tail moments of ONE object, the spectrum —
the (r,g) truncation grammar verbatim (rank-D truncation of (1−1)^C; the deficit is the
Vandermonde/binomial tail; the first failure is the cheapest (D+1)-fold stack = the corona
onset, here a (D+1)-fold resonance — THM-598/604's "depth-D resolution" hypotheses are
exactly "no cheap champion stacks").

**(RIGIDITY)** If x is covering with max multiplicity 2 (the least possible, since arc mass
13/7 > 1), then (R0)+(R1) force **the UNIQUE spectrum (μ₁, μ₂) = (1/7, 6/7)**. The tight AP
{1,…,13} at x = 1/7 realizes it (as do 1/13, 1/14 — 5 instances found in the sweep, all with
the exact profile). This is the LRC twin of THM-870's Kakeya optimum: there, points in ≤ 2
needles made Bonferroni-2 tight (the witness = the K₆ pair-incidence); here, multiplicity ≤ 2
makes the cover minimal and pins the spectrum. **Tight objects are minimal-multiplicity
exact covers, in both worlds.**

**(MID-BAND LAW)** For the AP cluster at x = 1/q: maxC = min(13, ⌊q/7⌋ + [7∤q]) (verified
q = 14..200). Hence: covering band q ∈ {7,…,14} has maxC = 2 (rigid spectra, B₂-exact);
**B5's exactness zone is maxC ≤ 5 ⟺ q ≤ 35 plus the covering band; for q > 91 the cluster
stacks 13-deep** (spectrum degenerates toward μ₀ = 6/7·(1 − o(1)), the maximally-lonely
single-cluster band). The density route's architecture — finite mid-band checks + explicit
far-tail transfers (THM-710/Lemma A) — is FORCED by the spectrum: the far tail is precisely
where B_D leaks and precisely where the cluster is one stacked block with an explicit profile.

**(VEIN)** S₂ = pair-overlap sum: the AP at generic x sits far below the independent value
C(13,2)/49 (equal spacing = maximal pair correlation; random clusters approach independence,
TV-to-Binom(13,1/7) quantifies resonance). THM-863's ρ(a,b) ≥ 1/78 floors these summands;
the (1,12) minimizer is the d = 2 champion pair — the corona onset one level down.

## Addendum (same session): the Fejes Tóth floor — the second-moment wall, exact

**(W)** Covering x forces S₂ ≥ 6/7 (convexity of C(k,2): minimal covering support {1,2}), so
`S₂(x) < 6/7 ⟹ x is a witness` — a criterion from the DIFFERENCE SET alone:
S₂(x) = Σ_{i<j} max(0, 1/7 − ‖(e_i−e_j)x‖).
**(FT)** But the criterion is VACUOUS, and its vacuity is a theorem: S₂ is a pair energy with
convex decreasing kernel, so by Fejes Tóth the regular 13-gon minimizes it:
min = 13·(1/7 − 1/13) = **6/7 exactly, for ALL 13-point configurations** (4000 exact random
trials; the tight AP attains it at x = 1/7). The universal pair-energy minimum EQUALS the
covering budget. Hence **pair data can never certify loneliness** — the repo's second-moment
wall ("unsigned OffLine bounds provably fail", the transform blindness memory), now pointwise,
exact, and attributed: it is Fejes Tóth exactness. The minimal certifying moment order is 3.
**(M3)** S₂ − overlap excess = Σ C(k−1,2)μ_k, so for covering x: **S₂ − 6/7 = the ≥3-fold
mass** — one number measuring distance-from-rigidity both configurationally (FT deficit) and
spectrally (champion-stack mass). The d = 3 pincer should be built on (M3), not on raw S₃.

## Next steps (named)
1. Characterize ALL covering x with maxC = 2 across admissible clusters (are they exactly the
   clock witnesses of the sieve? the 5 found are all q ∈ {7, 13, 14}).
2. The deficit identities are decide-shaped exact rational statements — Lean batch candidate
   alongside LRCDiscreteBonferroni.
3. d = 3 pincer: combine (S) at d = 3 with the 3-speed floor c₃ = 2/105 (THM-863 T3) to squeeze
   the spectrum between second- and third-moment constraints on covering sets.
