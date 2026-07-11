---
id: THM-690
title: The 13-torsion test-point theorem — for EVERY co-offset set E with |E| = k ≤ 12 and every a coprime to 13, m_E(a/13) ≥ 2/13 − 1/7 = 1/91 (the ≤ 12 centers cannot fill the 13-net: an empty class forces a gap ≥ 2/13 — SHARP, attained by consecutive E at k = 12); and for every P with pmax ≤ 12, a/13 lies in the interior of G_P with clearance 1/182 (13 is prime, so pa mod 13 ∈ {1..12} and frac(pa/13) ∈ [1/13, 12/13] ⊂ (1/14, 13/14)); hence μ∞(P,E) > 0 for EVERY two-scale class with pmax ≤ 12 — ALL E, ALL k, the E-quantification ELIMINATED; consecutive-maximality (THM-689(B)) is SUPERSEDED everywhere except the corner 13 ∈ P, which is exhaustively P-censused (794 sets × adversarial E: zero degenerate, min μ∞ = 0.0441)
status: PROVED (the two-part proof below is elementary counting + prime arithmetic + continuity; the E-side margin verified sharp — worst over nasty batteries = exactly 1/91 at consecutive k = 12). The corner (13 ∈ P): per-class decidable, exhaustively censused over ALL P's containing 13 (sizes 1..5 = 794 sets, k = 8..12) against adversarial E batteries (consecutive/detuned/AP2/random, 3176 evaluations): ZERO zeros, min μ∞ = 0.044079 at P = {1,3,4,13}, E = {0..8}. Wholesale closure of the corner = the downgraded remaining question (consecutive-maximality now needed ONLY there, if at all — per-class checks cover practice).
source: klein-2026-07-10-S239 (HYP-5920; superseding the consecutive-maximality lead of THM-689)
depends_on:
  - THM-687 / THM-688   # the limit measures
related:
  - THM-689 (k ≤ 7 rigidity; the measure criterion this supersedes at pmax ≤ 12)
  - klein-S207 LRCRulerPoints (the inversion: a family is never lonely on its own ruler; the test ruler is the COMPLEMENT's)
  - kps-S127 (covering ⟹ q ≥ 15; the bounded small-modulus obligation — complementary modular structure)
---

# THM-690 — the 13-torsion test-point theorem

## Statement and proof

Let E ∋ 0 be a set of k ≤ 12 distinct non-negative co-offsets, P ⊆ [1,13]
with pmax ≤ 12, and a ∈ {1,…,12}. Then:

**(E-side.)** At α = a/13 the centers {e·a/13 mod 1 : e ∈ E} lie on the
13-net {i/13} and occupy at most k ≤ 12 of its 13 classes. Some class is
empty; the circular gap across an empty class is ≥ 2/13; hence

> m_E(a/13) ≥ 2/13 − 1/7 = **1/91**,

uniformly in E — no hypothesis on the shape or size of the co-offsets.
SHARP: consecutive E at k = 12 attains exactly 1/91 (verified).

**(P-side.)** 13 is prime and p ≤ 12, so pa ≢ 0 (mod 13):
frac(p·a/13) = (pa mod 13)/13 ∈ [1/13, 12/13] ⊂ (1/14, 13/14), with
clearance min(1/13 − 1/14, 13/14 − 12/13) = **1/182**. So a/13 lies in the
INTERIOR of G_P — for every a.

**(Conclusion.)** m_E is continuous (piecewise linear) and positive at an
interior point of G_P, so it is positive on a neighborhood inside G_P:

> **μ∞(P,E) = ∫_{G_P} m_E > 0 for every two-scale class with pmax ≤ 12 —
> every E, every k ≤ 12.** ∎

(Quantitative version: the neighborhood radius min(1/(182·pmax),
1/(182·2Σ_E e)) gives an explicit per-class floor; positivity is what the
wall pipeline consumes.)

## Why 13, structurally

The test modulus must satisfy q > k (so ≤ k centers cannot fill the q-net)
and 2/q > 1/7 (so an empty class certifies a gap), i.e. k < q < 14. With
k ≤ 12 the window contains exactly q = 13: **the test ruler is n − 1 for
the n = 14 problem.** And the P-side works precisely because 13 ∉ P — the
family's own 13-runner would sit at 0 at every 13-torsion time. This is
klein-S207's ruler-points theorem inverted: a family is never lonely on its
own ruler; the universal witness ruler is the COMPLEMENT's ruler.

## The corner 13 ∈ P

The only two-scale classes the theorem misses. Exhaustive census: all 794
P's containing 13 (|P| = 1..5, k = 8..12) × adversarial E batteries
(consecutive, detuned, AP₂, random — 3176 exact evaluations): **zero
classes with μ∞ = 0; minimum 0.044079 at P = {1,3,4,13}, E = {0,…,8}** (the
consecutive apex shape, as at every previous extremum). Wholesale closure
of the corner is the one place consecutive-maximality (THM-689(B)) would
still serve; per-class exact decidability covers it in practice, and the
corner P's are the tight-AP shadow (they contain 13 — the knife-edge
families), consistent with everything the project knows about where the
hardness lives.

## Status of the two-scale dead-zone question after this theorem

[k ≤ 7: THM-689(A) rigidity, unconditional] ∪ [pmax ≤ 12: THIS THEOREM,
unconditional, all E] ∪ [corner 13 ∈ P, k ≥ 8: P-exhaustive census, all
positive, per-class decidable]. The E-quantification — the genuinely
infinite side — is eliminated everywhere except the corner, where P is
finite.

## Verification & files

`04-computation/lrc14_test_point_theorem_klein_S239.py` (+ `.out`): the
E-side sharpness battery (worst = exactly 1/91), the P-side arithmetic, the
3176-evaluation corner census.

## Formalization (S243)

`LRCTestPointCore.lean` (kernel-pure, root-wired): `net_value_strictly_in_band`
(the P-side band bound), `qstar_p_nonzero` + `residue_in_range` (the
coprimality nonvanishing), `pigeonhole_missed_class` (the E-side), and THE
FATTENING LEMMA `qstar_cert` — the test point generates SafeIvStrict
certificates on [(4732a−q)/4732q, (4732a+q)/4732q], feeding
LRCMeasureTransfer's transfer pipeline directly.
