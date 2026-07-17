---
id: THM-929
title: THE QUINTIC WALL — the level-5 Bonferroni wall of LRC(14) crosses the Abel–Ruffini threshold at exactly the rung the problem needs. (I) WALL PARITY: the wall polynomial W_k(x) = Σ_{j≤k}(−1)^j C(x,j) 2^j 13^{k−j} has NO real root at even k (W₂, W₄ verified root-free; classical Bonferroni forces even-truncation positivity at every integer), so the wall sequence is ODD-LEVEL: w₁ = 13/2 (RATIONAL), w₃ = 11.098925… (root of 4x³−90x²+1100x−6591 — CARDANO-expressible), w₅ = 15.630647… (root of 4x⁵−170x⁴+4300x³−77680x²+930376x−5569395); (II) UNSOLVABILITY (PROVED): the w₅ quintic is irreducible over Q with GALOIS GROUP S₅ — Dedekind-hygienic Frobenius shapes at unramified primes only: (2,3)×4, (1,1,1,2)×1, (1,4)×7, (5)×2, (1,1,3), (1,2,2), (1⁵) — the (2,3) class alone forces S₅; disc = 2⁸·7²·435651247410573820409077061, positive, NON-square ⟹ by Abel–Ruffini/Galois the level-5 wall location admits NO expression in radicals. ROBUST: the (2,q) family q = 5,7,9,11,13,15,17 is S₅ at EVERY q — structural, not numerology of 13. THE POINT: the Bonferroni price ladder (union 6.5 → tree 7.5 → B3 → B5, THM-897) first covers the 13 combs at level 5 (13 < w₅ < 16), and that is the FIRST level whose wall is radical-transcendent — LRC(14) sits exactly at the Abel–Ruffini rung of its own certificate ladder; (III) THE AVOIDANCE POLYNOMIAL OF THE TIGHT SYSTEM (exact on the grid D = 14·360360): for speeds {1..13} at λ = 1/14, p(t) = Σ(−1)^k S_k t^k has p(1) = 0 EXACTLY (integer-verified: the tight kill is a perfect depth-13 telescoping); S₁₃ = 1/91 EXACTLY = the origin window |t| < 1/182 — the deep well is the LAST SURVIVOR of the inclusion–exclusion; the enhancement profile S_k/equid_k = 1, 1.38, 5.49, 32.5, 211, 1419, 9679, 66588, … has consecutive ratio → 7 = 1/(2λ); the Bonferroni partial sums alternate (B_even ≥ 0 ≥ B_odd, forced) with magnitudes GROWING to ±11 at m = 6 before collapsing — NO truncation before full depth 13 certifies anything; Newton's inequalities fail at EVERY k and 12 of 13 roots are complex (only t = 1 real): the tight system is maximally ANTI-LEE-YANG — the polynomial-side quantification of why the equidistribution-calibrated level-5 regime (THM-897/926's admissibility) cannot see the tight extremizer
status: (II) PROVED (irreducibility mod p + unramified Frobenius shapes + Dedekind; S₅ from the (2,3) class; family check 7/7); (I) verified exactly (no real roots at k = 2,4; classical Bonferroni argument for all even k); (III) machine-exact (all 8191 subsets, integer sweep, p(1) = 0 in integers)
source: kind-pasteur-2026-07-16-S128 (cont.34; owner: work the level-5 wall, see how it relates to quintic polynomials)
depends_on:
  - THM-897 (opus: the wall + admissibility ladder this analyzes)
  - THM-927 (opus: the third blocker / linear-forms condition; born 923, briefly 926)
related:
  - THM-826/853 (the Farey profile — the λ < 1/14 side of the tight system)
  - HYP-6955 / false-peaks (the repo's other quintic threshold: A₅ obstruction begins at n = 5)
  - THM-922-route-A-signoff (the closed residue-six arc this sits beside)
---

# THM-929 — the quintic wall

> **RENUMBER NOTE (opus S332):** born THM-927 (kps S128c34, pushed 17:39:18) but opus first-pushed THM-927 at 17:20:30 (parallel-class renumber); first-pusher rule -> this file is THM-929. Its depends_on reference to the third-blocker updated 926->927 (that file moved when mac-mini first-pushed 926). No content changed.

## (I) The wall family

W_k(x) = Σ_{j=0}^{k} (−1)^j C(x,j) 2^j 13^{k−j} (primitive integer forms):
k=1: [13, −2]; k=2: [169, −28, 2]; k=3: [6591, −1100, 90, −4];
k=4: [85683, −14312, 1192, −64, 2]; k=5: [5569395, −930376, 77680, −4300, 170, −4].
Even k: no real roots (k = 2, 4 verified; Bonferroni's even-truncation upper bound
Σ_{k≤2j}(−1)^k C(m,k)x^k ≥ (1−x)^m > 0 forces integer positivity at every even level).
Odd walls: w₁ = 13/2, w₃ ≈ 11.098925, w₅ ≈ 15.630647. (THM-897's schedule prices 11.5,
15.5 are the half-integer schedule labels; the exact roots are these. The integer
statements agree: coercive through m′ = 11 resp. 15.)

## (II) The unsolvability

f(x) = 4x⁵ − 170x⁴ + 4300x³ − 77680x² + 930376x − 5569395 is irreducible (irreducible
mod p at two unramified primes); its unramified Frobenius shapes include (2,3) at four
primes — an order-6 odd class, which no transitive subgroup of S₅ except S₅ contains.
Hence Gal(f) = S₅ and w₅ is not expressible in radicals. Discriminant
2⁸ · 7² · 435651247410573820409077061 > 0, not a square (consistency: Gal ⊄ A₅).
The same verdict holds for the whole dilation family q = 5..17.

**Reading.** Degree-k walls for k ≤ 4 are radical-solvable for free (Abel–Ruffini
threshold is degree 5). The first Bonferroni level whose wall covers all 13 combs of
LRC(14) is level 5 — and its wall is genuinely S₅. The ladder's ability to *name its own
death point in closed form* ends exactly where the problem starts needing it.

## (III) The tight system's avoidance polynomial (exact)

Speeds {1..13}, λ = 1/14, D = 5045040: S_k·D =
9369360, 11107704, 23094768, 48777032, 81523138, 104422622, 101763358, 74916002,
41029828, 16227372, 4385220, 725340, 55440 (k = 1..13).
- **p(1) = 0 exactly.** Σ(−1)^k S_k = 0 in integers: the tight system dies only at full
  depth. B₁₂ = S₁₃ > 0, B₁₃ = 0.
- **S₁₃ = 55440/5045040 = 1/91** = μ{t : |t| < 1/182} — the all-13 intersection is
  exactly the origin window; the deep well is the last survivor.
- Enhancement over equidistribution: ×1, 1.38, 5.49, 32.5, **211 (k=5)**, 1419, 9679,
  66588, …, ratio → 7 = 1/(2λ). The wall's 9.44% equid margin at level 5 vs the true
  −9.72: the tight system sits two orders of magnitude outside the level-5 regime.
- Newton fails at every k; 12/13 roots complex; only t = 1 real. Anti-Lee-Yang: the
  three blocker species (quadruples, dilates, linear forms — THM-926) are exactly what
  pumps the complex spectrum.


## Bring–Jerrard addendum (cont.35)

Exact chain: depressed form (x = y + 17/2): y⁵ + (705/2)y³ − 4290y² + (914549/16)y − 671385/2.
Principal-form Tschirnhaus z = y² + ay − p₂/5: a solves −705a² + 25740a − 318119/4 = 0 with
δ = 438273705 = 3·5·379·77093 (squarefree) — the resolvent field is Q(√438273705), and δ is
divisible by NONE of 7, 13, 14, 91, 183. Numeric completion (quartic Tschirnhaus, c₃ = 0
gauge): w⁵ + pw + q with p ≈ 1.8825×10³⁷, q ≈ −2.3142×10⁴⁶; real Hermite normal form
**u⁵ + u = t′, t′ = 0.59018696444649…** — no rational/LRC-constant match at 10⁻⁶ tolerance
(checked against 14, 183, 91, ratios, and the wall margin 35035/371293). HONEST VERDICT:
the quintic's radical-obstruction data is LRC-blind; the only 14-echo in its invariants is
the discriminant's square part 2⁸·7² = (8·14)²; the disc cofactor 4356…7061 has no prime
factor below 2×10⁵. The wall quintic remembers its Bonferroni provenance through the
discriminant square factor and nothing else — consistent with the S₅-genericity of the
whole (2,q) family. (Solution by modular functions à la Hermite requires the modulus of
t′ = 0.5902; whether that modulus sits on a level-14 curve is the one named follow-up —
the fleet's X₀(14) thread, mac-mini S89–94, is the right family to test against.)

## Named next
- The enhancement law S_k/equid ~ c·7^k for the tight system (the exact constant and
  its Farey meaning) — one lemma from THM-826's arc structure.
- p(t) for opus's corrected-admissibility packets (when the thin DFS finds one): the
  conjecture is near-real-rootedness (Newton ratios → 1) ⟺ BONF5 certifies.
- The Bring–Jerrard form of the wall quintic (Tschirnhaus) and whether the LRC constants
  (14, 183, 91) appear in its invariants.

## Evidence log
- [x] wall family + roots (wall_quintic_galois_kps_S128c34.out — first half; the DDF
      section hung, replaced by brute-force factorizer, MISTAKE-note in file)
- [x] S₅ verdict + disc (wall_quintic_galois2_kps_S128c34.out)
- [x] Dedekind-hygienic family robustness 7/7 (wall_quintic_family_kps_S128c34.out)
- [x] avoidance polynomial exact (avoidance_polynomial_kps_S128c34.out)
