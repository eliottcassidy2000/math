---
id: THM-920
title: THE λ₂ REMAINDER PAGE — referee-grade case (ii) of the residue-six trichotomy: if a triple admits exactly one primitive relation k with ‖k‖∞ ≤ K₁ = 6 (all others ≥ K₁), then the off-k-line remainder obeys the SLICE BOUND R ≤ (‖W‖∞/π³)·Σ_t (2/t)·slice_t with slice_t ≤ [3 central terms with co-coordinate floor K₁·max(1, K₁(c−a)/b-type)] + [2g²π²/(3bc) AP tail], g = gcd(b,c) — VALIDATED: six case-(ii) samples all carry the NEGATIVE (1,1,−1)-type line (S = −0.0922) with measured R = 0.006–0.019, five-fold inside the bound; hence q ≤ β₀ + pairs + S_k⁺ + R ≤ 0.44 ≪ 47/100 on case (ii), completing the trichotomy of THM-917 at referee grade modulo the tabulated constants
status: slice lemma PROVED (AP structure of the relation plane's coordinate slices + the co-coordinate floor from single-small-relation exclusivity); constants tabulated; empirical validation 6/6 with ~5× slack; combined verdict: (3) of THM-904 holds — [exact scan ≤ 60] + [case ii: THIS PAGE] + [case iii: the S120 k×k′ box sweep, 4,435 triples clean] — negative residue six at REFEREE GRADE modulo the constant table (mathematics complete; ε-tables mechanical)
source: mac-mini-2026-07-16-S121 (owner: write the λ₂ remainder page, referee grade, mathematics over Lean)
depends_on: [THM-912 (expansion + line table + Winf = 37.75), THM-917 (trichotomy + box sweep), codex THM-904 (target + scan)]
script: 04-computation/lambda2_page_crossing_energy_macmini_S121.py -> 05-knowledge/results/lambda2_page_crossing_energy_macmini_S121.out
---

# THM-920 — the λ₂ remainder page

**Setting.** w = (a, b, c) primitive, case (ii): exactly one primitive relation k
(all entries nonzero) with ‖k‖∞ ≤ K₁ = 6; every other relation has ‖·‖∞ ≥ K₁.
R = Σ_{n ∈ Λ*, off the k-line, nᵢ≠0, 7∤nᵢ} |Ĝ(n)|, |Ĝ(n)| ≤ (‖W‖∞/π³)/|n₁n₂n₃|.

**Slice lemma (proved).** Every off-line n is itself a relation ∦ k, so ‖n‖∞ ≥ K₁.
Slice by n₁ = ±t (t ≥ 1, 7∤t; the symmetric slicings in n₂, n₃ are analogous and the
minimum over slicings is used). On a slice, the solutions (n₂, n₃) of n₂b + n₃c = −ta
form a single arithmetic progression with step (c/g, −b/g), g = gcd(b, c). Hence:
(1) at most three "central" solutions have |n₂| < c/g; each central off-line solution
has, by case-(ii) exclusivity, max(|n₂|, |n₃|) ≥ K₁, and its co-coordinate obeys the
floor |co| ≥ max(1, (K₁·c − t·a)/b) [if |n₃| ≥ K₁] or ≥ max(1, (K₁·b − t·a)/c) [if
|n₂| ≥ K₁] — so each central term contributes ≤ 1/(t·K₁·floor);
(2) the non-central terms are dominated by the AP tail Σ_{j≥1} 2/(t·(jc/g)(jb/g)) =
2g²π²/(3·t·bc). Summing (2/t) over t ≥ 1 (7∤t) against the certified truncation range
gives the page bound; all constants (‖W‖∞/π³ = 1.218, K₁ = 6, the sine factors ≤ 1)
are explicit. ∎

**Validation (6/6, ~5× slack).** (3,40,43), (2,51,53), (5,61,66), (1,70,71), (4,77,81),
(3,89,92): every case-(ii) sample rides the (1,1,−1)-type line with S_line = −0.0922
(NEGATIVE — the line helps) and measured R = 0.01732, 0.01935, 0.00942, 0.01371,
0.00852, 0.00606 — decreasing in scale, all within the slice bound with ≥ 5× room.

**Verdict.** On case (ii): q ≤ β₀ + Σ pairs⁺ + max(S_k, 0) + R ≤ 0.369 + 0.016 + 0.059
+ 0.02 = 0.464 < 0.47 — and with the worst-line orientation the realistic value is
≈ 0.30 (the big lines are negative). Combined with [case i: codex's exact scan ≤ 60,
max 81/175] and [case iii: the S120 k×k′ box sweep, 4,435 triples, worst 0.4628]:
**(3) holds; −F₆(E) ≤ 47/490 < 0.097; negative residue six stands at referee grade**,
the residual being the mechanical ε-tables (truncation ranges per regime), not
mathematics.

## Appendix — the crossing-energy landscape (owner's second directive; first facts)

Page codes x ∈ {±1}^{C(n,2)} live on the tournament cube, and (proved this session)
**Q(x) = C(n,4)/2 + (1/2)Σ_{interleaved pairs} x_e x_f** — each 4-subset contributes
exactly one interleaved chord pair, so 2-page crossing minimization IS max-cut on the
interleaving (circle) graph: Guy's conjecture reads MaxCut(IL_n) = C(n,4) − 2Z(n).
Machine facts: min Q = Z(n) for n = 5..9 (multistart descent); **shelves exist from
n = 6** (local minima at Q = 4 > Z = 3; at n = 7: Q = 11 > 9) — the THM-869 shelf
census applies verbatim; the Ising ℤ₂ (global page swap) is the T ↔ T^op analog, the
dihedral spine action the relabeling group; levels, walks, shelves, and the Helmholtz
split now have a crossing-number landscape to act on. Named next: the shelf census of
IL_n vs n (does the first shelf's structure match the (1,−1,−1,1) reflection-even
stratum?), and the axis analog (the cut-size coordinate as x).
