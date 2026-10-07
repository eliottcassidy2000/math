---
id: THM-4558
title: "The Heegner rungs Q(sqrt-3, sqrt-(4N-1)) with N = 2 mod 3 are 3-chromatic, the N = 3n rungs (n odd Loeschian) are 4-chromatic, the Heegner compositum is 4-chromatic and the Polymath field Q(sqrt-3, sqrt-11, sqrt-15) is 5-chromatic: residue colourings at primes fixed by complex conjugation (a KNOWN technique: Woodall 1973, Fischer 1990, Moorhouse 2010, Madore 2015, MildlyMeticulous 2026, decalion89 2026) applied to the repo's Hadwiger-Nelson field tower, refuting the HYP-2276/2277/2278 Heegner roadmap"
status: "The residue-colouring lemma (Lemma R) is KNOWN: decalion89 notes/local_colourings.md (dated 2026-09-24), Proposition A, credited there to Madore arXiv:1509.07023, Prop. 3.2; MildlyMeticulous hn-2adic-obstruction (2026-07), Theorem A/A' at places over 2; the coset/additive-colouring mechanism goes back to Woodall 1973, K. G. Fischer 1990 (Discrete Math. 82, Thm 1) and Moorhouse 2010. KNOWN values: chi(Q(sqrt2)^2) = 2, chi(Q(sqrt3)^2) = 3 (Madore 2015); chi(Q(sqrt3, sqrt11)^2) = chi(Q(sqrt-3, sqrt-11)) = 4 (K. G. Fischer 1994, known to us only through the zbMATH review quoted by decalion89; hn-2adic Cor. B; decalion89 Thm 1); the Heegner-compositum bound and Corollary 1 (implied by hn-2adic Cor. B / B'). NEW here (PROVED, elementary): the N = 2 mod 3 rungs are 3-chromatic (prime over 3); the N = 3n (n odd Loeschian) rungs are 4-chromatic (generalized spindle, checked exactly for N = 3, 9, 21, 27 by the audit; upper bound = hn-2adic); the explicit Polymath-field value; Corollary 2; the refutation of the repo's HYPs (MISTAKE-576). FINITE-EXACT: kappa(q) for q <= 11 and 16 (three independent SAT computations). INDEPENDENTLY AUDITED 2026-10-06 (audit A: every decomposition re-derived, PASS WITH CORRECTIONS, applied; the missing attributions above came from that audit; MISTAKE-579)."
session: mac-mini-2026-10-06-oaimath2 (the lemma and table were found by the session's reader of openai/math #158; re-checked here and by the audit)
source: 05-knowledge/results/oai2_openai_math_second_reading_20261006.md
scripts:
  - 04-computation/experiments/oai2_20261006_fields_ramsey_checks.py (+ .out, ALL CHECKS PASSED)
related:
  - THM-418 (lattice dichotomy; this is the field-plane version)
  - THM-412, THM-431, THM-440 (Eisenstein / Moser-field unit-distance work)
  - THM-4552 (G_7 = Cay(F_49, mu_8), so chi(G_7) = kappa(7) = 4)
  - HYP-2276, HYP-2277, HYP-2278 (historical index; the Heegner roadmap, refuted here; MISTAKE-576)
  - D. R. Woodall, Distances realized by sets covering the plane, JCTA 14 (1973) (Q^2 mod 2)
  - K. G. Fischer, Discrete Math. 82 (1990) 181-195, Thm 1; K. G. Fischer, A planar geometric graph of chromatic number four, Congr. Numer. 104 (1994) 73-79 (via zbMATH review)
  - G. E. Moorhouse, draft (2010), Lemmas 4.2, 8.2; Axenovich-Choi-Lastrina-McKay-Smith-Stanton, Graphs Combin. 2014, Thm 2.3
  - D. A. Madore, The Hadwiger-Nelson problem over certain fields, arXiv:1509.07023 (2015)
  - github.com/MildlyMeticulous/hn-2adic-obstruction (2026-07): Thm A/A', Cor. B/B', Thm C/C1
  - github.com/decalion89/chromatic-number-of-the-plane, notes/local_colourings.md (2026-09-24): Prop. A, Thms 1-2, chi(G_q) table (kappa(13) = 6, kappa(17) in [5,6], kappa(19) = 5), Prop. B, field screen
  - G. Exoo, D. Ismailescu, The chromatic number of the plane is at least 5: a new proof, DCG 64 (2020) 216-226, arXiv:1805.00157 (a 5-chromatic graph in Q(sqrt3, sqrt11, sqrt247)^2)
  - Polymath16 (2018): homomorphic 4-colouring of the ring Z[omega_1, omega_3] (P. Gibbs, a four-element group, per Goucher; D. Speyer, reduction mod 2, per decalion89); 5-colourings of Z[omega_1, omega_3, omega_4(, omega_7)]; several of Heule's 5-chromatic graphs lie in Z[omega_1, omega_3, omega_4] (Goucher, cp4space 2018)
---

# THM-4558 — residue colourings of the Hadwiger–Nelson field tower

## Setting

`L ⊂ C` is a number field with `c(L) = L` (c = complex conjugation) and `L ⊄ R`; `L⁺ = L ∩ R`, so `[L : L⁺] = 2`.
* The **field plane** `G(L)` has vertex set `L` and an edge `x ~ y` when `(x − y)·c(x − y) = 1`. For a real field `K`, the plane `K²` is `G(K(i))`.
* The Moser field `M = Q(√−3, √−11)` contains the Moser spindle. The rotations `ω_t` (real part `1 − 1/(2t)`, modulus 1) have minimal polynomial `t z² − (2t − 1) z + t` and generate `Q(√−(4t−1))`.
* The repo's Hadwiger–Nelson field tower (HYP-2276/2277) is `L_N = Q(√−3, √−(4N−1))`: a norm-`N` Eisenstein rhombus rotated so that its tips meet at unit distance.

For a prime power `q`, put **`κ(q) = χ(Cay(F_(q²), μ_(q+1)))`**, where `μ_(q+1)` is the kernel of the norm `F_(q²)^× → F_q^×`. This is decalion89's `G_q`.

## Lemma R (KNOWN; restated with proof)

Let `𝔓` be a prime of `L` with `c(𝔓) = 𝔓`.
1. Every unit vector `u` (`u·c(u) = 1`) is a `𝔓`-unit: `v_𝔓(u) + v_𝔓(c u) = 0` and `v_𝔓(c u) = v_(c𝔓)(u) = v_𝔓(u)`.
2. `c` induces an automorphism `c̄` of the residue field `k_𝔓`. The residues of unit vectors lie in `U_𝔓 = {z : z·c̄(z) = 1}`.
3. **`χ(G(L)) ≤ χ(Cay(k_𝔓, U_𝔓))`.**
   * The components of `G(L)` are cosets `x₀ + Γ`, where `Γ` is the additive group generated by unit vectors. `Γ` lies in the local ring at `𝔓`.
   * Colour `x ∈ x₀ + Γ` by `φ(red_𝔓(x − x₀))`, with `φ` a proper colouring of the Cayley graph.
   * An edge changes the residue by `red(u) ∈ U_𝔓`, which is nonzero.

`𝔓` cannot split over `L⁺`, since then `c𝔓 ≠ 𝔓`. So:

| `𝔓` over `L⁺` | `c̄` | `U_𝔓` | bound |
|---|---|---|---|
| ramified | identity | `{±1}` | `χ ≤ 2` if `char k_𝔓 = 2`, `χ ≤ 3` otherwise |
| inert, residue field `F_q` of `L⁺` | Frobenius `z ↦ z^q` | `μ_(q+1)` | `χ ≤ κ(q)` |

**κ(q).**

| q | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 | 13 | 16 | 17 | 19 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| κ(q) | 4 | 3 | 4 | 4 | 4 | 4 | 3 | 5 | 5–6 (6 per decalion89) | 4 | 5–6 | 5 (decalion89) |

* **Provenance.** The values for `q ≤ 11` and `q = 16` are FINITE-EXACT: three independent SAT computations (the reader, this session, audit A). For `q = 13, 17` there is no 4-colouring, and audit A found explicit 6-colourings; 5-colourability is undecided after 25–30-minute runs here and in the audit.
  * decalion89 reports `κ(13) = 6` from a completed cube-and-conquer run (not re-checked here).
* `κ(3) = 3`: `Cay(F_9, μ_4) = K_3 □ K_3`.
* **Hoffman** (NUMERICAL, margin `≥ 4·10⁻⁵`): `κ(p) ≥ 5` for primes `23 ≤ p ≤ 139` and `≥ 6` for `53 ≤ p ≤ 139` (audit A). decalion89, Prop. B, proves `κ(q) ≥ 6` for every prime `q ≥ 53` (Hoffman + Weil).
* **Knight torus.** The `7 × 7` knight set is the circle `x² + y² = 5 = (1+2i)μ_8` (THM-4552). So `G_7 ≅ Cay(F_49, μ_8)` and **`χ(G_7) = κ(7) = 4`**.

## Values

| plane | χ | proof | status |
|---|---|---|---|
| `Q(√−d)`, any `d` | `≤ 3`; `≤ 2` if `d ≡ 1, 2 (mod 4)` | a ramified prime (an odd prime dividing `d`, or 2) | field version of THM-418; also follows from Madore / hn-2adic |
| `Q(√2)² = Q(ζ_8)` | 2 | 2 is ramified in `Q(ζ_8)/Q(√2)` | KNOWN (Madore 2015, Prop. 3.9) |
| `Q(√3)² = Q(ζ_12)` | 3 | `(√3)` is inert in `Q(ζ_12)/Q(√3)`, `κ(3) = 3` | KNOWN (Madore 2015) |
| `M = Q(√−3, √−11)` and `M(i) = Q(√3, √11)²` | 4 | each prime over 2 (`g = 2`) is fixed by `c`, unramified over `L⁺`, residue `F_4`; lower bound the Moser spindle | KNOWN (Fischer 1994; hn-2adic Cor. B; decalion89 Thm 1) |
| **`L_N = Q(√−3, √−(4N−1))`, `N ≡ 2 (mod 3)`** | **3** | 3 is inert in `Q(√−(4N−1))` and ramified in `Q(√−3)`; the prime over 3 (`e = f = 2`, `g = 1`) is inert over `L⁺`, so `κ(3) = 3` | **NEW (PROVED)** |
| **`L_N`, `N = 3n`, `n` odd and Loeschian** | **4** | upper: 2 is inert in both quadratic subfields, residue `F_4` (= hn-2adic Thm C1). Lower: two copies of a rigid triangulated Eisenstein patch through the pivot 0 and a point `v ∈ √−3·Z[ω]` with `\|v\|² = N`, the second rotated by `ω_N` (checked exactly for `N = 3, 9, 21, 27` by audit A) | **NEW (PROVED)** |
| Heegner compositum `Q(√−3, √−11, √−19, √−43, √−67, √−163)` | 4 | all six `d ≡ 3 (mod 8)`, so 2 is inert in each `Q(√−d)` and `Frob_2 = c`; residue `F_4`. It is the only `c`-stable prime below 400 | bound implied by hn-2adic Cor. B′; stated here explicitly |
| Polymath field `Q(√−3, √−11, √−15)` | 5 | upper: each prime over 11 (`e = f = 2`, `g = 2`) is fixed by `c` and inert over `L⁺ = Q(√33, √5)`, `κ(11) = 5`. Lower: several of Heule's 5-chromatic graphs lie in `Z[ω₁, ω₃, ω₄]` (cited) | stated here explicitly (cf. decalion89 Thm 2 for `Q(√−3, √−11, √−247)`) |

* The `N ≡ 2 (mod 3)` row covers **every** lucky-Euler / Heegner rung `N = 2, 5, 11, 17, 41` (`4N − 1 = 7, 19, 43, 67, 163`).
* These `N` are not Eisenstein norms. So the field contains no rhombus-type spindle, and no unit-distance graph with `χ ≥ 4` at all.

**Corollaries.**
1. **(KNOWN; implied by hn-2adic Cor. B′.)** Any finite unit-distance graph whose edge vectors lie in `Q(ζ_6, ω_t : t odd)` is 4-colourable. For odd `t`, `4t − 1 ≡ 3 (mod 8)`, so every finite subcompositum has `Frob_2 = c`.
   * So no 5-chromatic graph has all its edge vectors in this field.
   * The known ones leave it via `ω₄` (de Grey, Heule), `ω₁₆` (de Grey), or `(119 + 3√−247)/128` (Exoo–Ismailescu).
2. **(NEW remark.)** If `G(L)` contains a graph with `χ ≥ 4`, then `L/L⁺` is unramified at every finite prime. So `L` lies in the narrow Hilbert class field of `L⁺` and `h⁺(L⁺)` is even. For example, `M` is the genus field of `Q(√33)`.
3. **The Heegner roadmap is false** (MISTAKE-576).
   * HYP-2277 ("the χ = 4 junction field ranges over exactly the class-number-one fields"; "each chromatic step adjoins a class-number-one rotation field; χ = 5 ↦ √−19"): `Q(√−3, √−(4N−1))` is 3-chromatic for every rung except Moser's, and `Q(√−3, √−19)` in particular.
   * HYP-2278 (4) ("χ = 2 + #independent Heegner rotations"): the compositum of the six rotation fields with `d ≡ 3 (mod 8)` is 4-chromatic, not 8.
   * HYP-2276's conjectural lower bound ("χ ≥ the number of pairwise-incommensurate imaginary-quadratic rotations forceable into one unit-distance graph") fails for the same reason.
   * **Class number does not govern the 5-chromatic step.** 5-chromatic graphs are known in `Q(√−3, √−11, √−15)` (Heule; `h(−15) = 2`) and in `Q(√−3, √−11, √−247)` (Exoo–Ismailescu; `h(−247) = 6`; exactly 5-chromatic by Lemma R at 11). decalion89 also reports one in `Q(√−3, √−7, √−11)`.
   * What bounds χ from above in these fields is how small primes decompose in `L/L⁺`.

## Not claimed

* Nothing here bounds `χ(R²)`. Every finite unit-distance graph has a realization in some conjugation-stable number field, but no single prime works for all fields.
* openai/math #158 (unrefereed; its Lean development was not built here) claims `χ(R²) ≥ 6`. If true, its compactness graph `H` can be realized with real-algebraic coordinates. Lemma R excludes realizations in the Moser, Heegner and Polymath fields and in the whole odd-`t` compositum (CONDITIONAL on #158).
* Fields with no conjugation-stable prime of small residue field are untouched. The field screen of decalion89 (§5) and the results note list the first open cases.
