---
id: THM-4592
title: "The eleven-squares Fibonacci torus is Collatz's period-5 cycle set: reading the parity word of a periodic point of the standard Collatz map (x -> 3x+1, x/2) in base phi gives kappa_n : {C-periodic points of period dividing n} -> R_n = Z[phi]/(phi^n - 1) with kappa_n(Cx) = phi kappa_n(x), a bijection for odd n (FINITE-EXACT, n <= 22; for even n exactly 0, -1, -2 collapse); for n = 5 the eleven cells of the Fibonacci torus of the eleven-squares paper are 0, the 1/13 cycle (labels = the quadratic residues H = {1,3,4,5,9} mod 11, the paper's gold cells) and the -5 cycle (labels -H, the teal cells); the paper's group G_5 = <x+1, phi x> on F_11 is BS(1,3) = <x+1, 3x> mod 11, the relation group of THM-4581"
status: >
  PROVED (elementary): equivariance, the labels at n = 5, and G_5 = BS(1,3) mod 11 (since <3> = <4> = H in F_11^x).
  FINITE-EXACT: bijectivity of kappa_n for n <= 22 (golden reader), and the n = 10 identification of the 121-cell cover.
  DICTIONARY: the paper's statements about R_5, R_10 and G_5 become statements about rational Collatz cycles.
  The eleven-squares paper (Pingyou Ltd) is a trusted source per owner directive. Found by the session's golden reader
  (inheriting collatz_golden_carriers_20261004 section 4, which had the denominators 1, 2, 11, 76); re-verified here.
session: mac-mini-2026-10-07-golden
source: 05-knowledge/results/golden_collatz_resonance_20261007.md
scripts:
  - 04-computation/experiments/golden_20261007_readers/golden/ (golden_cycle_dictionary.py, pair_regularity.py + .out)
related:
  - THM-4528 (Collatz parity in base phi: the golden beta map; Theta), THM-4581 (BS(1,3) is the pair-relation group), THM-4553 (the p = 7 Borel case)
  - The eleven-squares paper, "The Fibonacci Geometry of Eleven Squares", Theorems C and D, Proposition 5.4
---

# THM-4592 — the Fibonacci torus of eleven squares is Collatz's period-5 cycle set

## Setting

* `C` is the standard Collatz map on `Z_2 ∩ Q`: `C(x) = 3x+1` for odd `x`, and `x/2` for even `x`.
* Its parity words never contain `11`, because `3x+1` is even. These are exactly the words of base-φ normal form (Bergman, Zeckendorf).
* A point of period dividing `n` has a cyclic parity word `w = (w_0, …, w_(n−1))`. It is the rational fixed point of the composed affine map, `x_w = b_w/(1 − 3^k 2^(−(n−k)))`, where `k` is the number of odd steps.
* `O = Z[φ]` with `φ² = φ + 1`, and `R_n = O/(φ^n − 1)`. The eleven-squares paper has `|R_n| = L_n − 1 − (−1)^n`.

## Statements

1. **Golden reading.** `κ_n(x) = Σ_j w_j φ^(n−1−j) mod (φ^n − 1)` satisfies `κ_n(C x) = φ·κ_n(x)`. The map `C` is a shift of the cyclic word, and the wrap-around digit dies modulo `φ^n − 1`.
2. **Counting.** The cyclic words without `11` of length `n` number `L_n` (Lucas). For odd `n`, `κ_n` is a bijection onto `R_n`. For even `n`, the words `0…0`, `1010…` and `0101…` all read as 0; these are the points 0, −1 and −2, the cycle `−1 → −2`. FINITE-EXACT for `n ≤ 22`.
3. **`n = 5`: the eleven cells.** `R_5 ≅ F_11` via `φ ↦ 4`. The eleven period-5 points and their labels:

   | point | label |
   |---|---|
   | 0 | 0 |
   | 1/13 cycle: `1/13, 2/13, 4/13, 8/13, 16/13` | `3, 9, 5, 4, 1` = `H`, the quadratic residues (gold) |
   | −5 cycle: `−5, −14, −7, −20, −10` | `8, 10, 7, 6, 2` = `−H`, the non-residues (teal) |

   * `C` acts on labels by multiplication by `φ = 4`, which has order 5.
   * The paper's cell exchanges `S` (`x ↦ x+1`) and `P` (`x ↦ 5x = φ²x`; its Prop. 6.6, Cor. 6.5) generate `G_5`. `C` acts as the paper's `Q` (`×φ = ×4`), and `⟨4⟩ = ⟨5⟩ = H`.
  * Which 5-cycle is gold depends on the sign convention of `κ`: replacing `κ` by `−κ` swaps `H` and `−H`. The canonical content is "one 5-cycle ↔ `H`, the other ↔ `−H`".
4. **`G_5` is `BS(1,3)` mod 11.** `G_5 = ⟨x+1, φx⟩ = ⟨x+1, 3x⟩ = {x ↦ ax + b : a ∈ H}`, because `3 = 4^4` in `F_11`.
   * This is the Borel subgroup of `PSL(2,11)`, acting regularly on the 55 unordered pairs. The Paley tournament `QR_11` is its arc-regular orbit.
   * It is also the reduction mod 11 of `BS(1,3) = ⟨x+1, 3x⟩`, the group whose relations THM-4581 shows merge almost surely.
5. **`n = 10` (FINITE-EXACT, golden reader).** The paper's 121-cell cover is `R_10 = O/(11) ≅ F_11²`.
   * The eleven β = 0 cells are the images of the points of period dividing 5; the zero cell also receives −1 and −2.
   * A period-5 point `x` sits in cell `(2κ_5(x), 0)`. So under the paper's projection `(α, β) ↦ α`, it lies over base cell `2κ_5(x)`, which has the opposite colour (`2 ∉ H`). For example, 1/13 (gold) lies over teal cell 6.
   * `C^5` acts as `(α, β) ↦ (α, −β)` on the 110 primitive period-10 rational points. This is PROVED from equivariance, since `4^5 ≡ 1` and `8^5 ≡ −1 (mod 11)`, and checked by applying `C` five times (audit B). So the paper's 55 exchanged pairs are these points up to `C^5`.

## Proof of 1 and 3

* Shifting `w` to `w' = (w_1, …, w_(n−1), w_0)` gives `κ(w') = φκ(w) − w_0(φ^n − 1)`.
* For `n = 5`: `φ^5 − 1 = φ³(4 − φ)`, so `O/(φ^5 − 1) = O/(4 − φ) ≅ F_11`, with `φ ↦ 4` (paper, Proposition 5.2).
* The labels are computed directly (script `golden_cycle_dictionary.py`; re-checked independently in this session).
* Each point is the rational solution of the word's affine fixed-point equation, and its 2-adic parities match its word. ∎

## Remarks

* **Why 11.** `L_5 = 11`. The period-5 golden cyclic words are one constant word and two necklaces of length 5. The −5 cycle exists on `Z` because its word `10100` has `2^3 − 3^2 = −1` (THM-4591). The 1/13 cycle has denominator `2^4 − 3 = 13`.
* **The −17 cycle.** Its standard-map period is 18 (7 odd, 11 even steps). It lives in `R_18 = O/(φ^18 − 1) = O/(76) ≅ O/(4) × O/(19)`, and 19 splits in `O`.
  * The 18/19 meeting is generic Fermat: every prime `p ≡ ±1 (mod 5)` divides `φ^(p−1) − 1`, so `O/(p)` is a quotient of `R_(p−1)` (e.g. `R_10 = O/(11)`).
  * The only specific input is that the −17 cycle has `C`-period 18. NUMEROLOGY (audit B).
* **The golden pair chain (golden reader).**
  * On `Z_2[φ]` with `x ↦ x/2` or `(φx + ρ)/2`, THM-4581's chain table holds verbatim (60,000 steps, 0 mismatches). Every branch contracts, so the Lyapunov weight degenerates to `|f|^θ`, and merging is almost sure (sketch level). The tail is `T^(−1/2)` with `√T q ≈ 3.8`.
  * Conversely, for maps `(mx+1)/2` the balance `½((m/2)^θ s^(−1) + 2^(−θ) s) = 1` has a solution with `0 < s < 1` for some `0 < θ < 1` iff `m < 4`. At `m = 4` the only solution is `s = 2^θ > 1`. The condition `m < 4` is negative Terras drift.
  * So `m = 3` is the only nontrivial odd multiplier where THM-4581's mechanism exists. This is consistent with the conjectured (unproved) divergence of typical 5x+1 orbits.
  * The tail constant `√T q ≈ 3.8` rests on 2000 paths, about ±15%.

**Audit (2026-10-07, independent audit B).**
* Nothing refuted. Re-verified with independent code:
  * equivariance and bijectivity (`n ≤ 22` from the points, `n ≤ 29` from the words);
  * the `n = 5` labels;
  * `G_5` = Stab(∞) in `PSL(2,11)`, regular on the 55 pairs, with Paley arc orbit;
  * `C^5` at `n = 10`.
* Corrected above (MISTAKE-584): notation `C`, the paper's `P` = ×5, the factor-2 cell twist, the `0 < s < 1` qualifier, and the `R_18` retyping.
