# Audit D report: THM-4654, the Mahler–Collatz note, and the Chowla exponent ledger

2026-10-14. Independent adversarial audit run by a separate agent. This file condenses its returned report; the scripts and outputs are in this folder. Corrections are recorded as MISTAKE-612 and applied.

## Verdicts

| File | Verdict | Notes |
|---|---|---|
| THM-4654 | CORRECT WITH FIXES | Statements 1–7 are true and their proofs valid; no counterexample to statement 7. Faults are in the consequences and the Reading. |
| Mahler note | CORRECT WITH FIXES | Inherits M1–M3; misses Dubickas 2009. |
| Chowla ledger | CORRECT WITH FIXES | Part I arithmetic right; Part II (iii) wrong (M4); §3.3 overtyped (M5). |

## Major findings

- **M1. The dichotomy (d) was mis-stated.**
  - (b) gives only `exp(−C(log H)^s)` bounds.
  - (c) concludes a finite irrationality exponent, which is the generic case.
  - A transcendence input would exclude nothing.
  - Correct statement: either `limsup N(ℓ)/ℓ = ∞`, or `z_A` has a finite exponent. This is not by itself a restriction.
- **M2. The memory lemma was vacuous; the memory is 2-adic.**
  - Two points of `(1/v)Z` with a common `T`-code prefix of length `ℓ` are congruent mod `2^ℓ` (Terras/Everett), so they are equal or at distance at least `2^ℓ/v`. The lemma's hypothesis forced equality.
  - Real contraction is needed only to enter the window. Inside it, the next `⌈log₂(2.1v)⌉` bits fix the state.
  - Statement 7 stays valid and improves; (c) improves to `μ <= log₂3·(K+1)` in `v`.
- **M3. "(c) covers Sturmian, automatic and morphic codes" was false.**
  - General morphic words can have complexity `Θ(ℓ log log ℓ)`, `Θ(ℓ log ℓ)` or `Θ(ℓ²)` (Pansiot; Devyatov).
  - Sturmian codes are impossible: Dubickas, Glasgow Math. J. 51 (2009), Thm 3, gives `liminf N(ℓ)/ℓ >= log 2/log(3/2) ≈ 1.71` for the Mahler code of every positive integer.
- **M4. The ledger's ρ = 1/20 claim was wrong.** `ρ` is used throughout §§14–16. Lemma 14.6 gives `ρ < 0.124975`; Lemma 16.3 gives `ρ > 8ε`. So `0.0008 < ρ < 0.125`, not "only `ρ > 4·10^(−4)`".
- **M5. Ledger §3.3 was overtyped.** STRUCTURAL should be ANALOGY. The paper's prohibited words come from a divisibility condition and are deleted at a density cost; THM-4654 deletes nothing.

## Minor findings

- **Constants.**
  - `log₂3/c = 3.36026`: use 3.361.
  - `1/c = 2.12009`.
  - The derivable exponent was `3.3603(K+1)`.
- **Statement 6.**
  - The additive constant is optimal: for `w = 0^(p−1)1` the gap tends to `log₂(3^p/(3^p − 2^p))`.
  - The `A_t` gap is 0.562 at `t = 0`.
  - The proof's "M is a 2-adic isometry on each parity class" is false; cite instead that Mahler cylinders are residue classes (THM-2228 §1).
- **Attribution.**
  - Statement 1 is THM-4072's follower graph in slack coordinates.
  - Statement 3's count is THM-3848's renewal law (constant `3/Σ j d_j (2/3)^j = 1.55104519`).
  - Cite Rényi–Parry and Flatto 1992.
- **Notation.** `n₁` was used both for the first one and for the entry time.
- **Ledger.**
  - `A >= 1000` does not come from `D/Y`.
  - The `τ` range has a cap: `τ ∈ (1.000224, 1.18744)`.
  - The 6/5 paragraph: for weight `k^ω` one needs `φ(m) > 2k/(k−1)`.
  - `(1,1,1)` is a witness for every `p`; it is the unique avoider only at `p = 2`.
  - `log log D/log D` is used quantitatively in (12.17).
  - Omitted 20s: Lemma 3.1, Theorem 5.4, **Lemma 6.5's `exp(−L^(21/20))`**, and §15's `t₀/20`.
  - §3.1: "literally" overstates.
  - §3.2 is generic.
  - §3.3: `1 <= k <= 4`.
  - §3.4: cite Tao's Remarks 1.4 and 1.16.

## Independently verified

**Statements 1–4.**
- Four acceptors agree on 524,287 words.
- Exhaustive search over 396,551 safe words of length `<= 28`.
- `N(ℓ)` and `k_max` to `ℓ = 400`.
- `sup T_w(0) = 11/10 − 4.6·10^(−29)`.

**Statements 5–7.**
- Entry window: 214,575 rationals and 1.46 million orbit points.
- Statement 6: `10⁶` runs.
- Statement 7 for all positive `A` at once:
  - 592,561 cases with `v <= 151`;
  - 1.45 million adversarial cases (`|u|` up to `10^15 v`, zero runs, 1,953 periodic safe cycles);
  - 0 violations, minimum slack 4.75 bits.
- The original 4,980,199-pair scan was replicated.

**Ledger.** All numbers recomputed against the PDF (pp. 1–85 read).

## Prior art

- Flatto, "Z-numbers and β-transformations", Contemp. Math. 135 (1992).
- Boyland–de Carvalho–Hall, TAMS 2016 (arXiv:1308.4437): frequency sets.
- Dubickas 2009.
- Andrieu–Eliahou–Vivion, arXiv:2510.11723.
- Eliahou–Verger-Gaugry, arXiv:2504.13716.
- Lagarias's survey (arXiv:2111.02635), for Terras/Everett.
- Adamczewski–Bugeaud; Adamczewski–Cassaigne.
- Pansiot; Devyatov.
- Tao, arXiv:1909.03562.
- No prior art found for statements 4–7 as stated.

## Not checked

- Flatto 1992 in full text.
- Dubickas's theorem verbatim (seen only as quoted by Andrieu–Eliahou–Vivion).
