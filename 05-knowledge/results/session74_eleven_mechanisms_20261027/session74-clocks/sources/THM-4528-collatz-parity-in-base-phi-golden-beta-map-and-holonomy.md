---
id: THM-4528
title: "Collatz parity words are golden-mean words; their base-phi reading Theta semi-conjugates n -> n/2, 3n+1 to the golden beta-map x -> phi x mod 1 and sends the trivial cycle to {cos 72, cos 60, cos 36 degrees}; Collatz <=> (1 - z^3) F_n(z) is a polynomial for every n >= 1 <=> 7 Phi_T(n) in Z (z = 2, 2-adic) <=> 2 Theta(n) in Z[phi] (z = 1/phi), the last because the golden beta-map has exactly two cycles in (1/2)Z[phi]; no divergent orbit <=> every parity series F_n is D-finite (Szego)."
status: >
  PROVED (elementary) + FINITE-EXACT (the lattice census of (1/2)Z[phi]: 10 points, 2 cycles; all identities checked
  exactly for n <= 20000). INDEPENDENTLY AUDITED (2026-10-01, blind: G1-G3 SOUND; wording corrections applied,
  MISTAKE-556: semi-conjugacy, the exceptional point x = -2, the Haar push-forward is not stationary).
  Notation: T(n) = n/2, 3n+1; b_j(n) = T^j(n) mod 2; F_n(z) = sum_j b_j(n) z^j; Phi_T(n) = F_n(2) in Z_2;
  Theta(n) = sum_j b_j(n) phi^-(j+1).
  (G1) parity words contain no '11'; Theta(T x) = phi Theta(x) mod 1 for all x in Z_2 (as an equality in [0,1) it fails
       only at x = -2, word (01)^inf);
       Theta(1), Theta(4), Theta(2) = phi/2, 1/(2 phi), 1/2.
  (G2) for n >= 1: orbit reaches 1 <=> (1 - z^3) F_n in Z[z] <=> 7 Phi_T(n) in Z <=> 2 Theta(n) in Z[phi].
  (G3) for n >= 1: orbit bounded <=> F_n rational <=> F_n D-finite.
  Non-consequences: none of this proves Collatz; Theta is 2-adically continuous and the conditions are tail conditions.
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), twenty-first note, answering the owner's base-phi prompt ("11 = 100; think of 11 as bit shift with memory")
depends_on:
  - Parry (beta-expansions; admissible golden-mean words), Szego 1922 (power series with finitely many coefficients), Lagarias 1985 (bijectivity of the T1 parity map)
related:
  - 01-canon/theorems/THM-4527-syracuse-in-owner-labels-microcosm-and-difference-spectrum.md (the 2-adic restatement 7 Phi_T(n) in Z)
note: 05-knowledge/results/collatz_golden_holonomy_20261001.md
scripts: 04-computation/experiments/collatz_golden_holonomy_20261001.py
output: 04-computation/experiments/collatz_golden_holonomy_20261001.out (ALL CHECKS PASSED)
---

# THM-4528 — Collatz parity in base φ

**Status: PROVED (elementary); INDEPENDENTLY AUDITED (2026-10-01; wording corrections applied).**
Full statement and proofs: [`05-knowledge/results/collatz_golden_holonomy_20261001.md`](../../05-knowledge/results/collatz_golden_holonomy_20261001.md).

## Statement

Let `T(n) = n/2, 3n+1`, `b_j(n) = T^j(n) mod 2`, `F_n(z) = Σ_j b_j(n) z^j`, and `Θ(n) = Σ_j b_j(n) φ^-(j+1)`.

1. **Golden words.** Parity words never contain `11`, so they are base-φ normal forms. Θ semi-conjugates `T` to the
   golden β-map `x ↦ φx mod 1`. This holds as a congruence mod 1 for all `x`; as an equality it fails only at
   `x = −2`. The trivial cycle goes to
   `{φ/2, 1/(2φ), 1/2} = {cos 36°, cos 72°, cos 60°}`.
2. **One identity, two evaluations.** For `n ≥ 1` the following are equivalent:
   - the orbit of `n` reaches 1;
   - `(1 − z³)F_n(z) ∈ Z[z]`;
   - `7·Φ_T(n) ∈ Z`, at `z = 2` read 2-adically;
   - `2·Θ(n) ∈ Z[φ]`, at `z = 1/φ`.

   The golden converse uses that the β-map has exactly two cycles in `½Z[φ] ∩ [0,1)`: `{0}` and the pentagon
   3-cycle.
3. **Holonomy.** For `n ≥ 1`, the orbit is bounded ⟺ `F_n` is rational ⟺ `F_n` is D-finite. The last equivalence is
   Szegő's theorem: a series with finitely many distinct coefficients is rational or has the unit circle as a natural
   boundary.
