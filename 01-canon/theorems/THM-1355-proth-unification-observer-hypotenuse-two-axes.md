---
id: THM-1355
title: The n·2^x+1 unification — the repo's two master "+1" constants (the observer 2n+1 and the hypotenuse 2^x+1) are the two axes of ONE two-parameter family, the PROTH numbers; the whole family (x≥1) is odd, the Rédei/observer parity signature; its distinguished slices are the Fermat numbers (n=1, x=2^m), the observer/odd numbers (x=1), and the Cullen numbers (x=n).
status: >
  VERIFIED (arithmetic, exact) — f(x,n) = n·2^x+1: the x=1 slice is 2n+1
  (odd numbers), the n=1 slice is 2^x+1, the diagonal x=n is the Cullen numbers
  n·2^n+1, the n=1,x=2^m entries are the Fermat numbers 2^(2^m)+1; f is odd for
  every x≥1; f is exactly the Proth-number form k·2^x+1. INTERPRETIVE (a
  dictionary/organizing identity, honestly labeled — NOT a deep new theorem;
  the Proth form itself is textbook) — the two axes are the repo's two
  independently load-bearing constants: 2n+1 = the observer principle
  (A_k ↔ n=2k+1; Rédei-odd; the LRC worry modulus 2n−1) and 2^x+1 = the
  hypotenuse H=1+2^d (single-flip, THM-250/284), the Fermat rigidity rungs
  (THM-871), the Cayley-Dickson tower rungs n=2^j+1 (THM-448/868), and the
  2-adic vacuum digit of H (THM-466). The odd-parity of the whole family is the
  Rédei signature, and on the Jacobian side the odd-degree conjecture
  (boxeph-S146) is its mirror.
source: death-star-2026-07-20-S59s (HYP-8165; owner's n·2^x+1 table + triangle). Corroborated by three repo-mining sweeps (this session).
depends_on: []
related:
  - THM-466 (2-adic digits of H; the +1 = α₀ = vacuum term)
  - THM-871 (Fermat rungs 3,5,17,257 = rigidity), THM-448/THM-868 (CD/Mersenne tower)
  - THM-250/THM-284 (single-flip H = 1+2^{k-1}); THM-001/THM-002 (Rédei odd; OCF)
  - the observer principle (opus-2026-06-30; jacobian-dixmier-…-S59m §4; dc1-…-S59r)
  - boxeph-S146 (odd-degree conjecture); HYP-6865 (harmonic = staircase resistance)
  - opus-S317 / mac-mini-S109 / klein-S313 (the polygonal/polyhedral figurate triangles)
scripts:
  - 04-computation/triangle_and_n2x1_table_deathstar_S59s.py -> 05-knowledge/results/triangle_and_n2x1_table_deathstar_S59s.out
  - 04-computation/proth_and_three_figurate_triangles_deathstar_S59s.py -> 05-knowledge/results/proth_and_three_figurate_triangles_deathstar_S59s.out
---

# THM-1355 — the Proth table: observer × hypotenuse

## 1. The family (VERIFIED)

f(x, n) = **n·2^x + 1**, with the owner's boundary reads f(0,n) = n, f(x,0) = 1.
This is exactly the **Proth number** form k·2^x+1. Its structure:

| slice | value | repo meaning |
|---|---|---|
| **x = 1** | 2n+1 (odd numbers) | the **observer** axis: A_k ↔ n=2k+1; Rédei-odd; LRC modulus 2n−1 |
| **n = 1** | 2^x+1 | the **hypotenuse** axis: H = 1+2^d (THM-250/284); Fermat rungs (THM-871); CD tower (THM-448/868) |
| **n=1, x=2^m** | 2^(2^m)+1 | the **Fermat numbers** 3, 5, 17, 257, 65537 |
| **x = n** | n·2^n+1 | the **Cullen numbers** 3, 9, 25, 65, 161 |

**Every entry with x ≥ 1 is ODD** (n·2^x is even), so the entire table (off the
x=0 edge) lives in the odd numbers — the Rédei/observer parity signature.

## 2. Why this is a genuine unification, not a coincidence

The two axes were, before this, two *separately* central "+1" constants of the
project, each with its own theorem stack:

- **2^x+1** is the anti-diagonal of the tiling triangle: H = 1 + 2^d, where d is
  the distance of a flipped tile from the hypotenuse (THM-250: one skip-k flip
  gives H = 2^{k-1}+1; apex gives 1+2^{n-2}). Its arithmetic refinements are the
  Fermat rigidity rungs (THM-871: n−1 a 2-power ⟹ Aut = ℤ_n exactly, at
  3,5,17,257) and the Cayley-Dickson tower n = 2^j+1 (THM-448/868). Its 2-adic
  meaning is deepest: THM-466 makes the leading **+1 = α₀ = 1**, the empty
  odd-cycle collection, the *vacuum term* of H = I(Ω,2).
- **2n+1** is the observer axis: the observer principle (A_k has 2k generators;
  its tournament has 2k+1 vertices — "the +1 is the observer"), Rédei's
  H ≡ 1 (mod 2), and the LRC worry modulus 2n−1 = 2N+1 (THM-401; at n=14 it is
  27 = 3³, and the gate law is keyed on the factorization of 2N+1, opus-S123).

The single formula n·2^x+1 makes these the x=1 and n=1 slices of one object.
The observer's +1 and the hypotenuse's +1 are the same +1 — the constant term
that THM-466 identifies as the vacuum, that the observer principle identifies as
the marked vertex, and that the S59p–S59r Jacobian arc tracked as the conserved
unit (u = 1+xy). The table is the two-dimensional home of that +1.

## 3. Odd degree and tournaments (the owner's closing question)

The family is odd for x ≥ 1, and oddness is the tournament parity invariant:
Rédei's theorem says #Hamiltonian paths = 1 + 2·(#odd-cycle collections) is
*always odd* — a parity-protected quantity. The Jacobian-side mirror
(boxeph-S146's **odd-degree conjecture**, corroborated by the verified JC
counterexample's fiber = 1 + 2 = 3, THM-1305/opus-S418): Keller counterexample
cover-degrees are odd, and composition preserves oddness — "Rédei-odd ↔
odd-degree." So the odd numbers 2n+1 that fill this table are exactly the values
a tournament count, or a Keller cover-degree, is *allowed* to take. Odd degree
relates to tournaments as the shared parity floor: the OCF's 2-adic bottom digit
α₀ = 1 is why H is odd, why the escape is positive, and why n·2^x+1 never lands
on an even number.

## 4. Honesty

The Proth form is elementary; §1 is verified arithmetic and §2–§3 are an
organizing dictionary tying it to the theorem stacks above, not a new deep
theorem. The value is consolidation: the project's two most-cited "+1" constants
are one family, and its parity is the OCF vacuum digit. The broader weave — the
owner's triangle as the *third* figurate triangle, the Moser break, the
cubic-Fibonacci, the harmonic staircase-resistance — is in the companion
reflection.
