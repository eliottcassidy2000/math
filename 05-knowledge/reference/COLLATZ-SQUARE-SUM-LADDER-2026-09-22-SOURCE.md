# User-supplied "Delta=4 ladder / Tao transfer / square-sum horizons" conversation: preserved source

**Status: partly REFUTED, partly restating audited results; source preserved for provenance.**
Do not import any claim below as established. The audits with exact verdicts are the
wave-six lane notes `collatz_mod6_20260922_w6_{summand_closure_square_filter,square_sum_hamiltonicity,tao_minus_sheet_ladder,sign_specific_probes,lean_paste_audit_w6}.md`
under `../results/`, summarized in
[collatz_mod6_20260917_synthesis](../results/collatz_mod6_20260917_synthesis.md).
Earlier attachments by the same author are preserved in the other
`COLLATZ-*-SOURCE.md` files in this directory.

Pasted on 2026-09-22 (condensed by the session lead to plain text; mathematical
content unchanged). The user's own framing preceding the paste: consider the summand
graph (edges x->z and y->z iff x+y=z) built one natural-number vertex at a time, the
"3 distinct independent chains" recalled as {1,4,6}, how this relates to arranging
1..n (n from 15 to 25) so that consecutive pairs sum to squares, how these repeat
fractally in the macrocosm, and out-of-the-box proof strategies for the major open
problems.

---

1. The Delta = 4 shift and the torsion horizon. The prime roots 3, 7, 11, 17 form an
exact arithmetic progression with common difference 4 (3 -> 7 -> 11 -> 17). 3 and 7 are
the first two non-trivial prime targets of the accelerated 3n+1 and 3n-1 maps; 7 and 11
are the endpoints of the audited 7 -> 11 step, which maps onto the primitive
Pythagorean triple (77, 36, 85); 17 is the base root of the hyper-extended loop on the
minus sheet, driven by (2^11 - 3^7)(-17) = (-139)(-17) = 2363. The constant difference
4 is the spatial manifestation of the binary scale matrix; the modulus 4 separates
numbers with a single division by 2 from those with multiple divisions, so this ladder
regulates the flow velocity of the ordered carry B_L.

2. Transferability of Tao's almost-all theorem to the minus sheet. Tao's theorem proves
that Collatz^N(x) decays below any f(x) -> infinity on a set of logarithmic density 1.
Because Tao's proof relies on mapping trajectories as a renewal process on a dyadic
martingale, and the local prefix-descent counts modulo 2^J are identical for 3n+1 and
3n-1, Tao's theorem transfers completely to the minus sheet. Since the three known
disjoint negative loops (-1, -5, -17) exist, an almost-all decay theorem does not imply
that all numbers collapse to a single root; a strictly sign-specific algebraic invariant
is needed showing that the positive sheet has zero capacity for non-trivial Diophantine
cancellations.

3. Square-sum graph horizons and the triple-braid deficit. 15 is the first integer such
that {1..15} can be arranged in a line with every consecutive pair summing to a square.
Building the graph G_sq where X ~ Y iff X+Y = k^2: 1 <= N <= 13 highly disconnected;
N = 14 first unification horizon (3 independent components join); N = 15 first
Hamiltonian path, braiding the square horizons 9, 16, 25; N = 16, 17 connected;
N in [18, 22] parity desert (connectivity stalls); N = 23 the 23 valve restores global
connectivity via the Delta = 4 pivot; N = 24 stalled; N >= 25 asymptotic conjectured
connectivity. The triple-braid deficit: in the N = 15 path the number 4 is the edge to
avoid (a restricted binary vertex 2^2); the squares 9, 16, 25 must interleave their
bitwise weights.

4. Integration into the minimal Lean 4 descent framework. The gaps 18-22 mirror the
capacity horizons of the ordered carries; 23 is the universal valve (mod 256 guard).

```lean
import Mathlib.Data.Nat.Basic
import Mathlib.Combinatorics.SimpleGraph.Basic
def square_sum_graph (n : ℕ) : SimpleGraph (Fin n) := {
  Adj := λ x y => ∃ k : ℕ, (x.val + 1) + (y.val + 1) = k^2,
  symm := sorry,
  loopless := sorry }
theorem square_sum_fourteen_connectivity : (square_sum_graph 14).Connected := by sorry
theorem delta_four_prime_step_invariance :
  ∀ p ∈ ({3, 7, 11, 17} : Set ℕ), (p + 4 ∈ ({7, 11, 15, 21} : Set ℕ)) := by trivial
structure SignSpecificCertificate (n : ℕ) where
  (L : ℕ) (K_L : ℕ) (B_L : ℕ)
  (is_positive_sheet : True)
  (descent_inequality : 3^L * n + B_L < 2^K_L * n)
```

Proposed next: the adjacency matrix of the N = 23 square-sum graph; the exact martingale
decay limits per sheet; octonion multiplication paths across the Delta = 4 ladder.

Second block (restating the wave-five audit and extending it): the Paley 7-tournament has
14 cyclic triples partitioning into two cyclic Fano planes dev{0,1,3} and dev{0,1,5}, the
first carrying the octonion rule e_r e_{r+1} = e_{r+3}; for p = 3 mod 4 the cyclic
triples form a 2-(p,3,(p+1)/4) design; a block matrix with Re(lambda) < 0 for every
trajectory violates Perron-Frobenius, the consistent structure being block triangular
with spectrum x^2 chi_A chi_B; the entropy of the 8/pi^2 odd coprime mask is 0.700 bits
and no lossless code compresses below 1 bit per bit. "Any proof route depending on
sheet-blind invariants is dead": T_+(-n) = -T_-(n), prefix-descent counts mod 2^J agree,
the minus sheet has three cyclic basins and a fourth cycle needs L > 2738, the gate q=1,
cylinders, stopping densities and bounded-strip exclusions are sheet-blind, and a
pairwise contradiction search finds none. Avenues for the sign-specific ingredient: (3.1)
from 2^68 + Legendre + Hardy-Wright 171 any nontrivial positive cycle has more than
14,878,203,146 odd terms; the "196 geometric horizon" 1+3+5+7+11+13+17+19+23+29+31+37
= 196 = 14^2 where positive carries build constructively and negative carries shear
subtractively into the cancellations (2^11-3^7)(-17) = 2363; (3.2) gate-testing the
surviving minus-sheet clocks by mapping trajectories onto the octonion multiplication of
the dev{0,1,3} Fano plane: on the negative sheet the octonion monodromy permits closed
loops, on the positive sheet the sign change breaks the composition balance and the
paths bleed complexity into the trivial 1-loop past the 2^68 boundary.
