---
id: THM-4720
title: "k-valued tournaments by Burnside: a k-tournament assigns to each unordered pair a value in Z/k with the swap acting by x -> 1 - x (k = 2: tournaments; k = 3: oriented graphs; even k: (k/2)-arc-coloured tournaments; odd k: ((k-1)/2)-arc-coloured semicomplete digraphs with the unique tie state 1/2 mod k); the number of isomorphism classes is T_k(n) = sum over cycle types lambda of k^(o_+(lambda)) c^(o_-(lambda)) / z_lambda with o_+ = sum floor((L_i - 1)/2) + sum_(i<j) gcd(L_i, L_j), o_- = the number of even cycles and c = [k odd]; hence for even k only permutations with all cycles odd contribute (the tournament shape) while for odd k all do (the graph shape), and the rigidity defect T_k(n) n!/k^C(n,2) - 1 is 2 C(n,3) k^(-(2n-4)) (even k) versus C(n,2) k^(-(n-1)) (odd k) to leading order; the sequences for k = 4..8 are new (not in OEIS), k = 3 is A001174; the arity-3 analogue, 3-hypertournaments (a cyclic orientation on every triple), has 1, 1, 1, 3, 12, 1512, 6818336, 1787143624656, ... classes (new), and 3-coloured complete graphs are A004102"
status: >
  PROVED: the Burnside formula and the parity dichotomy (elementary orbit count on pairs), the identification of the states,
  the leading rigidity terms. FINITE-EXACT: T_k(n) for k = 2..8, n <= 14 (k_tournaments.out; k = 2, 3 agree with A000568,
  A001174); 3-hypertournaments to n = 10 by Burnside with the sign of the induced permutation on each triple orbit,
  validated by brute force at n = 4, 5; coloured graphs to n = 10. NUMERICAL/NULL: the Hamiltonian-path colour polynomial
  of 2-coloured tournaments has no structure mod 2 or mod 4 beyond Redei's parity (redei_colour_polynomial.out). OEIS
  searches on 2026-10-19 found none of the k = 4..8 or hypertournament sequences; a literature search found Klarner 1968
  (k-coloured graphs) and the 2023 'tournaments and even graphs are equinumerous' theorem but no coloured-tournament
  counts, so novelty is plausible but not asserted. Pre-audit.
session: mac-mini-2026-10-19-twothree
source: 05-knowledge/results/two_three_axis_20261019_numbers_fixed_points_relations.md
scripts:
  - 04-computation/experiments/two_three_20261019/k_tournaments.py (+ .out, .json), hyper_and_coloured.py (+ .out, .json), validate_hyper.py (+ .out), redei_colour_polynomial.py (+ .out)
related:
  - THM-001 (Redei), THM-1980 (the 2-adic edge of the Redei count), THM-4533 (Redei graphs), THM-283 (SC Burnside formula), THM-214 (Paley Hamiltonian cycles by Burnside)
  - the S48 Burnside/Polya enumerator (burnside_enum_v2.c: tournaments, oriented graphs, digraphs, hypergraphs)
  - HYP-9339 (the 2-3 axis synthesis)
---

# THM-4720 — k-valued tournaments by Burnside

## Setting

* A **k-tournament** on `[n]` is a map `f` from ordered pairs `(i, j)`, `i ≠ j`, to `Z/k` with `f(i, j) + f(j, i) = 1`. Equivalently each unordered pair carries one of `k` states, and reversing the pair applies the involution `τ(x) = 1 − x`. The `τ`-orbits of `Z/k` are the pairs `{x, 1 − x}` (an **arc colour** with its two orientations) and, when `k` is odd, the fixed point `x = 2^(−1) mod k` (the **tie**). So a `k`-tournament is a `(k/2)`-arc-coloured tournament for even `k`, and a `((k−1)/2)`-arc-coloured semicomplete digraph (ties = symmetric pairs) for odd `k`; `k = 2` gives tournaments, `k = 3` oriented graphs (tie = non-edge).
* `S_n` acts by relabelling; `T_k(n)` is the number of isomorphism classes. For a permutation of cycle type `λ = (L_1, …, L_r)`, an orbit of unordered pairs is **flipping** if the pair returns with its orientation reversed, **preserving** otherwise.

## Statements

1. **Orbit count (PROVED).** A cycle of length `L` contributes `⌊(L−1)/2⌋` preserving orbits of pairs inside it and, if `L` is even, one flipping orbit (the antipodal pairs `{i, i + L/2}`); two cycles of lengths `L_i, L_j` contribute `gcd(L_i, L_j)` preserving orbits. Hence `o_+(λ) = Σ_i ⌊(L_i − 1)/2⌋ + Σ_(i<j) gcd(L_i, L_j)` and `o_−(λ) = #{i : L_i even}`.

2. **Burnside formula and the parity dichotomy (PROVED).**

       T_k(n) = Σ_(λ ⊢ n) k^(o_+(λ)) c^(o_−(λ)) / z_λ,     c = 1 if k is odd, c = 0 if k is even,

   because a flipping orbit must carry a `τ`-fixed state, which exists iff `k` is odd (`2x = 1` solvable mod `k`). So for even `k` only permutations with all cycles odd contribute — the shape of the tournament count `A000568` — while for odd `k` every permutation contributes — the shape of the graph counts. The number 2 decides the shape through the solvability of `2x = 1`.

3. **Rigidity defect (PROVED to leading order).** `T_k(n) · n!/k^(C(n,2)) = 1 + 2 C(n,3) k^(−(2n−4)) + …` for even `k` (the first contributing non-identity permutations are the 3-cycles) and `= 1 + C(n,2) k^(−(n−1)) + …` for odd `k` (transpositions contribute). Checked: at `k = 4, n = 12` the defect is `4.002·10^(−10)` against the predicted `4.002·10^(−10)`; at `k = 3, n = 12` it is `3.731·10^(−4)` against `3.726·10^(−4)`.

4. **Sequences (FINITE-EXACT).** `T_k(n)`, `n = 1, …, 9`:
   * `k = 2`: 1, 1, 2, 4, 12, 56, 456, 6880, 191536 (A000568);
   * `k = 3`: 1, 2, 7, 42, 582, 21480, 2142288, 575016219, 415939243032 (A001174);
   * `k = 4`: 1, 2, 12, 176, 8784, 1492288, 872687552, 1787154671104, 13013584213369088 (2-arc-coloured tournaments; not in OEIS);
   * `k = 5`: 1, 3, 25, 695, 82880, 42598925, 94739454950, 924263496657375, 40104876653759021125 (not in OEIS);
   * `k = 6`: 1, 3, 38, 1956, 504108, 653050728, 4352574626712, 152305122284204256, 28423789739316734337072;
   * `k = 7`: 1, 4, 63, 5012, 2364495, 6599821760, 110842443237533, 11408784389272827144, 7307504981018314650408502;
   * `k = 8`: 1, 4, 88, 10944, 8948544, 48867301888, 1830034253721088, 479732468882948644864, 894286137761306641148252160;
   to `n = 14` in `k_tournaments.json`.

5. **Arity three (FINITE-EXACT).** The number of **3-hypertournaments** (one of the two cyclic orientations on every 3-subset) up to isomorphism is `Σ_λ 2^(h(λ))/z_λ`, where `h(λ)` is the number of triple-orbits of a permutation of type `λ` on which the return permutation is even, and the term is 0 if some orbit returns with an odd permutation: `1, 1, 1, 3, 12, 1512, 6818336, 1787143624656, 53303607568835303424, 366299602013040506470327589888` for `n = 1, …, 10` (brute force confirms 3 and 12 at `n = 4, 5`; not in OEIS). The 3-coloured complete graphs (pairs, no twist) are A004102: 1, 3, 10, 66, 792, 25506, …; 4-coloured: 1, 4, 20, 276, 10688, 1601952, ….

6. **Redei in colour (NULL, FINITE-EXACT).** For 2-coloured tournaments the Hamiltonian-path colour polynomial `P_T(x) = Σ_paths x^(#red arcs)` has `P_T(1)` odd (Redei) and an odd number of odd coefficients, but its reduction mod 2 shows no further law (87 distinct patterns among 300 random 8-vertex tournaments), and `P_T(−1) mod 4` is 1 or 3 about equally often.

## Proofs

**1.** Inside a cycle `(i, i+1, …, i+L−1)` the pair `{i, i+d}` moves to `{i+1, i+1+d}`; its orbit has length `L` for `d ≠ L/2` with the orientation preserved, and length `L/2` for `d = L/2` (`L` even) with the orientation reversed after `L/2` steps; the distinct `d` with `1 ≤ d ≤ L/2` give the count. For two cycles, the pairs `(i, j)` with `i` in the first and `j` in the second form `gcd(L_i, L_j)` orbits of length `lcm`, orientation preserved. ∎

**2.** Burnside: `T_k(n) = (1/n!) Σ_σ |Fix(σ)|`; `σ` fixes `f` iff `f` is constant on preserving orbits (any of the `k` values) and takes a `τ`-fixed value on flipping orbits. ∎

**3.** For even `k`, the identity contributes `k^(C(n,2))/n!`, and the next contributions come from 3-cycles: `2 C(n,3)` of them, each with `C(n,2) − o_+ = 2n − 4`. For odd `k`, transpositions contribute `C(n,2)` terms with `C(n,2) − o_+ = n − 1`. ∎

**4–5.** Direct computation (`k_tournaments.py`, `hyper_and_coloured.py`); for triples the return permutation of each orbit is computed explicitly and its parity read off. ∎

## Remarks

* **Binary versus ternary information.** Oriented graphs are the `Z/3`-tournaments: the third state is the tie `1/2 mod 3 = 2`. The arity-3 version (information on triples, two states) is the hypertournament; both are "numbers greater than 2 instead of binary", in the value and in the arity direction.
* **Relation to even graphs.** Tournaments and even graphs are equinumerous (a 2023 theorem by Cauchy–Frobenius counting); the `k`-tournament counts for even `k` are the natural candidates for a coloured version of that identity.
