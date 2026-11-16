---
id: THM-5272
title: "One head over a base. (A) Cut form: H(B + x) = H(B) + 2 sum over a in X, b not in X of W(a,b), W(a,b) = sum of H(B - V(P)) over paths P from a to b with an even number of vertices; second differences are -2(W(e,f) + W(f,e)) <= -2 H(B - e - f). (B) Bank rigidity: 2^k head patterns over a cube of contacts S giving distinct counts inside L consecutive odd numbers force H(B - e - f) <= 2(L - 1) for all e, f in S. (C) Commutator carry: the triadic step T: h -> 3h and the dyadic step D: h -> 2h - 1 satisfy D(T(h)) - T(D(h)) = 2; the rays D^i(3) and T^j(3) meet only at 9. (D) Arithmetic of 27, 63, 1473, 1729 in the signed descent. (E) Finite-exact: no complete 3-bank through base order 9; port spectra of all counts up to 104"
status: >
  PROVED: (A) from Claim A (KNOWN: THM-4876 (1)); (B) elementary from (A); (C) the identity, the uniqueness of
  2^a + 1 = 3^b, and the word-count bound c X^sigma (2^-sigma + 3^-sigma = 1, sigma = 0.787885) for the orbit of
  the semigroup <2n, 3n+1> on half-counts; (D) S^3(63) = 27 for the signed descent S, the dyadic rays
  55, 109, 217, 433, 865, 1729 and 47, 93, 185, 369, 737, 1473, 2945, the Collatz twins 62, 63 (27 + 1024 j and
  63 + 2304 j reach 182 + 6561 j), H(11101001000) = 1729 (literal check).
  FINITE-EXACT: (E) complete 2-banks on 1, 2, 5, 9, 15, 21 bases of orders 4..9 and no complete 3-bank on any
  base of order <= 9 (all isomorphism classes); 336 strictly submodular quadratic bijections of the 3-cube onto 8
  consecutive integers exist abstractly; longest run of consecutive counts on one base 24 (order 8), 29 (order 9),
  pure tape 4, 5, 6, 6, 7, 9, 9 at orders 8..14; port spectra Pi(d) for the 50 counts d <= 104 (both digits -1, +1
  at 26 of them); D(T(h)) and T(D(h)) both physical for h = 3, 5, 9, 15, 17, 19, 23, 25, 27, 31, 33;
  P(63) = {1, 11, 25, 33, 37}; 1473 is not a one-apex count.
  NOT CLAIMED: Chan-Pak (OPEN); a complete low bank at large counts; a selector for arbitrary targets.
  The radix-four physical rollover and its budget bound are physical_radix_carry_20261116, the non-embedding of
  27 in 63 and the square-sum treatment of 8, 9, 1472, 1729 are coalescent_27_63_20261116 (both other sessions).
  Not independently audited (owner directive 2026-11-13: ideas over audits).
session: mac-mini-2026-11-15-carry
date: 2026-11-15
source: 05-knowledge/results/chanpak_carry_banks_ports_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_carry_banks_ports_20261115/scripts/cut_form_check.py, bank.c, ext_runs.c,
    tape_block.c, quad_bank3.py, port_spectra.py, descent_numbers.py
related:
  - HYP-9822, THM-5262, HYP-9812, THM-4876, HYP-9469
  - 05-knowledge/results/physical_radix_carry_20261116/README.md
  - 05-knowledge/results/coalescent_27_63_20261116/README.md
  - 05-knowledge/results/six_tutte_collatz_20261116/COUNTS.md
---

# THM-5272 - one head over a base: cut form, bank rigidity, and the commutator carry

Source: [the note](../../05-knowledge/results/chanpak_carry_banks_ports_20261115/README.md).

`H` counts Hamiltonian paths; `B + x` is a tournament `B` with a new vertex `x` of out-set `X`.

**(A) Cut form.** `H(B + x) = H(B) + 2 Σ_{a ∈ X, b ∉ X} W(a, b)` with
`W(a, b) = Σ_P H(B - V(P))` over the directed paths `P` from `a` to `b` in `B` with an even number of vertices.
*Proof.* Claim A: `H(T) - H(T - v) = 2 Σ_C H(T - V(C))` over odd cycles through `v`; a cycle through `x` is
`x → a ⇝ b → x`. Consequently, for `e, f ∉ X`, the second difference of `X ↦ H(B + x)` in the directions `e, f` is
`-2(W(e,f) + W(f,e))`, and `W(e,f) + W(f,e) >= H(B - e - f) >= 1` because one of the arcs `e → f`, `f → e` is a
path with two vertices.

**(B) Bank rigidity.** Let `S` be a set of `k` vertices outside `X0`. If the `2^k` out-sets `X0 ∪ Y`, `Y ⊆ S`, give
pairwise distinct counts inside a window of `L` consecutive odd numbers, then `H(B - e - f) <= 2(L - 1)` for all
distinct `e, f ∈ S`; for `L = 2^k` (a complete bank) `H(B - e - f) <= 2^(k+1) - 2`.
*Proof.* Four numbers inside the window have second difference at most `4(L - 1)` in absolute value; apply (A).

**(C) Commutator carry.** `T(h) = 3h` (ordinal sum with a triangle) and `D(h) = 2h - 1` give `T(D(h)) = 6h - 3`,
`D(T(h)) = 6h - 1`. On half-counts they are `n → 3n + 1`, `n → 2n`; the orbit of 0 under the semigroup they
generate has at most `c X^σ` elements below `X`. `2^a + 1 = 3^b` only for `(a, b) = (1, 1), (3, 2)`, so
`D^i(3) = T^j(3)` only at `9`, and `(2^a + 1)(2^b + 1) = 2^c + 1` only for `a = b = 1`, `c = 3`.

**(D)** With `S` the signed descent of THM-5262: `S(63) = 47`, `S(47) = 35`, `S(35) = 27`, `S(27) = 21`;
`1729 = 54·2^5 + 1` and `1473 = 46·2^5 + 1` descend to 21 through 55 → 41 and 47 → 35 → 27 respectively.

**(E)** See the status block and `out/` of the note.
