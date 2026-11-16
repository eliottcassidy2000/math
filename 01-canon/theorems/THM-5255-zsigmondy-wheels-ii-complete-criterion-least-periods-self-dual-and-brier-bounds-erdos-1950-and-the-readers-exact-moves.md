---
id: THM-5255
title: "Zsigmondy wheels II (continuing THM-5254; W(N) = the primes of 2^N - 1, one class modulo ord_p(2) each; N = 2^a m, m odd; the odd wheel = the odd parts of the non-dyadic orders). (A) PROVED, complete criterion: W(N) covers Z/N iff a >= 6 (the dyadic classes cover Z/2^a, by Kraft on the binary tree with nu(2^i) = 1 for i <= 5 and nu(64) = 2) or the odd wheel covers Z/m; the fibre argument on Z/2^a x Z/m shows no non-nested arrangement ever helps, settling HYP-9823 (W1) in its corrected form. (B) PROVED: for a prime-power odd part q^e and a <= 5 the criterion is the Kraft density sum of q^(-j) over the odd wheel >= 1 (36 covers because 2/3 + 4/9 >= 1; 18, 54 fail). (C) PROVED/FINITE-EXACT: a(q) := min{a : W(2^a q) covers} satisfies a(q) <= 6 for every odd prime q, with a(3) = 3, a(5) = a(7) = 4, a(11) = a(13) = 5, a(17) = a(19) = a(23) = 6 (periods 24, 80, 112, 352, 416, 1088, 1216, 1472), from the exact factorisations of Phi_(2^i q)(2) (for q = 23 the factordb factorisation of Phi_736(2), verified locally); conjecturally a(q) = 6 for all q >= 17. (D) PROVED: for N = 2m a wheel is self-dual (serves both signs of one coefficient) iff its flipper wheel {m' : a prime has order exactly 2m'} covers Z/m, the plus sheet being the translate of the minus sheet by (m+1)/2; for any covering set containing 3 that serves both signs, sum over keepers (4 | ord) of 1/ord + 2 sum over flippers (ord = 2 mod 4, ord > 2) of 1/ord + (1/2) sum over odd orders of 1/ord >= 1 (parity-sheet / Brier bound); the bound is below 1 for every period <= 150, so no wheel of period <= 150 is self-dual, and W(24), W(36), W(64) cannot host Brier numbers (0.667, 0.722, 0.5) while the period-720 Filaseta-Finch-Kozek set scores 1.476. (E) Erdos's 1950 covering system {0 (2), 0 (3), 1 (4), 3 (8), 7 (12), 23 (24)} is class for class the first covering pattern of the Riesel wheel W(24). (F) FINITE-EXACT, synchronized binary reader (I, Z, g) = (-11, 82, 9): for n <= 22 about 2.2 % of the 2^n addresses have another address with exactly 2 more Hamiltonian paths, single-bit such partners occur only at the critical bit (n+2)/2 and realise +2, +4 or +8 according to the width, 0.13-0.34 % of addresses have a +2 partner keeping both exterior phases, and every phase cell reaches every residue modulo 9; a factor-1 shield (common block keeping the unit step) is a transitive block, hence a zero run multiplying the slope by a power of two, and cannot supply a new coordinate."
status: >
  PROVED (klein-2026-11-16-wheels2; elementary, not independently audited): (A), (B), (D), (E) with full proofs in the packet's
  WHEELS2.md; (C) PROVED for the bound a(q) <= 6 and for each listed value by exact factorisation (sympy; receipts nu_*.out),
  CONJECTURE beyond. FINITE-EXACT (F): exhaustive enumeration of all addresses for n = 12..22 (address_moves.py); the shield
  statement cites BLOCKS.md of variable_reader_coverage_20261116 (h_B = 1 forces q_B = 2^s). Checks: wheel_theorems_check.py
  compares (A), (B) and (D) with the exact search at every N <= 150 (no mismatch); classical_covers.py verifies (E) and the
  bound values. Pre-registration in the packet's PREREGISTERED.md (P-A1 split, P-A2 confirmed, P-F1 confirmed for 19 and 23; P-S1, P-S3 confirmed and P-S2 refuted for spine 14, see SPINE14.md).
  CITED: Zsigmondy/Bang, Erdos 1950, Sierpinski 1960, Selfridge, Riesel, Filaseta-Finch-Kozek (Brier numbers), THM-5254,
  HYP-9823, HYP-9842 (no third digit in the seed monoid), variable_reader_coverage_20261116 (PROOF.md, BLOCKS.md),
  live_exterior_reader_20261116 (EXTERIORS.md, SWITCH.md, CONTACT.md).
session: klein-2026-11-16-wheels2
source: 05-knowledge/results/wheels_spine14_shield_20261116/README.md
related:
  - 01-canon/theorems/THM-5254-zsigmondy-wheels-covering-periods-for-k-2-n-plus-minus-1-reduce-to-the-odd-part-least-periods-24-80-112-the-classical-covering-sets-are-full-wheels-and-the-riesel-primes-meet-24-only-through-von-staudt-clausen.md
  - 05-knowledge/hypotheses/HYP-9823-zsigmondy-wheel-conjectures-nested-dyadic-is-without-loss-least-periods-grow-like-sqrt-q-no-self-dual-wheel-and-the-base-b-analogue.md
  - 05-knowledge/results/logarithmic_repairs_covers_20261116/CLOCKS.md
  - 05-knowledge/results/live_exterior_reader_20261116/EXTERIORS.md
  - 05-knowledge/results/variable_reader_coverage_20261116/BLOCKS.md
---
# THM-5255 — Zsigmondy wheels II, and the synchronized reader's exact moves

**(A)** `Z/N ≅ Z/2^a × Z/m`; a class modulo `2^i m′` is `A × B` with `A` a class modulo `2^i` and `B` a class modulo `m′`;
for a dyadic prime `B = Z/m`. If some `u ∈ Z/2^a` lies in no dyadic `A`, its fibre must be covered by the `B`'s of a
sub-multiset of the odd wheel, impossible when the whole odd wheel cannot cover `Z/m`. So a cover forces either a dyadic
cover of `Z/2^a` or an odd-wheel cover of `Z/m`; both suffice (nested dyadic classes leave one `u`, over which every
non-dyadic class is placed with its prescribed `B`). Classes modulo `2^i` are subtrees of the binary tree, so they cover
iff `Σ 2^(−i_p) ≥ 1`, i.e. iff `a ≥ 6`. ∎

**(B)** Classes modulo `q^j` in `Z/q^e` are subtrees of the `q`-ary tree; a multiset of depths covers iff `Σ q^(−j) ≥ 1`
(nondecreasing depth placement; measure). Combine with (A). ∎

**(C)** `W(64q) ⊇ W(64)` covers; for `a ≤ 5`, (A) with `m = q` gives the criterion `Σ_{i≤a} ν(2^i q) ≥ q` (THM-5254 (C)).
Factorisations: `c(a, q)` for `a = 0..5` equals `1, 1, 2, 3, —` (`q = 3`), `1..5` (`q = 5`), `1, 2, 4, 5, 7` (`q = 7`),
`2, 3, 5, 7, 10, 12` (`q = 11`), `1, 2, 5, 7, 8, 13` (`q = 13`), `1, 2, 5, 7, 10, 13` (`q = 17`), `1, 2, 5, 8, 11, 17` (`q = 19`), `2, 3, 7, 8, 11, 14` (`q = 23`; `ν(736) = 3`, receipt `nu_23_factordb_verify.out`).

**(D)** Let 3 cover the even residues on the minus sheet, hence the odd residues on the plus sheet. For `N = 2m`, a prime of
order `2m′` on an odd class `r` covers `t ≡ (r − 1)/2 (mod m′)` among the odd residues `x = 2t + 1` on the minus sheet and
`t ≡ (r + m′)/2 ≡ (r − 1)/2 + (m′ + 1)/2 (mod m′)` among the even residues `x = 2t` on the plus sheet; `(m′ + 1)/2` is the
restriction of `(m + 1)/2 ∈ Z/m`, a common translation, so one sheet covers iff the other does. For the bound: a keeper's
class and antipode share a parity and help one sheet by `N/d`; a flipper on an odd class helps both sheets by `N/d`; an
odd-order class helps its one sheet by `N/(2d)`; the two sheets need `N/2` each. ∎

**(E)** Direct verification (receipt `classical_covers`); the pattern `(0, 1, 0, 7, 3, 23)` of CLOCKS.md in the modulus
order `2, 4, 3, 12, 8, 24` is Erdős's system.

**(F)** Packet SHIELD.md §2 (census) and §1 (shields, after BLOCKS.md).

**What is new.** (A)–(D), the least-period values, (E), the census (F). Known and used: THM-5254, Zsigmondy–Bang, the
classical covering sets, BLOCKS.md, CONTACT.md.
