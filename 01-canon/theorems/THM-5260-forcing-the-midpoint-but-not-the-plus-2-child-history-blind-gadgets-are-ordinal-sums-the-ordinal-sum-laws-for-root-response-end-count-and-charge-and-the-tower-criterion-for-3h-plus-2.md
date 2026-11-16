---
id: THM-5260
title: "(A) History-blind extensions are ordinal sums: substituting an arbitrary tournament T for a vertex tau of a gadget Q gives a count H(Q[T]) that is a function of H(T) alone iff tau lies on no directed cycle of Q, and then H(Q[T]) = H(Q_1) H(T) H(Q_2) (Q = Q_1 -> tau -> Q_2); in particular no uniform gadget sends H to 3H + 2 (the midpoint 3H is forced by C_3 -> T, the +2 child must read the realiser). (B) Ordinal-sum laws for the root response B_t, the end count e_t and the charge a_t = B_t - e_t - h + 1: for U = A -> B and t in B, B_t(U) = h(A) B_t(B), e_t(U) = h(A) e_t(B), a_t(U) = h(A)(a_t(B) - 1) + 1; for t in A, B_t(U) = B_t(A) P_B(2), e_t(U) = 0, a_t(U) = B_t(A) P_B(2) - h(A) h(B) + 1. Hence a transitive prefix changes nothing and a transitive suffix of length k sends B_t to 2^k B_t and e_t to 0, so the digit becomes 2^k B_t - h and the charge 2^k B_t - h + 1. (C) Tower criterion: for a count h all of whose realisers are paddings of strong realisers S (every prime; every h < 125 with one nontrivial strong component), 3h + 2 is the count of a height-two tower over some realiser of h at some order iff some strong realiser has a vertex with 2^k B_t(S) = h + 1 for some k >= 0; digit 2 needs B_t = h + 2 exactly, while digit 0 is always available at some order because a padded sink has B_s = h (MISTAKE-776 corrects the first statement); a rooted doubling child 2h - 1 exists at some order iff some strong realiser has a_t = 0 or 2^k B_t = h - 1 with k >= 1. Also B_t - e_t = 2 E_t (the two charge formulas agree)."
status: >
  PROVED (klein-2026-11-16-force; (A) forward direction: tau a singleton strong component makes Q[T] an ordinal sum and H
  multiplies (THM-1370); converse: FINITE-EXACT for all 110 pointed gadgets of order <= 4 (scripts/module_blind.py), with
  the structural reason that every cycle of Q[T] meets T in one run so H(Q[T]) is a positive combination of run sums
  that equal-count tournaments do not share; (B) elementary, 286 random checks (scripts/padding_check.py); (C) from (B)
  and THM-4865's tower h + 2B_t; the identity B_t - e_t = 2E_t checked on every rooted tournament of orders 4-9
  (scripts/digits.c); not independently audited). CITED: THM-002, THM-1370, THM-4865, THM-4876, wave 45 of
  chanpak_waves_20261113 (charges), THM-5272 (doubling failures), wave 50 (sterile forks).
session: klein-2026-11-16-force
date: 2026-11-16
source: 05-knowledge/results/forcing_digits_analogy_additive_20261116/FORCING.md
scripts:
  - 05-knowledge/results/forcing_digits_analogy_additive_20261116/scripts/module_blind.py (+ receipts/module_blind.out)
  - 05-knowledge/results/forcing_digits_analogy_additive_20261116/scripts/padding_check.py (+ receipts/padding_check.out)
  - 05-knowledge/results/forcing_digits_analogy_additive_20261116/scripts/strongB.c, digits.c (+ receipts/strongB_*.out, digits_*.out)
related:
  - 01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 05-knowledge/results/forcing_digits_analogy_additive_20261116/ANALOGY.md
---
# THM-5260 — forcing: the midpoint yes, the +2 only by reading the realiser

**(A)** If `τ` is on no cycle of `Q`, it is a singleton strong component, `Q = Q_1 → τ → Q_2`, and `Q[T] = Q_1 → T →
Q_2` with `H = H(Q_1)H(T)H(Q_2)`. If `τ` is on a cycle `C`, the odd-cycle formula expands `H(Q[T])` over packings of
`Q[T]`; a cycle of `Q[T]` enters `T` from a vertex that beats all of `T` and leaves to one that loses to all of `T`, so
it meets `T` in a single run, a directed path `P` of `T`; the mixed cycles through `C` contribute
`2 Σ_P H(…)` over paths `P` of a fixed parity, which differ between the 3-cycle and the 3-cycle-with-sink (both `H = 3`).
Exhaustive at gadget order ≤ 4. ∎

**(B)** A path starting at `t ∈ B` never enters `A` (no arc from `B` to `A`), so in `B_t(U)` the complement of `A`
must contain all of `A`, giving `h(A)·B_t(B)`; a Hamiltonian path of `U` is a path of `A` followed by one of `B`, so
`e_t(U) = h(A) e_t(B)` for `t ∈ B` and `0` for `t ∈ A`; a path starting at `t ∈ A` is its `A`-part followed by any
Hamiltonian path of the rest of `B`, giving `B_t(A)·P_B(2)`. For `B = TT_k`: `h = 1`, `P_B(2) = 2^k`. ∎

**(C)** The digit of an inherited vertex `t ∈ S` after padding is `2^k B_t − h` and the tower gives `3h + 2σ`; the charge is
`2^k B_t − h + 1` (`k ≥ 1`) or `a_t` (`k = 0`). The padding vertices themselves have digits `h(2^j − 1)` (the `j`-th vertex from
the end of the suffix; `0` at the sink) and at least `h` (a prefix source), and charges `h(2^j − 1) + 1`, `1 − h` (sink) or at
least `h + 1` — never 1, 2 or charge 0. ∎ (Erratum MISTAKE-776: the digit-0 clause originally omitted the sink.)
