---
id: THM-4764
title: "Modular cancellation depth of odd regular gap actions, and the fifty-vertex realization: for a regular action of a p-group of odd order p^e the subset-orbit polynomial P_A(t) reduced mod p vanishes at t=-1 to order p^e - p^J + p^(J-1), J the largest j with p not dividing the number c_j of cyclic subgroups of order p^j (cyclic: p^(e-1); abelian non-cyclic: p^e - p + 1); for a prime l not dividing |A| the order is min{|A|/o : l does not divide e_o}; in general the order mod l is at least d = min cycle count with equality iff l does not divide the leading coefficient kappa, and a candidate-order lower bound holds for every regular action. The tournament Cay(Z5^2,S)[TT2] on 50 vertices has Aut = Z5^2 regular on its 25 gaps, one invariant toggle, depth 5 over Q (kappa = 3000), 21 mod 5, 25 mod 2 and mod 3; its cyclic twin Circ(Z25,{1..12})[TT2] has depths 1, 5, 25, 1."
status: >
  PROVED (klein, 2026-10-25; elementary, complete proofs in the source packet). Theorem B (p-groups), B' (primes
  prime to |A|), B'' (general candidate-order lower bound) and Lemma B0 are proved for regular actions of odd
  order. The corollary J = 1 for NON-ABELIAN non-cyclic p-groups is CITED (classical divisibility of cyclic-subgroup
  counts; Miller / P. Hall 1933 / Berkovich I §1, secondary source) and verified on seven non-abelian groups.
  Equality in Theorem B'' for mixed order is OPEN (HYP-9374), verified on 22 regular actions.
  FINITE-EXACT: both 50-vertex tournaments (automorphisms enumerated from arcs with no group assumed; 2^25
  orbit census in C agrees with Burnside). Not yet independently audited by another agent.
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/clone_depth_modular_twist_average_20261025/README.md
scripts:
  - 05-knowledge/results/clone_depth_modular_twist_average_20261025/check_clone_depth.py (+ clone-normal.json / clone-optimized.json)
  - 05-knowledge/results/clone_depth_modular_twist_average_20261025/census_orbits.c (+ census_z5sq.out, census_z25.out)
related:
  - 05-knowledge/results/session62_agent_theme_bridges_20261025/clone/README.md (root order over Q, toggle rank, 24-vertex C3^2 realization)
  - 05-knowledge/results/ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md (clone fibre = cut orbits)
  - 01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md (module paths)
  - 01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md (2-adic depth of A000568)
  - 05-knowledge/hypotheses/HYP-9374-mixed-order-clone-depth-equals-the-candidate-minimum.md
  - 01-canon/theorems/THM-4765-quadratic-twists-of-a-square-root-cover-average-to-the-branch-complement.md
---

# THM-4764 — Modular cancellation depth of odd regular gap actions

## Setting

`A` finite of odd order acting on a finite set `G`; `P_A(t) = sum_{[C] ∈ 2^G/A} t^|C| ∈ Z[t]`;
`s = 1+t`; `d = ord_{t=-1} P_A = min_g c(g)` and `κ = [x^d] P_A(-1+x) > 0` (Session62). For a prime `ℓ`
the depth mod `ℓ` is `d_ℓ = ord_{s=0}(P_A mod ℓ)`. In the clone fibre of a tournament `T` (CLONE_INVOLUTION,
THM-4463 §2) `A = Aut(T)` acts on the canonical module gaps and `P_A` is the generating polynomial of the
weighted-quotient isomorphism classes inflating to `T`; the signed binomial moments
`sum_[C] (-1)^(n-|C|) C(|C|,j)` are therefore divisible by `ℓ` for `j < d_ℓ`.

## Statements

1. **Lemma B0.** `d_ℓ ≥ d`, with equality iff `ℓ ∤ κ`.
2. **Theorem B′.** Regular action, `ℓ ∤ |A|`: `d_ℓ = min{ |A|/o : ℓ ∤ e_o }`, `e_o` = number of elements of order `o`.
3. **Theorem B.** Regular action of a `p`-group of odd order `p^e`, `c_j` = number of cyclic subgroups of order `p^j`:
   `c_1 ≡ 1 (mod p)` and `d_p = p^e - p^J + p^(J-1)`, `J = max{ j ≥ 1 : p ∤ c_j }`. Cyclic: `d_p = p^(e-1)`.
   Abelian non-cyclic: `d_p = p^e - p + 1`. Non-abelian non-cyclic: `d_p = p^e - p + 1` given the cited count.
4. **Theorem B″.** Regular action, `|A| = p^e m`, `p ∤ m`: with `E_{o'}` the number of elements whose `p'`-part has
   order `o'` and `e_{j,o'}` the number of elements of order `p^j o'`, one has `p^e | E_{o'}`, `p^(j-1) | e_{j,o'}` and
   `d_p ≥ min({|A|/o' : p ∤ E_{o'}/p^e} ∪ {|A|/o' - p^j + p^(j-1) : p ∤ e_{j,o'}/p^(j-1)})`.
5. **Realization.** `T = Cay(Z5^2, S)[TT2]`, `S = {(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(3,1)}`:
   50 vertices, `|Aut| = 25`, homogeneous pairs = the 25 clone blocks, gap action = regular `Z5^2`, toggle rank 1,
   `P = ((1+t)^25 + 24(1+t^5)^5)/25`, `|Ω| = 1,342,208`, `(d_Q, κ) = (5, 3000)`, `(d_2, d_3, d_5, d_7) = (25, 25, 21, 5)`.
   The cyclic twin `Circ(Z25,{1..12})[TT2]` has the same order, `|Aut|`, gap count and toggle rank, with
   `P = (1/25) sum_{d|25} φ(d)(1+t^d)^{25/d}`, `|Ω| = 1,342,184`, `(d_Q, κ) = (1, 20)`, `(d_2, d_3, d_5) = (25, 1, 5)`.

## Proofs

Complete proofs are in the source packet, §2. The mechanism of Theorem B: `1 + t^{p^j} = s^{p^j} + p·s·g_j(s)`
with `g_j ≡ u_j s^{p^(j-1)-1} + …` mod `p` (Kummer); in the binomial expansion of `(s^{p^j} + p s g_j)^{p^{e-j}}`
divided by `p^e` only `m = 0, 1` survive mod `p`, the `m = 0` terms sum to `s^{p^e}`, and the `m = 1` term of
index `j` is `c_j (p-1) s^{p^e-p^j+1} g_j`, of exact order `p^e - p^j + p^(j-1)` when `p ∤ c_j`; distinct `j`
have distinct orders. Theorem B″ repeats this with `s' = 1 + t^{o'}` in place of `s` and uses Frobenius'
theorem for `p^e | E_{o'}`. Theorem B′ and Lemma B0 are one-line valuation arguments.

## Mechanism, boundary, and failure anatomy

- Over `Q` the depth is a **minimum** (one element), the toggle rank a **whole-group** invariant (Session62);
  modulo `p` the depth is a third thing: it reads the `p`-adic size of the class multiplicities and the
  cyclic-subgroup counts. The three numbers `(d_Q, d_p, r)` distinguish `Z5^2` from `Z25` on identical
  `(|V|, |Aut|, |G|, r)`.
- First failed implication for a naive extension: in `C21` modulo 3 the order-21 and order-7 contributions to
  `[s^3]` are `760` and `98`, which cancel (`858 ≡ 0`), so the `p`-group proof step "distinct orders do not
  cancel" fails in mixed order; the candidate-order bound repairs it and is sharp in every tested case.
- Hostiles: even actions (transposition `1+t+t^2`, no zero) and complex scalar actions (`1+t^3`) are outside
  scope; the odd-order and real/permutation hypotheses of Session62 are retained.
- `p^2 - p + 1 = Φ_6(p)` for `e = 2` is a numeral coincidence; `e ≥ 3` gives `p^e - p + 1`.

## Non-consequences

No faster enumeration of A000568, no black-path transport, no statement about tournaments whose
automorphism group acts non-regularly on gaps beyond Lemma B0, and no geometric or LRC transfer.
