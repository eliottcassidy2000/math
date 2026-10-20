---
id: THM-4726
title: "HYP-9172 (b), what a pairing proof of the level-branch integrality must carry: for every n, THM-479's level branch is N_lev(n) = sum over S_n-orbits of pairs (C, x), C a switching class of tournaments and x a nontrivial 2-element of Aut C, of rho(C_{Aut C}(x)), the odd-order proportion of the centraliser (the 2-part/odd-part factorisation g = g_2 g_{2'} of even-order automorphisms); when Aut C has a normal 2-complement this is sum over the nontrivial elements x of a Sylow 2-subgroup P of 1/|P|, so N_lev = sum_C (|P_C| - 1)/|P_C| and Theorem B of THM-4531 (n = 2 mod 4, |P| = 2) is the case where every term is 1/2 and the decorated flip (C, g) -> (C + M_g, g) pairs the terms. For n = 0 (mod 4) the decorated flip on pairs (C, g), g an involution class, has fixed points (2 of 2 at n = 4; 9 of 31 at n = 8) and the weights 1/|Syl_2 C(x)| take the values 1/2, 1/4, 1/8 (19, 14, 16 pairs at n = 8; 19/2 + 14/4 + 16/8 = 15 = N_lev(8)), so any pairing must mix Sylow types and weights; the Sylow-type strata contribute 19/2 (C_2), 3/2 (C_4), 1 (V_4: Aut = V_4 with 3/4 and Aut = A_4 with 1/4), 7/4 (C_8), 5/4 (D_8: Aut = D_8 with 7/8 and Aut = PSL(2,7) with 3/8)"
status: >
  PROVED (pre-audit): the decomposition N_lev = sum_{(C,[x])} rho(C_{Aut C}(x)) (elementary: g -> (g_2, g_{2'}) is a
  bijection from even-order elements to commuting pairs (nontrivial 2-element, odd-order element)); the normal-2-complement
  specialisation; the statement that Theorem B is the |P| = 2 case. FINITE-EXACT: the n = 4, 6, 8 census of even-Aut
  classes with |Aut|, Sylow-2 order, numbers of involution classes and of 2-element classes, 1 - rho, the formula value,
  and the decorated-flip fixed-point flags (branch_carrier.out: N_lev = 1, 2, 15; the formula agrees class by class).
  NOT PROVED: a fixed-point-free pairing or decorated carrier for n = 0 (mod 4) (the integrality itself is proved
  arithmetically in the 2026-10-20 switching_branch_integrality packet, CITED). The theorem fixes what such a carrier must
  do: weights 2^{-k} with k = v_2(|Syl_2 C(x)|), unequal across Sylow types, and the V_4 stratum already pairs as
  3/4 + 1/4.
  CITED: THM-479 (branch split; N_odd = sum rho(Aut C)), THM-4531 (Theorem B; Lemma N1), Babai-Cameron 2000 (Theorem
  5.1: 2-subgroups of class automorphism groups are semiregular; involutions are fixed-point-free), the switching
  branch integrality packet 05-knowledge/results/switching_branch_integrality_20261020/ (Frobenius root counting), the
  lane module procgen_tbij_20261001_nat.py (class_orbits, group_elements).
source: mac-mini-2026-10-20-jcweb (owner's prompt: work HYP-9172 (b))
depends_on:
  - 01-canon/theorems/THM-479-level-holonomy-and-branch-split.md
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
related:
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
  - 01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md
note: 05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md
scripts: 04-computation/experiments/jcweb_20261020/branch_carrier.py
outputs: 04-computation/experiments/jcweb_20261020/branch_carrier.out
---

# THM-4726 — the level branch as a weighted count over 2-element decorations

**PROVED (pre-audit) for the decomposition; FINITE-EXACT for the census.** Note: [jcweb_20261020_uniform_rules_branch_integrality_jacobian_web](../../05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md).

## 1. The decomposition

THM-479 splits the number of switching classes of tournaments on `n` vertices as `A049313(n) = N_odd(n) + N_lev(n)` with `N_odd = sum_{[C]} rho(Aut C)`, `rho(G)` the proportion of odd-order elements of `G`, hence `N_lev = sum_{[C]} #(even-order elements of Aut C)/|Aut C|`.

**Proposition.** For every finite group `G`, `#(even-order elements)/|G| = sum_{[x]} rho(C_G(x))`, the sum over the `G`-conjugacy classes of nontrivial 2-elements `x`. *Proof.* `g -> (g_2, g_{2'})` is a bijection from the even-order elements onto the pairs (nontrivial 2-element `x`, odd-order element `y`) with `xy = yx`; for fixed `x` the second entry ranges over the odd-order elements of `C_G(x)`, and the class of `x` has `[G : C_G(x)]` members. ∎

Hence `N_lev(n) = sum rho(C_{Aut C}(x))` over the `S_n`-orbits of pairs `(C, x)`. If `Aut C` has a normal 2-complement (Burnside: in particular when its Sylow 2-subgroup `P` is cyclic, which Babai–Cameron make the case for `n = 2 mod 4`), the odd-order elements of `C(x)` are its 2-complement and `rho(C(x)) = 1/|P|` for every `x in P` (all `x` of the class meet `P`), so the class `C` contributes `(|P| - 1)/|P|`. For `n = 2 (mod 4)`, `|P| = 2`, each even-Aut class contributes `1/2` through its unique involution class, and THM-4531's Theorem B pairs these terms by `(C, g) -> (C + M_g, g)` — the decorated form of the flip.

## 2. What happens at n = 0 (mod 4)

`branch_carrier.out` (lane module `class_orbits`, dreadnaut automorphism groups):

| n | even-Aut classes | decorated pairs (C, involution class) | flip-fixed | weights of the pairs (C, 2-element class) | N_lev |
|---|---|---|---|---|---|
| 4 | 2 (|Aut| = 4, 12) | 2 | 2 | 1/4 × 3 (C_4-type), 1/4 (A_4) | 1 |
| 6 | 4 (Sylow 2) | 4 | 0 | 1/2 × 4 | 2 |
| 8 | 27 | 31 | 9 | 1/2 × 19, 1/4 × 14, 1/8 × 16 | 15 |

At `n = 8` the Sylow-type strata contribute `C_2: 19/2`, `C_4: 2 × 3/4`, `V_4: 3/4 (Aut = V_4) + 1/4 (Aut = A_4) = 1`, `C_8: 2 × 7/8`, `D_8: 7/8 (Aut = D_8) + 3/8 (Aut = PSL(2,7))`. The flip `(C, g) -> (C + M_g, g)` has nine fixed decorated pairs (three with `|Aut| = 2`, two `C_4`-type classes, two of the three involution classes of the `V_4` class, two of the three of the `D_8` class), so it does not pair the carrier; a pairing proof must mix the strata `C_2 ∪ C_4` (`19/2 + 3/2 = 11`, the number of flip-transpositions) and `C_8 ∪ D_8 ∪ V_4` (`7/4 + 5/4 + 1 = 4`). The arithmetic integrality (2026-10-20 packet) is a theorem; the decorated carrier for it is still open and the table above is its specification.

## 3. Relation to the affine criterion

By THM-4724 Theorem A with `U` the cut space, `|Fix_C(<g>)| = sum_{X in Z^{<g>}} sgn_X(g)` for every `g`, so `N_lev(n) = (1/n!) sum_{g of even order} sum_{X Euler, gX = X} sgn_X(g)`: the level branch is the signed count of Euler graphs over the even-order permutations. The signed sum is nonzero only for level `g` (all cycle lengths of the same 2-adic valuation, THM-479), which recovers the branch formula and shows that a pairing on the Euler-graph side would have to be an involution on `{(g, X) : gX = X, sgn_X(g) = +1}` minus one on the `-1` part — the same mixed-weight problem in another carrier.
