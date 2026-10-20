---
id: THM-4726
title: "HYP-9172 (b), what a pairing proof of the level-branch integrality must carry: for every n, THM-479's level branch is N_lev(n) = sum over S_n-orbits of pairs (C, x), C a switching class of tournaments and x a nontrivial 2-element of Aut C, of rho(C_{Aut C}(x)), the odd-order proportion of the centraliser (the 2-part/odd-part factorisation g = g_2 g_{2'} of even-order automorphisms); when Aut C has a normal 2-complement this is sum over the nontrivial elements x of a Sylow 2-subgroup P of 1/|P|, so each such C contributes (|P| - 1)/|P| and Theorem B of THM-4531 (n = 2 mod 4, n >= 6, |P| = 2) is the case where every term is 1/2 and the decorated flip (C, g) -> (C + M_g, g) pairs the terms. For n = 0 (mod 4) the decorated flip on pairs (C, g), g an involution class, has fixed points (2 of 2 at n = 4; 7 of 31 at n = 8) and the exact weights rho(C_{Aut C}(x)) take the reciprocal values 1/2, 1/4, 1/8 in the finite n = 8 census (19, 14, 16 pairs; 19/2 + 14/4 + 16/8 = 15 = N_lev(8)), so a denominator-clearing integral-packet interpretation cannot keep every Sylow stratum separate; the Sylow-type strata contribute 19/2 (C_2), 3/2 (C_4), 1 (V_4: Aut = V_4 with 3/4 and Aut = A_4 with 1/4), 7/4 (C_8), 5/4 (D_8: Aut = D_8 with 7/8 and the order-168 group with 3/8)"
status: >
  PROVED (weighted identity independently audited; finite flip correction dated 2026-10-20): the decomposition N_lev = sum_{(C,[x])} rho(C_{Aut C}(x)) (elementary: g -> (g_2, g_{2'}) is a
  bijection from even-order elements to commuting pairs (nontrivial 2-element, odd-order element)); the normal-2-complement
  specialisation; the statement that Theorem B is the |P| = 2 case. FINITE-EXACT: the n = 4, 6, 8 census of even-Aut
  classes with |Aut|, Sylow-2 order, numbers of involution classes and of 2-element classes, 1 - rho, the formula value,
  and the corrected decorated-flip fixed-point flags (independent weighted_branch_carrier_audit_20261020 packet: N_lev = 1, 2, 15; true fixed counts 2, 0, 7). The original branch_carrier.out fixed count at n = 8 is superseded: its 9 counts unmarked type coincidences, not decorated fixed points.
  NOT PROVED: an explicit denominator-clearing decorated integral-packet construction for n = 0 (mod 4) (the integrality itself is proved
  arithmetically in the 2026-10-20 switching_branch_integrality packet, CITED). The exact carrier weights are rho(C_Aut(C)(x)); these are reciprocal Sylow-2 orders when the centralizers have normal 2-complements, as in the stated finite census. The V_4 stratum has total 3/4 + 1/4 = 1. A weight-preserving involution alone need not prove integrality for weights 1/4 or 1/8; integral packets must clear the full denominator.
  CITED: THM-479 (branch split; N_odd = sum rho(Aut C)), THM-4531 (Theorem B; Lemma N1), Babai-Cameron 2000 (Theorem
  5.1: 2-subgroups of class automorphism groups are semiregular; involutions are fixed-point-free), the switching
  branch integrality packet 05-knowledge/results/switching_branch_integrality_20261020/ (Frobenius root counting), the
  lane module procgen_tbij_20261001_nat.py (class_orbits, group_elements).
source: "mac-mini-2026-10-20-jcweb (owner's prompt: work HYP-9172 (b))"
depends_on:
  - 01-canon/theorems/THM-479-level-holonomy-and-branch-split.md
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
related:
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
  - 01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md
note: 05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md
scripts: 04-computation/experiments/jcweb_20261020/branch_carrier.py
outputs: 04-computation/experiments/jcweb_20261020/branch_carrier.out (historical; n=8 decorated fixed count superseded)
audit: 05-knowledge/results/weighted_branch_carrier_audit_20261020/WEIGHTED_BRANCH_AUDIT.md
---

# THM-4726 — the level branch as a weighted count over 2-element decorations

**PROVED (weighted identity independently audited; finite flip correction dated 2026-10-20) for the decomposition; FINITE-EXACT for the census.** Note: [jcweb_20261020_uniform_rules_branch_integrality_jacobian_web](../../05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md).

## 1. The decomposition

THM-479 splits the number of switching classes of tournaments on `n` vertices as `A049313(n) = N_odd(n) + N_lev(n)` with `N_odd = sum_{[C]} rho(Aut C)`, `rho(G)` the proportion of odd-order elements of `G`, hence `N_lev = sum_{[C]} #(even-order elements of Aut C)/|Aut C|`.

**Proposition.** For every finite group `G`, `#(even-order elements)/|G| = sum_{[x]} rho(C_G(x))`, the sum over the `G`-conjugacy classes of nontrivial 2-elements `x`. *Proof.* `g -> (g_2, g_{2'})` is a bijection from the even-order elements onto the pairs (nontrivial 2-element `x`, odd-order element `y`) with `xy = yx`; for fixed `x` the second entry ranges over the odd-order elements of `C_G(x)`, and the class of `x` has `[G : C_G(x)]` members. ∎

Hence `N_lev(n) = sum rho(C_{Aut C}(x))` over the `S_n`-orbits of pairs `(C, x)`. If `H = Aut C` has a normal 2-complement `K` and Sylow 2-subgroup `P`, then for `x in P`,

    rho(C_H(x)) = 1/|C_P(x)|.

Indeed `C_H(x)` maps onto `C_P(x)` under `H -> H/K = P`, with odd-order kernel `C_K(x)`. Conjugacy of elements of `P` in `H` is the same as their conjugacy in `P`. Summing over nonidentity `P`-conjugacy classes gives `sum 1/|C_P(x)| = sum |[x]_P|/|P| = (|P|-1)/|P|`. Equivalently one can sum `1/|P|` over the individual nonidentity elements of `P`. The individual class weight equals `1/|P|` only when `x` is central in `P`; the previous unrestricted claim was false. For `H=P=D_8`, the four nonidentity 2-element conjugacy classes have weights `1/4,1/4,1/4,1/8`, not four equal weights.

A cyclic Sylow 2-subgroup gives a normal 2-complement by Burnside. For `n = 2 (mod 4)`, `n>=6`, semiregularity forces `|P|=2`, each even-Aut class contributes `1/2` through its unique involution class, and THM-4531's Theorem B pairs these terms by `(C,g) -> (C+M_g,g)`. At `n=2` the matching is a cut, the flip is fixed, and `N_lev(2)=1/2`.

## 2. What happens at n = 0 (mod 4)

The [independent exact audit](../../05-knowledge/results/weighted_branch_carrier_audit_20261020/WEIGHTED_BRANCH_AUDIT.md) enumerates the normalized switching quotient and every permutation, then independently checks the marked flips by adjacency matrices. The original `branch_carrier.out` agrees on weights and branch totals but overcounts decorated fixed points:

| n | even-Aut classes | decorated pairs (C, involution class) | flip-fixed | weights of the pairs (C, 2-element class) | N_lev |
|---|---|---|---|---|---|
| 4 | 2 (|Aut| = 4, 12) | 2 | 2 | 1/4 × 3 (C_4-type), 1/4 (A_4) | 1 |
| 6 | 4 (Sylow 2) | 4 | 0 | 1/2 × 4 | 2 |
| 8 | 27 | 31 | 7 | 1/2 × 19, 1/4 × 14, 1/8 × 16 | 15 |

At `n = 8` the Sylow-type strata contribute `C_2: 19/2`, `C_4: 2 x 3/4`, `V_4: 3/4 (Aut = V_4) + 1/4 (Aut = A_4) = 1`, `C_8: 2 x 7/8`, `D_8: 7/8 (Aut = D_8) + 3/8 (the order-168 group)`. These totals reproduce.

The decorated flip `(C,g) -> (C+M_g,g)` has **seven** fixed pairs: three with `|Aut|=2`, two C4-type classes, and two of the three involution classes of the V4 class. Its other 24 decorations form **12 transpositions**. The original code compared only the unmarked canonical types, not isomorphisms commuting with `g`. In the D8 class its two reflection decorations are exchanged, although both underlying types are fixed; this caused the erroneous count nine. The audit supplies explicit state 4722, the two marked involutions, noncommuting coarse isomorphisms, and exhaustive checks of all 384 commuting relabellings for each marking.

The marked stabilizer `C_Aut(C)(g)` is unchanged by this flip, so its weight is preserved. At n=8 the 31 involution decorations have mass 12, the seven fixed ones have mass 5/2, and the noninvolution 2-element decorations contribute the remaining mass 3. This does not produce the required integer-packet carrier. In particular a fixed-point-free weight-preserving involution alone would be insufficient when a 2-cycle contributes `2 x 1/4` or `2 x 1/8`.

The nonintegral individual Sylow-stratum masses preclude an integral-packet interpretation that keeps every stratum separate. They do not force the particular unions `C_2 union C_4` and `C_8 union D_8 union V_4`, or identify their arithmetic totals with numbers of decorated-flip transpositions. Numerical integrality is proved; an explicit decorated or weighted integral-packet construction at `n=0 mod4` remains open.

## 3. Relation to the affine criterion

By THM-4724 Theorem A with `U` the cut space, `|Fix_C(<g>)| = sum_{X in Z^{<g>}} sgn_X(g)` for every `g`, so `N_lev(n) = (1/n!) sum_{g of even order} sum_{X Euler, gX = X} sgn_X(g)`: the level branch is the signed count of Euler graphs over the even-order permutations. The signed sum is nonzero only for level `g` (all cycle lengths of the same 2-adic valuation, THM-479), which recovers the branch formula and gives another exact signed carrier. A proposed cancellation construction must retain the sign, the permutation, and the `1/n!` normalization. The identity alone neither supplies a sign-cancelling involution nor forces every combinatorial proof to use an involution.
