---
id: THM-4723
title: "Two reductions of HYP-9340 (a). (i) Sylow-profile form: for an odd prime p, tau_T(P) = tau_E(P) holds for every p-subgroup P of S_n iff the tournaments and the even graphs on n vertices have the same distribution of Sylow p-subgroups up to S_n-conjugacy (Moebius inversion on the poset of p-subgroup classes). (ii) Defect form for the first open case H = Z_3 of type (3,1^{n-3}): tau_T(H) = sum over graph types X of b(X), where b(X) is the alternating count of pair-even markings of X by disjoint independent twin triples with the contracted sign (marked identity + blow-up + odd-order inclusion-exclusion); so HYP-9340(H) says sum_X (b(X) - e(X)) = 0, e(X) = [X even with a twin triple]. The defect b - e is supported on even graphs with a twin class of size c, floor(c/3) even, or with two swappable twin classes (b = 0, e = 1), and on odd graphs with two adjacent swappable twin triples (b = 1, e = 0); it is balanced for n <= 8 (0, 2, 4, 12 defective types)"
status: >
  CONTEXT (2026-10-19/20): HYP-9340 (a) for general odd-order H is REFUTED by the session-10 packet's C15
  counterexample (H = <(012)(34567)> <= S_10: 32 tournament types versus 31 even-graph types; the missing graph is
  K_{5,5}, globally odd but N(H)-even), re-verified here independently (c15_n10_check.out: 32 and 31 from the 512
  invariant objects); and HYP-9340 (a) for p-groups H is PROVED there (SYLOW_ACTION_ADDENDUM.md: a Walsh matching
  preserving the S_n-class of a Sylow p-subgroup, via a signed p-group orbit-localisation lemma), which by Theorem P
  below is exactly the equality of Sylow-p profiles. So the 3-cycle case the owner asked for is settled by that
  packet; Theorem B gives a second, self-contained route whose last step (the defect cancellation) is not carried out.
  PROVED (pre-audit): Theorem P (Sylow-profile equivalence); Theorem B (sum_X b(X) = tau_T(Z_3), with the contracted
  sign s'); Corollary (HYP-9340 for the 3-cycle type iff the defects cancel); the single-class evaluation
  b(X) = [floor(c/3) odd] for an even graph whose only big twin class has size c >= 3 and full symmetric action.
  FINITE-EXACT: sum_X b(X) = tau_T = 4, 11, 44, 280 and the defect lists for n = 5..8 (defect_z3.out; the first
  implementation with a formula for the sign difference was wrong - the contracted sign is the correct one).
  NOT PROVED here: the cancellation itself (now a corollary of the session-10 p-group theorem; HYP-9341 (b)); the pairing suggested by the data (split a twin class as
  3 + 3 + rest and join the two triples completely) is a conjecture.
  CITED: THM-4721 (Theorems G', R, S, Proposition Q), THM-4531 (Lemma N1), Sylow.
source: mac-mini-2026-10-19-hyp9172
depends_on:
  - 01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md
related:
  - 05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md
  - 05-knowledge/hypotheses/HYP-9341-what-remains-open-after-thm-4722.md
note: 05-knowledge/results/hyp9172_20261019_sylow_hall_violator.md
scripts: 04-computation/experiments/even_graphs_20261019/defect_z3.py
---

# THM-4723 — two reductions of HYP-9340 (a)

## 1. The Sylow-profile form

Fix an odd prime `p`. For a type `x` (tournament or even graph) let `σ(x)` be the `S_n`-conjugacy class of a Sylow `p`-subgroup of `Aut x`, and `a_T(Q)`, `a_E(Q)` the numbers of tournament resp. even-graph types with `σ = [Q]`.

**Theorem P.** `tau_T(P) = tau_E(P)` for every `p`-subgroup `P <= S_n` iff `a_T(Q) = a_E(Q)` for every `p`-subgroup class `[Q]`.

*Proof.* A `p`-subgroup of `Aut x` lies in a Sylow `p`-subgroup of `Aut x`, and these are conjugate, so `Aut x` contains a conjugate of `P` iff `σ(x)` contains a conjugate of `P`. Hence `tau(P) = Σ_{[Q]} a(Q) [P <=_c Q]` on both sides. The matrix `([P <=_c Q])` over `p`-subgroup classes is unitriangular when the classes are ordered by group order, so the two vectors `(tau(P))_P` and `(a(Q))_Q` determine each other. ∎

So HYP-9340 for `p`-subgroups says that tournaments and even graphs have the same Sylow-`p` statistics. THM-4721 Theorem S is the top class; the conjecture is the statement for all lower classes, and by Theorem P it is a statement about how often each `p`-subgroup is the full `p`-part of an automorphism group.

## 2. The defect form for the 3-cycle type

Let `H = Z_3` of type `(3, 1^{n-3})`. A tournament has an automorphism of this type iff it has a 3-cycle module; a graph iff it has an independent twin triple (three pairwise non-adjacent vertices with the same neighbourhood). Distinct 3-cycle modules of a tournament are disjoint (two overlapping ones would force a vertex of the second module to beat or lose to the whole first module, contradicting the 3-cycle); twin triples of a graph overlap inside twin classes of size `>= 4`.

For a graph `X̂` on `n` vertices and a set `Ĵ` of `j` pairwise disjoint independent twin triples, contract each triple to a vertex: a graph `X` on `n - 2j` vertices with a marked `j`-set `J`; every automorphism `a` of `X̂` preserving `Ĵ` induces an automorphism `α` of `(X, J)`, and `s'_Ĵ(a) := sgn_X(α)` is a character of `Stab(Ĵ)`. Put

    b(X̂) = Σ_{j >= 1} (-1)^{j+1} #{Aut(X̂)-orbits of such Ĵ with s'_Ĵ trivial on Stab(Ĵ)},    e(X̂) = [X̂ even and has a twin triple].

**Theorem B.** `Σ_{X̂} b(X̂) = tau_T(H)`, the sum over all graph types on `n` vertices.

*Proof.* Let `A_j` be the number of types of tournaments on `n` vertices with `j` marked 3-cycle modules and `B_j` the number of types of pairs `(X̂, Ĵ)` with `s'` trivial on the stabiliser. Blowing a marked vertex up into a 3-cycle module (resp. an independent twin triple) is a bijection of types between tournaments (graphs) on `n - 2j` vertices with a marked `j`-set and tournaments (graphs) on `n` vertices with `j` marked modules (disjoint twin triples), and under it `Stab(X̂, Ĵ) -> Stab(X, J)` is onto with `s' = sgn_X`. THM-4721 Theorem G' for the group `S_{n-2j}` acting on pairs (object, marked `j`-set) gives `A_j = B_j` (the `g`-fixed pairs are `(fixed j-sets) × (fixed objects)` on both sides, the sign enters only through the object). On the tournament side `Σ_j (-1)^{j+1} A_j = Σ_T (1/|Aut T|) Σ_{a} Σ_{∅≠S ⊆ M_T, aS = S} (-1)^{|S|+1}`; every `a`-orbit on the set `M_T` of modules has odd size (`|Aut T|` is odd), so the inner sum is `1 - Π(1 + (-1)^{odd}) = 1` when `M_T ≠ ∅` and `0` otherwise, giving `tau_T(H)`. On the graph side the signed Cauchy–Frobenius count of orbits with trivial character gives `Σ_j (-1)^{j+1} B_j = Σ_{X̂} b(X̂)`. ∎

**Corollary.** HYP-9340 holds for `H` on `n` vertices iff `Σ_{X̂} (b(X̂) - e(X̂)) = 0`.

**Where the defect lives.** If `X̂` is even and has no automorphism swapping two adjacent twin triples, then `s' = 1` and `b = Σ_j (-1)^{j+1} (number of orbits of j-matchings)`; for a single twin class of size `c` carrying the full symmetric group this is `Σ_{j=1}^{floor(c/3)} (-1)^{j+1} = [floor(c/3) odd]`, so `b - e = -[floor(c/3) even]`: the empty graph on 6, 7, 8 vertices, `K_{1,6}`, `K_{1,7}`, `K_{2,6}` are defective (`b = 0`), as are `2K_{1,3}` and `K_{4,4}` (two swappable classes). Odd graphs with two adjacent swappable twin triples (`K_{3,3}` and its extensions by isolated, universal or pendant vertices, `K_{4,4}` minus an edge, the double star) have `b = 1`. Census (`defect_z3.out`): `n = 5, 6, 7, 8`: defective types `0, 2, 4, 12`, exactly half even and half odd, totals `Σ b = Σ e = 4, 11, 44, 280`. The data suggest an involution "split a twin class as `3 + 3 + rest` and join the two triples completely" (`∅_6 ↔ K_{3,3}`, `K_{1,6} ↔ K_{1,3,3}`, `K_{2,6} ↔ K_{2,3,3}`, `2K_{1,3} ↔` double star, `K_{4,4} ↔ K_{4,4} - e`); proving that it is a sign-reversing involution on defective types for every `n` would prove HYP-9340 for the 3-cycle type (HYP-9341 (b)).
