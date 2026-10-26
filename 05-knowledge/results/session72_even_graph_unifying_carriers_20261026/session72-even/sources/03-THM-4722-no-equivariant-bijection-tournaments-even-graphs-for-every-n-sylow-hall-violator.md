---
id: THM-4722
title: "HYP-9172 (a), tournament half, for every n >= 5: no S_n-equivariant map from tournaments to graphs is bijective from tournament types onto even-graph types. Proof: the Sylow refinement (THM-4721 Theorem S) gives tau_T(P_3) = tau_E(P_3) and tau_T(P_5) = tau_E(P_5) for Sylow 3- and 5-subgroups P_3, P_5 of S_n; no odd-order subgroup of S_n contains conjugates of both (a 3-cycle and a 5-cycle sharing a point generate a group of even order, by Jordan); the empty graph has both; so the stabiliser-containment graph of THM-4531 Lemma N1 has the Hall violator S = {T : Aut T contains a conjugate of P_3 or of P_5} with |N(S)| <= |S| - 1. The same argument proves the switching-class half for every odd n >= 5, via THM-4531 P1 and an odd-n Euler parity lemma (the invariant cuts carry no self-paired orbits, V = Z + B)"
status: >
  PROVED (pre-audit): Lemma J (Jordan), Theorem H (tournament half, all n >= 5), Lemma E' (odd-n Euler parity lemma),
  Lemma Q (odd-n realisability of parity vectors), Theorem R' (restricted principle for switching classes, odd n),
  Theorem S' (Sylow refinement for switching classes, odd n), Theorem H' (switching half, all odd n >= 5).
  FINITE-EXACT: deficit >= 1 certificates n = 5..9 from the census (hall_sylow_certificate.out); |<c,d>| = 60, 360, 2520
  for the three overlap patterns; tau_C(P) = tau_EE(P) for every odd-order class at n = 5, 7 and for fifteen classes at n = 9
  (switching_sylow.out, switching_sylow_n9.out: e.g. 51 = 51 for Z3 of type (3,1^6), 1 = 1 for Z3 wr Z3) and at n = 6
  (switching_even_n.out; even n, outside Theorem H'); the Euler parity lemma holds for all 24 + 284 pairs at
  n = 5, 7 and fails for 36 + 36 pairs at n = 6, 8 (euler_parity_lemma.out), as the proof predicts.
  CITED: THM-4721 (Theorems R, S); THM-4531 (Lemma N1, Theorem G, Theorem O3, P1: for odd n every switching class
  contains exactly one tournament with d+ = d- (mod 4) at every vertex, an isomorphism of S_n-sets); Jordan's theorem
  (a primitive permutation group containing a 3-cycle contains the alternating group); Sylow; Frattini;
  Royle-Praeger-Glasby-Freedman-Devillers 2023 (the open question answered here in the equivariant sense);
  Mallows-Sloane 1975.
  CONCURRENT: the session-10 packet 05-knowledge/results/uniform_equivariance_obstruction_20261020/ (commit
  25752f463a, dated 2026-10-20 by its clock, pushed while this session was running) proves HYP-9172 (a) for every
  n >= 5 for BOTH halves at every n by a different mechanism (a signed Walsh determinant forcing matched stabilisers to
  preserve odd-prime valuations, then primes 3 and 5); THM-4722 is an independent proof of the tournament half and of
  the odd-n switching half via the Sylow refinement, and claims no priority. Its even-n switching gap is closed there.
  OPEN here: the switching-class half for even n by THIS route (certificates 5 <= n <= 100 in THM-4531; the refined identity holds at
  n = 6 for every odd class although the Euler parity lemma fails there, so the even case needs a different proof;
  HYP-9341).
  Notation as in THM-4721: tau_T(H), tau_E(H) = numbers of isomorphism types of tournaments / even graphs whose
  automorphism group contains a conjugate of H; tau_C, tau_EE the same for switching classes / even Euler graphs.
source: mac-mini-2026-10-19-hyp9172 (owner's prompt: prove or disprove HYP-9172 thinking in odd primes and their powers)
depends_on:
  - 01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
related:
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
  - 05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md
  - 05-knowledge/hypotheses/HYP-9341-what-remains-open-after-thm-4722.md
note: 05-knowledge/results/hyp9172_20261019_sylow_hall_violator.md
scripts: 04-computation/experiments/even_graphs_20261019/{hall_sylow_certificate,switching_sylow,switching_even_n,euler_parity_lemma}.py
---

# THM-4722 — no equivariant bijection between tournaments and even graphs, for every n >= 5

**PROVED (pre-audit).** Note: [hyp9172_20261019_sylow_hall_violator](../../05-knowledge/results/hyp9172_20261019_sylow_hall_violator.md).

## 1. The group-theoretic lemma

**Lemma J.** Let `n >= 5` and let `P_3`, `P_5` be Sylow 3- and 5-subgroups of `S_n`. No subgroup of odd order contains a conjugate of `P_3` and a conjugate of `P_5`.

*Proof.* Let `G` have odd order and contain `P_3' ~ P_3`, `P_5' ~ P_5`. `P_3'` contains a 3-cycle on each of its 3-blocks, so its 3-cycles cover `n - (n mod 3) >= n - 2` points; the 5-cycles of `P_5'` cover `n - (n mod 5) >= n - 4` points; for `n >= 7` the two sets of moved points share a point, and for `n = 5, 6` they do as well (`P_5'` moves 5 points, `P_3'` moves 3 resp. 6). Take a 3-cycle `c` of `P_3'` and a 5-cycle `d` of `P_5'` through a common point and let `Ω = supp(c) ∪ supp(d)`, `|Ω| ∈ {5, 6, 7}`. `<c, d>` is transitive on `Ω` and contains the 3-cycle `c`. If `|Ω|` is 5 or 7, a transitive group of prime degree is primitive; if `|Ω| = 6`, the stabiliser of the point fixed by `d` contains `d`, transitive on the other five, so the group is 2-transitive, hence primitive. By Jordan's theorem `<c, d>` contains `A_Ω`, of even order, so `|G|` is even. ∎ (Computation: `<c,d>` is `A_7`, `A_6`, `A_5` for overlaps 1, 2, 3.)

## 2. The tournament half

**Theorem H.** For every `n >= 5` there is no `S_n`-equivariant map from labelled tournaments to labelled graphs on `[n]` that induces a bijection from tournament types onto even-graph types.

*Proof.* By THM-4531 Lemma N1 such a map exists iff the bipartite graph on types, `T ~ X` iff `Aut T` is contained in a conjugate of `Aut X`, has a perfect matching (the two sides have the same size by Royle et al.). Let `S` be the set of tournament types whose automorphism group contains a conjugate of `P_3` or of `P_5`. By Lemma J no type is counted twice (`Aut T` has odd order), so `|S| = tau_T(P_3) + tau_T(P_5)`. If `T ~ X` with `T ∈ S` then `Aut X` contains a conjugate of `P_3` or of `P_5`, so `N(S) ⊆ E_3 ∪ E_5`, `E_p` = even types with a conjugate of `P_p` in their automorphism group. The empty graph lies in `E_3 ∩ E_5`, so `|N(S)| <= tau_E(P_3) + tau_E(P_5) - 1`. THM-4721 Theorem S gives `tau_E(P_p) = tau_T(P_p)` for `p = 3, 5`. Hence `|N(S)| <= |S| - 1`, Hall's condition fails, and no perfect matching exists. ∎

For `n = 5..9` the census gives `|S| = 5, 3, 8, 28, 49` against `|N(S)| <= 4, 2, 7, 27, 48`; the exhaustive deficits of THM-4531 (1, 2, 5, 12, 44) are larger because `S` is only one violator.

This answers the question of Royle, Praeger, Glasby, Freedman and Devillers in the equivariant sense for every `n >= 5`: the counting identity cannot come from a construction that commutes with relabelling.

## 3. The switching-class half for odd n

Fix `n` odd. `V = F_2^{E(K_n)}`, `B` = cut space, `Z = B^⊥` = Euler graphs; `V = Z ⊕ B` as `S_n`-modules (odd `n`). Switching classes `𝒞 = Tour/B`. THM-4531 P1: the map "the unique member with `d+ = d- (mod 4)` everywhere" is an isomorphism of `S_n`-sets `𝒞 ≅ 𝒯_4`, the tournaments with all scores `≡ (n-1)/2 (mod 2)`; in particular stabilisers of classes have odd order.

**Lemma E' (odd-n Euler parity lemma).** Let `P <= S_n` have odd order, `h ∈ N(P)` of even order, `G = <P, h>`. Some `G`-invariant Euler graph `X` has `r(h, X) = 1`. *Proof.* THM-4721 Lemma P gives `Y ∈ V^G` with `r(h, Y) = 1`. Write `Y = Y_Z + Y_B` with `Y_Z ∈ Z^G`, `Y_B ∈ B^G` (the decomposition is `G`-equivariant). A `G`-invariant cut is `δ(A)` with `gA ∈ {A, A^c}`, and `|A| ≠ |A^c|` forces `gA = A`; an edge `{a, b}` of `δ(A)` with `a ∈ A`, `b ∉ A` is never swapped by a power of `h`, so `δ(A)` has no self-paired `h`-orbit and `r(h, δ(A)) = 0` (Royle et al. Lemma 4.3). Hence `r(h, Y_Z) = r(h, Y) = 1`. ∎

**Lemma Q (realisability).** For odd `n` the degree-parity map `V -> F_2^n` restricts to a bijection from `B` onto the even-weight vectors (`δ(A) ↦ 1_A` or `1_{A^c}`, whichever has even size), `S_n`-equivariantly; so for every subgroup `G`, every `G`-invariant even-weight parity vector is the parity vector of a `G`-invariant graph.

**Theorem R'.** For odd `n`, odd-order `P` and `Γ <= N(P)`: `#(Γ-orbits on Fix_𝒞(P)) = #(Γ-orbits on the Γ-even Euler graphs in Fix(P))`. *Proof.* Cauchy–Frobenius on both sides; fix `h ∈ Γ`, `G = <P, h>`. Via P1, `Fix_𝒞(G) = Fix_{𝒯_4}(G)`. If `|G|` is odd, choose a `G`-invariant tournament `T_0`; `T_0` reversed along `X ∈ V^G` has all scores `≡ s` iff the parity vector of `X` equals a fixed `G`-invariant even-weight vector, so by Lemma Q the solutions form a coset of `Z^G` and `|Fix_{𝒯_4}(G)| = |Z^G| = Σ_{X ∈ Z^G} sgn_X(h)` (`h` of odd order reverses nothing). If `|G|` is even, `Fix_{𝒯_4}(G) = ∅` (stabilisers in `𝒯_4` have odd order) and `Σ_{X ∈ Z^G} (-1)^{r(h,X)} = 0` by Lemma E'. ∎

**Theorem S'.** For odd `n`, an odd prime `p` and a Sylow `p`-subgroup `P` of `S_n`: `tau_C(P) = tau_EE(P)`. *Proof.* As THM-4721 Theorem S: Sylow gives multiplicity one on both sides (stabilisers of classes have odd order and contain `P` as a Sylow subgroup; so do the automorphism groups of even Euler graphs containing `P`), Frattini makes `N(P)`-even equal to even on `Fix(P)`, and Theorem R' with `Γ = N(P)` finishes. ∎

**Theorem H'.** For every odd `n >= 5` there is no `S_n`-equivariant map from switching classes to Euler graphs that is bijective from class types onto even-Euler-graph types. *Proof.* As Theorem H with Theorem S' in place of Theorem S: stabilisers of classes have odd order, so Lemma J keeps the two families disjoint; the empty graph is an even Euler graph in both `E_3` and `E_5`. ∎

For even `n` the refined identity still holds at `n = 6` for every odd-order class, but Lemma E' fails there (36 of 132 pairs), classes can have even-order stabilisers, and the argument above does not apply; THM-4531's certificates cover `5 <= n <= 100`.
