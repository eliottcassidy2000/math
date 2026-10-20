---
id: HYP-9341
title: "What remains open after THM-4722/THM-4723: (a) the switching-class half of HYP-9172 for even n (refined identity true at n = 6 for every odd class, Euler parity lemma false there); (b) HYP-9340 beyond Sylow, reduced for the 3-cycle type to the cancellation of an explicit defect, conjecturally by the involution 'split a twin class as 3+3+rest and join the two triples'; (c) a uniform rule for the prime-circulant bijection: no quantifier-free rule with <= 3 affine atoms passes more than 4 of the primes 5..23, and any rule S -> D_c(S) sends the K_max-invariant tournaments to the empty graph iff c lies in K_max, which varies with p"
status: >
  (a) OPEN; FINITE-EXACT at n = 6 (switching_even_n.out: tau_C = tau_EE for Z5, Z3 (3,3), Z3 (3,1,1,1), Z3xZ3);
  THM-4531 certificates to n = 100. (b) OPEN; PROVED reduction (THM-4723 Theorem B); the defect lists for n <= 8 are
  balanced and the proposed pairing matches every listed case (defect_z3.out). (c) OPEN; FINITE-EXACT negative search
  (circulant_rule_search.out: 10 atoms a x + b, all boolean combinations of <= 3 atoms containing x in S, primes
  5..23: best rules pass {7, 11, 19, 23} only, all of them variants of D_4).
source: mac-mini-2026-10-19-hyp9172
related:
  - 01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md
  - 01-canon/theorems/THM-4723-sylow-profile-form-of-hyp-9340-and-the-defect-reduction-for-the-three-cycle-type.md
  - 05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
note: 05-knowledge/results/hyp9172_20261019_sylow_hall_violator.md
---

# HYP-9341 — what remains open after THM-4722 and THM-4723

## (a) Switching classes at even n

THM-4722 proves the switching half of HYP-9172 for odd n through THM-4531 P1 (class ≅ its mod-4-balanced member) and the odd-n Euler parity lemma (`V = Z ⊕ B`, invariant cuts carry no self-paired orbit). For even n both fail: classes can have even-order stabilisers and the lemma has 36 counterexamples at n = 6 and n = 8. Yet `tau_C(P) = tau_EE(P)` holds at n = 6 for every odd-order class. Conjecture: the refined identity holds for switching classes at every n and every odd-order P, and the Hall violator {Sylow-3, Sylow-5} (with the empty graph) still works once one shows that no switching class on an even number of vertices has a stabiliser containing both Sylow subgroups.

## (b) HYP-9340 beyond Sylow

By THM-4723, for `H = Z_3` of type `(3,1^{n-3})` the conjecture is `Σ_X (b(X) - e(X)) = 0`. Conjectured involution on defective graph types: an even defective graph has a twin class `C` with `floor(|C|/3)` even or two swappable twin classes; split `C` as `3 + 3 + rest` (resp. take the two classes), join the two triples completely, keeping everything else — the result is odd (the swap of the two triples reverses an odd number of edges) with `b = 1`; and conversely. All 18 defective types for n <= 8 pair up this way. A proof must handle twin classes of size >= 9 (where `floor(c/3)` is odd and no defect occurs) and nested configurations. For general odd H the analogous reduction needs "H-blow-ups" (Theorem S' of THM-4722 handles Sylow H; THM-4723 Theorem P reformulates the p-subgroup case as equality of Sylow profiles).

## (c) A uniform rule for prime circulants

THM-4721 Corollary C gives an AGL(1,p)-equivariant bijection between circulant tournaments and even circulant graphs of order p for every prime p, by Hall on the lattice of odd multiplier subgroups. Obstruction to uniformity: the tournaments invariant under the largest odd multiplier subgroup `K_max` (QR_p for p ≡ 3 mod 4) must go to the empty graph, while a multiplier derivative `D_c(S) = {x : [x ∈ S] ≠ [cx ∈ S]}` sends them to the empty graph iff `c ∈ K_max`, a condition on p (for `c = 4`: iff 4 is in the odd part of `Z_p^*`, which fails at p = 5, 13, 17 and holds at 7, 11, 19, 23 — exactly the passing primes of the search). A uniform rule would have to be a formula whose behaviour on `K_max`-invariant sets is "empty" for every p; whether a first-order formula in the ring language (with `S` as a predicate) can do this and be bijective is open. The natural guess is that no such formula exists because the odd part of `p - 1` is not first-order definable uniformly in p.
