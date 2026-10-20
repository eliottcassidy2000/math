---
id: HYP-9341
title: "Current boundary after THM-4722/THM-4723: uniform equivariance closure, complete mixed-prime census through twelve, weighted twin defects, and unbounded Boolean degree"
status: >
  (a) PROVED: the HYP-9172(a) no-equivariant-type-bijection assertion is closed for both halves at every n>=5
  by the independently reviewed uniform-obstruction packet; this does not prove every stronger odd-subgroup
  switching refinement proposed in the historical route below.
  (b) PROVED: all p-group equalities, exact fusion/sign-escape defect identity, least failure degree ten,
  weighted false-twin formula and sharp order-twelve defect boundary. FINITE-EXACT: the complete mixed-prime
  census through n=12. The all-order unit-weight split/join proposal is REFUTED; a decorated structural
  cancellation remains OPEN, while the global Z3 defect sum already vanishes by the p-group theorem.
  (c) PROVED: all-prime nonlinear equivariant circulant type construction; Boolean-affine existence iff
  p=3 mod4 or p=5; unbounded necessary Boolean degree; exact minima two at p=13 and three at p=17.
  First-order ring definability of some other prime-field bijection remains OPEN. The older bounded search
  and initial questions below are preserved as historical lineage, not as the current status.
  Incoming THM-4724 has author status PROVED (pre-audit), not RESERVED; its new claims are
  preserved below as attributed incoming work, not independently audited dependencies of this checkpoint.
incoming_status_20261020_pre_audit: >
  CONTEXT: the concurrent session-10 packet (uniform_equivariance_obstruction_20261020, prime_local_walsh_bijection_20261020,
  c15_symmetry_counterexample_20261020) proves HYP-9172 (a) for both halves at every n >= 5, proves HYP-9340 (a) for
  all p-groups, refutes it for mixed-prime H (C15 at n = 10, verified here: c15_n10_check.out), and gives an explicit
  per-prime circulant construction ("derivative path"); the items below are what remains after that packet.
  (a) CLOSED by the Sylow route too (THM-4724, 2026-10-20: the affine fixed-point criterion gives the restricted principle
  and the Sylow refinement for switching classes at every n, and Lemma L supplies the Hall violator at every even n >= 6);
  the packet's valuation route is independently audited there; FINITE-EXACT at n = 6 (switching_even_n.out: tau_C = tau_EE for Z5, Z3 (3,3), Z3 (3,1,1,1), Z3xZ3);
  THM-4531 certificates to n = 100. (b) OPEN; PROVED reduction (THM-4723 Theorem B); the defect lists for n <= 8 are
  balanced and the proposed pairing matches every listed case (defect_z3.out). (c) PARTLY CLOSED: THM-4724 Theorem W proves that no derivative rule D_c (c != 1) is a type bijection for all primes
  (c = -1 mod the odd part of Z_p^* at infinitely many p sends the symmetric tournaments to the odd complete graph) and gives
  the forcing analysis for general dilation rules; the general non-existence remains OPEN. Earlier: FINITE-EXACT negative search
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

## Current reconciliation, 20 October 2026

The [ongoing checkpoint map](../results/overnight_prime_arithmetic_20261020/README.md) owns the current continuation. The [uniform obstruction and completed audit](../results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md) settle the no-equivariant-type-bijection assertion in (a), not every stronger refinement contemplated below. The [complete census](../results/mixed_prime_census_extension_20261020/README.md) settles the least mixed-prime failure degree and gives tau_T-tau_E=F_E+Z-F_T, retaining the actual permutation embedding. The [weighted twin-quotient formulas](../results/twin_quotient_defect_weights_20261020/README.md) refute an all-order unit-weight split/join cancellation, with the first magnitude-two defects and the first negative odd defect at order twelve; decorated structural cancellation remains open. The [affine](../results/affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md) and [Boolean-degree](../results/circulant_boolean_degree_20261020/BOOLEAN_DEGREE_BOUNDARY.md) theorems separate affine impossibility from the already existing all-prime nonlinear algorithm.

The original sections below are retained as **historical questions and finite-search provenance**. Their “open”, “first found”, “pending audit”, proposed pairing and forced-empty language must be read with the current status above and [HYP-9340](HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md), which also corrects the labelled/type and largest-odd-subgroup distinctions.

## (a) Switching classes at even n

THM-4722 proves the switching half of HYP-9172 for odd n through THM-4531 P1 (class ≅ its mod-4-balanced member) and the odd-n Euler parity lemma (`V = Z ⊕ B`, invariant cuts carry no self-paired orbit). For even n both fail: classes can have even-order stabilisers and the lemma has 36 counterexamples at n = 6 and n = 8. Yet `tau_C(P) = tau_EE(P)` holds at n = 6 for every odd-order class. Conjecture: the refined identity holds for switching classes at every n and every odd-order P, and the Hall violator {Sylow-3, Sylow-5} (with the empty graph) still works once one shows that no switching class on an even number of vertices has a stabiliser containing both Sylow subgroups.

## (b) HYP-9340 beyond Sylow (p-groups proved and mixed primes refuted by the packet; what is left)

The refutation shows the exact failure mechanism predicted by THM-4721 Proposition Q: for mixed-prime H the N(H)-even odd graphs (K_{5,5} for C15) and the multiplicity defects need not balance. Open: a formula for tau_T(H) - tau_E(H) for general odd H in terms of these two defects, and the smallest mixed-prime H with a nonzero difference (C15 at n = 10 is the first found; Z_15 at n = 8, 9 balances).

By THM-4723, for `H = Z_3` of type `(3,1^{n-3})` the conjecture is `Σ_X (b(X) - e(X)) = 0`. Conjectured involution on defective graph types: an even defective graph has a twin class `C` with `floor(|C|/3)` even or two swappable twin classes; split `C` as `3 + 3 + rest` (resp. take the two classes), join the two triples completely, keeping everything else — the result is odd (the swap of the two triples reverses an odd number of edges) with `b = 1`; and conversely. All 18 defective types for n <= 8 pair up this way. A proof must handle twin classes of size >= 9 (where `floor(c/3)` is odd and no defect occurs) and nested configurations. For general odd H the analogous reduction needs "H-blow-ups" (Theorem S' of THM-4722 handles Sylow H; THM-4723 Theorem P reformulates the p-subgroup case as equality of Sylow profiles).

## (c) A uniform rule for prime circulants

THM-4721 Corollary C gives an AGL(1,p)-equivariant bijection between circulant tournaments and even circulant graphs of order p for every prime p, by Hall on the lattice of odd multiplier subgroups. Obstruction to uniformity: the tournaments invariant under the largest odd multiplier subgroup `K_max` (QR_p for p ≡ 3 mod 4) must go to the empty graph, while a multiplier derivative `D_c(S) = {x : [x ∈ S] ≠ [cx ∈ S]}` sends them to the empty graph iff `c ∈ K_max`, a condition on p (for `c = 4`: iff 4 is in the odd part of `Z_p^*`, which fails at p = 5, 13, 17 and holds at 7, 11, 19, 23 — exactly the passing primes of the search). A uniform rule would have to be a formula whose behaviour on `K_max`-invariant sets is "empty" for every p; whether a first-order formula in the ring language (with `S` as a predicate) can do this and be bijective is open. The natural guess is that no such formula exists because the odd part of `p - 1` is not first-order definable uniformly in p.


**Incoming-work scope:** The following update and the `incoming_status_20261020_pre_audit` metadata preserve the concurrent author's text. [THM-4724](../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md) is labelled **PROVED (pre-audit)** by its author. This checkpoint does not independently audit or endorse its new theorems, and does not use them as proved dependencies. The reviewed reconciliation above controls the older general defect and affine/Boolean-degree questions.

## Update 2026-10-20 (THM-4724, mac-mini-2026-10-20-oddprimes)

- (a) is closed by the Sylow route: [THM-4724](../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md) proves the per-h identity `|Fix_{T/U}(<P,h>)| = Σ_{X ∈ (U-perp)^{<P,h>}} sgn_X(h)` for every submodule `U` at every `n` (affine fixed-point criterion; checked for `U = 0` and `U = cuts` on 1832 pairs, n = 5..8), hence the refined identity `tau_C(P) = tau_EE(P)` for Sylow `P` at even `n`, and Lemma L (no class stabiliser at even n contains a p-cycle and a q-cycle for two primes above n/2, resp. the 5-cycle/7-cycle and the Sylow 3-subgroup at n = 6, 10) gives the Hall violator. The Euler parity lemma was never needed: the criterion decides existence of fixed classes directly.
- (b) is settled by the packet (p-groups proved via the Sylow action theorem, mixed primes refuted by C15); THM-4724 independently recounts C15 (32 vs 31) and confirms equal Sylow-p profiles for p = 3, 5, 7, n <= 8. The defect involution for the 3-cycle type remains a nice open bijective exercise, no longer needed for the count.
- (c) the derivative family is dead for every `c ≠ 1` (Theorem W). New angle: the power-sum detector (THM-4724 Proposition D): a circulant tournament of prime order has a nontrivial multiplier symmetry iff all power sums `σ_k(S)`, `gcd(k, p-1) = 1`, vanish, provided the odd part of `p-1` is a prime power or `p ≡ 3 (mod 4)` with two odd primes; it fails at `p = 61` and `p = 211`. In the regimes where it holds, the least unit index with `σ_k ≠ 0` gives a canonical representative of every rigid orbit, uniformly in `p`, so the residue of (c) is a uniform bijection between normalised rigid tournament sets and normalised rigid even symmetric sets. OPEN.
