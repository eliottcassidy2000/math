---
id: HYP-9342
title: "What a uniform rule for the prime-circulant bijection (circulant tournaments versus even circulant graphs of prime order) can still be after THM-4725: (a) at p = 3 (mod 4) the Paley rule beta is canonical (Gamma-equivariant, one quantifier in the multiplicative structure) and the question is whether it is the ONLY Gamma-equivariant type bijection; (b) at p = 1 (mod 4) no Gamma-equivariant rule exists and no dilation rule works at infinitely many primes, so a uniform rule must use the additive structure of F_p (convolutions, representation numbers, Gauss-sum data): is there one, and is there a p = 1 (mod 4) analogue of the Paley reference (a canonical tournament set invariant under the fourth powers, e.g. built from the biquadratic character)? (c) the residue of HYP-9341 (c): no multiplier-equivariant rule definable in (Z_p^*, ., S) with constants is uniform"
status: >
  (a) PARTLY SETTLED (FINITE-EXACT, additive_rules.out): the Gamma stabiliser-containment graph has exactly ONE perfect
  matching at p = 7 and p = 11, so the Paley rule is the unique Gamma-equivariant type bijection there; at p = 19 there are
  96 perfect matchings (the Paley rule is one of them), so uniqueness fails as an abstract equivariant bijection and
  'canonical' must mean definable. Hall's condition for Gamma-equivariant injections holds at every p <= 31 and the
  Gamma-orbit counts agree at p = 3 (mod 4) (power_equivariance*.out).
  (b) OPEN. THM-4725 Theorem D excludes every Gamma-equivariant rule at p = 1 (mod 4), p >= 13 (Gamma-orbits 4 < 6, 3 < 10,
  65 < 125, 673 < 1335 at p = 13, 17, 29, 37); Theorem U excludes every dilation rule at infinitely many p = 9 (mod 16).
  A rule using + is neither excluded nor known. First tests (additive_rules.out): R1(S) = {x : |S ∩ (x + S)| odd} is
  symmetric and multiplier-equivariant; at p = 13 and p = 29 it is injective (hence bijective) from tournament types onto
  symmetric-set types but some images are odd circulants; at p = 17 all images are even but two types collide; the
  representation-number rule {x : #(a + b = x) odd} and the 3-term-progression rule are not even symmetric. For p = 5 (mod 8) the tournaments invariant under the fourth powers form
  one multiplier type and must go to the empty graph (THM-4724 Lemma F); a biquadratic reference tournament exists
  (the fourth-power classes {1, i, -1, -i} split as a tournament set only if i and -i are separated, which the biquadratic
  character does for p = 5 mod 8 but then the set is not multiplier-canonical).
  (c) PROVED by THM-4725 Theorem U for rules whose atoms are dilations by fixed integers; parameter-free rules are Gamma-
  equivariant and fall under Theorem D.
source: mac-mini-2026-10-20-jcweb
related:
  - 01-canon/theorems/THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md
  - 01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md
  - 05-knowledge/hypotheses/HYP-9341-what-remains-open-after-thm-4722.md
note: 05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md
---

# HYP-9342 — what a uniform rule can still be

THM-4725 settles the two natural classes: rules with fixed dilation constants (dead at infinitely many primes by the `K_max` collapse) and rules natural for the multiplicative structure (exist exactly at `p = 3 (mod 4)`, where the Paley tournament is a canonical reference and the rule is `beta(S) =` the even-class member of `{backward classes of S relative to QR_p, its complement}`). What is left is the additive structure.

1. **Uniqueness of the Paley rule at `p = 3 (mod 4)`.** Is `beta` the only `Gamma`-equivariant type bijection, or at least the only one definable in `(Z_p^*, ·, S)`? Test: enumerate all `Gamma`-equivariant type bijections for `p = 7, 11, 19, 23` (Lemma N1 matchings) and compare with `beta`.
2. **An additive rule at `p = 1 (mod 4)`.** Every `Gamma`-equivariant rule fails there because the flip `x -> x chi(x)` is a `Gamma`-odd symmetry of the even circulants with an odd number of non-residue classes. A rule using `+` breaks the power symmetry. Candidates to test for `p = 13, 17, 29`: `S -> {x : |S ∩ (x + S)| odd}`, `S -> {x : |S ∩ (x - S)| odd}`, the parity of the number of 3-term progressions through `x` inside `S`, and Gauss-sum-weighted versions. Each must send the fourth-power-invariant tournaments to the empty graph (THM-4724 Lemma F for `p = 5 (mod 8)`).
3. **The multiplier-canonical reference problem.** At `p = 3 (mod 4)` the reference is `QR_p`, the unique tournament set invariant under the index-two subgroup. At `p = 1 (mod 4)` no tournament set is invariant under a subgroup of index two; the largest odd subgroup `K_max` has index `2^e >= 4` and its invariant tournament sets form at least one type (two for `e >= 3`). A reference would have to be a canonical choice among them, which the power symmetry forbids; this is the structural reason behind Theorem D.
