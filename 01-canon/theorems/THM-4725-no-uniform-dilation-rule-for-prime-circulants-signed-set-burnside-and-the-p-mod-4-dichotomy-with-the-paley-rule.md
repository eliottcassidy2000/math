---
id: THM-4725
title: "Uniform rules for the prime-circulant bijection (tournaments versus even circulant graphs of prime order), resolved: (U) no dilation rule S -> {x : phi(1_S(c_1 x), ..., 1_S(c_k x)) = 1} with fixed integers c_i is a type bijection for all primes: at the infinitely many primes p = 9 (mod 16) at which every c_i is an eighth-power residue (Chebotarev; the exceptional congruence obstruction vanishes once v_2(p-1) >= 3) the rule is constant on the two K_max(p)-invariant tournament types; (B) the signed-set Burnside principle: for any finite group G acting on a set with a free commuting involution, G-orbits of sections equal G-orbits of G-even subsets (flipped-orbit sign), which contains Royle et al., THM-4721 R, and the circulant case with the evenness criterion as the sign of the multiplier -1; (D) the dichotomy for Gamma = AGL(1, Z_(p-1)) = multipliers and power maps x -> a x^j (the symmetries available to a parameter-free multiplicative rule): for p = 3 (mod 4) an even circulant is one with an even number of distance classes, every even circulant is Gamma-even, and the PALEY RULE beta(S) = the member of {backward classes of S relative to QR_p, its complement} with an even number of classes is an explicit canonical Gamma-equivariant type bijection; for p = 1 (mod 4), p >= 13, the flip x -> x chi(x) is a Gamma-odd symmetry of every even circulant with an odd number of non-residue classes, so Gamma-orbits of even circulants outnumber Gamma-orbits of tournaments and no Gamma-equivariant type bijection exists"
status: >
  PROVED (pre-audit): Theorem U (dilation collapse; Chebotarev CITED, the compatibility lemma for v_2(p-1) >= 3 proved
  by a norm computation in Q(zeta_(2^e))); Theorem B (signed-set Burnside, same proof as THM-4721 Theorem R);
  Proposition C (for p = 3 mod 4: even iff an even number of classes; sgn_S(a x^j) = (-1)^{(|S|/2)[a non-residue]});
  Theorem D (the Paley rule is a Gamma-equivariant type bijection for p = 3 mod 4; no Gamma-equivariant type bijection
  for p = 1 mod 4, p >= 13).
  FINITE-EXACT: Gamma-orbit counts for 5 <= p <= 37 (power_equivariance.out, power_equivariance_fast.out): tournaments =
  Gamma-even circulants = the signed-set Burnside sum at every p; all even circulants are Gamma-even at p = 7, 11, 19, 23,
  31 and not at p = 13 (16 of 28 sets), 17 (25 of 117), 29 (4096 of 8128), 37 (65536 of 130816), where the Gamma-orbit
  counts are 4 < 6, 3 < 10, 65 < 125, 673 < 1335; the Hall test for Gamma-equivariant injections passes at every p (Lemma N1), so the obstruction at
  p = 1 mod 4 is purely the orbit count; the collapse primes for c in {2, 3, 5, 6, 7, 10}, {2, 3}, {2, 3, 5}
  (dilation_collapse.out: 73, 89, ... for c = 2; 601 for {2, 3}; 17881 for {2, 3, 5}); THM-4724's dc_family.out.
  CITED: THM-4721 (Corollary C, Theorem R, the evenness criterion), THM-4724 (Theorem W, forcing analysis, Proposition
  D), THM-4531 (Lemma N1), Royle et al. 2023 (Lemma 4.3), Chebotarev's density theorem, Kummer theory
  (Q(zeta_(2^e))(c^(1/2^e)) and its quadratic subextensions).
  OPEN: rules that use the additive structure of F_p (convolution/representation-number rules) at p = 1 (mod 4): neither
  Theorem U nor Theorem D applies to them (HYP-9342).
  Notation: p an odd prime, Z_p^* the multipliers, QR_p the squares, K_max(p) the odd part of Z_p^* (the 2^e-th powers,
  2^e || p-1), classes = the pairs {x, -x}; a circulant tournament of order p is a set S with S u -S = Z_p^* (one element
  of each class); a circulant graph is a symmetric set (a set of classes); types = multiplier orbits (Turner); Gamma = the
  group of maps x -> a x^j, a in Z_p^*, gcd(j, p-1) = 1, acting on subsets of Z_p^* (it commutes with negation).
source: mac-mini-2026-10-20-jcweb (owner's prompt: prove the general non-existence of uniform dilation rules)
depends_on:
  - 01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md
  - 01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md
related:
  - 05-knowledge/hypotheses/HYP-9342-what-a-uniform-rule-can-be-after-thm-4725.md
  - 05-knowledge/hypotheses/HYP-9341-what-remains-open-after-thm-4722.md
note: 05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md
scripts: 04-computation/experiments/jcweb_20261020/{power_equivariance,power_equivariance_fast,dilation_collapse}.py
outputs: 04-computation/experiments/jcweb_20261020/*.out
---

# THM-4725 — uniform rules for the prime-circulant bijection: collapse, signed sets, and the p mod 4 dichotomy

**PROVED (pre-audit).** Note: [jcweb_20261020_uniform_rules_branch_integrality_jacobian_web](../../05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md).

## 1. Dilation rules collapse (Theorem U)

A *dilation rule* is `f_phi(S) = {x in Z_p^* : phi(1_S(c_1 x), ..., 1_S(c_k x)) = 1}` with fixed nonzero integers `c_i` and a fixed Boolean function `phi`; it commutes with multipliers for every `p`. (The derivatives `D_c` of HYP-9341 (c) and THM-4724 Theorem W are the case `k = 2`, `phi = XOR`.)

**Lemma U1.** Let `e >= 3`. For every finite set of nonzero integers `c_1, ..., c_k` there are infinitely many primes `p` with `v_2(p-1) = e` such that every `c_i` is a `2^e`-th power residue modulo `p`, i.e. `c_i in K_max(p)`.

*Proof.* Let `L = Q(zeta_(2^(e+1)), c_1^(1/2^e), ..., c_k^(1/2^e))`. A prime `p` unramified in `L` has `v_2(p-1) = e` iff its Frobenius acts on `zeta_(2^(e+1))` as `zeta -> zeta^(1+2^e)`, and has all `c_i` as `2^e`-th powers iff its Frobenius fixes `Q(zeta_(2^e), c_i^(1/2^e))` pointwise. By Chebotarev it suffices to produce one `sigma in Gal(L/Q)` with both properties, i.e. to show that `1 + 2^e` (which fixes `zeta_(2^e) = zeta_(2^(e+1))^2`) extends from `Q(zeta_(2^(e+1)))` across `M = Q(zeta_(2^e), c_i^(1/2^e))`; the obstruction is `zeta_(2^(e+1)) in M`. By Kummer theory over `Q(zeta_(2^e))` the quadratic subextensions of `M` are `Q(zeta_(2^e))(sqrt(c))`, `c in <c_1, ..., c_k>`, so the obstruction is `zeta_(2^e) = c b^2` with `c in Q^*`, `b in Q(zeta_(2^e))`. Take the norm to `Q(zeta_(2^(e-1)))`: `N(zeta_(2^e)) = zeta_(2^e) * zeta_(2^e)^(1+2^(e-1)) = -zeta_(2^(e-1))`, while `N(c b^2) = c^2 N(b)^2` is a square. For `e - 1 >= 2`, `-zeta_(2^(e-1)) = zeta_(2^(e-1))^(1 + 2^(e-2))` has an odd exponent, and `zeta_(2^(e-1))^odd` is a square in `Q(zeta_(2^(e-1)))` only if `zeta_(2^(e-1))` is, i.e. only if `zeta_(2^e) in Q(zeta_(2^(e-1)))`, which is false (for `e = 3`: `-i` is not a square in `Q(i)`, as `(a+bi)^2 = -i` forces `a^2 = b^2 = 1/2`). ∎ (For `e = 2` the lemma fails, e.g. `c = 2`: `i = 2((1+i)/2)^2`, and indeed `2` is never a fourth power modulo `p = 5 (mod 8)`; this is the `-4 = (1+i)^4` phenomenon of THM-4724 and the reason the earlier search found `D_4` failing exactly at `5, 13, 17`.)

**Theorem U.** No dilation rule is a type bijection for all primes; more precisely, for every rule there are infinitely many primes `p = 9 (mod 16)` at which `f_phi` takes a single value on the two `K_max(p)`-invariant tournament types.

*Proof.* Take `p` from Lemma U1 with `e = 3`. If `S` is `K_max(p)`-invariant then `1_S(c_i x) = 1_S(x)` for every `i`, so `f_phi(S) = {x : phi(b, ..., b) = 1 with b = 1_S(x)}`, which is `∅`, `Z_p^*`, `S` or `-S`; the last two are not symmetric sets, so `f_phi(S) in {∅, Z_p^*}` and is independent of `S`. The `K_max(p)`-invariant tournament sets are the antisymmetric functions on `Z_p^*/K_max(p) = Z_8` (`s(t+4) = 1 - s(t)`), sixteen of them, on which the multiplier group acts through `Z_8` by translation; no proper translation fixes such a function (`d = 4` contradicts antisymmetry, `d in {1, 2}` forces it), so there are exactly two types. A type bijection is injective on types. ∎

The same collapse kills any rule all of whose constants lie in `K_max(p)` at infinitely many `p`; rules whose output depends on `S` only through finitely many dilates are therefore never uniform. (Rules using the additive structure of `F_p` are not touched by this argument.)

## 2. The signed-set Burnside principle (Theorem B)

Let a finite group `G` act on a finite set `Ω~` together with a free involution `ι` commuting with `G`, and let `Ω = Ω~/ι` be the set of pairs. A *section* chooses one element of each pair (a "tournament"); a *subset* `X ⊆ Ω` is a "graph". For `g in G` and a `g`-invariant `X`, an orbit of `g` on `X` is *flipped* if `g^m x = ι x` for the first return `m`; put `sgn_X(g) = (-1)^{#flipped g-orbits in X}`. Fixing a reference section `R`, `r(g, X) := #{pairs in X whose R-member is sent by g to a non-R-member} mod 2` equals the number of flipped `g`-orbits in `X` mod 2 (telescope the flip bits along each orbit), so `sgn_X` is a character of the stabiliser `G_X`. Call `X` *G-even* if `sgn_X` is trivial on `G_X`.

**Theorem B.** `#(G-orbits on sections) = #(G-orbits on G-even subsets)`.

*Proof.* `|Fix_sections(g)| = 2^{c(g)}` if no `g`-orbit on `Ω` is flipped and `0` otherwise (`c(g)` = number of `g`-orbits on `Ω`); `sum_{X in Fix(g)} sgn_X(g) = prod_O (1 + chi_O(g))` over the `g`-orbits `O`, `chi_O(g) = -1` exactly for the flipped ones, which is the same number; Cauchy–Frobenius on both sides. ∎

Instances: `Ω~` = ordered pairs of `[n]`, `ι` = swap, `G = S_n` (Royle et al.; THM-4721 R for `G <= N(P)` on `Fix(P)`); `Ω~ = Z_p^*`, `ι = -1`, `G` = multipliers: a multiplier `a` of even order flips every one of its class-orbits and has `|S|/ord(a)` of them inside `S`, so `sgn_S(a) = (-1)^{|S|/ord(a)}` — the evenness criterion of THM-4721 is the flipped-orbit sign, and Corollary C there is Theorem B for the multiplier group; `G = Gamma` (adding the power maps) gives `#(Gamma-orbits on circulant tournaments) = #(Gamma-orbits on Gamma-even circulants)`, checked for `p <= 31` against the Burnside sum `(1/|Gamma|) sum_psi 2^{c(psi)}[psi flips no class]`.

## 3. p = 3 (mod 4): the Paley rule (Proposition C, Theorem D, first half)

**Proposition C.** Let `p = 3 (mod 4)`, so `QR_p` is a section (a tournament set: `-1` is a non-residue) of odd order `(p-1)/2`. For `psi(x) = a x^j in Gamma` and any `psi`-invariant symmetric `S`: `sgn_S(psi) = (-1)^{(|S|/2) [a is a non-residue]}`. Hence: (i) a circulant graph of order `p` is even iff it has an even number of classes; (ii) every even circulant is `Gamma`-even.

*Proof.* Use `R = QR_p` in the telescoping formula: for `x in QR_p`, `psi(x) = a x^j` lies in `QR_p` iff `a` does (`j` odd), so every pair of `S` contributes the same flip bit `[a not in QR_p]`, and `r(psi, S) = (|S|/2)[a not in QR_p]`. (i) The multiplier `-1` (a non-residue) stabilises every symmetric `S`, so even forces `|S|/2` even; conversely every even-order multiplier is a non-residue (`QR_p` has odd order) and the condition `2 ord(a) | |S|` of THM-4721 reduces to `|S|/2` even because `ord(a)/2` is odd. (ii) is immediate from the formula. ∎

**Theorem D (p = 3 mod 4).** For a tournament set `S` let `X_S` be the set of classes on which `S` disagrees with `QR_p` (the backward classes relative to the Paley tournament), and let `beta(S)` be whichever of `X_S`, `X_S^c` has an even number of classes (exactly one does, as the total `(p-1)/2` is odd). Then `beta` is a `Gamma`-equivariant map from circulant tournaments onto even circulant graphs that is bijective on types.

*Proof.* `X_{aS} = a X_S` for `a in QR_p` and `(a X_S)^c` for `a` a non-residue, and `X_{S^j} = X_S^j` (`QR_p^j = QR_p`); the even-class member is unaffected by complementation, so `beta` commutes with `Gamma`. The image is even by Proposition C (i). On labelled sets `S -> X_S` is a bijection onto symmetric sets and `X -> {X, X^c}` is two-to-one, so `beta(S) = beta(S')` iff `S' in {S, -S}`; hence `beta(S') = a beta(S) = beta(aS)` forces `S' in {aS, -aS}`, i.e. `beta` is injective on multiplier types, and the two type sets have the same size (THM-4721 Corollary C). ∎

`beta` is the construction of THM-4721 Theorem G' ("reverse a Γ-invariant tournament along the graph") for the odd-order multiplier subgroup `QR_p`, extended to the whole multiplier group by the complement. It is definable without parameters in `(Z_p^*, ·, S)` with one quantifier (`x in QR_p` iff `∃y x = y^2`): a uniform natural rule exists for every `p = 3 (mod 4)`. For `p = 7` it sends `QR_7` to the empty graph and the rigid circulant to the complement of `C_7`, as THM-4721 required.

## 4. p = 1 (mod 4): no power-equivariant rule (Theorem D, second half)

**Theorem D (p = 1 mod 4, p >= 13).** `#(Gamma-orbits on even circulants) > #(Gamma-orbits on circulant tournaments)`; hence no `Gamma`-equivariant map from circulant tournaments to even circulant graphs is bijective on `Gamma`-orbits, and in particular no parameter-free rule in the multiplicative structure (which is automatically `Gamma`-equivariant) induces the type bijection.

*Proof.* The flip `psi_0(x) = x^{(p+1)/2} = x chi(x)` lies in `Gamma` (`(p+1)/2` is odd and coprime to `p-1`) and fixes every symmetric set, acting as the identity on residue classes and as `x -> -x` on non-residue classes; so for a symmetric `S`, `psi_0` flips exactly the non-residue classes of `S` and `sgn_S(psi_0) = (-1)^{#non-residue classes of S}`. Take `S = {±x, ±y}` with `x` a residue, `y` a non-residue, `y != ±ix` (possible for `p >= 13`): its multiplier stabiliser is `{±1}` and `|S|/2 = 2`, so `S` is even, but `sgn_S(psi_0) = -1`, so `S` is `Gamma`-odd. `Gamma`-oddness is `Gamma`-invariant, so the `Gamma`-orbits on even circulants strictly contain those on `Gamma`-even circulants, whose number is the number of `Gamma`-orbits on tournaments by Theorem B. An equivariant map induces a map on orbit sets, and a bijection needs equal counts. ∎

Together: **a uniform rule that is natural for the multiplicative structure exists exactly at the primes `p = 3 (mod 4)`**, and no rule with fixed dilation constants works at all primes. The remaining freedom is the additive structure of `F_p` (HYP-9342).
