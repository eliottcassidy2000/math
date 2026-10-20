---
id: THM-4727
title: "Mixed-weight torus reduction of Keller maps and the dimension profile of the Jacobian conjecture. (R) Ramification exponent: a polynomial map F of A^n equivariant for a C*-action with one positive weight (source weights (1,-w_2,...,-w_n), target weights (1,-w'_2,...,-w'_n), orbit coordinates x_1, X_1, orbit scaling l = X_1/x_1 an invariant) satisfies det JF = x_1^{W-W'} l^{1-W'} det J(Fbar) on the quotient planes, W = sum w_i; so a Keller map forces W = W' and det J(Fbar) = c l^{W-1}: the quotient map is Keller for the degenerate volume form l^{1-W} dv (a 'twisted planar Jacobian equation' {V,T} = c l^2 for weights (1,-1,-2)); F and Fbar have the same generic degree, Fbar is etale off {l = 0}, and over {l = 0} each torus orbit is collapsed into the fixed hyperplane X_1 = 0 by x_1 -> (A x_1^{-w_2}, B x_1^{-w_3}), gcd(w_2,w_3)-to-one generically and w_3-to-one on {A = 0} (the Alpoge collision). (E) Every tangent-sweep Keller map of Gao-Gallagher (THM-4625; THM-1300 is the first member) is (1,-1,-2)-equivariant with LINEAR orbit scaling l = gamma_0 + alpha v + beta t, so the whole known dimension-three Keller counterexample corpus is the W = 3, deg l = 1 stratum of (R): quotient Jacobian c l^2 (THM-3543's square for the first member). (D) Dimension two is the W - 1 = 0 case for weights (1,-1), and for weights (1,-w) the quotient is a line on which the integrality of Y = phi(t)/(x^w l^w) forces l constant and phi linear: hyperbolic planar equivariant Keller maps are linear (THM-1345 recovered), while the first dimension with a two-dimensional mixed quotient is three. (M) The same anatomy carries the Markus-Yamabe counterexample: the Cima-van den Essen-Gasull-Hubbers-Manosas field X = (-x + z d^2, -y - d^2, -z), d = x + yz, is equivariant for weights (1, 2, -1) with quotient plane (xz, yz^2), its Keller companion x + (X + x) = (x + z d^2, y - d^2, z) is a triangular automorphism, and its escaping solution (18 e^t, -12 e^{2t}, e^{-t}) is a torus orbit over the EQUILIBRIUM (18, -12) of the quotient field (a', b') = (-2a + (a+b)^2, -3b - (a+b)^2): Markus-Yamabe (true for n <= 2, false for n >= 3) fails by a torus orbit over a 2-dimensional quotient exactly as the Jacobian conjecture does. (P) The truth set of JC is an initial segment {1, ..., N} (padding), N = 2 is open and N >= 3 is excluded; HYP-9343 records the profile atlas (Markus-Yamabe, tame generation, cancellation, Hilbert 14 / LND kernels, Dehn's invariant, Banach-Tarski, rank-one versus higher-rank automorphism groups) and what a planar counterexample would have to be (the l = 1 case of (R): a fold that has gone to infinity)"
status: >
  PROVED (pre-audit): Lemma R (ramification exponent; block-triangular Jacobian in orbit coordinates; verified
  symbolically for eight weight configurations, including W != W', ramification_exponent.out); Corollary R' (Keller forces
  W = W' and det J Fbar = c l^{W-1}; Fbar etale off {l = 0}, equal generic degrees, the orbit collapse over {l = 0}); Proposition E (the
  Gao-Gallagher family is (1,-1,-2)-equivariant with linear l; read off THM-4625's displayed form); Proposition D (planar
  hyperbolic rigidity by the integrality-degree argument; THM-1345 CITED as the original); Proposition M (Markus-Yamabe
  anatomy: equivariance, nilpotent JX + I, the automorphism companion, the quotient field and its equilibrium; all
  verified symbolically, torus_quotient.out).
  FINITE-EXACT: the Alpoge quotient map (BC, AC^2) with Jacobian 2(t + 3v - 2)^2, fibre over (0,0) = {l = 0} u {(0,0)},
  descended collision (torus_quotient.out, agreeing with THM-3543's (6)-(9)); the planar search finding no nonlinear
  weighted-homogeneous Keller pairs for weights (1,-1), (1,-2), (2,-3), (1,-3) within the degree box.
  CITED: THM-1300 (the map), THM-1305 (equivariant normal form), THM-1345 (planar equivariant JC), THM-1370 (elliptic JC
  all dimensions), THM-3543/3549 (the quotient of the fixed map, its square Jacobian, compression and correction no-gos),
  THM-4625 (tangent-sweep family; Gao arXiv:2608.00222), THM-4620/4626 (dual fields; JC_2 iff dimension-three surjective
  dual fields), Cima-van den Essen-Gasull-Hubbers-Manosas, Adv. Math. 131 (1997) (the Markus-Yamabe counterexample),
  Fessler 1995 and Gutierrez 1995 (planar Markus-Yamabe), Jung 1942 / van der Kulk 1953, Shestakov-Umirbaev 2004,
  Bisi-Furter-Lamy 2014 (tame automorphisms of A^3 acting on a 2-dimensional complex), Dehn 1901, Banach-Tarski 1924.
  NOT PROVED: any classification of the solutions of the twisted planar equation {lA, l^2 B} = c l^2 beyond the known
  families; any statement about non-equivariant Keller maps; JC_2.
  Notation: weights act by lambda.(x_1, ..., x_n) = (lambda x_1, lambda^{-w_2} x_2, ...); invariants v_i = x_1^{w_i} x_i;
  Fbar the induced map on Spec of the invariant rings (polynomial planes when the invariant ring is free, as for
  (1,-1,-2)); Keller = det JF a nonzero constant.
source: mac-mini-2026-10-20-jcweb (owner's prompt: classification of counterexamples in higher dimensions, why dimension two would be true and all higher false, a structural analogue from another field, the possibility of a planar counterexample)
depends_on:
  - 01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md
  - 01-canon/theorems/THM-3543-torus-quotient-ramification-square-no-go.md
  - 01-canon/theorems/THM-4625-tangent-sweep-keller-maps-fibrewise-punctured-lines-hyperelliptic-fibres-and-symmetric-monodromy.md
related:
  - 01-canon/theorems/THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow.md
  - 01-canon/theorems/THM-1370-elliptic-jacobian-conjecture-all-dimensions.md
  - 01-canon/theorems/THM-3549-torus-quotient-correction-no-go.md
  - 01-canon/theorems/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md
  - 05-knowledge/hypotheses/HYP-9343-dimension-profile-of-the-jacobian-conjecture-and-its-neighbours.md
note: 05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md
scripts: 04-computation/experiments/jcweb_20261020/{torus_quotient,ramification_exponent,nagata_weights}.py
outputs: 04-computation/experiments/jcweb_20261020/{torus_quotient,ramification_exponent,nagata_weights}.out
---

# THM-4727 — mixed-weight torus reduction of Keller maps and the dimension profile

**PROVED (pre-audit).** Note: [jcweb_20261020_uniform_rules_branch_integrality_jacobian_web](../../05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md). The repository's state: the Jacobian conjecture is false in every dimension `n >= 3` (THM-1300, July 2026; padding) and open in dimension two. This theorem isolates the mechanism of every known counterexample as a torus phenomenon over a two-dimensional quotient, shows that the same mechanism is empty in dimension two, and matches it with the Markus–Yamabe conjecture.

## 1. The ramification exponent (Lemma R)

Let `C^*` act on the source `A^n` with weights `(1, -w_2, ..., -w_n)`, `w_i >= 1`, and on the target with weights `(1, -w'_2, ..., -w'_n)`, and let `F` be a polynomial map with `F(lambda . x) = lambda . F(x)`. On `x_1 != 0` use the orbit coordinates `(x_1, v_2, ..., v_n)`, `v_i = x_1^{w_i} x_i` (the `v_i` generate the invariants), and likewise `(X_1, V)` on the target. Equivariance says `X_1 = x_1 l(v)` for an invariant `l` (the orbit scaling) and `V = Fbar(v)` (the quotient map). The Jacobians of the coordinate changes are `x_1^W` and `X_1^{W'}` (`W = sum w_i`, `W' = sum w'_i`), and the map `(x_1, v) -> (x_1 l(v), Fbar(v))` has block-triangular Jacobian `l(v) det J Fbar(v)`. Hence, as rational functions,

    det JF = x_1^{W - W'} l(v)^{1 - W'} det J Fbar(v).                                          (R)

**Corollary R'.** If `F` is Keller (`det JF = c != 0`), then `W = W'` (otherwise a power of `x_1` survives in a function of `v`) and

    det J Fbar = c l^{W - 1}:

the quotient map is a Keller map for the degenerate volume form `l^{1-W} dv_2 ∧ ... ∧ dv_n`, etale exactly off the divisor `{l = 0}` (where its Jacobian vanishes to order `W - 1`). Off `{l = 0}` every torus orbit is mapped bijectively onto a torus orbit (`X_1 = x_1 l` is linear in `x_1`), so `F` and `Fbar` have the same generic degree (three for THM-1300: THM-1315's generic 3:1 and THM-3549's cubic fibre equation are the same number), and the non-injectivity of `F` away from `{l = 0}` is exactly that of `Fbar`. Over `{l = 0}` each orbit is collapsed into the fixed hyperplane `X_1 = 0`: with `F = (x_1 l, A_2/x_1^{w_2}, ..., A_n/x_1^{w_n})` the orbit map there is `x_1 -> (A_i x_1^{-w_i})_i`, which is `gcd(w_2, ..., w_n)`-to-one where all `A_i` are nonzero and `gcd` of the remaining weights where some vanish — for weights `(1,-1,-2)` two-to-one on the curve `{l = 0} ∩ {A = 0}`, which is where THM-1300's `lambda -> lambda^2` fold lives, the third preimage coming from the fixed plane `x_1 = 0`, which also maps into `X_1 = 0`. Checked symbolically for eight weight configurations (`ramification_exponent.out`).

For weights `(1,-1,-2)` the invariant rings are polynomial planes (THM-3543 §1), `W = 3`, and with `F = (B/x^2, A/x, x l)` (`A, B, l` polynomials in `v = xy`, `t = x^2 z`, integrality `A in (v,t)`, `B in (v^2, t)`) the quotient map is `(V, T) = (lA, l^2 B)` and the Keller condition is the **twisted planar Jacobian equation**

    {lA, l^2 B} = c l^2,   equivalently   l {A, B} + A {l, B} + 2B {A, l} = c.                     (T)

The ordinary planar Jacobian equation is the case `l = 1`.

## 2. The known counterexamples are the deg l = 1 stratum (Proposition E)

THM-4625 displays the Gao–Gallagher tangent-sweep family as `F = ((q(w) + gamma w)/(gamma x)^2, (p(w) + 2 gamma)/(gamma x), gamma x)` with `gamma = gamma_0 + alpha xy + beta x^2 z` and `w = gamma(1 + xy)`: this is exactly `(B/x^2, A/x, x l)` with `l = gamma_0 + alpha v + beta t` linear, `A = (p(w) + 2gamma)/gamma`, `B = (q(w) + gamma w)/gamma^2` (Gao's side conditions are the integrality of `A` and `B`). So every known dimension-three Keller counterexample, THM-1300 (`p = -3w^2 + 4w`, `l = 2 - 3v - t`) and Gao's degree-four map included, is `(1,-1,-2)`-equivariant with linear orbit scaling, and by Corollary R' its quotient Jacobian is `c l^2`; THM-3543's `2(2 - 3v - t)^2` is the first instance. The triple collision of THM-1300 descends to the fibre `{l = 0} ∪ {(0,0)}` of the quotient map over the origin (`torus_quotient.out`, agreeing with THM-3543 (8)–(10)): one preimage on the fixed plane `x = 0`, the orbit over `(-3/2, 13/2) in {l = 0} ∩ {A = 0}` folded two-to-one by `lambda -> lambda^2`; generic targets have three preimages, one per sheet of the quotient map, in accordance with Corollary R'.

Classifying the `(1,-1,-2)`-equivariant Keller maps is therefore the planar problem (T) with its integrality conditions; the Gao–Gallagher family is the solution set found so far, with `l` linear. Whether (T) has solutions with `deg l >= 2`, or non-sweep solutions with `l` linear, is open (HYP-9343 (c)).

## 3. Dimension two is empty (Proposition D)

For weights `(1, -w)` on `A^2` the quotient is the line `t = x^w y`, `W = w`, and a Keller map has the form `F = (x l(t), phi(t)/(x^w l(t)^w))` with `phi' = c l^{w-1}`. Integrality of the second coordinate forces `phi in (t)` and `l^w | phi`; writing `phi = l^w psi` gives `w l' psi + l psi' = c`, whose top-degree coefficient `(w deg l + deg psi) lead(l) lead(psi)` is nonzero unless `deg l = 0`; then `psi` is linear and `F` is linear. (For `w = 1` the exponent `W - 1` is already zero: the quotient Jacobian is `phi'` itself and no fold is possible.) This is THM-1345's hyperbolic case re-derived from (R); the search in `torus_quotient.out` finds no nonlinear weighted-homogeneous Keller pair for weights `(1,-1), (1,-2), (2,-3), (1,-3)`. Together with THM-1370 (positive weights: automorphisms in every dimension) the equivariant Jacobian conjecture is true in dimension two for every torus, and false in dimension three for the first mixed type whose quotient is a plane. The mechanism needs two things that first coexist in dimension three: a mixed-sign action (so that `l` can vanish) and a quotient of dimension at least two (so that the twisted equation (T) is not an ODE).

## 4. Markus–Yamabe has the same anatomy (Proposition M)

The Markus–Yamabe conjecture (a polynomial vector field whose Jacobian has eigenvalues of negative real part everywhere has a globally attracting equilibrium) is true for `n <= 2` (Feßler, Gutierrez 1995) and false for `n >= 3`: the field `X = (-x + z d^2, -y - d^2, -z)`, `d = x + yz`, has `JX = -I + N` with `N^3 = 0` and the solution `(18 e^t, -12 e^{2t}, e^{-t})`. `torus_quotient.out` verifies: `X` is equivariant for weights `(1, 2, -1)`; the invariant ring is the plane `k[a, b]`, `a = xz`, `b = yz^2`; the Keller companion `H = (x + z d^2, y - d^2, z)` (`det JH = 1`) is a triangular automorphism (it preserves `d` and `z`), so the Markus–Yamabe failure is not a Keller collision; and the escaping solution projects to the equilibrium `(a, b) = (18, -12)` of the quotient field `(a', b') = (-2a + (a+b)^2, -3b - (a+b)^2)`: the solution is the torus orbit over that equilibrium, moving along the orbit with `lambda = e^t`. In dimension two a mixed action `(1,-1)` has the line `t = xy` as quotient; an equivariant field `(x p(t), y q(t))` escapes along an orbit only over a zero of `p + q` with `p != 0`, and the planar theorem says the eigenvalue hypothesis forbids it. So the two conjectures with the profile "true for `n <= 2`, false for `n >= 3`" fail, in their first known counterexamples, by the same device: a torus orbit over a two-dimensional quotient — a fold of the orbit map (Keller) or a run along the orbit over a quotient equilibrium (Markus–Yamabe).

The Nagata automorphism `(x - 2yw - zw^2, y + zw, z)`, `w = xz + y^2`, the wild generator of dimension three (Shestakov–Umirbaev), is equivariant for the mixed weights `(3, 1, -1)` with quotient plane `(yz, xz^3)`, on which it descends to the planar automorphism `(a, b) -> (a^2 + a + b, b + a^2 - (a^2 + a + b)^2)` (`nagata_weights.out`). So the three dimension-three phenomena with the profile "trivial in dimension two, nontrivial from three" — the Keller collision, the Markus–Yamabe escape, the wild automorphism — are all mixed-weight torus phenomena over a quotient plane.

## 5. The dimension profile (Proposition P) and the planar question

If `F` is a Keller counterexample on `A^m` then `F × id` is one on `A^n` for every `n >= m`, so the set of dimensions in which the Jacobian conjecture holds is an initial segment `{1, ..., N}`; THM-1300 gives `N <= 2`, and `N = 2` is the open case. The reductions that raise dimension (Bass–Connell–Wright, de Bondt–van den Essen; THM-1430) are compatible with this and say nothing about `N = 2` versus `N = 1`. HYP-9343 records the profile atlas and the analogies. For the planar case, (R) says what a counterexample would have to be: an `l = 1` solution of the planar Jacobian equation, i.e. a map whose fold divisor has left the finite plane; the known constraints (Moh: degree `>= 101`; THM-3550: no prime degree, pencil height `>= 8`; THM-3544: target pencil degree `>= 6`; THM-4626: equivalently a dimension-three Keller map with a surjective dual field; Meng–Yang: equivalently HC(4), CITED in CORE-PAPERS) are all constraints on where that fold can hide at infinity.
