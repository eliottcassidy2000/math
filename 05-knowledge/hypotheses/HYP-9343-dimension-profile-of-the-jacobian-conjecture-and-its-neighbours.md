---
id: HYP-9343
title: "The dimension profile of the Jacobian conjecture and its neighbours. (a) PROFILE: the truth set of JC is an initial segment {1..N} with N <= 2 (THM-1300 + padding); conjecture N = 2 (JC_2 true), and the structural reason proposed here: every mechanism producing a Keller collision needs a torus orbit over a quotient of dimension >= 2 (THM-4727), which first exists in dimension three. (b) ATLAS of statements with the same profile 'true for n <= 2, false for n >= 3' and the mechanism in each: Markus-Yamabe (torus orbit over a quotient equilibrium, THM-4727 M); tame generation of Aut(A^n) (Jung-van der Kulk amalgam in dimension two, Nagata wild in dimension three: rank one versus rank two, Aut(A^2) acts on a tree, Tame(A^3) on a 2-complex); Dehn's third problem (Bolyai-Gerwien in the plane, the Dehn invariant length (x) angle in space); Banach-Tarski (Isom(R^2) amenable, SO(3) contains free groups); Poincare-Bendixson (planar flows have no recurrence beyond cycles); the two-dimensional cancellation theorem (Fujita-Miyanishi-Sugie) against the open/failed higher cases. (c) CLASSIFICATION: all (1,-1,-2)-equivariant Keller maps of A^3 are the solutions of the twisted planar Jacobian equation {lA, l^2 B} = c l^2 with A in (v,t), B in (v^2,t); conjecturally every solution with l linear is a Gao-Gallagher tangent-sweep map, and solutions with deg l >= 2 exist. (d) PLANAR COUNTEREXAMPLE: it would be the l = 1 case, a fold at infinity; any planar counterexample is torus-rigid (admits no C*-equivariance, THM-1345/THM-4727 D) and has the known degree constraints"
status: >
  (a) OPEN (JC_2). The profile statement is PROVED (padding); the 'reason' is a conjecture supported by THM-4727 (all
  known dimension-three counterexamples are torus-driven) and by THM-1345/THM-1370 (every equivariant Keller map in
  dimension two is an automorphism).
  (b) DICTIONARY: each row is a CITED theorem with the profile; the common mechanism ('a second independent mixed
  quantity appears in dimension three': quotient dimension two, Dehn's tensor, free subgroups, rank two) is an analogy,
  typed as such; the Markus-Yamabe row is verified structurally (THM-4727 M), the others are not claimed to share a
  proof.
  (c) OPEN. The reduction is PROVED (THM-4727 R, E). Numerical search of the twisted equation in bounded degree not yet
  run (the natural next computation; it is bilinear in (A, B) for fixed l).
  (d) OPEN; the constraints are CITED (Moh 1983 degree >= 101 via the sub-125 classification, THM-3550, THM-3544,
  THM-4626, Meng-Yang HC(4)).
source: mac-mini-2026-10-20-jcweb (owner's prompt)
related:
  - 01-canon/theorems/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md
  - 01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md
  - 01-canon/theorems/THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow.md
  - 01-canon/theorems/THM-1370-elliptic-jacobian-conjecture-all-dimensions.md
  - 01-canon/theorems/THM-3543-torus-quotient-ramification-square-no-go.md
  - 01-canon/theorems/THM-4625-tangent-sweep-keller-maps-fibrewise-punctured-lines-hyperelliptic-fibres-and-symmetric-monodromy.md
  - 05-knowledge/hypotheses/HYP-9335-farey-mediant-mandelbrot-dictionary-collatz-angles-and-cross-conjecture-rows.md
note: 05-knowledge/results/jcweb_20261020_uniform_rules_branch_integrality_jacobian_web.md
---

# HYP-9343 — the dimension profile of the Jacobian conjecture and its neighbours

## (a) The profile

A Keller counterexample on `A^m` pads to one on `A^n`, `n >= m`. So JC holds exactly in dimensions `1, ..., N` for some `N in {1, 2}` (THM-1300 excludes `N >= 3`). The owner's question "why would it be true in two and false in all greater" is the question `N = 2`. The proposed reason: a collision of a Keller map is driven, in every known case, by a `C^*`-orbit over a point of a quotient of dimension at least two (THM-4727: a quotient map of generic degree at least two whose only non-etale locus is the zero divisor of the orbit scaling `l`; the fold surface over it covers the fixed hyperplane together with the fixed plane, and the orbit fold `x -> -x` appears where the `x^{-1}`-coordinate vanishes on it), and in dimension two every `C^*`-quotient is a line, on which the twisted Jacobian equation is an ODE with only linear solutions (THM-4727 D, THM-1345, THM-1370). The conjecture is that nothing else can drive a collision: every Keller counterexample is, after composing with automorphisms on both sides, torus-driven in this sense. Hostile control: this is a statement about all Keller maps, and the torus-free planar strata of the repository (THM-2063..., THM-3016/3025, the row programme THM-4380-4438) are where it will be tested.

## (b) The atlas

| statement | true | false / open | mechanism in the first failure |
|---|---|---|---|
| Jacobian conjecture JC(n) | n = 1 (and n = 2 open) | n >= 3 (THM-1300) | a generically 3:1 quotient map on a 2-dimensional quotient plane, etale off the line {l = 0}; the fold surface over that line and the fixed plane together cover the fixed hyperplane 2 + 1 (THM-4727) |
| Markus-Yamabe | n <= 2 (Feßler, Gutierrez, Glutsyuk 1995) | n >= 3 (CEGHM 1997) | torus orbit run over an equilibrium of a 2-dimensional quotient (THM-4727 M) |
| tame generation of Aut(A^n) | n <= 2 (Jung 1942, van der Kulk 1953: amalgamated product) | n = 3 (Shestakov-Umirbaev 2004: Nagata's map is wild) | rank one (tree) versus rank two (the Bisi-Furter-Lamy 2-complex); the Nagata map `(x - 2y w - z w^2, y + z w, z)`, `w = xz + y^2`, is itself equivariant for the mixed weights `(3,1,-1)` with a two-dimensional quotient `(yz, xz^3)` (nagata_weights.out), in the same way the Markus-Yamabe field is built on `x + yz` with weights `(1,2,-1)` |
| Zariski cancellation | n <= 2 (Fujita, Miyanishi-Sugie) | n = 3 in char p (Gupta 2014); open over C | exotic affine structures need a second coordinate direction to hide in |
| Hilbert 14 / LND kernels finitely generated | n <= 2 (Zariski) / LND kernels n <= 3 (Miyanishi) | n >= 5 (Nagata; Daigle-Freudenburg) | — |
| Dehn's third problem | n = 2 (Bolyai-Gerwien) | n = 3 (Dehn 1901) | the Dehn invariant: a tensor of two incommensurable quantities (edge length, dihedral angle) |
| Banach-Tarski | n <= 2 (Isom amenable) | n >= 3 | free subgroups of SO(3) |
| Poincare-Bendixson | planar flows | n >= 3 (Lorenz) | a second direction for recurrence |

Common shape, stated as an analogy: in dimension two the "mixed" data form a one-dimensional object (a single invariant `xy`, a single angle, an abelian isometry group, a tree), and in dimension three a second independent mixed quantity appears (a second invariant, length ⊗ angle, a free group, a rank-two complex), which is exactly the room the counterexamples use. The rank-one/rank-two row is the one most likely to become a theorem: `Aut(A^2)` acts on the Bass–Serre tree of the Jung–van der Kulk amalgam (HYP-9335 (d) already records the Farey-tree reading), `Tame(A^3)` acts on a two-dimensional complex, and the Keller counterexample of THM-1300 sits at the tame/wild boundary (opus-S413 reflection); a precise statement would relate the degree-reduction theory of Shestakov–Umirbaev to the twisted equation (T).

## (c) The twisted planar Jacobian problem

THM-4727 (T): `l {A, B} + A {l, B} + 2 B {A, l} = c`, `A in (v, t)`, `B in (v^2, t)`. Known solutions: the Gao–Gallagher family (`l` linear). Questions: (1) is every solution with `l` linear a tangent-sweep map? (2) are there solutions with `deg l >= 2` (quotient Jacobian `c l^2` with a conic or cubic fold)? (3) for the other mixed weight types in dimension three (`(1,-1,-1)`: W = 2, exponent 1, invariant ring `k[xy, xz]`; `(1,-2,-3)`; `(2,-1,-1)` with two positive weights not covered by Lemma R's normalisation), what are the twisted equations and do they have non-injective solutions? (4) the equivariant automorphism group acts on solutions; normal forms for `l`. The equation is linear in `B` for fixed `(l, A)` and linear in `A` for fixed `(l, B)`, so an alternating linear search in bounded degree is feasible.

## (d) What a planar counterexample would be

The `l = 1` case of (T) is the planar Jacobian equation; a planar counterexample is torus-rigid (no `C^*`-equivariance, THM-1345), of degree at least 101 (Moh's classification below 125, as recorded in THM-3549's corollary), of non-prime degree with pencil height at least 8 (THM-3550), with target pencil of total degree at least 6 (THM-3544), and equivalently a dimension-three Keller map `(P, Q, z)`-like with a surjective (locally nilpotent) dual field (THM-4626) or a counterexample to HC(4) (Meng–Yang, CITED in CORE-PAPERS). In the language of THM-4727 its fold has gone to infinity: the Jelonek set is the only place a non-injectivity can live, and every known constraint is a constraint on that set.
