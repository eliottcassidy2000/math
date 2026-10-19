# Automata, arithmetic lifting, rational cycles, and reconstruction geometry

Research session of 19 October 2026. This map covers the cellular-automata/Gilbreath request and its broader rational-dynamics and geometry extension. It is a navigation and status document, not a claim that a major open conjecture has been settled. Historical novelty has not been established for the individual ingredients or assembled theorems.

## Read first

The strongest directly targeted advance is a denominator-based rational-cycle enumeration theorem: for x²+c, no rational cycle of length greater than three exists when the reduced denominator of c is at most 1,000,000, with no numerator cutoff. A proved finite algorithm enumerates cycles of length at least three across all numerators at each fixed point denominator; the global all-denominator conjecture remains open.

For the Gilbreath-adjacent size-six automaton, all finite extinction delays are even. A proposed logarithmic excess and ratio-to-one prediction is refuted by explicit regular-language contributions. The prime-number Gilbreath conjecture remains open. These are exact finite-strip and alphabet statements, not an iid model substituted for the primes.

The broader work identifies actual maps and information lost under them: integral reduction loses lifting obstructions; boundary restriction loses polynomial classes; tangent matroids lose conductor information; Hilbert ranks can lose individual branch conductors; critical-group quotients can lose a proposed parity observable; moments can recover an integer set only when enough boundary data are retained.

## 1. Automata: reset, absorption, and the nonlinear reachability obstacle

[Finite entry and exact mixing](../gilbreath_binary_entry_20261019/BINARY_ENTRY_AND_EXACT_MIXING.md).
For width F, a length-F binary boundary word universally lowers a nonbinary maximum exactly when its first bit is one. Repetition gives a sharp constant-word bound; random input gives an explicit negative-binomial domination. Once binary, the F-step Pascal control matrix gives exact uniform mixing under fair input. Fixed forcing remains bijective on binary states, so this is not synchronization to one state.

[Alphabet-three size-six analysis](../gilbreath_alphabet_three_size_six_20261019.md) and [all finite delays even](../gilbreath_even_delay_20261019/GILBREATH_EVEN_DELAY_THEOREM.md).
A shielding lemma reduces the random boundary problem to an exact zero-run language. Seven all-length regular-language certificates already force a limiting lower ratio at least 285677/262144, contradicting the former ratio-one prediction. A seven-state pair-language invariant forbids every odd finite delay.

[Precise paired-halving stopping point](../gilbreath_paired_halving_frontier_20261019/GILBREATH_PAIRED_HALVING_STOPPING_NOTE.md).
An exact phase criterion characterizes when one boundary bit annihilates all remaining 2s. The all-length reachable-state version would imply paired halving, including the necessary age-zero-parent case. Correct parity trace and the proved pair invariant alone are insufficient: explicit hostile states share the complete parity trace of genuine reachable states. Exact infinite-sequence tests through front length 14 support two stronger hypotheses, but do not prove them.

Next gate: control the high bit of reachable age vectors. The paired-halving law and sharp dyadic horizon remain open.

## 2. Tournament statistics as an integral divided-power module

[Integral source-free basis and Hasse contractions](../tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md).
Every integer tournament statistic has a unique integral expansion in source-free induced-subtournament counts. This unimodular basis diagonalizes subset deletion, giving cokernel torsion and prime-adic obstruction depths from binomial valuations. Ordinary repeated deletion contains a factorial; it cannot be identified with divided-power deletion after reduction modulo a prime.

[Double-path selector and forest floor](../tournament_double_path_selector_20261019/DOUBLE_PATH_SELECTOR.md).
The eight-vertex double-path tree has 64 universally odd labelled orientations and 64 nonconstant orientations sharing one parity function. Every oriented forest has no nonempty primitive coefficients in orders one through six; order seven is sharp. This does not assert new forbidden Hamiltonian-path counts.

Next gate: exploit the explicit primitive coordinates for other orientation families, retaining same-order cores rather than misreading the expansion as graph reconstruction.

## 3. Wronskians, Frobenius carries, and unique compatible lifts

[Wronskian and scalar Frobenius packet](../wronskian_frobenius_20261019/).
For scalar Laurent polynomials, compatibility A(X)^p=A(X^p) modulo p^m has a complete support criterion and a unique coefficient lift. The lift is multiplicative but generally not additive. Hasse jets retain different valuation budgets from ordinary derivatives; sorted row/column budgets sharpen the determinant estimate.

Finite-logarithm evaluation at X=1 can vanish at a Wieferich prime while a polynomial coefficient obstruction remains nonzero. Scalar evaluation therefore cannot certify an entire Frobenius lift. No finiteness claim for Wieferich or Wilson primes follows.

The supplied Barreto paper's no-Landau–Siegel-zero conclusion remains a source claim, not a result certified by this session.

## 4. Ternary period relations: torsion, reconstruction, and an exceptional prime

[Ternary Frobenius reconstruction packet](../ternary_frobenius_reconstruction_20261019/).
The Ilin–Ono characteristic-zero setting admits new finite-characteristic questions. Universal interior families occur in degrees 3p+3 and 2p+12 under their stated prime restrictions. In odd degree n, the supported integral p-primary torsion is exactly floor(log_p n) copies of Z/p for p≥5. Every supported modular class is obstructed at its first support-preserving lift.

The tempting universal minimum-degree prediction is false: p=29 has an interior degree-42 class, below both proposed families. Its first lift is obstructed by an exact augmented-rank calculation. Full cyclic corrections and support-preserving corrections are different problems; the universal odd-supported theorem does not resolve every full-carrier lift.

[Integral lifting synthesis](INTEGRAL_LIFTING_CONNECTIONS.md) gives a common matrix obstruction map and an exact prime-power carrier. The torsion-order jump at odd degree n is p precisely when n is a power of p≥5. Its Dirichlet series is the logarithmic zeta derivative with the 2- and 3-power terms removed. This supplies an explicit structural connection, not a new estimate for primes or zeta zeros.

Next gate: classify interior exceptions and obtain uniform dual obstructions for full cyclic corrections.

## 5. Rational quadratic cycles: denominator signatures and local shells

[Rational denominator-cycle packet](../rational_quadratic_denominator_cycles_20261019/).
Bad-prime pair-sum valuations turn periodic orbit enumeration into finitely many signed signature seeds. Given the factorization of b, enumeration of cycles of length at least three across all reduced parameters c=C/b² is fixed-parameter tractable in the number of prime divisors of b, with polynomial dependence on log b at fixed support size. Exact sign reconstruction and orbit closure are essential: locally admissible triangle data can escape.

The complete reduced-parameter-denominator≤1,000,000 theorem combines cited denominator obstructions with 101 remaining point-denominator cases b≤1000; all 62 cycles found of length at least three have length three. Numerators are unrestricted. This is not a literature-wide computational record claim.

For c=−29/16, a 43-adic resonance theorem proves that the three residue balls around the known rational cycle contain no other Q_43-periodic orbit. Apparent period-six cycles in finite quotients lie in shells whose periods grow under lifting. New finite period-six cycles can appear in deeper shells, so simply counting one lift layer would be misleading. Two additional critical-cycle congruence families give infinite local exclusions.

Next gate: a uniform global obstruction across unbounded denominator support. Local valuation signatures alone do not provide it.

## 6. The parabolic three-cycle and the real seventh cyclotomic field

[Parabolic/cyclotomic supplement](../parabolic_cyclotomic_cycles_20261019/PARABOLIC_CYCLOTOMIC_THREE_CYCLE.md).
The rational three-cycle parameterization has an exact square formula for its distance below c=−7/4. Its boundary cubic has discriminant 49 and roots in the real seventh cyclotomic field. An explicit affine map intertwines the relevant three-point permutations, reversing the Galois generator; it is not a global polynomial conjugacy and does not preserve multipliers.

There are infinitely many rational attracting three-cycles approaching the parabolic parameter. Thus no uniform real multiplier gap can prove the rational period bound. Conversely, the rational denominator theorem excludes every rational neutral cycle beyond the two stated neutral fixed points. This excludes neutral cycles, not hyperbolic higher-period cycles.

Next gate: retain cubic splitting data and denominator information when using real parameter or multiplier coordinates.

## 7. Geometry from the supplied S6 paper

[Broad geometry packet](../broad_geometry_reconstruction_20261019/).
The supplied alpo.ge/s6.pdf is the previously encountered 108-page source. Its global complex-S6 headline is not established here.

Three scoped calculations are available:

1. Seifert gauge classes form Z plus a finite gcd torsion factor; an exact local arithmetic formula counts admissible marked classes at a fixed clutching invariant.
2. Rank-two nilpotent monodromy gives explicit exterior Smith forms. Under the stated free/proper translation-fan hypotheses, identical middle homology can coexist with unbounded central cusp Euler characteristic. Those quotient hypotheses are retained.
3. The conductor's smoothing line bundle depends on holonomy, not only degree. The displayed toric transition maps preserve the smoothing coordinate, providing a trivial frame under the local quotient hypotheses. This does not prove a global proper/Hausdorff compact quotient.

## 8. Independently developed connections beyond the prompt themes

These connections were pursued because they supplied usable intermediate structure:

- **Interpolation and conductor extrema.** [Sharp endpoint classification](../conductor_extrema_20261019/conductor_extremal_classification_session7.md). For r line branches with tangent matroid U(3,r), the seminormal defect ranges between explicit sharp endpoints. Rational power-curve and conic configurations attain them. Maximum equality for r≥6 characterizes smooth conics. The minimum grows as r^(3/2), the maximum as r². Even identical full Hilbert functions can have different branch-conductor vectors.
- **Critical groups and parity repair.** A proposed cut-parity observable descends to the critical group exactly when a mod-two Laplacian condition holds. For bipartite graphs this requires even weighted degrees. Minimal edge-copy repairs are T-joins; the exact core has enumerated minimal repairs. The changed lattice is not automatically an admissible physical dynamics.
- **Finite differences and Padé reconstruction.** An integer set with at most k runs is uniquely determined by its first 2k moments. Its boundary polynomial quotient follows from finite differences and a rational generating function. Integer endpoints, interlacing, and the zero endpoint must be checked. Discrete sharpness is only proved here for k=1,2.
- **Bockstein lifting and Kummer carries.** The same explicit integral obstruction framework explains why a modular relation may fail to lift, while factorial normalization changes the depth. The carrier and allowed corrections are specified for each application.
- **Cyclotomic parameter ramification.** The seventh cyclotomic cubic gives an exact boundary map for rational marked three-cycles, connecting Galois permutations, square-root bifurcation, and bad-prime repulsion without equating them.
- **Exterior monodromy versus fan index.** A Smith invariant visible in fan size is absent from middle exterior torsion. This supplies an explicit example of an invariant losing precisely the data needed for a geometric inference.

## 9. Hypothesis board and recommended continuation

Proved within stated scope: finite-strip entry/mixing; all finite size-six delays even; integral source-free tournament expansion; unique scalar Frobenius lifts; odd-supported norm torsion; denominator-bounded rational period exclusion and enumeration algorithm; local resonance exclusion; conductor extrema; bounded-run moment uniqueness.

Refuted: logarithmic size-six excess and its ratio-one consequence; the proposed minimum-degree formula at p=29; conductor determination from tangent matroid; branch-conductor determination from full Hilbert ranks; uniform real multiplier separation near the parabolic three-cycle.

Still open: prime Gilbreath; all-length paired halving and sharp age horizon; universal rational quadratic period bound; intermediate conductor spectra; uniform full-carrier ternary lift obstruction. Collatz, planar Jacobian, Hadwiger–Nelson, Goldbach/twin-prime adjacent targets, and Wilson-prime completeness are not resolved by these connections.

The most concrete next proof targets are the reachable high-bit automaton invariant, a uniform dual full-cyclic lifting obstruction, and an arithmetic compatibility constraint that eliminates denominator-signature seeds before orbit closure. The conductor intermediate spectrum is a smaller independent classification problem with exact endpoint control already in hand.

## Sources and publication

Primary supplied sources: [Ilin–Ono](https://arxiv.org/pdf/2610.19309), [the S6 PDF](https://alpo.ge/s6.pdf), and the user-provided Barreto text. Rational-cycle background: [Eliahou–Fares](https://arxiv.org/pdf/1903.08865) and [Morton–Raianu](https://arxiv.org/pdf/2105.07435). Individual packets include precise source and replay information.

All results are ordinary commits on main. The [publication ledger](PUBLICATION_LEDGER.md) records the verified checkpoints. Scripts, exact receipts, limitations, and history-preserving corrections are retained for future work.
