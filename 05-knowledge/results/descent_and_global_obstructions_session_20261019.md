# Descent, global periods, and inverse data: extended research session

19 October 2026. This session revisits the repository's aliquot, Collatz, tournament, prime-generation, poset, and planar inverse-map work to find usable connections across their existing proof boundaries. The results are scoped theorems and explicit reductions. None settles general Collatz convergence, Hadwiger–Nelson, the planar Jacobian conjecture, the general Gold/balance conjectures, or a prime-infinitude conjecture. Historical priority is not asserted.

The [preceding parabolic session](parabolic_scaling_session_20261019.md) supplies the near-cusp clocks, guarded arithmetic operations, and source-sensitive Collatz compactification. The present session moves from those interfaces into exact changing-core descent, stable global periods, and inverse-sheet geometry.

## 1. An old aliquot frontier closes, then strengthens to every cycle length

The [all-exponent three-cycle exclusion](dyadic_aliquot_three_cycle_exclusion_20261019/DYADIC_THREE_CYCLE_EXCLUSION.md) closes the earlier A3/D3 dyadic frontier using an interval trap, a residue condition, and exact valuations. It applies to aliquot cycles whose every member has the form 2^a p, with a>=1 and p an odd prime.

The [global cofactor classification](dyadic_cofactor_global_descent_20261019/COFACTOR_MAP_GLOBAL_CLASSIFICATION.md) goes further. For positive even n=cp, with c its exact power of two and p odd, define

    F(n)=(c-1)p+2c-1.

Every orbit reaches a fixed point c(2c-1) or one of the two explicit cycles

    8 -> 22 -> 14 -> 10 -> 8,
    16 -> 46 -> 26 -> 16.

The proof finds a smaller actual iterate within six steps outside fixed points and five exceptional starting values. Its cases cover all dyadic exponents. Six steps are genuinely needed at n=48.

On the native odd-prime domain, F is exactly the aliquot map. Both nontrivial F cycles leave that domain at a power of two. Consequently no nontrivial aliquot cycle of ANY length can remain entirely in this class; every aliquot trajectory that remains there forever reaches an even perfect number. The odd-prime condition is retained at every step, so this is not a general aliquot-convergence theorem.

The [extensions](dyadic_cofactor_global_descent_20261019/extensions/) retain information beyond eventual convergence:

- If M(n) is the largest value on the F orbit, then limsup (M(n)-3n)/n^(2/3)=3. The exact family n=d²(d-1), d dyadic and at least 8, has M(n)=3d³-2. The exponent 2/3 comes from a cubic input scale and quadratic correction.
- Immediate predecessors of the fixed point 2^A(2^(A+1)-1) correspond exactly to divisors of A. Immediate predecessors of 2^b correspond exactly to proper divisors of b. The inverse identity F(cp)-1=(c-1)(p+2) makes these actual Mersenne-divisor trees, not numerical resemblance.

The [cofactor deformation](cofactor_deformation_classification_20261019/README.md) tests which part of the mechanism survives a changed intercept. For F_t(cp)=(c-1)p+2tc-1, growing inputs whose p+1 has odd part at least 2t+1 descend within two steps. The remaining unbounded exponent lanes are completely resolved at t=2: every orbit reaches a fixed point or one of six listed cycles. This map is explicitly auxiliary, not an aliquot rule.

For every t, a separate exact theorem classifies ALL two-cycles by a Mersenne-gcd divisibility condition. Their cores satisfy c<=2t-2 and d<=8(t-1)^2; the quadratic bound is attained for infinitely many parameters. General t>=3 orbit convergence remains a research hypothesis. This supplies a controlled comparison where the cycle catalogue changes while a proved descent mechanism survives.

## 2. Unbounded Collatz interleaving now has a repetition exclusion

The [interleaved repetition theorem](interleaved_repetition_barrier_20261019/INTERLEAVED_REPETITION_BARRIER.md) permits arbitrarily many expansion/reset episodes in the specified expansion-first bank language. A macroscopic consecutive power of an expansion word, anywhere in a first-crossing word, forces strict descent at the endpoint itself.

The proof transports actual source height through the earlier context instead of deleting that context. A normalized carry bound e/S<N/6, an exact rational logarithm bracket covering the entire short range N<=2^52, and a cited explicit logarithmic-form estimate cover all lengths with a proved overlap. A finite census does not supply the universal quantifier.

The [matched-block continuation](interleaved_repetition_barrier_20261019/extensions/MATCHED_BLOCK_EXTENSION.md) broadens the mechanism to two equal long parity factors at arbitrary bit positions. A non-descending endpoint then forces the two actual orbit values to coincide. A noncontracting separator or an observed mismatch in their later aligned parity bits excludes that coincidence and proves descent.

A quantified low-factor-complexity hypothesis initially forced an actual internal periodic orbit at any remaining non-descending endpoint. The next [image-guard checkpoint](interleaved_repetition_barrier_20261019/image-guards/IMAGE_GUARD_DESCENT.md) removes that alternative for an unbounded upper-convergent family by retaining the past as well as the future: V outputs are 82 modulo 243, whereas VF outputs are 1 modulo 243. They cannot be an exact return even when every long matched future has identical continuations and a contracting separator. An exact rational-log certificate specifies a word of more than 10^42 bits without materializing it.

The [all-greedy closure](interleaved_repetition_barrier_20261019/all-greedy-closure/ALL_GREEDY_ENDPOINT_DESCENT.md) then settles the ENTIRE specified greedy V/F family. The decisive link is reverse-injective coding: V and VF have disjoint output images, so any repeated actual checkpoint would uniquely decode backwards to a forbidden return to the original source. All completed-iteration values are therefore distinct. The mechanical reset-symbol sequence has at most k+1 length-k factors; a finite pigeonhole argument gives a long repeated future, and the matched-block theorem now forces strict endpoint descent.

This also works for every fixed nonempty expansion template U over {K,V}, with greedy F resets after each copy and at least two copies. The arbitrary-template one-copy case is explicitly outside that corollary. This closes the endpoint goal behind the prior typed-factor target by proving an anchored-or-different-image disjunction; the literal stronger combinatorial hypothesis is not needed or silently asserted.

Other bank families can still have high factor complexity, overlapping terminal images, or no suitable source-anchored coding. Their actual endpoint equations, the full realized limit boundary, and arbitrary integer entry remain open. The proof does not promote a formal symbolic path or a short program to a lawful non-descending source.

## 3. Local stress plus a global period strengthens the finite-field graph bound

The [balanced-digit stress](balanced_digit_section_defect_20261019/README.md) starts with the norm Q(x,y)=x²+xy+y² over F_29. A section from precision one to precision two has edge errors satisfying a twelve-term mixed stress with coefficients 1,3,9. Under equality in the disjoint-cycle count, balanced-digit uniqueness forces binary straight seams. Eight offset-independent periods reject every orientation. This raises the full section-defect lower bound from 435 to 440.

The [sparse-error continuation](section_defect_stability_20261019/README.md) rules out the next possible error count on the three-direction subsystem. One extra failed edge gives a dipole. Cube incidence, a short-column binary-normalization lemma, and an explicit local motion remove that dipole while preserving the full carry coset, reducing to the forbidden seam configuration. The resulting bounds are 89 errors on the triple and 445 on the full graph.

The exact 88-error hostiles satisfy all the compact stress checks and straight-cycle sums, yet retain a nonzero global period. This identifies why local compatibility alone cannot establish a section. Neither 445 nor 89 is claimed sharp. These are prescribed-section defect bounds, not new Euclidean chromatic bounds.

The [two-extra-error analysis](section_defect_two_error_geometry_20261019/README.md) further proves that an 89-error section cannot put all three exceptional failures on one straight cycle. It certifies the complete 1,134-geometry universe and leaves only two distinct double-failure cycles. A small-rank local-motion criterion excludes separated configurations. A path-complement addendum removes its apparent same-direction sparsity exception while preserving the exact dipole and mixed cube. The explicit disjoint-support hostile has a six-dimensional coupled kernel, with two additional connecting modes: disjoint cube supports alone do not justify the four-dimensional correction assumption.

The numerical lower bound remains 445. The next obstruction is binary normalization after including connecting local-motion potentials and overlapping support geometries; no remaining pattern is claimed to attain 89.

## 4. A saturated incomparability barrier transfers one Gold strategy to arbitrary exteriors

The [saturated-barrier theorem](saturated_barrier_gold_20261019/SATURATED_BARRIER_GOLD.md) finds a connected dominating set of saturated incomparability vertices in a twelve-point seed. This forces every exterior vertex to sit uniformly below or above the barrier, with only specified boundary ports left available.

An exact chain compression then proves that every seven-thin poset containing this seed induced has the SAME strict-Gold comparison tree. Its sharp fixed-tree margin is 50695/453252, greater than 1/9. There is no chosen-critical-root assumption and no exterior-size cap.

A separate deletion expansion permits two ports to have unbounded attachment degree. The same tree has sharp margin infimum 386/6245, greater than 3/50, throughout that larger family. The unbounded result follows from a positive finite generator cone and an affine-tail identity, not extrapolation from large finite examples.

The general Gold and 1/3–2/3 problems remain open. The reusable next mechanism is to find another induced saturated barrier or extend the positive deletion expansion to several unbounded ports while retaining one common comparison strategy.

## 5. Prime generation needs both its divisor distribution and the right phase information

The [split-sum prime packet](split_sum_prime_phase_20261019/PRIME_BRANCHING_PHASE_AND_FERMAT_BALANCE.md) separates Chua's greedy split-sum process from the original Euclid–Mullin rule. With the published Booker–Pomerance support theorem, the smallest missing odd prime q is eligible exactly when the quadratic character of -N modulo q is positive. A nonresidue next prime reverses the blocked phase. Persistent omission is therefore exactly an eventual selected-prime trap in an index-two subgroup. At q congruent to 3 modulo 4 this is a literal Paley-tournament neighborhood condition.

The support theorem and phase-flip mechanism are explicitly attributed to the literature. The resulting greedy recurrence target is still open; freely choosing a useful branch does not prove the greedy process takes it.

Exact binary subset-product balance has a separate complete cyclic classification: the group order must be a power of two, and every logarithmic valuation layer must occur. At a Fermat prime q>3, a minimum balanced source cannot generate that prime; one additional nonresidue factor is needed. At 65537 the missing highest layer gives an exact finite optimization: 327941 is the smallest possible largest prime factor of a balanced squarefree source. An explicit 16-factor source attains the sharp generating threshold.

The [recovery companion](split_sum_prime_phase_20261019/KERNEL_PHASE_RECOVERY.md) refines the apparent information loss. The reflection centers of a binary subset-sum kernel form exactly the total-sum coset of its translation stabilizer. At a Fermat modulus q>3 a nonuniform kernel already determines the quadratic phase; uniform multiplicity one also forces a specific blocked phase. Uniform even multiplicity is the abstract ambiguity regime. At odd-order signed quotients the center is uniquely recoverable, but the quadratic character has already been lost in quotienting by ±1. Even a recovered phase need not identify the candidate root class without the appropriate support premise.

There is also a direct [perfect-core twin-edge bridge](perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md) back to actual aliquot dynamics. Perfect cores are exactly the guarded affine branches that translate the prime cofactor by +2. At core 6, twin-prime pairs p,p+2 with p>=5 correspond bijectively to actual aliquot edges 6p->6(p+2) with both marked cofactors prime. Yet every individual fixed-6 prime lane exits after at most two steps, by reduction modulo 3. Infinitely many legal edges across different sources and an infinite single trajectory are different quantifiers; the bijection supplies no new prime-count estimate.

## 6. Inverse trace moments measure the genuine planar boundary

The [inverse-sheet theorem](inverse_trace_pole_boundary_20261019/TRACE_POLE_BOUNDARY.md) uses the finite normalization of a polynomial Keller map. At every affine target divisor at least one finite inverse sheet survives. This saves one moment from the generic Newton bound: in a primitive monic transverse coordinate, the first d-1 inverse power traces detect exactly the nonproper-value divisors of a degree-d map.

If the first d-2 moments are regular at a boundary divisor, all d-1 escaping sheets have the same pole order and leading coefficients forming a root-of-unity packet. In even degree, polynomial regularity of those moments forces the inverse discriminant to be a square and the monodromy into the alternating group. For degree 4 this excludes the S4 lane under that explicit trace hypothesis; it leaves A4.

No theorem here asserts those traces are polynomial for an arbitrary Keller map. Rational cotangent and polynomial non-Keller examples show why finite-sheet survival, local volume preservation, and zero residues do not alone remove the remaining packet. The next obligation is global polynomial source geometry, retaining divisor classes and weighted boundary incidence rather than only an unweighted picture.

The [tetrahedral weighted-boundary continuation](tetrahedral_weighted_boundary_20261019/A4_WEIGHTED_BOUNDARY.md) supplies an explicit global polynomial non-Keller map with exactly the surviving A4 pole pattern. Its Jacobian is 64e, so it is étale off the source line e=0. The retained finite sheet, cubic pole packet, vanishing first two traces, regular commuting dual derivations, and unweighted boundary forest all coexist.

Its finite normalization is the cone b²=ac with divisor class group Z/2. The inherited free-boundary-class criterion therefore rules out EVERY affine-plane Keller chart for this exact field extension, independently of coordinates. A single exceptional (-2)-curve explains the weight-sensitive obstruction that an unweighted forest misses. More generally, finite domination of a normalization by a factorial affine variety excludes a nontrivial Keller realization.

The packet also classifies the A4 quartic cover with globally standard cusp branch support. It explicitly credits Chau's stronger 1999 theorem when the entire nonproper-value set is homeomorphic to C. This is an explanatory normalization calculation and a scoped source-realization exclusion, not a resolution of the general A4 lane.

## 7. How the connections are being used

Three operations recur with concrete, different proofs:

1. A local expression becomes a global result only after its actual domain is retained: odd-prime cofactors, lawful parity sources, realizable poset lift counts, or genuine affine inverse sheets.
2. An intermediate representation can retain more information than first expected. A divisor kernel can recover its total phase via reflection symmetry; an inverse trace list can detect every boundary pole because one sheet survives.
3. Local constraints can miss a global compatibility condition. The graph examples retain a nonzero period despite local stress agreement, and planar boundary forests require divisor/intersection data in addition to their incidence graph.

These observations are guides to the next proofs, not a universal theorem identifying the underlying problems. The packets give their exact maps, hypotheses, proofs, controls, and still-open targets so later work can extend the strongest valid statement rather than restart from numerical coincidences.
