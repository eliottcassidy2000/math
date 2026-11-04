# Prime-five closure and the original-module collision frontier

## 1. Frozen inputs and physical convention

The full named response F_T(a,c) counts Hamilton paths after adding actual ports x->y, where a_v=[v->x] and c_v=[v->y]. Physical names are never permuted between source and target. This integration uses these separately proved and reviewed results:

1. `../prime_five_all_ambient/PROOF.md`: every prime five-name quotient converse with arbitrary nonempty named children is separated in every fixed ambient. Source child orders and profiles may be unequal. Target child interiors may be restored to any identical corresponding complete ordered-cover profiles. The finite mod4, mod8 and special mod16 branches cover both total parities.
2. `../cyclic_ambient_synthesis/PROOF.md`: the corresponding arbitrary-child, arbitrary-fixed-ambient result for the directed triangle, together with its highest-node consequence.
3. `../module_decomposition/PROOF.md` and `../highest_quotient_reduction/PROOF.md`: simultaneous edge-XOR of modules of the ORIGINAL tournament preserves a common recursive partition on the original strong-module tree. Scalar Hamilton equality implies hereditary equality of every ordered-cover profile. Highest changed nodes are disjoint; all changed edges lie inside them. A highest changed linear quotient is detectable in every ambient. A changed prime quotient is its converse.
4. `../monotone_dual_queries/PROOF.md`: baseline quotient-converse cover-profile contrasts for z and Jz=(1-c,1-a) are exact negatives, with the same fixed exterior port incidences.

A strong module is a module overlapping no other module, not a synonym for a strongly connected tournament. All tree statements concern the original strong-module tree, not a claim that it is the target's canonical tree. A selected module contributes an edge toggle exactly when it contains both ends; repetitions are reduced modulo two.

## 2. Precise source-relative prime-five family consequence

Let the source core be T=Q[A_0,...,A_4], where Q is prime and every A_i is a nonempty tournament on a fixed physical name set. Fix any rooted ambient R. Consider every family of original modules contained in this core, with the whole core selected oddly after parity reduction. Proper selected modules may overlap or cross. Do not select modules containing exterior names.

There is a child-uniform two-port query, depending on the source and this fixed ambient but not on the family target, that separates the source from EVERY such target. Both ports imitate the core at each exterior name. A fixed two-query set {z,Jz} suffices before choosing the sign.

Indeed, every proper original module is contained in a child of the prime root. The root quotient is therefore conversed, and hereditary original-module monotonicity gives q_h(B_i)>=q_h(A_i) for each target child B_i. Internal converse from the selected ancestor preserves q. Choose any baseline separator z supplied by input 1 and pair it with Jz, complementing and exchanging the CORE bits only. Their actual baseline ambient contrasts are negatives by input 4. At least one is strictly negative in the source-minus-baseline direction. Fix that query. Uniform port incidences keep every child a module of the actual augmented graph. Increasing its complete q profile cannot decrease the target Hamilton count, because occupied-run substitution has nonnegative coefficients. Hence the already negative source-minus-baseline contrast stays negative for all family targets simultaneously.

The sign-selected query can depend on the fixed ambient. A single member working with the same sign in all ambients is not asserted. This is detection relative to the specified source, not pairwise injectivity between arbitrary family outputs. It does not cover root-even families merely because the source has a prime-five root. The argument uses profile monotonicity and duality, not an exposed-anchor theorem requiring laminarity.

## 3. Highest changed prime-five nodes cannot occur in a full collision

Let U be obtained from a named T by simultaneous XOR of any family of ORIGINAL modules. Suppose T differs from U but F_T=F_U. Since F_T(0,0)=H(T), scalar Hamilton counts agree. Input 3 gives hereditary equality q_h(T[S])=q_h(U[S]) for every fixed old-name subset S and every h. It also gives a nonempty disjoint set of highest changed original strong nodes. All such quotients must be prime and conversed, because a highest changed linear quotient is already excluded.

Suppose one highest changed node M has a prime quotient of order five on original children A_0,...,A_4. Let V replace only T[M] by U[M], leaving all exterior names and arcs as in T. Each child has its original q profile by hereditary equality, and M has unchanged exterior incidences. Input 1 therefore provides a physical child-uniform query on T and V that distinguishes them in this fixed ambient. Both ports imitate M at exterior names.

For every other highest changed node N, the source and target complete q profiles agree. The two ports see N uniformly: M and N are disjoint old modules, all their cross arcs are uniform, and the ports have M's incidences. Thus N remains a module in the actual port extension. Replacing its source interior by its target interior preserves the queried Hamilton count exactly, by the occupied-run substitution identity. Replace all other highest nodes successively. Their disjointness and fixed boundaries preserve these hypotheses, and every changed edge lies inside one of them. The final result is U, with the same physical query value as V. This contradicts F_T=F_U.

This proof does not require the selected family to be laminar. It uses the common original-tree decomposition and scalar-equality restoration; it does not assume intermediate selected modules remain modules.

## 4. Certified residual and scope

Input 2 excludes highest changed prime quotients of order three by the same restoration argument. There is no prime quotient of order four, and a two-child quotient is linear. Combining those exclusions with Section 3 proves:

**Corollary.** In any still-possible nontrivial full named-response collision obtained by simultaneous XOR of original modules, EVERY highest changed quotient is prime of order at least six.

In particular, if all prime quotients in the original strong-module tree have order three or five, every original-module XOR target different from the source is detected by the full named response.

This is a necessary frontier, not an existence assertion for a collision at order six, an exclusion of larger prime quotients, or unrestricted named reconstruction. It is not extended to arbitrary sequences of inversions of later/current modules that were not original modules. Any later order-six exclusion requires its own separately reviewed input and integration; it is not part of this checkpoint.

## 5. Reproducibility

The accompanying standard-library checker binds every frozen input, checks the small no-prime-four assertion directly, and constructs three physical examples with two disjoint highest changed nodes: a prime-five component and a triangle in an ordered join. It finds a literal query at the five-node component, verifies exact restoration of the other conversed component, and retains all source/intermediate/target rows, physical masks and actual counts. These focused controls illustrate the restoration argument; the all-order conclusion is the proof above, not an extrapolation from them.
