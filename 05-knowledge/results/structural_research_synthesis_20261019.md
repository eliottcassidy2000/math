# Source flows Gold strategies and exact projection

Snapshot as of 19 October 2026, 11:04 UTC. This synthesis explains the mechanisms behind the latest Hamiltonian-flow and seven-thin-poset results. The strongest conclusions depend on keeping track of the time differential, the precise generating family and topology, and the ambient information retained by a projection.

## Algebraic motion requires control of time

Work over a characteristic-zero field K in the original coordinates x,t. For H=c+t²R(x,t), the derivation δ_H={−,H}, with {A,B}=A_xB_t−A_tB_x, raises t-adic order. Hence exp(λδ_H) exists in K[x][[t]] for every scalar λ. All classifications below concern these actual source-coordinate series at nonzero scalar times.

The [two-exponent classification](planar_two_exponent_flow_classification_20261019.md) is complete for

    H=x^k t^m(1+x²t),   m≥2, k≥0, a=k−2m.

Both coordinates are algebraic exactly when a=0, when a>0 divides m, or when a=−b with b≥2 a proper divisor of m+1. Outside this locus each coordinate is transcendental. Within it the pair is rational exactly when a=0 or (m,k)=(3,4); these rational maps are nonpolynomial. The remaining algebraic cases are nonrational as pairs. Individual coordinate degrees are not classified throughout this larger family.

The [m=2 slice](planar_monomial_carrier_flow_trichotomy_20261019.md) is particularly sharp: k=4 gives a rational nonpolynomial map; k=5 and k=6 give individually nonrational algebraic coordinates, with joint field degrees four and five; every other k gives individually transcendental coordinates.

The mechanism is stronger than rational parametrization of an invariant curve. Put u=x²t and E=H. For a≠0, over an algebraic closure C of the initial-value field, the moving invariant curve and its time differential become

    w^a=v^m(1+v),   dτ=−w dv/(aηE),   η^a=E.

Except at the resonant weight a=−1, weighted dilation promotes algebraicity at one nonzero scalar time to algebraicity in formal time. Trace down a finite algebraic cover then forces w dv to have a rational primitive on the invariant curve. Cyclic-character projection reduces exactness to

    a v(1+v)R′+[m+(m+1)v]R=a v(1+v),   R∈C(v).

Pole orders yield the divisor conditions. Explicit polynomial time primitives prove sufficiency. At a=−1, a separate logarithmic-norm obstruction handles the resonance. For example, k=3,m=2 has a rational invariant parameter u, yet its time differential du/[u²(1+u)] has nonzero logarithmic residues and every nonzero-time coordinate is transcendental. Rational invariant geometry alone does not determine algebraic motion.

The [quintic analysis](planar_quintic_s5_petersen_resolvent_20261019.md) identifies S5 as the k=6 geometric and arithmetic Galois group. Each coordinate generates the quintic field and cannot be expressed by radicals. The ten root-pair sums carry the Petersen graph's faithful S5 action; adjoining the square root of the discriminant leaves A5, abstractly the icosahedral rotation group. These are precise permutation-group connections, without a claim of spatial symmetry of the flow. Its exact E→∞ specialization is V⁵−z⁵, with pair-sum resolvent Y¹⁰+11z⁵Y⁵−z¹⁰. The normalized pair sums therefore have fifth powers −φ⁵ or φ⁻⁵, where φ=(1+√5)/2 and φ⁵−φ⁻⁵=11. The associated C5 rotation describes this root configuration, not the specialized Galois group.

## Polynomiality changes with the generating family and topology

For a single arbitrary H=c+t²R, [polynomiality has a simple criterion](planar_single_hamiltonian_polynomiality_20261019.md):

    exp(λδ_H) is polynomial ⇔ H_x=0,

and then the map is (x+λH′(t),t). A polynomial realization fixes t=0 pointwise and has Jacobian one. [Gwoździewicz's fixed-line lemma](https://arxiv.org/pdf/alg-geom/9305008) makes it an automorphism. Its irreducible second coordinate is divisible by t, so the source jet forces T=t. If H_x≠0, the first nonzero term of T−t is −λH_x, a contradiction. Finite commuting words in functions f_i(H) reduce to one effective Hamiltonian; when H_x≠0 they are polynomial exactly when the effective polynomial is constant, hence exactly when the word is the identity.

Different-invariant words require their own argument. For the [affine-x family](planar_affine_x_rational_word_rigidity_20261019.md)

    H_i=a_ixt²+f_i(t),   a_i≠0, f_i∈t²K[t],

every finite word is rational and is polynomial only when it is the identity. In canonical coordinates z=1/t,w=xt², a word becomes (z+L,w+G(z)). Every translated principal part contributes opposite coefficients at its two endpoint poles. The sum of coefficients at each pole order is therefore zero. Polynomiality first forces L=0, then permits a pole only at zero; conservation eliminates that last pole.

The [ordinary source t-adic closure](planar_source_closure_formal_shears_20261019.md) is much larger. Even the flows of just xt² and (x+1)t² generate a dense subgroup of the exact closure

    (x,t) ↦ (x(1+Lt)²+g(t), t/(1+Lt)),
    L∈K, g∈tK[[t]].

Its polynomial members are all shears (x+p(t),t), p∈tK[t]. Return loops give corrections (1+bt)^−2−1, and finite differences prescribe successive jets. The finite-word principal-part sums are discontinuous in this topology, so their obstruction need not survive a limit. This is a source closure, with no identification with a joint-boundary closure.

The [broader finite-word theorem](planar_full_source_group_polynomial_intersection_20261019.md) enlarges the generators to every x-dependent H=c+t²R. It realizes every source polynomial shear using three individually nonpolynomial rational factors. For g∈tK[t] nonzero, choose f′=g, f∈t²K[t], and set

    J=x³t²,   H_g=f(t(1−x²t)³).

Apply the J flow at −1, the H_g flow at +1, and the J flow at +1, in that chronological order. Their exact composite is (x+g(t),t). Together with the fixed-line upper bound, this identifies the broader group's polynomial intersection as precisely the shear group. Finite affine-x words, their source closure, and this larger generating class are three distinct questions.

## A sufficient Gold structure can fail while Gold succeeds

A seven-thin poset has at most seven incomparable neighbors per vertex; a balanced comparison has probability in [1/3,2/3]. Under Zaguia's definition, a good pair is an ordered incomparable pair in the poset or its dual with D(a)⊆D(b), chain U(b)\U(a), and p(a<b)≤1/2. Such pairs supply a Gold strategy. The universal height-two existence theorem does not extend to height three.

The [ten-point height-three witness](height_three_seven_thin_no_good_pair_20261019.md) is rigid and has 11,971 linear extensions. Every one of its 48 structurally eligible orientations has probability at least 6,547/11,971>1/2, so it has no good pair. Nevertheless, comparing vertices 3 and 5, then 2 and 4 regardless of the answer, gives extension cells

    (3281,3281,1899,3510).

If E is the initial extension count and t₁,t₂ are the remaining counts after the two answers, the four margins E−t₁−t₂ are (2128,2128,4663,3052), all positive. Thus failure of the sufficient structure coexists with a directly certified strict Gold strategy. Ordinal sums give arbitrarily large rigid seven-thin examples with no good pair and the same normalized positive margin; their heights grow. The unrestricted seven-thin Gold question remains open.

## Projection is exact when maximum span prevents hidden vertices

In a hypothetical seven-thin counterexample to balance, the strict >2/3 majority directions form a total order M. Choose an incomparable root pair x,y of maximum index span, and retain its full core

    Q={x,y}∪(D(y)\D(x))∪(U(x)\U(y))∪(I(x)∩I(y)).

Every exterior vertex is below both roots or above both. If another core pair ties the root's span after projection to Q, [the exact-transfer lemma](canonical_maxspan_transfer_20261019/canonical_maxspan_transfer_20261019.md) proves that no exterior vertex was deleted between its endpoints. Every common incomparable neighbor of a maximum-span pair lies between its endpoints. Therefore its core common-neighbor set is exactly its ambient common-neighbor set.

Choosing a maximum-span root with maximum common-neighbor count now permits a canonical-root screen. This restriction must remain explicit. Transfer to earlier any-root exclusions additionally bounds the secondary pair's exclusive neighbors using saturation and exterior degree slack. Common-neighbor exactness does not automatically make both exclusive counts exact, and projection does not preserve comparison probabilities.

The latter require the boundary-generator lift: actual exterior extension weights are nonnegative combinations of core generators. A weighted inequality valid on every retained generator transfers to every actual exterior. A feasible relaxed mixture need not describe an exterior or a counterexample. This mechanism underlies the completed [(2,4,0) and (3,3,0) exclusions](mixed24_mixed33_c0_exact_exclusion_20261019/mixed24_mixed33_c0_exact_exclusion_20261019.md), covering 56,265 residual orders by exact certificates without an ambient-size cutoff. Short spans explicitly use Peczarski's established six-thin theorem.

## Remaining problems

- Classify rational and algebraic source flows for general polynomial R, and individual coordinate degrees across the two-exponent algebraic locus. The fixed-line criterion settles polynomiality of one flow, not these finer questions or entry into a required source chart.
- Close the remaining seven-thin morphologies with complete coverage. The completed any-root exclusions leave 26 parameter triples in that catalogue, including (2,2,4), (2,3,3), and (2,4,1). Canonical-root exclusions have narrower scope and do not reduce this any-root count. At this snapshot the canonical (2,4,1) branch is awaiting its final fresh independent replay and acceptance; it is not counted as a completed exclusion here. Canonical reductions and partial scans alone do not close a morphology.
- For no-good-pair height-three posets, study the actual decision-tree quantity G(P)=max_tree min_path(E−t₁−t₂). The second comparison may depend on the first answer, and determined or repeated comparisons are allowed. With this convention, G(P)≥0 is weak Gold. The witness has G(P)≥2128; a universal argument needs more than good-pair existence.

These results isolate exact local mechanisms. They establish neither the general seven-thin balance or Gold conjectures nor the planar Jacobian conjecture.
