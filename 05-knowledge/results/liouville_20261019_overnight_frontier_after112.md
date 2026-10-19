# Overnight proof frontier: exact results and falsifiable next steps

Snapshot: 19 October 2026, 09:49 UTC. This is a new successor to the [earlier overnight frontier](liouville_20261019_overnight_frontier.md); it does not overwrite that record. The [structural synthesis](preserved_structure_and_proof_transfer_20261019.md) explains the proof mechanisms. This document reports completed proofs and explicitly distinguishes unfinished work. It makes no claim that publication has completed or that historical novelty is established.

## 1. Square-sum cycles: the prescribed seed is settled exactly

The [exact-threshold theorem](seeded_optimal_threshold_6242_20261019/OPTIMAL_THRESHOLD_THEOREM.md) proves:

> A square-sum Hamiltonian cycle on [1,N] contains 482,3362,6242 consecutively if and only if N≥6242.

The finite certificate covers every size 6242 through 128452 by protected insertions and reversals; the [eventual construction](seeded_eventual_cycle_theorem_20261019/SEEDED_INTERVAL_THEOREM.md) handles every N≥128453 by finite universal transition tables. The two components meet without a gap. Ordinary interval Hamiltonicity is inherited background; the extra conclusion is the fixed consecutive seed.

The [all-size connected-path theorem](square_sum_all_size_connected_paths_20261019/ALL_SIZE_CONNECTED_PATHS.md) additionally gives, simultaneously for every N≥6,294,221,

    (2401+t, 7203−t, 14406+t),
    t in [−1200,1200] excluding [−24,24] and [74,122].

There are 2,303 pairwise vertex-disjoint triples, with adjacent sums 98² and 147². They coexist with the seed; t=961 also forces its continuation to 15367. Disjointness is among the triples, not between all triples and the seed. The sufficient threshold is not claimed optimal.

**Next certificate target.** Supply a compatible all-size B-pair family with the longer prefix (1,3,6,19), then use the [finite-cut theorem](square_sum_finite_cut_retention_20261019/FINITE_CUT_RETENTION.md) to compute the surviving copies. Before any search, fix the permitted input invariants, states, and cuts. A proved parity/coverage obstruction in a required state falsifies that proposed transition system; an unsuccessful search does not. Prefix (1,3,6,10) is an unsuitable shortcut for the existing B-pair endpoint specification, because its odd-length member already terminates at 10.

## 2. Graceful reconfiguration: upper tails succeed; arbitrary frozen colors fail

For every acyclic choice of the fixed colors r+1,…,m on labels 0,…,m, the [entire completion fiber](graceful_upper_tail_two_color_connectivity_20261019/UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md) is nonempty and connected by moves changing at most two free colors. The proved diameter bound is 2^(r−1) for r≥1, and zero for r=0. The proof is a rank-staircase induction with closure-preserving double exchanges, not enumeration of small trees.

The [13-vertex counterexample](graceful_fixed_color_counterexample_20261019/ARITHMETIC_FIXED_COLOR_COUNTEREXAMPLE.md) disproves the analogous statement for an arbitrary frozen color set. Free colors 6,9,11 give precisely five completions in two components. Consequently a global connecting path can have to disturb an already-correct edge. Abstract tree shape is not preserved, so the Graceful Tree Conjecture is untouched.

**Sharpened hypothesis to test.** Every upper-tail fiber has width-two diameter at most r². This is unproved; it is a specific replacement for the exponential upper bound, not an asserted consequence of it. A fully enumerated fiber with exact diameter greater than r² would falsify it. The fixed-color counterexample does not, because it violates the upper-tail premise. A proof must permit temporary changes to correct edges.

## 3. Gold: good pairs give a strategy; a new structural class supplies them

The [attributed good-pair synthesis](good_pair_gold_synthesis_20261019.md) proves weak Gold when D(a)⊆D(b), U(b)∖U(a) is a chain, and the current probability of a<b is at most one half. The proof explicitly combines Zaguia's decreasing insertion bins with the half-crossing Gold calculation associated with Peczarski. The full-upset-chain and upset-difference-chain hypotheses are distinguished in the source note.

A very good pair permits the required orientation. The [full height-two theorem](seven_thin_height_two_very_good_gold_20261019.md) proves that every finite non-total seven-thin height-at-most-two poset contains such a pair, hence has weak Gold. Isolated vertices are allowed. This strengthens the [earlier c=0 theorem](seven_thin_height_two_c0_very_good_pair_20261019.md): no critical pair or common-neighbor restriction is required. An elementary separated-set-family argument proves the structural result; ordinary height-two balance was already known. It does not assume that higher-height orders can be replaced by an incidence matrix.

**Next structural hypothesis.** Every non-total seven-thin poset of height at most three contains a very good pair. This is unproved. A finite poset with no four-element chain, all incomparability degrees at most seven, and no very good pair in the poset or its dual is a complete falsifier. Failure would not disprove Gold; it would identify the need for a different strategy. Unlike the settled height-two case, both predecessor and successor data can vary at the middle level.

## 4. Seven-thin balance: the remaining maximum-span morphologies are explicit

Assume hypothetically that a seven-thin poset has no balanced pair. Its strict >2/3 majority tournament is a total order. Choose an incomparable root pair x before y with maximum index span, and put

    l=|D(y)∖D(x)|, u=|U(x)∖U(y)|,
    c=|I(x)∩I(y)|.

The roots are critical, and their degrees give l+c≤6 and u+c≤6. The [completed one-sided and min-exclusive-one reduction](mixed16_c0_exclusion_20261019/two_sided_two_exclusive_reduction_20261019.md) proves l,u≥2. By duality take l≤u.

Further completed exclusions are:

- mixed-(2,2): [c=0,1,2](mixed22_c02_exact_exclusion_20261019/mixed22_c02_exact_exclusion_20261019.md) and [c=3](mixed22_c3_exact_exclusion_20261019/mixed22_c3_exact_exclusion_20261019.md);
- mixed-(2,3): [c=0,1](mixed23_c01_exact_exclusion_20261019/mixed23_c01_exact_exclusion_20261019.md) and [c=2](mixed23_c2_exact_exclusion_20261019/mixed23_c2_exact_exclusion_20261019.md).

The proofs use exact exterior-defect screens and independently checked positive integer boundary-generator certificates, with no ambient exterior-size cutoff. Mixed-(2,2), c=3 checks 82,987 degree-valid cores and certifies 248,199 residual majority orders. Mixed-(2,3), c=2 checks 87,137 degree-valid cores and certifies 476,837 residual orders. These are completed universal morphology exclusions via exact finite reductions. The c≤2 mixed-(2,2) and c≤2 mixed-(2,3) results explicitly use Peczarski's six-thin theorem for short spans. Mixed-(2,2), c=3 does not: its common neighbors force span at least four.

At this snapshot the exact unclosed list, up to duality, is:

| Exclusive sides (l,u) | Common-neighbor counts still unclosed |
|---|---|
| (2,2) | 4 |
| (2,3) | 3 |
| (2,4), (3,4), (4,4) | 0,1,2 |
| (3,3) | 0,1,2,3 |
| (2,5), (3,5), (4,5), (5,5) | 0,1 |
| (2,6), (3,6), (4,6), (5,6), (6,6) | 0 |

These are 28 parameter triples, not 28 posets or 28 independent proofs. The full height-two theorem removes every height-at-most-two instance throughout this frontier. Any hypothetical remaining counterexample must therefore have height at least three; higher-height c=0 cases remain.

**Immediate boundary.** Mixed-(2,2), c=4 and mixed-(2,3), c=3 are still open in this work. A preliminary analytic census or a partial residual scan is not an exclusion theorem.

**Exact next certificate hypothesis.** For each remaining core and admissible projected majority order, after the documented saturation, directed-defect, and short-span screens, a positive weighted combination of forward comparisons is at most 2/3 on every retained boundary generator. A rational generator mixture making all those forward comparisons greater than 2/3 falsifies this stronger cone hypothesis. It does not falsify seven-thin balance: realizability and exterior comparison witnesses must then be investigated. The [actual mixed-(2,2) realization](mixed22_c02_exact_exclusion_20261019/mixed22_core_realization_and_span_dependency_20261019.md) already shows why even realizable core weights can require exterior witnesses.

## 5. Collatz: uniform reset heads have an exact compression boundary

The [finite-bank theorem](collatz_bank_reset_20261019.md) applies to any expansion-bank word followed by a reset-bank word, subject to full parity guards and a contracting total slope. Bounded expanding centers C and reset fixed points P give at most ceil(C+P) excluded initial source lifts. Adjacent companions share the same exact cutoff.

For the concrete banks {K,K⁴F} and {F,H}, the [modular refinement](collatz_modular_bank_head_20261019.md) gives uniform head 107, or 42,77,107,24 according to the parities of total length and odd count. Appending one further F or H after the first contracting boundary removes the head on the enlarged lawful cylinder. It is not licensed on every source of the shorter cylinder.

The [unbounded-center theorem](collatz_unbounded_return_centers_20261019.md) constructs actual expanding adjacent macroreturns with no descending lower prefix and unbounded compressed centers. Their slopes approach one from above. It refutes an unrestricted uniform-center assumption for all such macroreturns, while leaving the fixed-bank theorem intact.

**Exact next hypothesis.** Every contracting expansion-then-reset word in these particular banks has zero excluded head already, without the extra reset. In the notation T^N(x)=(Sx+H)/M, Δ=M−S>0 and lawful least residue ρ, the assertion is Δρ>H. One guarded word with Δρ≤H falsifies it. The 396 checked schedules all pass, but do not establish the universal assertion. Even a proof would remain a cylinder theorem, not arbitrary source entry or Collatz convergence.

## 6. Planar Hamiltonian repair: finite polynomial composition is closed

The [rigidity theorem](planar_carrier_factor_rigidity_20261019.md) excludes nonzero-time polynomial specialization of every nonconstant universal carrier c+p²(p³−y²)R. An invariant polynomial's full divisor support is permuted by a hypothetical Keller realization; the compulsory components t and 1+x²t force a finite-order iterate. The actual completed Hamiltonian jet forces infinite order, a contradiction.

The stronger [finite-composition theorem](planar_carrier_composition_rigidity_20261019.md) is also proved and independently reviewed: every polynomial realization of a finite composition of these universal flows is the identity, even with different invariants and arbitrary scalar times. Cancellation is allowed. Gwoździewicz's established injectivity-on-one-line theorem supplies automorphy; a faithful common-ring second-boundary evaluation then kills the only possible shear. A separate completed-cusp intersection argument gives the same final obstruction. No common Hamiltonian is assumed.

**Next precise hypothesis.** There is no nonidentity rational source realization of a finite composition of these same universal flows. This remains unproved. A rational pair exactly equal to the completed images of x,t, regular as the required source series, would falsify it if nonidentity; a finite jet match would not. Polynomial fixed-line and divisor arguments cannot simply be applied across poles. Other Hamiltonian classes and entry into the source-normal chart remain separate questions. None of these local results settles planar JC.
