# Session22: rational metagraph strata and route-free tail certificates

Research date: 2026-10-21. **PROVED:** the scoped structural and algebraic results below, using named inherited inputs. **FINITE-EXACT:** the documented executable controls and selected native scales. **OPEN:** uniform native 153N coverage, all deficit-six center coefficients, complete metagraph wiring, bridges and connectivity.

Incoming research was independently reconciled through commit `359ddd20b63b00e24bc64204e1e7e78ed5abd5e4`. Remaining-target statements are scoped to that reviewed input cutoff plus the session results below; later repository work may sharpen them.

This session produced new proofs and executable evidence after Session21's workspace loss. It does not claim to restore every lost file or audit. Reviewed packets were published incrementally rather than held only in local storage.

## 1. Every fixed fundamental black degree has a rational counting series

[Complete proof and independent controls](../fixed_degree_metagraph_rationality_20261021/README.md).

Use the fundamental converse-merged paired-black graph: degree is h=H/|Aut| at a non-self-converse pair and (h-b)/2 at a self-converse class, with loops counted twice. This is not the literal doubled-edge network.

Let p(z)=1-z-z^3 and A(z)=1/p(z). The inherited effective black-deficit bound makes every fixed-presentation strong-atom alphabet finite. Unique ordered SCC decomposition then gives, for positive odd h,

F_h(z) = sum over ordered nonrigid factorizations h=h1...hr of A(z)^(r+1) product P_hi(z),

where P_j is the finite polynomial of actual strong atoms with presentation count j. This proves rationality without pretending that the unknown P_j have been computed.

Combining this with the inherited exact converse-center divisor law proves that EVERY fixed-degree node-count series B_k is rational. For odd k its denominator divides

p(z)^(Omega(k)+1) p(z^2)^(L_k+1),

where Omega counts prime factors with multiplicity and L_k is the maximum Omega(x) over odd divisors x of 2k satisfying x(x-1)<=2k. For even positive k only the second factor is needed. Degree zero has the already known pure-blue series.

The consequences are stronger than a finite census:

- Every fixed-degree sequence has an eventual constant-coefficient recurrence, although unknown finite center/atom data still determine its numerator and onset.
- Even fixed degrees have at most polynomial times rho^(-n/2) growth, where rho is the positive root of p. Odd degrees with F_k nonzero have rho^(-n) growth times a polynomial. Existence of a self-converse degree-k node alone does not imply F_k nonzero.
- For k=3^r, the leading constant is explicit:

  [z^n]B_(3^r)(z) ~ rho^(4r-1) n^r rho^(-n) / [2 r! (1+3rho^2)^(r+1)].

- For every fixed K, nodes of degree at most K form a vanishing fraction of all merged nodes. This is not a connectivity result or a statement uniform in growing K.

The new contribution is the general fixed-degree synthesis and its consequences. The center law, rigid alphabet, trivalent special case and h=3 classification are explicitly inherited.

## 2. Deficit-six centers: combined degree at most 35, with seven unresolved coefficients

[Our odd-shore and score-run cutoff](../odd_shore_deficit_centers_20261021/README.md).
[Our two-shore exclusion and final reconciliation](../two_shore_deficit_exclusions_20261021/README.md).
[Incoming combined cutoff](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md).

For a marked Hamilton cycle write n=qm, where q is its exact rotation index and m is odd. The inherited cut-orbit bound gives q<=8 for deficit six.

Our structural contributions are:

1. A proper subset of one rotation shore can be deleted while retaining strong connectivity through the other shores and the marked-cycle matchings.
2. An eight-cycle in a shore of odd size at least nine has an external vertex with at least five arcs in one direction, forced by the exact odd cross-edge total. A half-supported power-of-two prefix gives too many black orbits.
3. Consecutive equal-score runs distinguish additional Hamilton-word orbits from cycle cuts. This excludes high-index score patterns and reduced the original order bound4096 to47.
4. A separate two-shore insertion construction proves Delta>=2m-2 for every self-converse q=2 frame with m>=5. It excludes orders10 and14 analytically. A complete128-frame order-ten computation provides additional independent controls, not the all-order proof.

During close-out, incoming work at f4353d7d composed our published odd-shore result with the older general closing-prefix theorem. That theorem permits choosing a Hamilton cycle with nonautomorphic successor at every odd order at least seven; even orders already have q>=2. Removing the q=1 branch gives the incoming bound35 and ten remaining orders. Our new q=2 exclusion removes14. A final incoming checkpoint at359ddd20 then supplied complete necessary-universe exclusions of12 and25, independently reviewed and replayed here. The combined frontier now has SEVEN possible unknown coefficients:

R6(z)=z^5+z^6+z^7 + sum_{n in U} r_n z^n,

U={15,18,20,21,24,28,35}.

The low coefficients are explicitly inherited. Each r_n is an unknown nonnegative integer; the list does not assert existence. Degree35 and the closing-cycle composition are credited to the incoming work; the intermediate nine-position intersection used our additional two-shore theorem, and the final seven-position intersection also uses the newer incoming exclusions.

A further necessary restriction comes from replacing score colors by the actual vertex orbits of Aut(T). Anti-automorphisms permute these intrinsic colors with matching cardinalities, so the same long-run obstruction forces vertex transitivity at orders25 and35. Our reconciliation and the newer incoming proof both retain this condition. The incoming [order-25 packet](../center_transitivity_20261021c/README.md) goes further: a transitive prime-square permutation group has a regular subgroup, so the two groups of order25 give8192 necessary Cayley frames. Seven literal black words with distinct full ordered-adjacency keys exclude every frame. This is not an isomorphism-class census, and the prime-square group argument is not exported to35.

The incoming [order-12 packet](../deficit_twelve_20261021c/README.md) authenticates all262144 necessary rotation frames with literal sufficient exclusion certificates. Existing consumers for12,14,25 pass in both normal and optimized execution; discovery programs were not rerun. These finite-complete results are credited to the incoming work. Order35 remains OPEN and must be vertex-transitive.

The exact minimal eventual recurrence order of the trivalent series remains18; the combined bound guarantees its homogeneous recurrence for n>=48. Neither actual degree35 nor earliest onset48 is claimed.

Our separately proved Cauchy–Davenport connectivity lemma gives kappa=(p-1)/2 for prime-order circulants, with a composite-order hostile. Its standalone center bound45 and20-position list remain valid intermediate results, but are superseded here: the older closing-cycle theorem already gives the stronger prime deficit lower bound p-1. The prime lemma must not be presented as the sharpest center frontier.

The shore mechanism also gives quadratic order bounds for bounded-deficit self-converse atoms and bounded-presentation non-self-converse atoms. These improve finite classification targets, not the missing classification itself.

## 3. Faster exact tournament enumeration and adjacent channels

[Degree-local arithmetic kernel](../degree_local_tournament_enumeration_20261021/README.md).
[Self-converse, merged and strong channels](../converse_tournament_count_channels_20261021/README.md).

Store a degree-d aggregate with the odd part O(d) of d!, instead of a common O(N). Every transfer coefficient

O(d+ell(k+e)) / [O(d) ell^k O(k)]

is an integer. It can be computed once per degree/multiplicity choice rather than divided separately in every response-state transition. The future-response equivalence relation and mathematical transition graph are unchanged.

Independent telemetry-matched shared-host medians show ordinary speedups of 1.83–2.16 at bounds 80–200 and 1.33 for the core series at bound100. These are finite implementation comparisons, not universal-fastest, bit-complexity or asymptotic runtime claims. The polynomial transfer-entry bounds describe the precomputation family, not full runtime.

The general weighted channel F_n(a,b) counts ordered a-tuples of tournaments on one vertex set with 2^b vertex colors. Its degree-local implementation supplies ordinary, self-converse, merged, strong, strong-self-converse and strong-merged counts. The half-order formulas C_(2m)=F_m(2,0), C_(2m+1)=F_m(2,1) and SCC formulas are inherited identities with fresh executable verification.

Incoming [cycle-count modules](../tournament_tail_kernels_20261021b/README.md) now provide a different structural-state reduction, with the actual cycle count retained as a polynomial payload. Independent reconciliation verified the normalization and terminal linear functional. In its paired samples, signed-core bound100 took0.2486s versus0.8517s for our scalar routine (about3.43x faster); ordinary bound200 scalar2.1967s versus hybrid2.1777s is treated as tied. These are incoming results, not a new benchmark run or a universal ranking. The current measured preference is scalar arithmetic for ordinary counts and the polynomial method for the tested signed-core workload.

A source explanation in THM-283 required correction: anti-automorphism consistency depends on endpoint swap as well as unordered-pair orbit length. A transposition's length-one pair orbit is allowed; a four-cycle's antipodal length-two orbit is not. The valid partition formula and numeric sequence are unchanged.

An explicit equivariant map also identifies a tuple (T1,...,Ta) with T1 decorated by a-1 simple graphs recording which pairs disagree with T1. This preserves actual relabelling structure. The resulting equality with an even-order self-converse count is not, by itself, a constructed canonical bijection to unmarked self-converse tournaments.

## 4. A finite prefix can certify every future arrival horizon

[Route-free gate, mixed strengthening and certificates](../routefree_suffix_gate_20261021/README.md).

At a finite cut let W,H be arrived old/new terminal vectors, Delta0 the current physical energy increment, and let the unfinished sources retain their indexed original charge vectors. Let D_o,D_n be their exact signed-zonotope energy maxima. Let L_n and L_o be the corresponding linear supports against W+H and H. Then for EVERY later finite cut, without assuming any unfinished source ever arrives,

Delta <= Delta0 + D_n + L_n + L_o + 2 sqrt(D_o D_n).

A separate exact 144-corner maximum replaces D_n+2sqrt(D_o D_n) by the no-larger joint quantity max[Q(S_n)+beta(S_o,S_n)]. Group-sign averaging proves the bound even for arbitrary missing sources, destinations and future clock signs. These signs are a proof device, not changes to the actual dynamics.

Fresh bounded-prefix certificates cover every horizon at the following individual scales:

| N | Finite cut | Unfinished nonzero new sources | All-future upper | Budget153N |
|---:|---:|---:|---:|---:|
|4096|105|386|617824|626688|
|8192|112|615|1211671|1253376|
|16385|121|961|2470224|2506905|
|32768|131|1320|4763921|5013504|

All earlier cuts are directly checked. Independent native replay covers4096/8192 and the N7 hostile; the two larger positive scales are producer-verified. These four scales do not represent an entire interval. The new proof can stop while many routes remain unfinished.

The method allocates NO negative pair credits. For an existing owned repartition its physical conclusion is exactly R'<=M-J', not an automatic R'<=M. The refunded term K is still added back to R'; no extracted pair is spent again. The actual source identities, mod15/mod35 compatibility and full past clocks remain authenticated.

### Incoming monotone block envelope and safe composition

The incoming [stopped-current-node envelope](../future_tail_blocks_20261021b/README.md) retains separate old/new vectors for sources that occupy the same unresolved node at the same native depth. They have a common future continuation. Its absolute local/pair envelope U_D is valid for all future cuts and nonincreasing under actual transport, merging and arrival. It yields stronger supplied certificates at some of the same scales, including4096 at cut105 with upper557164.

Independent reconciliation confirms that neither U nor our raw-source zonotope bound B dominates universally. At one fixed scale, taking the minimum of U_D and all valid earlier B certificates gives a monotone all-future upper bound. Prefix coverage is still separate; this is combining inequalities, not spending credit again.

A shortcut replacing optional mixed old/new blocks by only full signed placements is unsafe because their energy difference is indefinite. The incoming four-point zero-flow optimizer likewise improves compression, not physical energy: the refunded cross term changes together with removed zero norm. All these bounds concern the stated native-depth Nyquist energy, not every Fourier frequency or integer-clock cut.

### Final incoming geometric and excursion refinements

The newer [linked tail geometry](../tail_geometry_20261021c/README.md) keeps old/new blocks coupled, expresses their indefinite energy as a difference of positive quadratic norms, and applies an exact planar zonotope majorant. Its separately checked bounds at the SAME frozen positive prefixes are551844,1118881,4142214 for4096,8192,32768. This is improved certificate evaluation, not a new native scale or source replay. Pair-cluster bounds are separate complete upper bounds; their gains must not be added to angular gains as if they were independent credit.

The [moving-prefix ledger](../uniform_prefix_20261021c/README.md) supplies exact debit and associative excursion summaries. Its stronger moving potential rises by4215 between cuts17 and18 atN4096, while staying under budget. This refutes pointwise monotonicity, not coverage. A block must retain its maximum PREFIX debit, not only its net debit, to pay every intermediate horizon. The existence of uniformly payable blocks or uniformly successful stopping certificates remains OPEN.

## 5. A precise scale obstruction, not a false uniform induction

[Conditional propagation and native counterexamples](../tail_scale_certificate_obstruction_20261021/README.md).

When old sources have finished, six antipodal pending counts, six terminal supports and the current increment determine the sufficient certificate value. They do NOT determine future evolution. Two explicit sources with the same guards, charge and first-edge clock have different next-arrival status because their magnitudes and frontier membership differ.

There is a useful exact interface: at the same horizon, the next scale's old-source domain is the previous old-plus-new domain. Thus W_(64N)=W_N+H_N and pending-old axis counts add. The unknown outer population is not supplied by this identity.

The native affine lifts add six division-clock units, not six primitive arrival steps, and cover only 1/64 of primitive source addresses. A literal subset of4096 still-pending sources at N262144 and cut111 has signed energy7499929, greater than64 times the base cut105 value91876. Thus the naive depth+6 contraction is false. This does not show failure of153N at the enlarged scale.

A further logical obstruction matters: demanding that all old sources finish at some finite cut on every64-scale tower is equivalent to universal primitive entry into the fixed frontier. It cannot be silently assumed in a uniform induction. The general pending-old gate avoids that prerequisite.

## 6. Next obligations and stopping boundary

The strongest remaining coverage task is a uniform estimate on pending old/new zonotopes and terminal supports at suitable finite cuts, together with prefix coverage. The existing short-horizon cone and the new tail certificates do not automatically overlap uniformly.

For metagraphs, classify the seven remaining center coefficients and retain literal ports when moving from counts to adjacency. Fixed-degree rationality does not determine bridges or connectivity. For enumeration, the incoming polynomial payload now supplies a genuine structural-state reduction; the next task is to reduce its coefficient work while preserving the terminal functional and comparing fairly with the scalar method.

All computations above keep their stated finite universes, controls, source pins and independent-review boundaries. The session's stopping point is the reviewed package set, with publication status verified separately, not a claim that any of these open general problems is solved.
