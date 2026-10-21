# Incoming stopped-block envelope: proof review and exact reconciliation

2026-10-21. **PASS:** conditional block-envelope proof and monotonicity. **PASS:** exact composition with the published individual-source joint-zonotope gate. Native numerical scope is saved-receipt arithmetic, not an additional native source census. Uniform all-N coverage and arbitrary-source entry remain OPEN. Repository read-only; no Git or remote write was performed. Review artifacts are outside the repository.

## 1. Source, carrier, and exact future semantics

The incoming packet uses the same fixed primitive source carrier below 256N, native primitive map, 627-point first-stop frontier, clock sign (-1)^clock, old/new split at 4N, and native-depth Nyquist energy as our published route-free gate. An old source keeps its old-scale charge at scale64N. The quadratic observable is

    delta(w,h) = Q(w+h)-Q(w),
    kappa(a,b) = delta(a+b)-delta(a)-delta(b).

At cut D, each unresolved source has taken exactly D primitive edges without hitting the frontier. Grouping by its actual current node is therefore lawful: all group members share every additional edge, extra clock, future first-hit depth relative to D, and terminal, if a first hit occurs. Past clock signs may differ and are already included in the separate sums w and h. Subsequent transport multiplies BOTH sums by the same sign. A future cut includes the entire block with one shared sign or includes none of it. No terminal, future clock, or eventual entry is guessed.

This is stronger information than merely grouping equal vectors. Equal vectors at distinct nodes retain distinct blocks. The implementation's vector-type compression only accelerates the pair arithmetic and includes the binomial internal-pair term. Source multiplicity does not become one physical source. Grouping equal nodes from different native depths would be invalid without an additional synchronization argument; the implementation does not do that.

Both channels must vanish for free erasure. For b=((1,0),(-1,0)), the total charge is zero but delta(b)=-1, and its interaction with an old-only unit block is -2. Deleting such a mixed block changes the old baseline. Actual both-channel-zero fibres have identically zero contribution at every future native cut even if no constituent ever arrives. Erasing their numerical aggregate does not assert entry. The decoder (N,D,s), deterministic rule and frontier retain reconstructible literal source prefixes, including full clocks and compatible registers; aggregate rows alone are not a lossless source ledger or an all-frequency certificate.

## 2. Envelope and monotonicity

Let C be all arrived terminal blocks (zero ports may be included), B the distinct pending node blocks, Delta=sum delta(c),

    m_b = max(0, max_c |kappa(c,b)|),
    L_b = max(0, delta(b)+m_b),
    U = Delta + sum_b L_b + sum_{b<c} |kappa(b,c)|.

At any future cut, polarization gives Delta plus each included block's delta, its signed interaction with its selected arrived terminal, and exactly one pair interaction for each pair that lands at the same terminal. Included and absent local contributions are both at most L_b. Absolute pair costs cover the pair terms whether the pair is present, absent or split across terminals. Thus every finite future Delta_h<=U. The enlarged placements are proof comparisons, not native route claims.

Signed transport preserves delta and absolute kappa. For a merge a,b, m_(a+b)<=m_a+m_b and

    L_(a+b) <= L_a+L_b+|kappa(a,b)|.

All other pair costs decrease by triangle inequality. For signed arrival epsilon*a at terminal t, the change in Delta is delta(a)+epsilon*kappa(c_t,a)<=L_a. Every remaining m_b, and therefore L_b, increases by at most |kappa(a,b)|. These are exactly the disappearing pair costs. Both-channel-zero erasure changes nothing. Simultaneous transitions can be decomposed into these operations, so U_(D+1)<=U_D.

This proof is sound even with negative local delta. Positive parts cannot be dropped. The supplied producer preserves source signs, tests first arrival before subsequent steps, and removes only both-channel-zero current-node groups. No ROOT continuation or unresolved route resolver is used by its stopped-prefix routine.

## 3. Exact reconciliation with the joint-zonotope gate

On the same actual carrier and cut, let B_D denote our published bound

    B_D = Delta_D + G(F_old,F_new) + L_old + L_new,

where G is the full independent old/new signed-source zonotope maximum. The raw unfinished source multisets must remain authenticated, including sources whose synchronous zero aggregate is omitted from the block consumer. They can be reconstructed from source addresses; the compressed aggregate alone does not supply their original inventory.

Both U_D and B_D bound the identical physical Delta_h for every h>=D. Hence

    Delta_h <= min(U_D,B_D) for h>=D.

Neither bound universally dominates the other. Three separate new blocks with charges (1,0),(0,1),(-1,-1), and no arrived charge, have U=6 but B=4: independent pair absolute values lose a sign frustration. Conversely two opposite new sources already coalesced into a synchronous zero block have U=0 but the raw individual-source B=4: source-independent relaxation loses forced cancellation. These are algebraic examples, not claims about native occurrences.

For any authenticated cuts d<=D with an all-future valid bound B_d, define

    H_D = min_{d<=D, evaluated} B_d,
    K_D = min(U_D,H_D).

H_D and K_D are nonincreasing as actual prefixes advance and new valid certificates are added. Every h>=D is in every included certificate's validity domain, so Delta_h<=H_D and <=K_D. Our raw B_D need not be monotone: in the three-charge example B rises from4 to6 after one unit charge arrives, but the old bound4 remains valid for all later cuts. A successful tail certificate does not cover preceding cuts. For an all-horizon conclusion, authenticate those earlier cuts separately, or use an earlier certificate whose domain covers them. No future arrival is presumed by retaining a past certificate.

Do NOT directly apply the individual-source full-sign maximum argument to coupled mixed blocks. Delta is indefinite on four coordinates. For a single b=((1,0),(-1,0)), every full common-sign sum has delta=-1 but optional absence contributes0. The random-completion proof for the source zonotope relied on the nonnegative new-source norm added in expectation, which fails for an arbitrary coupled mixed block. Safe options are the proved U, the original individual-source B, or a separately proved optional-block optimizer. Taking min of the two existing valid bounds needs no new relaxation or change of source ownership.

## 4. Ownership, flow compression, and no energy gain

For any single authenticated ownership partition,

    Delta_h=A_h+I_h, A_h+M_h=153N, I_h=R'_h+J'_h.

Any direct bound C>=Delta_h implies R'_h<=C-A_h-J'_h. Once C<=153N this gives exactly R'_h<=M_h-J'_h, not the stronger sufficient R'<=M. Changing a partition moves R' and J' oppositely. Combining independently proved upper bounds uses their minimum, never the sum of overlapping negative credits.

The four-point model's points are (0,0),(0,1),(-1,-1),(1,0). Its twelve directed differences are the native palette; zero selected charge is equivalent to divergence t*(-3,1,1,1). Ordinary circulation alone misses nonzero t. For the present reconciliation it is enough to check that a mature selected zero subset has the same terminal, actually arrived signed source rows, integral capacities and no repeated source ID. It does not change H or physical Delta. If Z is its removed diagonal norm and Gamma=delta(W,H)-(A-Z), then Gamma-Z=delta(W,H)-A identically. Improving Z raises Gamma by the same amount. At a pending common node, removing an authenticated new-only zero subset leaves (w,h), hence U, unchanged. A mixed-total-zero subset cannot be erased unless both channels vanish. Cut-indexed optimal inventories may be reassigned but their credits cannot be accumulated across cuts. The non-nesting example in the incoming packet correctly explains why an optimal irreversible enlarged packet bank need not exist.

No new flow optimality census is claimed here. The seven cut inequalities, residual-potential interpretation, and adjacent-integer concavity certificate were inspected for consistency; the task's no-energy-gain assertion follows from the exact identity independently of optimizer quality.

## 5. Native scope and numerical differences

Incoming saved U rows:

- N4096,D105: 195 pending nonzero nodes, 373 sources represented in those nodes; Delta421276, local74675, pairs61213, U557164, budget626688.
- N8192,D110: 348 nodes, 689 represented sources; Delta793715, local153217, pairs191475, U1138407, budget1253376.
- N32768,D130: 760 nodes, 1401 represented sources; Delta2838302, local651015, pairs752165, U4241482, budget5013504.

Our published individual-source rows are N4096,D105,B617824 (386 unfinished nonzero new sources); N8192,D112,B1211671 (615); N16385,D121,B2470224 (961); N32768,D131,B4763921 (1320). At N4096 the same-cut difference is legitimate: the block consumer has erased synchronous zero groups, whereas the individual-source consumer retains their unfinished original sources. At the other two shared scales, cut depths differ; the listed inequalities are not a same-state dominance comparison. The incoming stopping cuts are earlier at N8192 and32768. Our N16385 is absent from the incoming production list, so neither packet subsumes the other's entire finite scale list. N32768 is an incoming additional scale relative to its older contiguous census, but not new relative to our already-published route-free scale list.

The incoming production receipt claims a full bounded-prefix replay at all three successful scales. This review independently recalculates its four envelopes (including N7) from literal pending pairs, checks current-node uniqueness and the saved earlier-cut predicate, and reviews the prefix implementation; it does not rerun its multi-million source checks. Our published independent native replay covers N4096 and8192, not the two larger producer scales. Do not blur producer, prior independent, and present arithmetic-only scope.

## 6. Reproduction and stopping boundary

Run from the workspace root:

    python3 session22-incoming-tail-review/audit.py --repo-root /path/to/math
    python3 -O session22-incoming-tail-review/audit.py --repo-root /path/to/math

The independent checker imports no producer. It exhausts all two-block tuples in {-1,0,1}^4 against optional common-sign placement at three ports, and tests signed arrivals, merging, transport, polarization, both comparator directions, raw-B nonmonotonicity, the invalid coupled-full-sign shortcut, and no-energy-gain identities. It recomputes the saved envelopes literally rather than using type grouping. All checks are explicit and remain active under -O. Source pins identify the exact files read; no finite control proves scale-uniform stopping.

The missing theorem remains a successful cutoff and prefix coverage for every intended N. Neither monotonicity, a running minimum, a finite native scale list, nor a lifted mature subset supplies that theorem or arbitrary-source Collatz entry.

## Freshness and publication-portability addendum

Incoming commit 5872b76a adds section8 to the stopped-block README. Its explicitly attributed comparison with the route-free packet and its nonrecycling selector are consistent with section3 of this review: retain authenticated past-cut all-future bounds, take a minimum with monotone U, and separately cover the earlier cuts. The same-cut difference617824-557164=60660 is correct; no general dominance is claimed. The source-multiset reconstruction guard is explicit. This is a reconciliation of parallel incoming mechanisms, not a claim that the selector originated in this review. Section8's scale-obstruction discussion retains the distinction between within-scale monotonicity and cross-scale contraction, and between integer clock increments and primitive depth. Its additional native obstruction figures are attributed to that separate packet, not independently replayed here.

The checker now requires an explicit --repo-root path and no longer assumes a sibling math checkout. It reads the same pinned receipt and writes only its JSON result to stdout. Normal and optimized replays pass with byte-identical results. The refreshed source pin supersedes only the prior stopped-block README pin; the prior value is retained in the pin history. No additional native sweep or mathematical expansion was performed.
