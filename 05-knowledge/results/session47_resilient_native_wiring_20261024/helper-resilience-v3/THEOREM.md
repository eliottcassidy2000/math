# Resilient native charge banks: 46 helpers, exact failure sets, and a ten-source alternative

Status: PROVED conditional all-size controller theorem; FINITE-EXACT selected order-eight supply. This is not an all-order native-supply theorem. The incoming 48-helper claim was not present as a frozen certificate in the available `compatible_charged_helpers_20261023` API. We construct an independent 48-helper bank and a stronger 46-helper bank from its authenticated 53-source inventory.

## 1. Exact deletion resilience is a transversal problem

Let H be a finite set of actual binary helper slots, with fixed full U and V records at each slot, label lambda(h) in a finite abelian group G, and exactly r helpers at V. Slot indices, not merely labels, distinguish capacity. For charge g define the r-uniform witness hypergraph

    E_g = {J subset H : |J|=r, sum(lambda(j), j in J)=g}.

Assume here that every listed record choice is compatible with the other records and immutable background. If shared constraints remove some choices, replace E_g with the genuinely feasible family for the current core record; do not use the unrestricted charge table.

For D subset H, forcing the helpers in D to remain U retains charge g exactly when some J in E_g is disjoint from D. Thus D destroys g exactly when D is a transversal of E_g, meaning it meets every witness. With tau(E_g) the least transversal size, the bank survives every set of at most f forced-U helpers if and only if

    min_g tau(E_g) > f.                                      (1)

For a missing charge, put tau(empty)=0; if a witness is the empty set, no transversal exists and tau=+infinity. Deletion here removes V eligibility only: the helper slot and its occupied U record still exist at every arrival. Formula (1), not witness multiplicity alone, is the exact static criterion. A common helper in all representations makes a single deletion fatal even when there are many representations. In particular, single-U robustness is exactly the condition that E_g is nonempty and the intersection of all its members is empty, for every g.

A particularly transparent sufficient certificate is f+1 pairwise disjoint members of every E_g. A forbidden set with at most f vertices cannot meet all f+1 disjoint witnesses. This sufficient packing test need not be necessary. The three edges of a triangle have empty total intersection but no two disjoint edges.

For 1<=r<=h and 0<=f<=h, counting gives necessary but usually weaker lower bounds: choose(h-f,r)>=|G|, and choose(h,r)>=(f+1)|G|, since each charge needs at least f+1 representations. For r=3, f=1 and |G|=2048, the second implies h>=31. Our 46-helper construction is not claimed globally smallest.

## 2. Restrictions on both binary choices

Let F and D be disjoint sets of helpers forced to V and U respectively. The possible charges at exact total helper occupancy r are exactly

    lambda(F) + S_(r-|F|)(H minus (F union D)),               (2)

where S_k uses k distinct actual slots. This follows by writing each permitted J uniquely as F union K. It retains exact occupancy and is valid in any abelian group; compatibility must still be checked separately.

In particular one forced-V helper at occupancy three needs complete pair-charge coverage in its complement. Merely counting those pairs would require h>=66 for |G|=2048. The full 72-new-source alphabet has the stronger native obstruction that its pair sums cover only 1468 charges and omit charge 2. For any specified forced helper v, total charge lambda(v)+2 is unavailable at occupancy three. This remains impossible even before considering native resource conflicts.

Conversely, single-U resilience at occupancy three gives complete occupancy-four coverage containing any prescribed one helper at V: subtract its label, use a residual triple avoiding that helper, and adjoin it. This changes the declared constant helper occupancy from three to four. It is a separate packet family, not a zero-cost operation inside the occupancy-three controller.

## 3. Causal lifting and what “at each arrival” means

Keep a fixed total charge q, a core move graph whose vertices are actual compatible constant-margin core records I and whose edges change at most c slots, and fixed helper occupancy r. All helper records remain the same full records throughout; a hidden kernel-label change at an unchanged bit is not free. Write t=min(r,h-r).

At destination arrival i, a restriction set D_i is revealed before the simultaneous packet event. It may depend on all previous choices. The controller is allowed to choose the new helper configuration after seeing D_i. Suppose every required core state I and helper debt g has an actually compatible witness J disjoint from D_i. Then choose such a J_i with charge q-lambda(I_i). No lookahead or common section across arrivals is needed. Since any two r-subsets of H differ in at most 2t slots, every lifted event changes at most

    c + 2t slots.                                           (3)

Both endpoint arrivals retain all required physical guards and capacities by hypothesis. If c=2, (3) gives eight for r=3 and ten for r=4. A helper-only repair costs at most 2t. A Johnson core path uses one event per paired core exchange; changing restrictions can also require helper-only events when the core waits.

This is an arrival-indexed, atomic replacement contract. It is not an instantaneous continuous failure model: D_i may contain a helper that was V at arrival i-1, provided the simultaneous event puts it at U at arrival i. It does not assert that every intermediate one-owner or one-resource transition is legal. If restrictions are cumulative, their union, not merely the number announced at one time, must satisfy the resilience bound.

These routes use capacity-one reversible borrowing: the same released alternative may be reoccupied at a later arrival, while no resource is multiply occupied at any one arrival. Under the stricter temporal rule that a released resource may NEVER be occupied again, any changed helper cannot return to its initial record; this controller and its exact-return guarantee do not apply.

The reveal order matters. If an adversary may choose a forbidden helper only after the controller irrevocably commits its r>0 selected helpers, it can choose one of those helpers and invalidate that arrival. Static coverage does not defeat this stronger adversary.

Prescribed endpoint helper records are retained only if they are feasible under their endpoint restrictions. In particular exact return to an initial subset J_* is impossible if the final restriction intersects J_*. If final restrictions admit J_*, use J_* at that last arrival and the same distance proof applies. Forcing a bank to return to “all U” is also incompatible with its fixed positive occupancy. The computation explicitly returns all literal initial records with an admissible final restriction.

## 4. Why static coverage still does not imply arbitrary online reachability

The budget in (3) is part of the theorem. At a smaller helper transition budget, a complete restriction-resistant table can have disconnected feasible states. Take six indexed labels (0,0,0,1,2,3) in (C2)^2 at occupancy two. Every charge fibre has empty total witness intersection, so the bank survives every single-U restriction. Its charge-one fibre is exactly

    {0,3}, {1,3}, {2,3}, {4,5},

where the numbers here are slot indices. The first three states connect by two-slot exchanges; the last is isolated under two-slot moves, since it differs from each by four slots. Forcing slot 3 to U leaves the last state as the only choice and traps a controller currently in the first component. Coverage and bounded transition reachability are different. The full four-helper-change budget avoids this obstruction.

Likewise arbitrary changing shared-owner constraints need genuinely compatible E_(I,g), not just labels. Full kernel records, normalized ports, integer clock costs, all aggregate guard rows and the unchanged background belong to those feasibility predicates. Equation (1) applies to the actual witness families and does not manufacture the missing native supply.

## 5. A certified 46-helper eight-source controller

All labels and actual records come from the existing complete 53-helper API at `compatible_charged_helpers_20261023/supply/selected_supply.json`, with its original targets U=1965 and V=2234 and group

    Q = C128 + (C2)^4, |Q|=2048.

An integer charge code is 16*(X mod 128) plus the four Y parity bits in its low four positions, with Y1 in the least significant bit. All numerical charges below use this exact encoding.

`robust_bank46.json` supplies these 46 actual new sources:

    368,556,650,676,753,760,776,789,811,827,862,1182,
    1248,1369,1391,1395,1400,1493,1574,1615,1624,1650,
    1702,1770,1882,1952,2066,2135,2283,2319,2384,2393,
    2397,2439,2483,2616,2638,2677,2728,2752,2948,2983,
    3105,3111,3139,3187.

All 15,180 distinct triples are enumerated. Every charge occurs, with multiplicities from 2 to 17, and every charge's entire witness family has empty intersection. Therefore this bank tolerates any one helper forced U at each destination arrival, even when the restricted identity changes, under the causal contract above. Coupled to the protected independent 28-core supply, every constant-margin core route of two-source exchanges lifts with eight-source events and at most 50 primary source-slot owner-clock changes per event.

The 46-bank is inclusion-minimal for this single-U property: deleting any one of its slots produces a bank for which some additional single-U restriction kills a charge. The receipt provides every removed/forced pair and missing charge. This is not a proof that every 45-source bank fails. An independently constructed 48-bank is also frozen and checked; the 46 result improves that construction, not an unavailable external certificate.

The old 40-helper bank still covers every unrestricted charge but is not single-U robust. At charge 10 its two complete witnesses are {776,1182,2218} and {1182,2728,2752}; forcing 1182 to U kills both. This illustrates why minimum multiplicity alone is not enough.

The 28+46 selected packet uses 74 actual slots, occupies 74 normalized resources / 296 closed words at each arrival, and reserves 148 alternative resources / 592 closed words. The 28+48 packet analogously occupies 76 resources / 304 words and reserves 152 / 608. These are selected packets at tournament order eight, not tournaments on 74 or 76 vertices or newly freed capacity in the fully occupied atlas. The immutable remainder must avoid the complete chosen reservation.

## 6. Sharp occupancy-three failure and an occupancy-four alternative

The entire authenticated 53-helper supply fails after two U restrictions. Charge 305 has exactly the four triples

    {776,1182,2638}, {1391,1948,3139},
    {1400,3111,3139}, {2439,2728,3139}.

Every triple meets {776,3139}. Hence forcing those two sources to U kills charge 305 in the full 53 supply and in every subset of it. Thus the certified 46,48 and53 banks have exact pure-U deletion radius one at occupancy three. This impossibility is scoped to these actual 53 supplied binary menus, not all possible source paths or every larger alphabet.

At occupancy four there is a substantial positive alternative: `quad_packings.json` gives twelve pairwise disjoint quadruples of each of the 2048 charges in the full 53 bank. Thus ANY eleven helpers may be forced U at each destination arrival, with a ten-source core-exchange controller. These packings are exact integer masks checked against the original native potentials; no solver objective or floating-point result is accepted as proof. The bound eleven is a certified lower bound on deletion tolerance, not a claim that twelve restrictions cause failure.

The uniform number of pairwise disjoint witnesses is exactly twelve. Thirteen disjoint quadruples would cover52 of 53 helpers and leave a single source v unused. They would imply

    13q = lambda(H) - lambda(v).

Multiplication by13 is invertible on Q. Since the 53 selected source labels are distinct, this necessary condition can hold at only 53 charges; the other1995 cannot have thirteen disjoint quadruples. Twelve is therefore a sharp uniform PACKING parameter. This does not establish the minimum transversal number or exact deletion radius: a hypergraph can have a much larger transversal than matching.

Yet even occupancy four cannot tolerate every simultaneous one-V and one-U restriction: forcing776 to V and3139 to U at total charge1368 requires a residual triple of charge305 avoiding both, which does not exist. The complete census of 2756 directed distinct(1V,1U) pairs finds2724 with full coverage and exactly 32 failures, each losing exactly one charge. Equivalently16 unordered pairs of U restrictions each kill exactly one triple charge. The receipt retains all controls. It does tolerate every single-V restriction alone by Section2. Pure-U redundancy and mixed binary restrictions are therefore sharply different resource requirements.

The full 28+53 packet has 81 slots, 162 alternative normalized resources /648 closed words, and occupies 81 resources /324 words at every arrival. Helper record clock costs have histogram {4:3,5:13,6:18,7:19}; a two-core/eight-helper event costs at most 8+56=64 primary owner-clock occurrences. This is separate from its ten-source support count. If inverse/converse copies are charged as separate work, retain the corresponding multiplicity.

## 7. An exact weighted certificate: thirteen restrictions, with only one unit left unresolved

The disjoint-packing number is not the robustness radius. `fractional_weights.json` supplies a stronger exact certificate for every charge q: each listed actual quadruple e has a positive integer weight n_e, and an integer denominator d satisfies

    sum_(e containing v) n_e <= d for every actual helper v,
    sum_e n_e > 13d.

If a set D of at most thirteen helpers met every listed quadruple, then counting the total weight through D would give

    sum_e n_e <= sum_(v in D) sum_(e containing v) n_e
                <= |D|d <= 13d,

a contradiction. Thus every such D leaves at least one actually supplied quadruple untouched. This is a fully integer certificate of transversal number at least 14, and it proves tolerance of ANY thirteen U restrictions at every revealed destination arrival, still with four V helpers, ten-source events and the same 64-clock bound. All original compatibility and causal timing hypotheses remain unchanged.

The certificate contains 108,309 positively weighted actual quadruples. Its least total-weight ratio is 13,145,902 /1,000,001 at charge 540, strictly larger than 13. Floating-point optimization proposed weights only; the verifier accepts only integer edge weights, exact native charges, integer vertex-load inequalities and the strict total inequality. It imports no optimization package. No optimizer's optimality status is used.

There is also an explicit 15-helper obstruction. At charge 353, force the following actual helpers to U:

    556,776,811,1493,1624,1948,2135,2218,
    2384,2397,2439,2483,2616,2752,3111.

The checker enumerates all 139 quadruple representations of charge 353 and verifies that every one meets this set. Therefore the full 53-bank occupancy-four pure-U deletion radius is either 13 or 14: it survives every 13-restriction set, and some 15-restriction set fails. Whether every 14-restriction set survives remains OPEN. No solver's report that a particular hitting set is optimal is substituted for a proof.

This also exposes a precise limitation of this fractional method. Since each witness has four vertices and there are 53 helpers, every fractional matching with vertex loads at most 1 has total weight at most 53/4=13.25. Such certificates can force the integer transversal to be at least 14, but can never certify at least 15. Closing the last unit therefore requires a genuinely integral argument, not just a better fractional solver. The actual native bank already separates a sharp uniform disjoint packing of 12 from a proven deletion radius of at least 13.

## 8. Native and electrical scope; reproducibility

The inherited supply verifier authenticates every selected full binary path domain, all automorphism kernels, normalized resource fibres and literal black weight two. In particular source 2439 is self-converse; multiple canonical-source lifts of its resource do not supply extra capacity. Our checker independently replays all 162 alternative records from the original matrices and owner map, their common clocks, chosen original X phases, normalized inverse/converse words, and actual converse paths.

For fixed target margin, the complete paired-code theorem identifies equality in Q with preservation of every original dyadic aggregate physical-row current. U and V have equal paired charge and both have X=146. The six primitive current tests use X modulo 256 and the five radical potentials modulo 4; the paired theorem, not numerical analogy, establishes completeness at all dyadic precisions. The factor-two black multiplicity and target degree-diagonal rows are retained. Exact integer X currents, odd-primary guards and every ordinary individual-port phase are NOT claimed. The ordinary harmonic 4e_U modulo 8 still distinguishes U from V.

For its route replay the seven authenticated helpers outside the 46-bank stay at U as a disjoint immutable background. Removing this unchanged background leaves exactly the 74-slot packet above.

The checker tests all 2048 total charges along a forward/backward 14-exchange core path, deliberately choosing destination restrictions that evict a currently active helper. It retains all actual physical rows at each arrival, fixed margins, eight-source/50-clock bounds and exact final full-record return under an admissible last restriction. An additional 2048 occupancy-four route replays choose thirteen forbidden helpers including all four currently active ones, find an actual surviving quadruple, and return under an admissible final restriction; every physical-row guard and the ten-source/64-clock budget is checked. Its exhaustive charge-family tests prove the universal finite claims; these selected route replays are additional implementation controls, not a substitute for the causal lifting proof.

Run `python check.py --results /path/to/05-knowledge/results` and the identical command with `python -O`. Explicit exceptions retain every check under optimization. Run the inherited `compatible_charged_helpers_20261023/supply/check.py` normally and with `-O` for the complete native-domain authentication. Frozen receipt and manifest hashes distinguish this exact weighted third checkpoint from the original nine- and twelve-packing checkpoints.

Unrestricted all-size wiring, universal jointly compatible charged supply, globally minimum helper counts and whether the occupancy-four deletion radius is 13 or 14 remain open. The all-size content proved here is the conditional full-record transversal criterion and causal lifting theorem, not native existence at every order.
