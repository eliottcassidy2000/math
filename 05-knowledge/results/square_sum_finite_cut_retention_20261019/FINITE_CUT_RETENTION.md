# Finite-cut transfer for retained prescribed paths

This is a structural transfer theorem for already certified affine gluing systems. It asserts retention and gives an explicit loss bound. It does not assert the existence of any missing gluing table, base family, or suitable host ports. No novelty claim is made.

## Certified setting

Let Q=q² be odd and h=(Q−1)/2. A compatible pair consists of Hamiltonian paths on consecutive integer intervals, with every common vertex at the same positional parity in both paths. A certified transition uses each offset c∈[−h,h] once, selecting one of the two input paths, replacing its vertex in one-based position j by Qv+(−1)^j c, cutting a specified finite set of internal edges, and gluing the resulting segments and bounded singleton vertices. The certificate proves exact coverage, square connector sums, connectivity, and (for pair outputs) common positional parity. Segment orientations are arbitrary. A final transition may instead produce one cycle.

Fix a simple square-sum path P=(p_0,…,p_(r−1)), r≥2. At level l, write

    P_(l,t)=(Q^l p_j+(−1)^j t)_j.

Suppose every relevant input pair contains all paths P_(l,t), t∈K_l, where K_l⊆[−h_l,h_l] and h_l=(Q^l−1)/2. The assumption is simultaneous containment in both input members. The paths may occur with opposite traversal orientations.

## Exact retention rule

For each t∈K_l, its input starting-vertex parity gives a sign σ_t∈{−1,1} shared across both members. In the block at offset c, the lifted copy is P_(l+1,u) with

    u=Qt+σ_t c.

It remains intact unless at least one certified cut removes an edge of that copy. Define D_l as any fixed set of integer parameters containing every such damaged u over all admissible input pairs, all possible transition states, and every output word. Then every output contains simultaneously all P_(l+1,u) with

    K_(l+1) = (union over t∈K_l of [Qt−h,Qt+h]) minus D_l.

Every interval [Qt−h,Qt+h] has Q integers; intervals for distinct integer t are disjoint. Consequently

    |K_(l+1)| ≥ Q|K_l| − |D_l|.

The statement applies equally to reversed segments and a final-cycle transition. It follows because cuts and external glues do not alter edges internal to any retained segment. Shared positional parity is what makes input-member selection harmless; no common orientation is required.

## A finite cut-budget bound

Suppose there are M allowed pair-output certificate states and each output word cuts at most b internal edges. Suppose every cut is specified by an offset c and an input edge whose membership in the finite prescribed path family can be determined. Since the P_(l,t) are vertex-disjoint, an input cut edge belongs to at most one of them. For that edge and that offset, at most two parameters Qt+c and Qt−c need be excluded if the starting parity is unknown. Unioning over both output words and all M states gives a sufficient set D_l of size at most4Mb. Known parities or repeated cut specifications can make the exact union much smaller. For single-output final states the corresponding bound is2Mb.

These are bounds for a fixed finite cut list, not for arbitrary cuts whose input endpoints may vary with unbounded input size. Such a variable-cut rule requires its own finite description or uniform damage-set proof.

## Vertex-disjointness and squares

For fixed p_j the labels Q^l p_j±t lie in [Q^l p_j−h_l,Q^l p_j+h_l]. Distinct p_j give disjoint intervals, and t is injective within each coordinate. Thus all P_(l,t) are vertex-disjoint whenever P is simple and K_l⊆[−h_l,h_l]. Adjacent sums equal Q^l(p_j+p_(j+1)), hence are squares. Actual interval inclusion is guaranteed by the certified expansion; direct bounds can also be used.

## Thresholds and isolated bases

For systems using output residues h,…,h+Q−1, a common-path input supply for all n≥T transfers to all output sizes n≥QT+h. After l layers the sufficient threshold is Q^l T+h_l. With only one compatible base pair of shorter length b, the exact l-fold descendant condition instead gives the band

    Q^l b+h_l ≤ N ≤ Q^l(b+1)+h_l−1,

provided the certified constructor uses exactly those l layers. This separates all-size supply from an isolated finite base witness.

## The three established instances

1. Matching transfer: P=(1,3), Q49. The first B layer can damage only the offset0 copy, giving K_1=[−24,24] minus{0}. Later retained edges avoid1–3 and2–7 cuts, so D_l is empty. This recovers the missing-central-band matching family.

2. Geometric-band connected paths: P is the twelve-vertex tail of the shared B54 prefix. It avoids1–3; all lifted copies then avoid every later cut. Thus every D_l is empty and K_l=[−h_l,h_l]. The input supply is an isolated base, giving geometric bands rather than an eventual all-size statement.

3. All-size triples: the complete weak-to-strong B table supplies P=(1,3,6) for every input size at least2621. The first protected A lift has fixed starting sign−1; its possible cut offsets0 and−2 therefore give exactly the sufficient damage set D_0={0,2}. The final lift has no damage to these copies. This yields2303 fixed paths with parameters [−1200,1200] excluding[−24,24] and[74,122], at threshold6294221.

The theorem abstracts what the complete certificates already justify. General finite prescribed-path embedding would additionally require a compatible all-size family containing that path or a proved mechanism that creates it. Finite cut budgeting alone does not supply that missing existence step.

## Proven instances and dependencies

- [Complete eventual lift construction](../seeded_eventual_cycle_theorem_20261019/SEEDED_INTERVAL_THEOREM.md) supplies the certified typed semantics and all original states (packet68).
- [Simultaneous matchings](../square_sum_simultaneous_matching_20261019/SIMULTANEOUS_MATCHING_THEOREM.md) is the one-edge instance (packet72).
- [Connected-path geometric bands](../square_sum_connected_path_bands_20261019/CONNECTED_PATH_TRANSFER.md) uses a single base pair (packet73).
- [All-size connected paths](../square_sum_all_size_connected_paths_20261019/ALL_SIZE_CONNECTED_PATHS.md) uses the complete weak-to-strong table (packet75).

This packet introduces no new search or finite census. Its universal conclusion is the proved recurrence under the stated hypotheses; the linked independently replayed certificates establish those hypotheses for the concrete examples. Elementary graph surgery and the existing affine lift are credited as tools, not claimed as new algorithms.
