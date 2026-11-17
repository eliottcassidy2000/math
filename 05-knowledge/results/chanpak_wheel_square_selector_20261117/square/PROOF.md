# A protected bridge in the parity-compatible square-sum ladder through 36

Let Q_n have vertex labels 1,...,n and an undirected edge for distinct x,y if x+y is a square. The archive's independently audited `collatz_mod6_20260922_w6_square_sum_hamiltonicity.md` reports path counts through 36, cycles beginning at 32, and Gerbicz's separate lift with a parity-compatible consecutive pair as the inductive data. Its 25-fold construction from the Q_35 chain is explicitly audited; it reports the 49-fold residue-complete construction from Gerbicz, not a new proof of it. The archive's protected-edge note `square_sum_affine_lift_protected_edge_obstruction_20261019.md` proves a different, conditional no-go for retaining an edge over one fixed entire residue range. Here the finite protected state is measured directly in Q_n.

## Exact residue-cut law

Given a path written in either allowed order with positions starting at 1, let A be the labels where position parity differs from the label's ordinary parity. Two paths on consecutive horizons are parity-compatible if all common labels have the same positional parity. The new vertex then necessarily occupies the parity that gained exactly one position, namely its own integer parity. Thus a whole consecutive compatible ladder is exactly a fixed anomaly set A throughout. Equivalently each Q_n must support a Hamilton path in the same colored-edge rule: keep only edges between opposite colors, color(x) = (x mod2) xor [x∈A], start the path in color 1, and restrict to [1,n]. Color classes must of course have the correct absolute sizes.

For normal vertices x outside A, an allowed normal-normal edge has odd square sum and hence x+y=1 mod8. Its endpoints lie in one of the four sheets {0,1},{2,7},{3,6},{4,5} modulo 8. Form two macrogroups U={0,1,4,5} and V={2,3,6,7}. Put each normal label in its native macrogroup. Put a flipped even label in its native macrogroup, and a flipped odd label in its opposite macrogroup.

**Claim.** In the colored-edge graph, the U–V cut consists exactly of the square-sum edges with both endpoints in A and opposite ordinary parity.

* Normal-normal: an odd square is 1 modulo 8, so the pair is in one original sheet and does not cross.
* Anomaly-normal: opposite colors mean equal ordinary parity, so the square sum is 0 or 4 modulo 8. For two even labels this gives the same native macrogroup; for two odd labels it gives different native macrogroups. The reassignment of anomalous odd labels makes both cases internal.
* Anomaly-anomaly: allowed colors are opposite iff ordinary parity is opposite. Their odd square sum puts their native classes in the same sheet; exactly the odd endpoint is now in the opposite macrogroup.

Hence at large n with both macrogroups populated, a compatible Hamilton path must use a fixed **finite menu** of bridges, and a singleton menu forces one particular edge in every path. A useful signed sidecar holds even for larger menus. Write `s(x)=2 color(x)−1` and `D_U=Σ_{x∈U} s(x)`. For every compatible path,

    2D_U = Σ_(path endpoints in U) s(x) + Σ_(used cross edges) s(its U endpoint).

Indeed each internal edge sums to zero signed degree mass. Thus simply counting bridge edges cannot substitute for compatible endpoint and sign data, the same distinction as retaining signed endpoint counts rather than merely their sum in a tournament.

## Exact 25–36 computation and a physical protected bridge

The unlabeled-direction path totals for n=25,...,36 are

    10, 12, 35, 52, 19, 20, 349, 392, 669, 4041, 17175, 12960,

matching the archive. Distinct anomaly-set counts are respectively

    2, 12, 20, 42, 8, 28, 101, 282, 223, 1526, 2615, 3262.

There is no coherent key across 25–36 or across 29–36. There is exactly one common key across 30–36:

    A={1,2,3,4,5,6,19,20,21,22,23,24,25,26,27,28,29,30}.

It occurs in 1,16,11,6,48,272,114 oriented paths respectively. The edge **19–30** crosses the two macrogroups and is used by *every* path carrying this key at *every* one of the seven horizons. The whole common forced-edge set is {7–18,9–27,19–30}; the first two do not cross macrogroups. Reproducible witness paths at each horizon and intersections are printed by `replay.cpp`.

The replay enumerates oriented paths, rejecting a partial path only when an unvisited vertex has no available neighbor, at least two have only one available neighbor (either would have to be the sole remaining terminal endpoint), or the unvisited induced graph cannot be traversed from the current vertex. It never assumes that an arbitrary degree-two vertex is interior. The first line's path count folds the two directions; anomaly keys and compatible counts intentionally retain orientations, since reversing an even-order path flips its positional color classes.

**Transfer discipline.** The finite phenomenon does not prove that one key survives indefinitely or that a square-sum bridge is equivalent to a tournament gadget. It identifies the constructive requirement: the analogue of the growing input is not merely an attained scalar, but a witnessed family all admitting a shared signed endpoint/bridge state. The Q_25–36 data specifically refute replacing this with one unrelated witness at each horizon.

## A literal all-height descendant of the protected contact

There is also a valid Hamilton path on Q_35 from 1 to 3 containing 19–30:

    1,8,28,21,15,34,2,7,18,31,33,16,9,27,22,14,35,29,
    20,5,11,25,24,12,4,32,17,19,30,6,10,26,23,13,3.

This is a different positional key from the seven-size ladder; do not identify them. It is in the inductive port class required by the audited 25-fold lift. That lift inserts all translated blocks T(c), including c=0, intact up to reversal, so it preserves every labeled old edge as its c=0 descendant. At depth m we obtain a Hamilton path of Q_{(71·25^m−1)/2} still from 1 to 3 and containing (19·25^m,30·25^m). The edge sum is 49·25^m=(7·5^m)^2, as required. More generally each occurrence has all 25 descendants, with internal opposite signs on the c shifts; adjacency survives block reversal.

This is an actual composable marked edge, not preservation of the original two labels across a full interval. `lift_edge.py` independently checks the new seed and replays the audited junction table through depth 3. In structural terms the viable replacement changes the scale and transports its port, rather than asking an unrelated absolute edge to fit after every enlargement.
