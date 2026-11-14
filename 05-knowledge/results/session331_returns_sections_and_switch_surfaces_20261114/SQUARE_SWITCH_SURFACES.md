# Switch surfaces: exact circuit counts and recursive three-state square-sum families

## What the edge identity does not tell us

Fix a spanning cycle in the square-sum graph on 1,…,N. A proposed two-switch cuts two cycle edges and reconnects their four endpoints by one of the other two pairings; require the replacement edges to exist in the square-sum graph. Take a list of such proposals with pairwise disjoint vertex support. Every subset S has a distinct literal spanning 2-factor, not necessarily a single cycle.

Associate to the list a symmetric matrix A over F2. Its diagonal entry δ_i is 1 exactly when switch i alone leaves one circuit. Off-diagonal A_ij is 1 exactly when the two removed edges of i alternate with the two removed edges of j around the original cycle. Then the circuit count for any S is exactly

    b(S) = 1 + nullity(A[S]).

This is a specialization of the classical Cohn–Lempel circuit-nullity mechanism and its [undirected extension](https://arxiv.org/abs/0903.4405), not a new general circuit-nullity theorem. The specialization to square-sum switches, physical witnesses, and recursive products below are the contributions here.

## Disk-and-band proof and the actual Gauss–Bonnet connection

Take a disk with the reference tour drawn on its boundary. For each selected switch attach a strip at the two short boundary intervals representing its two removed edges. The long sides of the strip perform one of the two reconnections. An untwisted band on the disk alone gives an annulus and two circuits; a half-twisted band gives a Möbius strip and one. Thus the twist is exactly δ_i. The entire resulting boundary traverses the switched 2-factor, and the surface is connected with Euler characteristic 1−|S|.

Basis loops over F2 go once through each band and return along a chord in the disk. Two such loops pair to 1 iff their chords cross, and a loop has self-pairing 1 iff its band is twisted. The pairing matrix is precisely A[S]. Capping all b boundary circles by disks gives a closed surface with nondegenerate mod-2 intersection pairing. The map on first homology has kernel spanned by the boundary classes with their one sum relation (by the relative homology exact sequence). Its dimension is b−1, proving the claimed nullity formula.

In the orientable case rank(A[S]) is twice the genus. In the nonorientable case it is the number of crosscaps. Thus the usual Gauss–Bonnet total 2πχ = 2π(1−|S|) by itself cannot distinguish different circuits: the same number of selected switches can spend the topological budget on handles or crosscaps instead of extra boundary. The intersection form retains that missing distinction. For two interlaced switches the three basic matrices (diagonal choices (1,1), (1,0), (0,0)) give respectively 2,1,1 circuits for the double selection; changing the order of the four cut locations alters that result without changing the local identities.

The repository's [fence density corner lemma](../procgen_fencelim_20260926_fence_density.md#3-the-corner-lemma) pushes boundaries inward and charges total turning by 2π(1−holes). The [chessboard ring note](../glued_chessboard_rings_scaffolds_20261006.md#6-theorem-3-curvature-shows-up-as-ring-growth-proved-finite-exact-check) distinguishes gluing and shell closure from a raw cone-angle count. Here boundary circuits play a similarly visible role, but the actual needed datum is the mod-2 pairing of switch cores.

## Exact reanchoring and addition versus products

If A[S] is invertible, switch S and use the new Hamilton cycle as anchor. Relative to it the very same switch palette (each switch now toggles backward if already used) has matrix equal to the principal pivot transform A*S. In blocks relative to S and its complement, write original A as [[P,Q],[Qᵀ,R]]. The new matrix is

    [[P⁻¹, P⁻¹Q], [QᵀP⁻¹, R+QᵀP⁻¹Q]].

Indeed in the relation Ax=y, exchanging x_i and y_i on S produces this matrix. A vector in the kernel of the new principal block on T is represented by x'_outside T=0, y'_inside T=0. Exchanging coordinates turns those conditions into x_outside(S symmetric-difference T)=0, y_inside(S symmetric-difference T)=0. Thus its nullity equals that of the old block on S symmetric-difference T, which is also the nullity determined by the new anchor's actual circuits. Nullities of singleton and pair principal blocks determine a symmetric binary matrix: the singleton determines its diagonal, and for a 2×2 matrix det = δ_iδ_j + A_ij. This proves the reanchoring rule. A single viable switch toggles all off-diagonal relations between its neighbors and also their diagonal bits.

There is a second exact reading of det(A[S]): it is the parity of matchings of the interlacement graph on S whose unmatched vertices all have δ=1. Terms with a permutation cycle of length at least 3 cancel in pairs by reversing a canonically chosen such cycle; every remaining transposition contributes one off-diagonal entry squared, and unmatched indices contribute their diagonal factors. Thus disjoint matching pairings produce multiplicative terms, but connectivity is their sum in F2. If the interlacement graph of a fixed palette is a path and all diagonals are 1, each selected consecutive run must have length not ≡2 mod 3. Indeed its determinant satisfies d_j=d_(j−1)+d_(j−2) over F2, with d_0=d_1=1. The exact generating function in path-palette length k is (1+z)/(1−z−z²−z³): decompose binary selection strings into zero gaps and positive runs with generating function (z+z³)/(1−z³). If instead all diagonals are zero, every selected run must have even length and the generating function is 1/(1−z−z²). These conditional path-interface laws do not assume a physical square-sum palette that realizes every path size.

## Same 99-vertex universe, same square identity, opposite geometry

The literal control paths in the replay both run from 1 to 3, with a common closing edge 3–1. In each there are two disjoint, interlaced switch proposals. All four proposals use exactly the identity

    3² + 11² = 7² + 9² = 130.

The identity comes from products, not just search: 130=2·5·13, and (1+i)(2+i)(3+2i)=−3+11i while (1+i)(2+i)(3−2i)=9+7i. Conjugating one Gaussian factor changes the pair of squares at fixed norm. The product identity supplies legal local edges but does not determine the topology of their placement.

* For path C, cut pairs (9,66) and (48,97), indexed by edge start from 0. Both band twists are 1, so the literal circuit counts for no, first, second and both are (1,1,1,2). One cannot simultaneously use the individually legal alternatives.
* For path R, cut pairs (0,56) and (39,81). Both band twists are 0, so the corresponding circuit counts are (1,2,2,1). Neither switch alone works, yet the combination repairs connectivity.

The same literal proposal occurs in both controls: replace {2–7,42–79} by {7–42,79–2}. It is viable in C, where the original tour direction is 7→2 and 42→79, but splitting in R, where the tour has 2→7 and 42→79. The S330 height rule gives eight physical four-cycles for this identity, at vertices (t,9−t,t+40,81−t), 1≤t≤8, with one birth at each height 73,…,80. Availability of this local motif predates its spanning use here.

In both controls, every changed edge except the untouched closing edge still joins opposite positions of the original 1-to-3 path. In C the cut-position gaps are odd for the reversal pairing; in R the gaps are even for the splitting pairing. Therefore both palettes preserve the full vertex parity labelling as an edge property, even when their resulting graphs are disconnected. This is stronger than comparing only the numerical identity.

## Physical unbounded products, not just abstract interlacement words

Use the [Gerbicz 25-block fixed-port return](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md#6-gerbiczs-25-fold-blow-up-s6-proved-by-a-finite-junction-check), also analyzed in [the preceding square reflection packet](../session329_activation_and_regeneration_20261114/SQUARE_SUM_SWITCH.md). For each offset c, label the copied vertices x by 25x+c or 25x−c according to their parity position on the unmodified anchor path. Each new switched edge joins opposite such parities, so its sum multiplies by 25 and remains a square. The chart's partition remains exact because it is defined from the anchor's parity, and the junction edges are unchanged. This applies even to choices that split off a cycle inside a block; we lift the actual edges using the anchor chart, rather than incorrectly assuming that every toggled option is itself a path.

At each iteration the copied interval of tour positions is disjoint from all other copied intervals. Therefore its two chords keep their internal interlacement, and chords in different leaf copies do not interlace. Old closing edges are never selected as cuts. At depth d the size is (199·25^d−1)/2 and there are k=25^d copies, so the signed interlacement matrix is a direct sum of k two-by-two blocks. Every switch choice gives a literal distinct spanning 2-factor. The numbers of circuits therefore have generating polynomial

    C-palette:  x(3+x)^k,       exactly 3^k Hamilton cycles,
    R-palette:  x(2+2x)^k,      exactly 2^k Hamilton cycles.

These are counts within these named constructed families, not global Hamilton-cycle counts. At a single lift one can also choose the C or R anchor by each matched offset pair ±c and independently at center c=0. The matched pair uses the same parity signature, so the lift remains bijective on the labels. With u C-blocks and v R-blocks it has polynomial x(3+x)^u(2+2x)^v, u+v=25, for that specified arrangement. Do not add different anchor arrangements without an injectivity check across arrangements.

The analogue for a tournament cutoff is exact but limited: a repeatable “good move” is not just a count or an equal-sum identity. Joining several must preserve the relevant interface on all simultaneous moves, or carry a computable rule for how the interface changes. Here the interface and its change are the signed interlacement matrix and its principal pivot. The resulting size progression is sparse; none of this implies a Chan–Pak cutoff.

## An exact directed bridge to tournament path moves

A directed successor-swap replaces a→b,c→d by a→d,c→b, provided the latter arcs exist. Along a directed tour it is the untwisted reconnection. On a Hamilton path add a formal closing edge, which is never allowed to be removed, then apply the same circuit law. Extra directed circuits are precisely those besides the one containing the formal edge. Consequently for a batch of disjoint directed successor-swaps all diagonals are zero. Such a batch retains a single Hamilton path with the same endpoints iff the interlacement graph has an odd number of perfect matchings. In particular every viable batch has even size.

For a literal tournament example, take the path 0→1→…→7 and include 0→5,4→1,2→7,6→3. No specified arc conflicts, so complete the remaining pairs arbitrarily as a tournament. The swaps on cut-index pairs (0,4),(2,6) individually split the path into one shorter directed path and one directed cycle. Their interlacement allows both together, producing 0→5→6→3→4→1→2→7. This solves simultaneous path legality, not the numerical effect of changing arcs on ALL Hamilton paths in the tournament.
