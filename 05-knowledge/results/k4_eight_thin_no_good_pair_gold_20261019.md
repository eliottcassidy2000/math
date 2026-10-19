# An eight-thin twelve-point height-two poset with no good pair, but with Gold

Date: 2026-10-19. Status: exact single-example certificate and elementary structural proof. No minimum-order or novelty claim is made.

## Construction and sharp structural threshold

Take a lower antichain of six vertices and an upper antichain of six vertices. Number upper vertices 1,...,6. The lower missing-incidence sets are

    S1={1,3,5}, S2={1,4,6}, S3={2,3,6}, S4={2,4,5},
    M1={1,2}, M2={3,4}.

A lower vertex lies below precisely the upper vertices outside its missing set. All other pairs within each level are incomparable. The four S rows have incomparability degree 5+3=8 and the two M rows have degree 5+2=7. Upper vertices 1,...,4 have three missing lower neighbors and degree 8; upper vertices 5,6 have two and degree 7. Thus the maximum incomparability degree is exactly8. Every lower vertex has three or four successors and every upper vertex has three or four predecessors, so there are no isolated vertices.

This poset has **no good pair**, hence no very good pair, but has an exact half-probability pair and therefore weak Gold. In conjunction with the [seven-thin height-two theorem](seven_thin_height_two_very_good_gold_20261019.md), the universal threshold7 is sharp for both guaranteed very-good-pair and guaranteed good-pair existence in non-total height-at-most-two posets. This does not make eight the threshold for Gold: the present example satisfies Gold.

## Conceptual K4 proof of no very good pair

Regard the six upper vertices as the edges of K4 on vertices s1,s2,s3,s4, labeled

    1=s1s2, 2=s3s4, 3=s1s3, 4=s2s4, 5=s1s4, 6=s2s3.

The first four missing rows are the four vertex stars. The last two are two of the three perfect matchings. Two stars share one edge, so both directed differences have size 2. The two matching rows are disjoint, so both differences have size 2. A star and a matching intersect once, so their directed differences have sizes 2 and 1. No two rows are near, where near means both directed differences have size at most 1.

For columns, a missing set consists of the two endpoint stars of an edge, together with its selected-matching row if it belongs to one of the two chosen matchings. Disjoint edges already have two different endpoint stars each way. Adjacent edges have one different endpoint star each way; at least one of them belongs to a selected matching, since the two unselected edges are themselves disjoint. Adjacent edges cannot belong to the same matching. Consequently at least one directed column difference has size at least 2. No two columns are near either.

Within a level, a very good pair is precisely a near pair of missing sets, possibly in the dual. Across levels, equal downsets cannot occur because lower downsets are empty while every upper downset is nonempty. Equal upsets cannot occur either, by the dual argument. Hence there is no very good pair anywhere.

## Exact no-good-pair certificate

We use the definition recorded in the [attributed good-pair synthesis](good_pair_gold_synthesis_20261019.md): in P or its dual, an ordered incomparable pair (a,b) is good if

    D(a)⊆D(b), U(b)\U(a) is a chain, and p(a<b)≤1/2.

This is Zaguia's probability-sensitive good-pair definition, [primary paper, Definition 1](https://arxiv.org/pdf/1610.00809v3). Absence of very good pairs alone would not prove absence of good pairs; the probabilities must also be checked.

Label lower vertices 0,...,5 in the displayed row order and upper vertices 6,...,11 in numerical column order. Their strict predecessor masks are

    (0,0,0,0,0,0,44,35,26,21,54,57).

The ideal recurrence F(S)=sum F(S∪{v}) over currently available v, with F(P)=1, gives exactly 984,960 linear extensions. Forward ideal counts multiplied by suffix counts give all exact ordered comparison numerators. Every incomparable ordered pair is tested in both P and its dual for the two structural good-pair conditions. Exactly 48 orientations satisfy them. Their reduced probabilities, including multiplicity across both orientations of the poset, are

    2377/4560: 16 orientations;
    717/760:   16 orientations;
    109/114:   16 orientations.

All three numbers strictly exceed one half. Thus no good pair exists. The receipt lists every eligible ordered pair, not just the aggregate histogram. The dynamic program uses exact integers and tests comparability and the chain condition directly from predecessor/successor masks; it does not assume that only same-level candidates need consideration.

## An explicit half pair and weak Gold

The K4 vertex permutation (s1 s2)(s3 s4) preserves both chosen matchings. It induces the poset automorphism, in zero-based labels,

    (0 1)(2 3)(8 9)(10 11),

fixing 4,5,6,7. It interchanges lower vertices 0 and 1, so applying it to linear extensions bijects the outcomes0<1 and 1<0. Each has exactly 492,480 extensions, half of 984,960. The verifier checks this automorphism and also independently recomputes the extension count after imposing 0<1.

Compare 0 and 1, then repeat the now-determined comparison. On either realized path t1=t2=E/2, giving E=t1+t2. This is a weak Gold strategy in the repeated/determined-comparison convention of the cited synthesis. It needs no good pair. No claim is made that absence of good pairs obstructs Gold or ordinary balance.

## Fano comparison and scope

For comparison, the fourteen-point height-two poset whose seven missing rows are the seven lines of the Fano plane has degree 6+3=9 at every vertex. Distinct rows intersect once, and distinct columns also intersect once, so both directed differences have size 2. It has no same-level very good pairs; nonempty predecessor/successor sets exclude cross-level ones. Thus the complement-incidence Fano construction is already nine-thin, but the K4 example above gives the sharper eight-thin obstruction to guaranteed very good pairs. This comparison does not claim absence of good pairs in the Fano example.

Run `python -O ../experiments/verify_k4_eight_thin_no_good_pair_20261019.py` from the results folder. Normal and optimized receipts accompany the proof. This is a single explicit example, not a large census, and no minimum twelve-point order claim is made. Packet 112 and its seven-thin theorem are unchanged. General incidence constructions, automorphism half-probability arguments and exact ideal counting are not claimed as original; the good-pair definition and Gold convention are explicitly attributed above.
