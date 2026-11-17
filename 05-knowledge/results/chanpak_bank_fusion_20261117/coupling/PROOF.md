# Migrating a unique endpoint through a variable contact catalog

This gives a literal unbounded supply of the **exact endpoint passport** needed by the reflected +2 seam, then a support bound for coupling it to the variable three-contact repair. The bound allows an arbitrary number of strong stages and some weak stages. It does not cover unrestricted arbitrary stages.

## 1. Strong rooted constructions with a unique ending path

For binary word u starting in 1 and ending in 0, make a tournament on a transitive spine with an extra vertex v, where bit 1 means v beats that vertex. This tournament is strong: the apex and the natural spine close a cycle. At the first spine vertex r, exactly one full path ends there. Indeed its sole in-neighbor is v, and the preceding spine vertices must appear in increasing order, whose last vertex beats v. Thus arbitrary internal bits are legal, not just the extreme-back-edge template.

More generally suppose the old tournament A has exactly one full path ending at marked r. Add a transitive spine B of length k and then a new marked vertex s. In the deletion of s put B→A. Let s-beats-B be any standard nonzero binary word u of length k. Let s beat every old vertex of A except r. The alternative u=0 means omit B entirely. This is a literal one-return connection between A and an arbitrary stage block S=B+s; its cut arcs are S→A except r→s.

The new marked end at s is unique: in its deletion B→A any full path ending at the necessary predecessor r goes through the unique B path and the unique path on A ending r. If A is strong of order at least three, the new tournament remains strong. For nonempty B, s reaches its first vertex and then the whole spine and old A; every old A vertex reaches r→s. If B is empty, s beats some old non-r vertex and the same argument works.

There is also a closed scalar/end transfer. Write h=H(A), R for optional two-slot covers whose first path is nonempty and ends at r, and S_0 for full paths starting at r. Put F(u) for the Hamilton count of B plus an unconnected apex with contact word u (F(0)=1), and RC_k for reverse complement in exactly k bits. Then

    h'   = h F(u)+(u+1)R−[u even]S_0,
    R'   = 2^k R+h(1+RC_k(u)),
    S_0' = h if u>0, else h−S_0,
    e'   = 1.

For the first line, in the old contact block A the selected end must be r, weighted sum R and signed even-minus-odd sum −S_0 by the signed endpoint identity. The local odd-cycle increment through s using only A is R−S_0, since any even selected run ending r begins elsewhere. The mixed B–A term is uR+[u odd]S_0, and the B-only increment contributes h(F(u)−1). For R', the first slot ending s is either singleton (h), enters from r (2^k R choices), or contains no A vertices and enters s from a bit-zero position j of B (2^j choices, times h). The full start at s enters the first vertex of B if present, giving h; without B it can start the old path at any vertex but r. These encompass all cases. Direct subset DP verifies all four coordinates in depth-three chains.

Any such old A:root r can be used as the left inside flank in the Farey construction. Prefix it by an arbitrary transitive binary word u and use any positive even ordinary flank w on the other side. The first owns exactly one marked full end and the second owns zero, so reflection remains exactly +2. Varying the inner migration words gives actual new strong contact blocks, not a variation of the two exterior addresses.

## 2. Entropy certificate for a large restricted alphabet

At each migration allow a singleton stage, or a strong stage S formed from an arbitrary spine whose word starts 1 and ends 0. All old odd-cycle packs in A and all in S survive the seam; hence with p=P₂,

    p(new) >= p(A) p(S).

This follows by retaining the pure internal terms in the positive Rédei–Berge expansion at two colors. It counts actual source complexity by multiplicative costs, without treating equivalent histories as distinct outputs.

For spine width k+2, the extreme inversion of the first and last bits alone guarantees a two-cover cost

    b_k=2^(k+3)+2(3^k+1).

Indeed an inversion at distance d in a spine of length ℓ contributes 2^(ℓ−d)(3^(d−1)+1); expand the odd cycle through the apex by choosing an even subset of its d−1 intermediate vertices. For k>=2, the first interior bit 1 forces one further long inversion, and the last interior bit 0 forces another, each of cost c_k=4(3^(k−1)+1); these two conditions are independent. Thus the reciprocal sum over width k+2 is at most

    2^(k−2) [1/b_k + 2/(b_k+c_k) + 1/(b_k+2c_k)].

Widths two/three contribute 1/12 and 2/32. For k>=9 the weaker sum is at most (1/2)(2/3)^k. At exponent s=99/100 the finite reciprocal sum through k=8 is inflated by at most 28/25, since every lower-bound denominator used is less than (28/25)^100; this only uses lower bounds on p. The tail is bounded by (51/100)(17/25)^9/(1−17/25), because 2^(−s)<51/100 and 2/3^s<17/25. The singleton adds less than 51/100. Their explicit rational sum is below 0.987084. Therefore

    Σ_(labeled histories) P₂(source)^(−99/100) < infinity

by geometric summation. A fixed initial rooted seed is enough; a free initial one-apex seed also works if desired, since its active-word cost grows as 3^ℓ and there are only O(ℓ·2^ℓ) choices including the root.

Weak stages can be admitted at bounded frequency. A nontrivial standard binary word ending in 1 is weak and consists of a shorter active strong stage plus a trailing ordinal singleton run, or is all ones. Each tail singleton multiplies p by 2, so the total weak weight is below

    (0.477084+0.51)·0.51/(1−0.51) < 1.03.

The first 0.477084 is the explicit upper bound for the strong contribution alone in the preceding rational sum. A fixed number of weak labels interspersed among unlimited strict labels preserves convergence. Even arbitrarily many weak labels are allowed provided at most 1/1000 of all stages are weak: weighting a weak letter additionally by t=1/1000 bounds the admissible length-j language weight by

    [1000^(1/1000)(0.987084+0.001·1.03)]^j,

a geometric rate strictly below 1 (1000^(1/1000)<1.007 by exact rational powers). This is an upper bound on physical output histories; collisions only improve it.

## 3. Applying the catalog to the repaired-count endpoints

For a fixed strong migrated source A rooted r, let a=H(A), p=P₂(A), L be its weighted starts at r, R_c the weighted ends outside r when the outer apex beats only r, and d its one-block apex increment. The usual marked interval law with transitive prefix and opposite even ordinary flank gives on odd u

    H=C_0+h(F(u)+F(w))+A_0 u+B_0 w+Puw,
    h=5a, P=256p, A_0=5R_c+142p, B_0=142a+256L.

All these are literal reads. Coefficients other than the product are O(P), with the absolute source bounded below by Puw from positive count contributions. In particular a≤p/2, L,R_c≤p, and d≤2p, the latter by choosing the even marked path in an optional two-slot cover.

Use the inherited three-contact pattern at dyadic G>=4, odd suffix r_0, complement s=G−1−r_0, and positive near-balance index z=x−y. Then with J=PG²,

    Delta = 3h y − (64J−h)z − B_(r_0),
    B_(r_0)=16J+8PG r_0+4PG+6h−4hG+8h rev(r_0)−4G(A_0−B_0).

This is the exact previous two-scale law, which applies to these asymmetric linear coefficients without needing the special R_m identities. Uniformly B_(r_0) is bounded above and below by constant positive multiples of J: the possibly negative linear corrections are at most 4G(h+327p), much smaller than 16J. Hence if |Delta|≤J, then y is Θ(J(z+1)/h) and the count on either side is at least a constant times P³G⁶(z+1)²/h². Reflection still adds exactly two at both endpoints by the unique-end passport.

For a fixed source and given G,r_0,z there are O(J/h) possible y's in this full displacement window. At most O(G) suffixes qualify and at most O(sqrt(X)h/(P^(3/2)G³)) z's can yield any endpoint ≤X. The product is O(sqrt(X)/sqrt(P)) per dyadic G. Both P=O(X) and G≤O(X^(1/6)) follow from the size bound and h≤P. There are O(log X) scales. Summing over the physically enlarged migration catalog using s=99/100 gives

    # endpoints <= X   = O(X^(99/100) log X).

Indeed sqrt(X)/sqrt(P) ≤ O(X^s P^(−s)) for P=O(X), and the history series in §2 converges. This includes all repairs of the stated form, both tracks and both endpoints, all congruence-admissible or inadmissible phase selections (we upper count them), and an arbitrary depth of the specified source enlargement. It is a support obstruction, not another finite census and not an obstruction to all targeted constructions. Unrestricted weak-stage frequency, multiple exposed roots, and a different placement grammar remain outside it.

Replay: `python replay.py` and `python -O replay.py`; exact rational convergence bounds and 125 literal subset-DP chains.
