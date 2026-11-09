# Independent audit: scoped loop minimality, placement, and diameter

## Verdict and exact minimality scope

PASS: twelve is the minimum positive length of a **unit-gain directed closed moving-core walk rooted at a physical core/cut of order at most9**, using the compatible singleton boundary cocycle. The12-step example supplies the upper bound. The mathematical reductions and independently rebuilt finite universe below supply the lower bound.

This is not a global minimality theorem over larger physical cores, a shortest unit-switch theorem, or a shortest consecutive-count construction. Nonclosed paths between different physical cuts, stationary pure-T switches, and different construction mechanisms are outside this loop-minimality statement. In particular length-zero switches do not contradict it. “Unit loop” here means a cyclic cocycle gain ±2, with the matching word-end boundary choices that make it repeatable; it does not mean an arbitrary one-off flank choice producing a two-count difference.

## Quantization, content, and length reductions

Let a nonempty physical core have h odd and q even. For a closed length-L move, let B=2^L and X,Y be its two singleton words. Closing A/G gives

    (B−1)A=q A_X−h A_Y,
    (B−1)G=q G_Y−h G_X.

With compatible boundaries, the gain is

    Δ=A G_Y−A_X G+hK,

where K is an integer containing the word-T and endpoint terms. Multiplying the first two equations respectively by G and A and subtracting yields

    q(A G_Y−A_X G)=h(A G_X−A_Y G).

Therefore h divides qΔ, so h/gcd(h,q) divides Δ.

Every positive-length closed augmented walk has p=A mod2 at EVERY vertex: each forward edge has A'=2A+hy−qx congruent to y=p' mod2. Its edge weight

    w=Ay+hp y−xG'−hx(G' mod2)

is consequently even, because h is odd, A+p is even, and G'+(G' mod2) is even. Thus the gain is divisible by2h/gcd(h,q). A unit gain ±2 forces h|q.

For d=gcd(h,q,A,G), d is odd. The forward equations keep both A and G divisible by d, since2 is invertible modulo d. Every edge weight is divisible by d; cyclic edges are also even. Thus every gain is divisible by2d. Together these two divisibilities give the valid obstruction lcm(2d,2h/gcd(h,q)). Once h|q, a unit gain requires gcd(h,A,G)=1.

Modulo h in that saturated case, the closing equations say h divides (2^L−1)A and (2^L−1)G. Bézout with gcd(h,A,G)=1 implies h divides2^L−1. For h>1, ord_h(2) therefore divides L. This uses the joint content gcd, not separate coprimality of A or G with h.

The boundary-state omission in the producer is justified but should be explicit: a physical seed with p different from A mod2 cannot close after any positive number of steps. The independent replay nevertheless tests BOTH p values.

### The h=1 exception

The modular-order loop must not be applied at modulus1. A physical tournament with exactly one Hamilton path is transitive. One proof: any Hamilton path in an induced subtournament extends by inserting a new vertex, and choosing one insertion for each path is injective under deletion. Hence Hamilton counts do not decrease when vertices are added. A directed triangle has three Hamilton paths, so an H=1 tournament contains no directed triangle and is transitive.

A transitive core together with singleton padding is again transitive. Its A is the binary integer of the entire selected cut word in the unique transitive order. Equal A with equal total length therefore forces the same word and identical full state, so every physical closed-loop gain is zero. This separately excludes all h=1 cores of every order. The independent order1..9 enumeration also finds exactly one such class at each order and q=2^n.

## Independently rebuilt finite universe

The audit does NOT read SATURATED3.jsonl through SATURATED9.jsonl to build its seeds. `minimal_seeds.cpp` reads the complete class lists, constructs every graph, independently computes all Hamilton endpoint/start counts by forward path extension, verifies saved H, computes q, filters h>1 and h|q, and enumerates all cut signatures. This includes every topology and every cut.

Completeness of the supplied class lists has independently checked pairwise canonical nonduplication, automorphism orders, and exact orbit mass at every order1 through9. MASS1_TO7.txt, MASS8.txt and MASS9.txt contain these receipts. The new filtering is therefore not conditional on an unverified saturated-row list.

Rebuilt totals match exactly:302644 distinct saturated seven-coordinate signatures summed over the orders, of which230568 have content one. The only surviving (h,q) banks with ord_h(2)<12, and their numbers of distinct physical A/G seeds, are:

    (33,660):176
    (33,1320):704
    (33,2640):2112
    (341,5456):492
    (341,10912):1968.

Every bank has order10. Thus among lengths1,...,11 only length10 can support a unit loop. There are5452 distinct relevant physical A/G seeds in total.

`minimal_walks.py` independently enumerates all four possible bit pairs at each state, solves the raw interface equations, rejects nonintegral or out-of-bound destinations, and exhausts every ten-step walk. It does not import the producer's parity-reduced edge routine. The strict finite-state square is complete for physical-to-physical paths by the previously audited forward/backward bound.

Results:10904 augmented seeds covering BOTH p states;216100 partial walks; **zero closed walks**, even before checking gain. The compatible-p half alone is108050 partial walks, matching the producer exactly. Normal and optimized receipts match. Consequently no candidate unit loop of length<12 exists in this precise physical core-order<=9 universe.

## Placement laws

PLACEMENT.md's fixed-context recurrence, modulo45 law, and lower bound pass both algebraic review and independent tests. For any fixed legal singleton prefix P of a bits and suffix R of b bits, the j=0 baseline has the form

    N_k=alpha·4096^k+(5870/91)k+gamma,

and hence

    N_(k+2)=4097N_(k+1)−4096N_k−264150.

This holds for each FIXED context. It is not a recurrence for choices of context that change with k. The coefficient of k remains unchanged because q,A,G are affine in4096^k and singleton composition gives T coefficient1.

The actual padded short12 interface has, modulo45, A=1,G=13 and h=q=S=E=0. Its bare baseline is9+5k modulo45. Applying the exact singleton-context law gives

    N_k ≡2^(b+1)−rev_b(R)−13P+7+5k (mod45).

Since gcd(2,45)=1, a window of at least45 consecutive odd values meets every residue modulo45. This is a residue statement, not scalar occupancy.

The lower bound has a direct injection. Let p be the first selected deletion vertex, and let M be the deletion graph with p removed. For every ordered covering pair of paths (U,V) of M, the vertex sequence U,v,r,p,V is a Hamilton path of the star. Empty U or V is permitted. All required arcs exist: U ends into v, v beats r, r beats selected p, and p beats every vertex of M. The image recovers U,V around the fixed vertices v,r,p, so it is injective. Singleton deletion removes a factor2 from the covering-pair count, giving

    H(star) >= q(M)=5760·4096^k·2^(a+b).

Thus a target N in this family can use only k<=floor(log_4096(N/5760)), when N>=5760. This is a restriction on this construction, not on all tournament representations of N. It neither proves nor refutes sufficient baseline placement after arbitrary legal contexts.

Independent `placement_audit.py` checks496 contexts with prefix/suffix lengths<=4, k<=6 and every position, totaling13888 scalar values. It also independently verifies the earlier h95 five-coordinate return modulo190 after three Y repetitions and residues149,169,119. All those h95 progression values are coprime to190, but no primality conclusion follows. Normal and optimized logs agree.

## Explicit linear diameter bound for fixed carrier parameters

The following statement is independently proved here, whether or not it is included in the producer's final text. Fix h,q and allow any bank of physical cores with those values. Give each core m singleton padding vertices. Within any common full-interface fiber,

    diameter(T) <= q+2h+m(q+2h−2).

The same bound holds for the associated star-count diameter. It is an O(m) statement with h,q fixed, not uniform over growing carrier parameters or nontrivial outer modules.

Proof: compare two core positions, cancel their common binary prefix/suffix using fixed-position decoding, and let d<=m be the intervening move length. Write the endpoint physical states as C_0,C_d. Every bookkeeping A_i,G_i is in (−h,q). Choose an initial bookkeeping bit p_0, put x_d=G_d mod2, and sum the edge cocycle W along the open path. A direct telescoping identity is

    W=A_0 G_Y−A_X G_d+h(T_Y−T_X)+hp_0 y_first−h x_last x_d.

For common exterior singleton contexts with last/first bits u,v (use coefficient0 when a side is absent), the actual T difference is

    T_0−T_d+W
    +(E_0−hp_0)y_first+x_last(hx_d−S_d)
    +u(S_0−hx_first)+v(hy_last−E_d).

Each of the four boundary corrections has absolute value<=h. Physical T counts only covering pairs with both paths nonempty, excluding the2h empty-sided covers. Thus0<=T<=q−2h and |T_0−T_d|<=q−2h. For x=y=1 the edge formula gives |w|<=q+2h−2; for exactly one of x,y equal to1 it gives |w|<=q+h−1; and for x=y=0 it is zero. Since h>=1, all three cases are bounded by q+2h−2. Summing gives the claimed bound. In particular it implies the producer's slightly weaker m(2q+h)+q+2h bound because q>=2h. Neither bound is conjectural. The same-position case follows immediately from its physical T and endpoint terms. No assumption that intermediate A/G states are physical is made.

`diameter_identity.py` independently verifies this open-path telescoping identity and its edge bound on all1200 prefixes of the old and new unit walks. Those tests are controls; the displayed telescoping proof gives the general result. Arbitrary non-singleton contexts instead multiply pure-T gaps by their Hamilton-count product, so that setting is not covered by the unscaled bound.

## Files

New independent artifacts: minimal_seeds.cpp, MINIMAL_SEEDS.txt/.stats, minimal_walks.py, MINIMAL_WALKS.json and normal/optimized logs; placement_audit.py, PLACEMENT_AUDIT.json and logs; diameter_identity.py and DIAMETER_IDENTITY.json; MASS1_TO7.txt. Earlier theorem/audit files are retained. No producer file was changed and no commit was made by this audit worker.
