# Place a one-arc +2 by changing the reader rather than the source

The old +2 grammar holds one source switch fixed and searches for a reader with a prescribed signed invariant. Here the inserted source is **identical on both sides**. The one changed edge is a reader contact, selected backwards from the desired shift. Its two coefficients change with the height; the construction does not use the old switching vector (10,-8).

## Physical reader and exact cut inversion

Fix any nonempty tournament B with H=h and two-cover count P₂=p. A transitive spine of length l uniformly beats B, B beats a new apex, and the apex beats just those spine vertices whose chronological binary word has a 1. Write u for that word's binary value, and start it with 1. All resulting tournaments are strong: the apex reaches the first spine vertex, then the whole spine and B; every vertex reaches the apex either directly or through B. Use

    U(v) = sum_(i>k) b_i b_k [2^(i-k-1)+[i=k+1]],

where the bits b_i are indexed from least significant. As in the earlier two-run reader, its Hamilton count is

    Phi(u;B)=u p+(1-U(u))h.

Indeed every visit to B must either end the path or return only via the apex; the two-run coefficient is the weighted start at marked spine vertices (the binary address u), and U counts the ordered covers whose left end and right start are both marked. This is a literal slot construction and places no new arcs inside B.

Let the switched bit have index j, write q=2^j, and put u=2q A+r, 0≤r<q, so bit j is initially zero. Write ρ=rev_j(r) (exactly j bits). Inserting a 1 in column j creates precisely the U-pairs with a mark at that column, giving

    U(u+q)−U(u)=A+[A odd]+ρ+[ρ odd].

For the count before the edit, the cross pairs spanning the blank column give

    U(u)=U(A)+U(r)+2Aρ.

Both formulas follow from the powers associated to the separation of the marked columns. In the first formula the two adjacent possible endpoints supply exactly the two parity corrections, and in the second there is no adjacent cross pair.

Now prescribe an even shift Δ. If M=(q p−Δ)/h is a positive even integer, choose **any even** ρ with 0≤ρ<q and A=M−ρ>0, and let r be its j-bit reversal. Then the one actual apex arc from column j changes the count by exactly Δ, and its starting count is

    N=2 h A²+2 Δ A+r p+h(1−U(A)−U(r)).                 (1)

The bit lengths and literal vertex labels do not change upon the edit, since A>0. Odd ρ and odd A can also be handled, but are not needed here. This solution lets a target prescribe the signed movement first while keeping a full free suffix in the absolute source; it is still not a promise to place an arbitrary N.

Already at j=0, with B transitive of length l≥2 (h=1,p=2^l) and Δ=2, taking A=p−2 gives the explicit all-height pair

    N=2p²−(9/2)p+3, N+2.

Here U(p−2)=p/2−2; the first pair is the prime pair 17,19. For the fixed nontransitive block with outgoing masks (2,4,9,3), (h,p)=(5,32), j=0 gives the actual pair 379,381, and every j≡0 mod 4 is eligible. A source can itself grow independently: take a copies in ordinal order of this four-vertex block and L−5a ordinal singletons, giving h=5^a,p=2^L when L≥5a. Then the precise shift-two criterion is L+j≡1 mod 4·5^(a−1) for a≥1. This follows from ord_(5^a)(2)=4·5^(a−1), obtained by iterating the elementary lifting of 2^4−1=15. The resulting reader still makes the final tournaments strong; the inner ordinal source need not be strong.

## Absolute inverse, density exponent, and residue tests for fixed B

Fix h,p and Δ=2. Eligible j are those with q p≡2 mod h and q p>2; evenness of M is automatic because h is odd and p even. They recur periodically or do not exist. Also p≥2h from the two empty/full slot splits. For all eligible j≥2, A>0 for every allowed even ρ. Thus each such height contains q/2 pairs.

For prescribed N in (1), there are O(1) candidate A at each height: 0≤U(v)≤v−1 for v>0 (by the cover interpretation: U is a subset of the v starting covers and excludes the full path with empty first part). Since q p=h(A+ρ)+2, p≥2h gives q≤A+2/h. Bound (1) by

    |N−2h A²| ≤ 8h A+7h.

For example rp<pq≤2hA+4, hq≤hA+2, and the remaining U terms are bounded in magnitude by h(1+A+q). Consequently a constant-width integer neighborhood of sqrt(N/(2h)) contains A, after which ρ=M−A and r are forced. Test these exact formulas and the physical parity conditions. There are O(log N) heights: from the original reader formula N≥u(p−h)≥q p. The inverse at a fixed source is polynomial in the input bit length; a varying physical source still has to be chosen.

These q/2 pairs are distinct at a fixed eligible height j≥2. Increase ρ by 2, so A decreases by 2. The leading drop in 2h A²+4A equals 8h(A−1)+8. The unfavorable r term costs at most p(q−1); the unfavorable U terms together cost at most h(A−1+q−2) for q≥4 (U(0)=0). Since h A≥q p−h q−2, the resulting lower bound is

    6pq−8hq+p−5h−6 >0,

using p≥2h and q≥4. Moreover N is between fixed positive multiples of q² p²/h, uniformly for those heights: A≥qp/(4h), u≥2qA and N≥u(p−h), while A≤qp/h and N≤u p+h. For a fixed B with infinitely many eligible heights, the progression of eligible j has bounded gaps. Summing the geometric number of suffixes over all heights then bounds the number of distinct starts ≤X above by O(sqrt(X)); any one suitable height gives the matching Omega(sqrt(X)). Thus this is a genuinely branching sqrt-scale family, not a universal spectrum.

There is full **mod-h steering** on the same physical +2 pairs once an eligible j is large enough. Observe that h and p are coprime from q p≡2 mod h, hence 3 does not divide h: the packing congruence P₂≡(-1)^|B| H mod 3 would otherwise force 3 to divide p. Reduce (1) modulo h:

    N ≡4M−4ρ+pr.

The coefficient c_k of the kth low-to-high bit of r is p·2^k−4·2^(j−1−k). Two adjacent types satisfy c_(k+1)−2c_k ≡3q·2^(-k) modulo h, which is a unit. As coefficients repeat when k increases by ord_h(2), a sufficiently long suffix has at least h−1 independently selectable bits of each of these two types, avoiding the leading zero that ensures ρ even. Taking any number of them gives all residues modulo h while preserving eligibility and the exact +2. The same proof accommodates any prescribed bounded fringe whose fixed bits contribute a translation. This is not a replacement for the absolute square-band test.

Replay `python replay.py` checks the literal paired graphs and their actual Hamilton counts, full reversal formulas, hundreds of exact inversions at four independently selected source profiles, and all 2048 monotonic and modular positions of one sizeable native band. The calculations use no external software.
