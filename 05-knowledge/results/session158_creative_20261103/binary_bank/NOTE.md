# An exact binary state bank for strong transitive-deletion tournaments

This is an exact membership test, not a global tournament-spectrum classification. It is a reformulation and consequence of the one-apex path-splitting argument in S157; no literature novelty claim is made.

For m>=1 define A(1)=2 and recursively

    A(2m)   = A(m)+m,
    A(2m+1) = A(m)+1 if m is even, and A(m) if m is odd.

Define F(m)=m+A(m). Then, for every integer T>=3, the following are equivalent:

1. There is a strong tournament with T Hamiltonian paths that becomes transitive after deletion of one vertex.
2. F(m)=T for some integer m satisfying ceil((T-1)/2) <= m <= T-2.

The witness is explicit: write m in binary and append 0; let the added apex beat exactly the transitive base vertices corresponding to 1s in the resulting word. There are at most floor(log2(T-2))+3 vertices. This does NOT test membership in the full tournament spectrum. For example, the code certifies that 13, 29, 35, and 41 have no strong transitive-deletion witness, with no assertion about general tournaments.

Proof. Let w start with 1. Write P(w) for the nonnegative integer represented by w. Define the accumulator A(w) to be 2 plus (i) the sum of 2^(i-j-1) over every earlier 1 at j and later 0 at i and (ii) the number of adjacent occurrences 01. The value 2 is an advance credit for the two exterior placements, not the path count of every prefix: for the one-bit prefix 1, the actual two-vertex tournament has only one Hamiltonian path. If a final 0 is present, those two exterior placements are indeed valid, and S157's split-at-the-apex count proves A(w) is the actual path count. The start 1 and end 0 are also necessary and sufficient for strong connectivity, because otherwise the first or last transitive vertex is respectively a source or sink, and when both hold every vertex reaches the apex and is reachable from it.

When a 0 is appended to w, all earlier 1s contribute their binary place values, summing to P(w). When a 1 is appended, it contributes 1 exactly if the preceding bit was 0. Hence the stated recursion for A(P(w)), with A(1)=2. For the strong word w=u0 its count is F(P(u)).

By induction 2<=A(m)<=m+1: the lower bound is immediate; for even 2m, A(2m)<=2m+1; for odd 2m+1, A(2m+1)<=m+2<=2m+2. Consequently T=m+A(m) forces m<=T-2 and m>=(T-1)/2. This proves the exhaustive finite interval and witness bound.

The recurrence supplies a bank of O(T) exact integer states through T-2; each state reads only state floor(m/2). Search takes O(T) integer additions/comparisons and O(T) space, with O(T log T) bit complexity using ordinary arithmetic. It gives certificates for nonmembership in the restricted family by listing F(m) over the proven interval. Multiplicative closure via ordinal sums then supplies exact tournament counts from products of bank entries, but generally loses strong connectivity; it does not give additive or exact universal coverage.

Reproduce with `python3 check.py` or `python3 -O check.py`. The script compares the bank to independently counted Hamilton paths using subset/endpoint dynamic programming for every binary word beginning 1 and ending 0 with length 2 through 8; checks the universal recurrence and proven bounds for all m<=100000; and verifies explicit restricted-family absences 13,29,35,41. It reports additional exact counts via products but does not claim completeness for general tournaments.

## A positive-density obstruction for this restricted family

Let S(X) be the number of distinct counts at most X attained by strong transitive-deletion tournaments. Then

    S(X) <= (47/120) X + O(sqrt(X)).

In particular, at least (13/120) X - O(sqrt(X)) positive odd integers up to X are absent from this family. These missing values need not be absent from the full tournament spectrum.

Proof. Each witness word w=u0 corresponds to one m=P(u). By the exact criterion m<=X-2. Moreover, if m is even, then A(m)=A(m/2)+m/2>=m/2+2. Thus no even m>2X/3 can contribute to S(X). This rules out X/6+O(1) candidates. Also, for m=4r+1 with r>=1, A(m)=A(r)+r+1>=r+3, so no such m>4X/5 contributes. These odd candidates are disjoint from the previous excluded even candidates and rule out X/20+O(1) more. Thus there are at most (47/60)X+O(1) witness states.

Reverse every edge of the tournament, then reverse the order of its transitive base. The result corresponds to w*=complement(reverse(w)); this is an involution on strong witness words and preserves the number of Hamilton paths. It pairs all non-fixed witness states with distinct states representing the same count. A fixed word has even length 2k, its second half is forced by the first, and its first bit is 1; there are 2^(k-1) such words of that length. The witness bound limits word length to floor(log2(X-2))+2, so there are O(sqrt(X)) fixed words in total. Therefore S(X) is at most half the number of witness states plus half the fixed points, establishing the stated inequality. The separate parity theorem for tournament path counts (or the bank recurrence) gives the final odd-integer comparison.
