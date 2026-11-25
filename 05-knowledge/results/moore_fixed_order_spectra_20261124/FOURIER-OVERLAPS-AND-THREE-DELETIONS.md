# Fourier overlaps shrink the order-11 Moore chart; deletion overlaps hide a +2 flip

This extends the [540-site voltage chart](ORDER11-ORBIT-PACKETS.md) and the earlier [three-deletion Möbius defect](../../../07-reflections/fkn-at-the-transitive-corner.md). The graph claim is conditional on the degree-57 Moore graph admitting an order-11 action. The Chan–Pak examples are unconditional and do not by themselves give an absolute selector.

## 1. Complete all missing pairs from 138–140 orbit rows

The action fixes a pentagon. Label its 295 moving orbits (25 ports in known groups of five, 270 remote), and use phases modulo 11. In the nontrivial Fourier channel ζ let $ B_{ij}=\sum_{d\in D_{ij}}\zeta^d $ be the adjacency, where D is the voltage set between two orbit labels. The fixed sites have no coordinates in this channel. It is Hermitian with eigenvalues 7 and −8. Thus the projection (P=(7I-B)/15) has rank (r=(2065+t)/15), with (t=c/5) and $c\in\{25,100,175\}$ from the [cyclic-orbit constraint](ORDER11-ORBIT-PACKETS.md). The trace computation there gives r=138, 139, 140.

Take a principal basis S of P, i.e. P_SS is invertible of order r. Record only the **off-diagonal** shift sets in pairs touching S. For i in S the projection identity gives
$$
 P_{ii}(1-P_{ii})=\sum_{j\ne i}|P_{ij}|^2.
$$
Its permitted diagonal values are ((7-a)/15) with a in ${0}\cup\{\zeta^k+\zeta^{-k}:1\le k\le5\}$. Any two different candidates could have the same product only if their a-values sum to −1. That is impossible: the only relation up to degree ten is a multiple of $1+z+\cdots+z^{10}$, whereas (1+a+a') has at most five nonzero coefficients. So the known off-diagonal voltages recover P_SS exactly, including its diagonal. Now complete all rows by the Gram identity
$$
                 P=P_{:,S}(P_{S,S})^{-1}P_{S,:}.
$$
Indeed the Gram vectors indexed by S span the range. From P recover B and each shift set: two binary shift polynomials with the same value at ζ have difference constant on all eleven coefficients. The only distinct binary subsets this could identify are empty and full, and the full set is impossible (all blocks have at most three shifts). In particular a recovered zero Fourier entry off the diagonal means an empty shift set. Fixed-to-moving incidences are determined by the port labels, and the fixed pentagon is prescribed.

One can split S into three panels of orbit rows. A pair payload lies in one or two panels; its exact union size is A+B+C−D−E−F+G (here G is zero for pair payload). The three panels need not know any outside-outside pair. Once their overlaps agree, the projection glues the whole voltage array. This is structurally different from simply writing down all three vertex-deletion charts of a tournament, whose labelled orientations recover every arc trivially. The replay recovers **all nine moving orbit rows from four** in the actual order-5 Hoffman–Singleton control, using all primitive fifth-root evaluations over F_11; that control keeps its diagonals, since the separate diagonal shortcut above uses order 11.

## 2. Certified contact budget: 298 sites with phase labels, 295 modulo phase origins

Write (m_{ij}=|D_{ij}|) for i≠j and $M=\sum_{i<j}m_{ij}=8395-c$. Use a random principal basis of the Hermitian projection, with probability det(P_SS). By Cauchy–Binet this is a probability law and
$$
\Pr$i\in S$=p_i,\qquad
\Pr$i,j\in S$=p_ip_j-|P_{ij}|^2,\qquad p_i=P_{ii}.
$$
For example, factor P=UU* with U*U=I and expand det(U*diag(z)U); its coefficients are squared maximal minors, and its first and second derivatives at z=1 give these formulas. See [Lyons, Determinantal probability measures](https://pmihes.centre-mersenne.org/articles/10.1007/s10240-003-0016-0/).

Let $I(S)=\sum_{i<j\in S}m_{ij}$ be the double-counted shifts and $L(S)=\sum_{i\in S}(57-\sum_{j\ne i}m_{ij})$ the degree not needing storage (one fixed neighbor per port, two internal neighbors per cyclic orbit). The touching panels use exactly (57r-I(S)-L(S)) off-diagonal shifts. For the expectation:

1. Lift p to a vector x constant on each moving orbit, zero on the fixed pentagon. The full adjacency obeys $A+8I\succeq(65/3250)J$, with $\sum x=11r$ and $\|x\|^2=11Q$. Here $Q=\sum p_i^2=((295-c)49+268t)/225$. For the five cyclic step values $a_s=\zeta^s+\zeta^{-s}$, we used $\sum a_s=-1,\sum a_s^2=9$. Thus the sum of the cyclic (p_i^2) is (268t/225). Remove the contribution of the diagonal internal graph edges to obtain
$$
\sum_{i<j}m_{ij}p_ip_j \ge \frac{11r^2}{100}-4Q-\frac{268t}{225}.
$$
2. $\sum m_{ij}|P_{ij}|^2\le \sum m_{ij}^3/225\le(M+4(1350-2c))/225$. Indeed $m^3\le m+4m(m-1)$ for 0≤m≤3, and the remote row test gives ten ordered repeated-shift units per independent orbit and six per cyclic orbit. Port pairs are simple. Count each unordered pair once.
3. Ports have (p_i=7/15); the cyclic p-values sum to (36t/15). Hence $\mathbf E L=(175+72t)/15$.

Add 1 and 3 and subtract 2. Some basis has I+L at least the ceiling of this expectation bound. Use a 2-bit branch index, a 295-bit basis mask, and the bounded-cardinality combinatorial rank of the occupied shift positions. There are (W=11r(589-r)/2) positions. All counts in the table compare exact integers:

| Cyclic orbits c | Rank r | Guaranteed I+L | Max stored shifts R | W | Bits with phase labels | Bits without phase-origin labels |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 138 | 1,806 | 6,060 | 342,309 | 44,223 | 43,206 |
| 100 | 139 | 1,887 | 6,036 | 344,025 | 44,128 | 43,111 |
| 175 | 140 | 1,969 | 6,011 | 345,730 | 44,026 | 43,009 |

For the last column, the support of the touching data is connected. To see this, P_SS is block diagonal across components of the touching support. For i outside S, P_iS has support in at most its own component; if that support is empty then P_ii=0 from Gram reconstruction, impossible. Thus the completed P and its voltage support could not connect distinct components. But the moving graph itself is connected: all remote vertices join by a path of length ≤2 avoiding the fixed pentagon, because no remote vertex has a fixed neighbor, and every port has remote neighbors. The quotient is connected as well.

Changing the origin of each moving orbit acts by translating its incident shift sets. A nonempty proper subset of F_11 has no nonzero translation stabilizer. Thus the only gauge transformations preserving the full touching support *and* its shift assignment are constant origin shifts; every data class has exactly (11^{294}) distinct assignments. Normalize shifts on a fixed spanning tree of the touching support (choose the least translate of each tree-edge shift set). This is equivalent to dividing the number of admissible raw records by (11^{294}). To make explicit ranking possible, first order types by the numbers of singles, doubles and triples, then rank the support and the allocation of multiplicities on the support. For each type with s≥294 occupied pairs, rank shift subsets with radix ${11\choose k}$ for nontree edges and ${11\choose k}/11$ on the 294 tree edges. The latter product exactly saves (11^{294}), independently of the multiplicities on the tree. Reserve slots even for disconnected supports of those types (they never encode a graph). The same capacity bound $\lfloor\sum_{j=0}^R{W\choose j}/11^{294}\rfloor$ still applies to these reserved indices, since for each support the product and the division factor are unchanged.

A marked tournament on 298 sites has 44,253 binary contacts, enough even for the phase-labelled version (this keeps input phase origins). One on **295** sites has 43,365 contacts, enough for the normalized version, which reconstructs the graph up to those independent origin choices. The c=175 stratum alone even fits 294 sites (43,071 contacts). Use a fixed lexicographic orientation with an arc reversed exactly at each bit equal to one. This is a reversible code, **not** an induced record of 295 original graph vertices, nor an unmarked tournament decoder. Its Hamilton count is odd and positive, at most (295!<2^{2003}); no exact degree-57 count is claimed.

The older [540 at order five](../../../07-reflections/tiling-hamiltonian-ratios-s95.md) is $\sum F(C)^2$, the ordered same-class pairs in a fixed-path tiling slice. The new compression does not use 540 as a parameter: what transfers is the lesson to count overlaps with their actual multiplicities. The Plücker weighting is also quadratic, but it measures retained Gram overlap, not isomorphism collisions. The second older 540, from a dyadic inverse, has a more direct connection below.

## 3. What the seven deletion counts miss in Chan–Pak

Let $U=\{a,b,c\}$ be three marked vertices, (T-S) an ordinary induced deletion, and $H(\varnothing)=1$. The odd-cycle collection identity of [Grinberg–Stanley](https://arxiv.org/abs/2307.05569) counts H as a sum over vertex-disjoint collections of nontrivial directed odd cycles, with a collection of q cycles having weight (2^q). Marked vertices not in a cycle of the collection are unconstrained. Inclusion–exclusion gives the *exact* recursion
$$
 H(T)=A+B+C-D-E-F+G+R_U(T),
$$
where A/B/C delete one mark, D/E/F delete two, G deletes all three, and (R_U) is the total weight of collections whose union contains all three marks. A modulo-four shadow only sees single odd cycles containing all three; the full weighted residual sees the rest. This makes the older deletion-lattice defect an explicit object for every choice of three marks.

There is a concrete one-arc seam **113 → 115** on eight vertices for which all seven deletion counts agree. Use mask 113248304 (bits for i<j in lexicographic order, bit 1 means i→j), flip the pair {1,3}, and mark {0,4,7}, disjoint from the edited pair. The single-deletion counts are **15,65,37**, the double-deletion counts **5,5,37**, and the triple-deletion count **5**, on both sides. Thus the seven-term sum is 75 on both sides and the residual is 38 versus 40. The missing information is not just an issue of labeling or deleting the edited arc. It can sit entirely in a three-way interaction. This is a statement about seven **scalar counts**, not matching every smaller induced tournament.

There is also a growing **two-mark** family where all three proper deletion readings are neutral and a single arc changes the exact count by two. Use mask 75588141, flip {2,5}, mark {3,4}, and replace vertex 1 by a transitive chain of length $m\ge1$. In this seed vertex 1 has either indegree one or outdegree one, and neither the edit nor taking marked deletions can increase that smaller degree. Any Hamilton word then has at most two nonempty chain runs: at most one run uses the unique outside entrance (or exit), and one more may start (or end) the word. For a fixed skeleton on the other vertices, one run has one filling; two have (2^m-2) fillings, by assigning each labelled chain vertex to a run, with both used and sorted. Thus (H_m=u+v(2^m-2)) for every face, determined by the two exact counts at m=1,2. The full readings are
$$
             54\,2^{m-1}-3\quad\longleftrightarrow\quad54\,2^{m-1}-1.
$$
The three proper face readings agree for all m because they agree at m=1,2. This is not a three-mark infinite-family claim.

The earlier [540-to-546 multiport inverse](../../../session204/multiport-inverse/REVIEW.md) proved that its fixed bank misses $(195\,2^n-9,336\,2^n+19)$ at every n≥1. Set m=n+3 in the new guarded construction: $216\,2^n-3$ and $216\,2^n-1$ are strictly inside **every one of those gaps**. This adds two particular absolute values per gap, not coverage of the interval or of all large odd integers.

Run `python3 fourier-overlap-check.py` for the exact bit bounds, rank/unrank unit checks, the actual Hoffman–Singleton Fourier control and the finite certificates for the deletion families.
