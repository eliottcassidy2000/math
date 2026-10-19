# Exact chain thresholds and bounded rank-window matching

19 October 2026. **PROVED general formulas and obstructions; FINITE-EXACT sharp repairs for chain lengths 2 through 7.**

This continues [the support-three Hall barrier](../gold_support_three_hall_20261019/GOLD_SUPPORT_THREE_HALL_BARRIER.md). A support bound counts moved labels in one exchange, equivalently unequal positions in the two full extension words. A matching or fractional matching must obey capacity one at every full extension. No historical-novelty claim is made.

The results concern specified first and second comparisons. They do not obstruct the existence of a different Gold strategy. In particular, the examples below contain proper non-total modules and have elementary alternative balanced comparisons. They are method barriers, not candidates against Gold or H7.

## 1. What the next support scale really costs

For every m>=2 there are strict-Gold comparison trees on posets of incomparability degree at most max(3,m) for which **every support-at-most-(m+1) matching or fractional matching falls short**. This is an interior obstruction: the first heavy probability is strictly less than 2/3, and every realized Gold path has a positive margin.

For m=7 the construction is seven-thin. Its smallest parameter choice used here has 32 vertices and

    E = 11260392,
    b = 7505784,
    b/E = 2187/3281 < 2/3,
    required matching = 3751176.

Every support-at-most-eight matching has capacity at most3750705, a shortfall of471. The actual heavy Gold margin is1716. A chain-only support-nine perfect matching, verified on just3432 chain words and lifted symbolically, proves that **nine is the exact minimum support bound for certifying this specified two-comparison tree**. No eleven-million-word enumeration is used.

More generally, for m=2,...,7 the exact threshold for the stated family is m+2, using six small chain certificates. For arbitrary m, the proved bounds are

    m+2 <= minimum sufficient support <= 2m.

The upper bound swaps the two full equal chains. The apparent sharp value m+2 beyond m=7 is not claimed.

## 2. Complete support-three capacity on two chains

Let C_m and C_n be disjoint chains with labels

    a1<...<am,     b1<...<bn,

and no relations between the chains. Let N=m+n. Extensions are binary shuffle words, with the labels of each chain appearing in order. Compare a1 with b1.

**Theorem 1.** The complete support-at-most-three exchange graph has maximum integral and fractional matching capacity exactly

    nu_3 = binomial(N−2,m−1).                                  (1)

This is already attained by support-two exchanges. Triple exchanges add no capacity.

**Proof for m,n>=2.** Apply the two-chain Hall cap from the companion note with (x,X,y,Y)=(a1,a2,b1,b2). No label is incomparable with all four: every remaining a-label is comparable with a1 and a2, and every remaining b-label with b1 and b2. Thus only the discordant cells contribute to the cap.

The event a1<b1 but b2<a2 is exactly the prefix a1 b1 b2. The opposite discordance is exactly b1 a1 a2. Their total number is

    binomial(N−3,m−1)+binomial(N−3,m−2)
      = binomial(N−2,m−1).

For the matching lower bound, pair a1 b1 w with b1 a1 w for every legal suffix w. There are exactly binomial(N−2,m−1) such pairs. Their endpoints are distinct. The same vertex-cover upper bound applies fractionally. If one chain has length one, the smaller orientation side has one word and the same initial swap gives a matching of size one, equal to (1). ∎

For equal chains, write C=binomial(2m,m). The full extension total is C, each orientation has C/2 words, and

    nu_3/C = m/[2(2m−1)].                                      (2)

Suppose an arbitrary first-comparison block with g extensions and h extensions in its heavy child is placed before the two-chain block by ordinal sum. Every legal support-at-most-three exchange reversing the chain-minima comparison fixes the first block: changing its word needs at least two moved labels, while reversing the chain orientation needs at least two more. Therefore the exact full-child capacity is h nu_3.

The second comparison is exactly balanced, so the chosen tree is actually Gold whenever p=h/g<=2/3. But the support-three method certifies it exactly when

    p <= (4m−2)/(7m−4).                                        (3)

For m=2,...,7 these thresholds are respectively3/5,10/17,7/12,18/31,11/19,26/45; they tend to4/7. This is an exact capacity statement for this class, not just a one-sided estimate.

A particularly small robust example takes the first block to have relations c<d,c<e with f isolated. It has g=8, and the first comparison f<d has h=5. Thus p=5/8 and the actual heavy margin is C/2, but support three fails for every m>=2. The resulting poset has degree at most max(3,m). For m=2 this gives an eight-vertex, three-thin strict-Gold tree whose support-three capacity is10 while12 is required.

## 3. Individual reachability stops one support level too early

**Lemma 2 (exact run cost).** In two disjoint chains, a word beginning with precisely r a-labels before b1 has minimum Hamming distance r+1 to any legal word with b1<a1.

**Proof.** The target must begin with b1. Each of the first r positions in the source contains ai at position i, whereas in the target ai cannot appear before position i+1. Those r positions and the old position of b1 must all change. Conversely, rotate the prefix

    a1 ... ar b1  ->  b1 a1 ... ar

and keep the suffix fixed. Both chain orders remain legal, and exactly r+1 positions change. ∎

For two equal length-m chains every word therefore has an opposite-orientation neighbor of support at most m+1. This does **not** supply a perfect matching.

Let

    U = a1 ... am b1 ... bm,
    V = a1 ... a_(m−1) b1 am b2 ... bm,
    W = b1 a1 ... am b2 ... bm.

**Lemma 3 (two-to-one Hall bottleneck).** For m>=2, the only opposite-orientation neighbor of either U or V under support at most m+1 is W.

**Proof.** In any b1-first target of U, the first m positions and the old b1 position must all change. These consume the entire m+1 budget, so the remaining suffix is fixed. The only legal target is W.

For V, its first m−1 a-label positions and its b1 position must change, consuming m moves. If am stays fixed at position m+1, the preceding labels force W. If am moves too, the entire budget has been consumed in the first m+1 positions, whose label set is unchanged; b1-first chain order again forces W, contradicting am having moved. ∎

Thus a matching has capacity at most C/2−1, even though there are no isolated vertices. This is the next collective obstruction after the individual run bound. It proves that perfect pairing requires support at least m+2.

The obstruction also survives an ordinal prefix. For every prefix extension L, consider (L,U) and (L,V). Their only support-at-most-(m+1) neighbor is (L,W): the required chain changes consume at least m moves, leaving at most one move in the prefix, and a permutation cannot change exactly one prefix label. Taking all h prefix words gives a Hall set of size2h with a neighborhood of sizeh. Even exchanges that are allowed to alter both blocks cannot evade this cap:

    nu_(m+1) <= h(C/2−1).                                     (4)

## 4. A degree-three root approaches two thirds from below

For t>=2, define K_t on two chains A1<...<At and B1<...<Bt, adding

    Ai < B_(i+2),    Bi < A_(i+2)       for 1<=i<=t−2.

Its extensions are exactly the shuffle words whose every prefix has A-count minus B-count between−2 and2. Add A1<B1 and At<Bt to obtain T_t. Thus extensions of T_t begin with A and end with B. Every vertex has at most three incomparable neighbors; adding the two extra relations cannot increase degrees.

**Lemma 4.** Put h=3^(t−2). Then

    e(T_t) = g = (3h+1)/2,
    e(T_t + B1<A2) = h.                                       (5)

For t>=3, the latter is the heavy branch. For t>=4 its probability lies strictly between3/5 and2/3.

**Proof.** After the first A and before the last B, the walk is at height+1. At odd times its possible heights are+1 and−1. Over each two-step segment the transition-count matrix is

    M = [[2,1],[1,2]].

For example, from+1 the words AB and BA return to+1 and BB reaches−1; AA is forbidden. The eigenvectors (1,1),(1,−1) have eigenvalues3 and1, so (M^(t−1))_(+,+)=(3^(t−1)+1)/2, proving the first count.

Because the word begins with A, the event B1<A2 is exactly that its second letter is B. The first two-step segment after the initial A then has choices BA and BB. Their continuation counts sum to a column sum of M^(t−2), namely3^(t−2). This proves (5). ∎

## 5. Strict Gold needs support m+2 in an entire seven-thin family

Take

    P_(t,m) = T_t ordinal-sum (C_m disjoint-union C_m).

Ask B1<A2 in the first block, and compare the minima of the two final chains second. Then

    E = g C,    b = h C,    each heavy second cell = h C/2.

The heavy Gold margin is

    E−3b/2 = C/2 > 0.                                       (6)

On the smaller first branch the same second comparison is also exactly balanced, and its margin is positive. Thus this is a fixed, strict-Gold two-comparison tree. The first heavy probability tends to2/3 from below as t grows; it never reaches the easy equality boundary.

The required matching capacity is

    2b−E = (h−1)C/2.

By Lemma3, every support-at-most-(m+1) matching or fractional matching has capacity at most h(C/2−1). Its shortfall is at least

    h−C/2.                                                   (7)

Consequently **h>C/2 is a sufficient strict obstruction**. The degree is at most max(3,m), independent of t. For every m>=3 this proves that a universal bounded-support theorem for every prescribed strict-Gold tree in degree-m posets would need support at least m+2.

For 2<=m<=7 the companion JSON contains an explicit support-at-most-(m+2) perfect matching on the equal-chain block. Keep every T_t word fixed and apply this matching to the final block. This gives hC/2 disjoint exchanges for every t, establishing a positive support-(m+2) certificate with heavy margin C/2. The same matching lifts through arbitrary ordinal padding, with all counts multiplied by the padding extension totals.

Combining this construction with (7), the exact support threshold is m+2 whenever h>C/2, for each m=2,...,7. The t choices4,5,6,7,8,9 respectively are enough. The corresponding base orders have12,16,20,24,28,32 vertices. For m=7, every t>=9 gives a seven-thin strict-Gold tree whose exact support threshold is9. This is an infinite family, both by increasing t and by ordinal padding.

The finite certificate has2352 matching edges in total across all six chain lengths;1716 belong to m=7. Verification needs only their distinct endpoints, opposite orientations, legal chain orders, and Hamming distances. It neither enumerates the full P_(t,m) extension set nor assumes a universal formula beyond m=7.

## 6. A positive bounded-window theorem for arbitrary thin posets

The preceding counterexamples distinguish local reversal from matching capacity. There is still a general positive construction that controls both support and overlap.

Let R be any finite poset of order n, and let x,y be incomparable. Write D(v),U(v),I(v) for strict predecessors, strict successors, and incomparable neighbors. Every extension position of v lies between |D(v)|+1 and n−|U(v)|. Define the deterministic rank interval

    l = min(|D(x)|,|D(y)|)+1,
    r = n−min(|U(x)|,|U(y)|),
    w = r−l+1.

Both x and y belong to this window in every extension. Moreover

    w <= |I(x) union I(y)|
      = d(x)+d(y)−|I(x) intersection I(y)|.                   (8)

To see this, the vertices comparable to both x and y are exactly their common strict predecessors and common strict successors: an element between them would make x,y comparable. Hence their remaining complement is I(x) union I(y). Also min|D|>=|D(x) intersection D(y)|, and similarly for U. Substituting proves (8). In a k-thin poset, w<=2k.

**Theorem 5 (disjoint rank-window matching).** Freeze all positions outside [l,r]. This partitions L(R) into disjoint fibers f. Let A_f and B_f be the counts of the two x/y orientations in one fiber. There is a support-at-most-w matching of size exactly

    nu_window = sum_f min(A_f,B_f)
              = (e(R)−sum_f |A_f−B_f|)/2                   (9)

inside these fibers. This is the optimum among matchings restricted to these fixed-outside fibers. In particular, if every nonzero fiber imbalance has the same sign, this matching saturates the globally smaller orientation side.

**Proof.** A nonempty fiber fixes its set S of w inside labels. Every linear extension of R[S] can be substituted at those positions: outside labels remain either before the whole window or after the whole window, so replacing the internal order cannot violate an outside relation. Conversely every word in the fiber restricts to an extension of R[S]. Thus the fiber is exactly L(R[S]). Any two of its words differ on at most w positions, so opposite orientations form a complete bipartite graph. Pair its smaller side arbitrarily with distinct members of the larger side. Different fibers have disjoint endpoints, so their capacities add without overlap. Summing proves (9). A common sign makes the sum of absolute imbalances equal the absolute global imbalance. ∎

For a first heavy branch of size b, the construction gives Gold whenever

    sum_f |A_f−B_f| <= 2E−3b.                                (10)

This is a sufficient condition, not a claim that every balanced global comparison satisfies it. Cancellation between differently biased fibers can make the global imbalance smaller than the sum in (10).

A useful finite representation avoids listing outside words. Let I and J be nested order ideals of R with |I|=l−1 and |J|=r. Put S=J\I and O=V\J. The multiplicity of the fiber type S for this ideal pair is e(R[I])e(R[O]). Thus (9) equals

    sum_(I subset J) e(R[I]) e(R[O])
      * min(e(R[S]+x<y),e(R[S]+y<x)).                         (11)

Every term is nonnegative, the decomposition is disjoint, and every inside poset has at most2k vertices in a k-thin ambient poset. Formula (11) does not remove the possible number of outside ideal pairs; it bounds the local exchange support, rather than claiming a polynomial-time algorithm.

For a minimal illustration of the lost sign, take C2 disjoint-union C2 and compare its minima. The canonical window is positions1 through3. Fixing the fourth label gives two fibers of size3, with opposite2-to-1 biases. Formula (9) gives capacity2, although the global sides both have size3. Support four repairs this. Thus common-sign or quantitative window control is a real hypothesis.

## 7. The pointwise interval bound and the remaining gap

There is also a simpler pointwise version. In any extension with x before y, every label strictly between them is incomparable with x or y; otherwise x<z<y would force x<y. The entire closed interval is therefore contained in I(x) union I(y). Reorder only that interval by any linear extension of its induced poset with y<x adjoined. Such an extension exists because x,y were incomparable. All outside relations survive the block replacement. Hence every extension has an opposite-orientation neighbor moving at most |I(x) union I(y)| labels, at most2k in a k-thin poset.

Neither this pointwise fact nor the rank-window construction proves that support2k always has enough matching capacity. The proved positive result is (9)–(11), including its common-sign saturation case. The family above forces support at least k+2 for selected strict-Gold trees when k>=3. The universal gap between such lower bounds and sufficient bounded-support capacities remains open.

## 8. Reproduction and status

Run `python verify_chain_support_and_windows.py` and its `python -O` variant beside `chain_support_sharp_certificates.json`. It verifies all six chain matchings, the degree-three root counts by independent ideal recurrence for a stated finite range, the m=7 consequence arithmetic, the run/Hall obstruction on small exact controls, the exact support-three formula on a bounded set of chain lengths, and the rank-window partition/ideal-sum identity on selected posets.

The proofs establish the all-m and all-t claims. The receipts establish the six claimed finite sharp repairs. The small-width source-count formula is exact symbolic arithmetic, not a census to a cutoff. No conjecture about seven-thin Gold or H7 is promoted.
