# The 30,030 boundary: an exact source and four forced chart changes

The earlier [prime-gap bridges](../session243_prime_gap_bridges_20261110/README.md) discuss 30,030=2·3·5·7·11·13 and its neighboring primes. In particular the symmetric prime offsets ±17 hide the closer 30,029 on the left: the crossing gap is 18, not 34. We use the same named integers as exact target probes, not as heuristic samples for density.

## 1. The raw 99-digit formula has a genuine geometric gap here

Even give the [S325 digit chart](../session325_absolute_cut_overlap_20261114/FULL_99_DIGIT.md) more freedom than its published phase bank: for every odd c=1,…,197, allow any phase r≥0 and any nonempty source U, with no padding requirement, and ask for

    N = 99·2^r Q(U)+cH(U).

The phase can be absorbed into adding r singleton SCCs, so it suffices to test the relaxed r=0 formula with any source. For N≤30047 this forces Q≤303. A complete profile catalogue at this budget consists of 34 nonempty pairs, plus (1,1) for the empty source. Its completeness is small for a structural reason: in the [Grinberg–Stanley odd-cycle packing formula](https://arxiv.org/abs/2307.05569), the empty packing contributes 2^n to Q and each directed triangle contributes 2^(n−1). A strong tournament of order seven has at least three triangles, since every vertex of a strong tournament lies in a directed triangle. For the latter fact, if no arc went from the out-neighbors of a vertex to its in-neighbors, its out-neighbors could not reach it. Three triangles give Q≥320. At order ≥8 one triangle suffices to give Q≥384. Thus strong letters at this budget have order at most six; singleton letters are allowed.

The replay enumerates all unlabeled tournaments through order six by all one-vertex extensions and degree-block canonical labeling, retains the strong letters, then multiplies pairs up to budget. This includes all tournaments at the stated Q cutoff, without assumptions on their total order, converse, or cuts.

| N | Feasible (H,Q,c) in the relaxed 99 chart |
|---:|---|
| 30013 | none |
| 30017 | (23,296,31) |
| 30029 | none |
| 30031 | none |
| 30047 | none |

The middle entry is also physical in the actual bank: digit 31 has phase zero, and the forward active word 11000 gives (H,Q)=(23,296). The literal identity is 99·296+31·23=30017. No arbitrary-target claim follows from the complete residue alphabet. In particular the three historical primes are not hidden in a different shift of these 99 digits.

## 2. Change the cut to match each prescribed target

For a literal left tournament I, let an apex beat a cut S in I and lose to the entire transitive middle block TT_m, with I→TT_m. There is no right block (or use one terminal singleton that is removed). With h=H(I), α and λ the star's two left-cut counts as in S325, the exact output is

    H = α·2^m + h − λ.

This is the original star formula with the right block empty; it also follows directly by splitting a Hamiltonian path at the apex. A candidate cut is found backwards: N−h+λ must be α times a positive power of two. Here are literal outgoing masks and cut masks, with every phase and interface checked directly:

| Target N | I outgoing masks | S | (h,α,λ) | m |
|---:|---|---:|---|---:|
| 30013 | (32,17,3,7,13,30,63) | 112 | (37,470,104) | 6 |
| 30017 | (96,17,3,7,77,30,46) | 80 | (115,468,50) | 6 |
| 30029 | (6,112,98,7,13,25,57) | 101 | (151,947,426) | 5 |
| 30031 | (2,36,17,7,11,29,63) | 92 | (29,472,206) | 6 |
| 30047 | (64,9,35,21,71,27,46) | 61 | (113,953,562) | 5 |

At middle length one, each template has a unique nonsingleton SCC containing the middle mark. Inflating that mark preserves the sole strong SCC; removing the other singleton SCCs keeps the count. The replay also evaluates each final core literally. These are five strong witnesses, including exact repairs at all four targets the relaxed 99 chart cannot reach. They show why choosing a cut from the actual target is substantive, not a different residue encoding of the same chart.

## 3. The numbers 99,100,198,200 are not interchangeable slot budgets

Rédei oddness makes H odd. Q is even, and for any tournament of order n≥2 it is divisible by four: pairing each subset with its complement in Q=Σ H(S)H(S^c) gives 2 times an even number of odd terms. A singleton has Q=2. So Q=99,198,30030 are impossible for any tournament (the last two are 2 mod4).

By contrast Q=100 and Q=200 are both physical, and the complete small-budget profile atlas proves each forces H=11. An order-five witness is (2,20,9,3,13) at (11,100); adding one singleton gives (11,200). In this case the appealing round Q budgets cannot even enter the coprime-to-99 residue decoder, because their forced H=11. The same clock doubling changes Q and not H. This is a concrete reason to treat 99→198 and 100→200 as arithmetic retimings but not as interchangeable moves in the physical profile lattice.

The same packing formula supplies a second exact obstruction: Q≡(−1)^n H mod 3 at order n. Each allowed odd cycle removes an even number from the power of the slot variable, hence all terms have the same parity of exponent; evaluate at 2≡−1 mod 3. In particular 3 divides one coordinate iff it divides the other. The round budget Q=30000 **is** physically realizable, for example on the strong order-nine masks (462,480,378,114,419,385,176,12,200), with H=1827. But no source at Q=30000 can have H coprime to 99, regardless of its size or arrangement. The obstructions at 30000 and 30030 therefore differ: at 30000 the slot budget exists but forces a nonunit H, while at 30030 the slot budget itself cannot exist.
