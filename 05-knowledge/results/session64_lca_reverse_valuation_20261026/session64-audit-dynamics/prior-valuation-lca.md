# A labelled p-adic collision carrier: trie, LCA and adjacent-LCP RMQ

Status: PROVED elementary statements with FINITE-EXACT independent controls. No novelty claim for tries, longest-common-prefix queries, range-minimum queries, or p-adic valuations. This packet provides an exact typed connection to the Session55 fixed-entry modular templates and Session53 guarded repetition cylinder. It does not construct an integer magic square of distinct squares, prove Collatz convergence, or supply universal coverage.

## 1. Inheritance and concept board

Closest proved mechanism: Session55's fixed-port CRT collision-mask theorem (`inputs/session55-crt.md`, sections 2–4). Canonical hostile: the crossed three-block partitions A and B on the same nine labels. Corrected near miss: a smooth root-vector lift tree is not the collision trie of nine entries, and neither is an integer-square certificate (`inputs/session55-lifts.md`, sections 2–4). Least-used sidecar: Session53's exact distance to the rational fixed point, recovered from the original guarded-pumping note and first-descent packet (`inputs/guarded-pumping.md`, `inputs/session53-base.md`), not guessed from a numerical analogy.

Board: named pair identity / residue-prefix depth / RMQ order / CRT partition intersection / root-height and native-path guards. Anchor is the exact carrier; niche is fixed-entry modular thresholds; wildcard is the guarded Collatz rational-path query. None transfers a global conclusion merely because its local carrier is exact.

## 2. Exact single-prime carrier

Let p be prime and x_0,...,x_(N-1) be a finite labelled list of integers, possibly negative or repeated. Set v_p(0)=infinity. For k>=0, let P_k be the partition of labels by equality of x_i modulo p^k. All labels share the root P_0, and P_(k+1) refines P_k.

**Theorem 1 (trie/LCA).** Write each integer in its infinite least-significant-first base-p expansion, with digits d_j(x)=floor(x/p^j) mod p in {0,...,p-1}. Then

    depth(LCA(i,j)) = LCP(d(x_i),d(x_j)) = v_p(x_i-x_j),
    i and j share a block of P_k iff depth(LCA(i,j)) >= k.

For equal integers the paths coincide forever, so the LCA-depth convention is infinity. Identical records must either be kept as several labels on one terminal or explicitly retain this identity relation; they must not be split artificially at a finite leaf. Diagonal queries i=j also return infinity.

Proof. Their first k digits agree exactly when their residues modulo p^k agree, exactly when p^k divides x_i-x_j. Taking the largest k proves the assertion for unequal records; equality gives every k. Euclidean division supplies the same expansion for negative integers, so no sign restriction is needed. The hierarchy is the partition family itself; digit values label branches but are unnecessary if only pairwise collision queries are wanted.

For N>=2 pairwise distinct integers the first k>=0 for which all residues are distinct is

    k_sep = 1 + max_(i<j) v_p(x_i-x_j).

For N=0 or 1 take k_sep=0. If distinct labels have equal integers, k_sep=infinity. Thus the commonly quoted maximum formula needs both the distinctness and nonempty-pair qualifications. This is a statement about these fixed records, not all modular configurations.

A finite sufficient precision is any H with p^H>max(x_i)-min(x_i) after duplicate values are collapsed. At that precision unequal values cannot coincide, since their nonzero difference has absolute value less than p^H. Every finite pair depth is <H. Empty or one-value sets permit H=0; duplicate identity is retained separately.

## 3. Exact ordered RMQ representation and implementation

Let y_0,...,y_(m-1) be the distinct values sorted lexicographically by their length-H least-significant-first digit strings. Retain the map from original labels to this order. Define adjacent depths a_t=v_p(y_(t+1)-y_t).

**Theorem 2 (adjacent LCP).** For 0<=r<s<m,

    v_p(y_s-y_r) = min_(r<=t<s) a_t.

Proof. The strings with a specified prefix form a contiguous interval in lexicographic order. The two endpoint strings share a length-k prefix iff every string in the interval shares it, iff every consecutive pair shares it. Thus their LCP is at least k exactly when every adjacent LCP is at least k. Both depths are finite and <H, proving equality. Consequently max pair depth equals max adjacent depth, so the all-distinct threshold is obtained from adjacent entries alone.

This is an exact compression to O(m) adjacent scalar depths plus an indispensable label-order sidecar. RMQ answers pair queries without expanding an N-by-N matrix. The scalar depths without their order and labels are not a lossless carrier. A finite trie with edge lengths can alternatively be recovered by grouping maximal consecutive runs above each attained depth; equal-depth branching is simultaneous, not an arbitrary binary refinement claiming new p-adic levels.

`verify.py` implements digit sorting, identity collapse, adjacent valuations and a sparse-table RMQ. With unit-cost integer arithmetic, digit-key construction takes O(NH), comparison sorting at most O(m H log m) digit work, adjacent valuation work O(mH), sparse-table preprocessing O(m log m) time/storage, and each named-pair depth query O(1) thereafter. Original labels require O(N) space; explicit keys need O(mH) temporary space. Big-integer costs and prime-validation cost are additional, not claimed constant. Listing all pair bits still costs O(N²), as it must for an explicit pair-bit output. This is a reusable carrier, not a claimed speedup over direct arithmetic for N=9.

Hostile: natural integer order 0,1,2 at p=2 has adjacent valuations 0,0, but the endpoints have valuation 1. Correct digit order is 0,2,1. Ordinary numeric sorting therefore invalidates the RMQ theorem.

## 4. CRT is a family of trees, not one common refinement chain

For fixed records and n=product p^(k_p), define C_p(k_p) as the set of named pairs with depth_p>=k_p. CRT gives exactly

    C_n = intersection_p C_p(k_p).

The result uses all named-pair relations, equivalently a labelled partition family or one labelled hierarchy per prime. It does not require literal storage of 36 coordinates for a nine-cell square; the RMQ carrier is one equivalent compressed representation. However, independent local maxima, partition shapes or leaf counts are insufficient.

The actual Session55 modulus-35 witness is

    0  9  1
    16 15 14
    29 21 30.

At p=5, level 1 is A={{0,4,8},{2,3,7},{1,5,6}}. At p=7 it is B={{2,4,6},{0,5,7},{1,3,8}}. Their masks are 0x61089888 and 0x8a054250 with zero intersection. The clusters {0,4,8} and {2,4,6} intersect in {4} while neither contains the other. Every pair of clusters in one rooted hierarchy is nested or disjoint. Therefore no single rooted hierarchy on the SAME nine leaves contains both partitions as cluster families. This is a precise no-go, not a claim that these data cannot be stored in a more general graph, product structure, or separate trees.

Each partition has shape (3,3,3), but using A at both primes leaves nine colliding pairs instead of zero. The missing coordinate is label alignment. A common refinement of A and B is discrete, but retaining only that refinement loses the local partitions and does not repair the carrier.

## 5. Fixed-entry Lo Shu towers and the root-tree distinction

From `inputs/session55-lifts.md`, take L=(8,1,6;3,5,7;4,9,2), A=3L-2 and B=8L-7. These are fixed integer records; their entries are not all ordinary integer squares.

Since L contains exactly 1,...,9, its positive differences range through 1,...,8. Thus

    max v_3(A_i-A_j) = 1 + max_(1<=d<=8) v_3(d) = 2,
    max v_2(B_i-B_j) = 3 + max_(1<=d<=8) v_2(d) = 6.

The first all-distinct residue layers are exactly 3 and 7, respectively. All A entries are 1 mod3 and all B entries are 1 mod8, so the inherited square-root arguments provide square residues at every corresponding prime-power level. This carrier recovers their fixed-template separation thresholds without enumerating all modular magic squares. The stronger assertion that no other template works below these thresholds uses Session55's independent opposite-pair and descent proof, not this trie theorem.

The p³-branching smooth lift tree of Session55 has whole parameter/root-vector solutions as vertices. The present collision trie has individual entry labels as leaves and residues as prefixes. They are different objects. Once nine entries are distinct at a base level, every compatible lift remains distinct, but the collision trie does not encode the number of root-vector lifts, their root choices, integer height, carries, or global integer squareness.

## 6. A marked rational path recovers exactly the inherited binary fuel

Use T(n)=n/2 if even and (3n+1)/2 if odd. The original Session53 first-descent packet, replayed directly here, gives the 59-bit word of 27, 37 odd steps, endpoint 23, and

    q=2^59, a=3^37, b=1100931843921811423,
    D=q-a=126176846412426125,
    F(n)=(a n+b)/q on the guarded 27 word class.

The positive rational fixed point z=b/D has odd denominator and is nonintegral, hence defines a marked infinite binary path in Z_2. Since D is odd,

    LCP(n,z)=v_2(n-z)=v_2(Dn-b).

The equality is in the 2-adic field; its right side is an ordinary integer valuation. For a finite query one can use z_H=b D^(-1) mod2^H: the visible common-prefix depth is min(H,v_2(Dn-b)). A depth equal to H is only a lower bound; it must NOT be called exact without more digits. If H>v_2(Dn-b), the finite LCA depth is exact.

The parity-cylinder theorem, not the trie alone, identifies the cylinder of k repeated blocks as Dn=b mod2^(59k). Consequently the inherited exact repetition budget becomes the marked-path query

    K(n)=floor(LCP(n,z)/59)=floor(v_2(Dn-b)/59).

For positive n this is finite because z is not an integer. Under a legal block,

    DF(n)-b=(a/q)(Dn-b), so K(F(n))=K(n)-1.

This is exact fuel for repetitions of this single word only. For 27, K=1 and F(27)=23 has K=0. It supplies neither a rank law after exhaustion nor a stopping-time bound, universal entry, or Collatz convergence. The map preserves the repeated-prefix predicate; it loses every other future itinerary unless the original dynamical guard and source are retained. Arbitrarily deep nonempty finite cylinders are not one integer with infinite fuel.

## 7. Closed finite controls and reproducibility

Run `python3 -B verify.py > normal.json`, `python3 -O -B verify.py > optimized.json`, then compare byte for byte. No load-bearing check uses Python assert.

- Exhaustive ordered records of lengths 0..4 over {-3,-1,0,1,3}, at primes 2,3,5,7, including negative entries, identities, empties and singletons.
- 100 fixed-seed length-25 records with values in [-10^12,10^12) at primes in {2,3,5,7,11}.
- Direct digit-prefix, pair valuation and residue-partition comparison at all depths through H+1: 3,224 generic record cases; 107,582 ordered pair checks including the two Lo Shu cases.
- 100 CRT exponent checks on the actual nine-cell witness, four prime pairs and each exponent 0..4.
- Independently derive 27's first-descent parity word and affine constants by native iteration; 60 marked-path and full repeated-parity replays, depths 1..12 and five explicit lifts per depth. Verify exhaustion, not just a lower bound on the word length.
- Hostiles: natural numeric ordering, duplicate labels, crossed-prime hierarchy, and the failure of 27 to support a second whole block.

The mathematical proofs carry the unbounded quantifiers. The finite controls are sanity checks, not extrapolation. Source snapshots are pinned in SOURCE_PINS.json; SHA256SUMS covers the closed packet.
