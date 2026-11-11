# Dyadic target control at the minimum-triangle boundary

Status: PROVED elementary synthesis conditional only on the explicitly inherited Moon interval-chain product law, with a full proof of every new assertion below. FINITE-EXACT computational controls are separately specified. Chan–Pak remains OPEN: the exact conjectural spectrum is all positive odd integers except 7 and 21. These two exclusions are known universal holes, not merely holes of this constructor; see `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. No literature-priority claim.

## 1. Inheritance, object and scope

The source is `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, at pinned `origin/main` commit `1ca9c5f83b55bcac7dacc0a6f5def38e02f245e7`. Its all-order Moon product law and elementary recognition theorem, not its finite census, are used. For parts p_i>=2, let IC(p) have vertices 0,...,N, N=sum p_i, all arcs increasing except one reversed endpoint arc in each successive interval of length p_i. Consecutive intervals share an endpoint. The inherited identities are

    H(IC(p)) = product_i (1+2^(p_i-1)),
    sorted scores = (1,1,2,...,n-3,n-2,n-2), n=N+1,
    cyclic triangles = n-2.

These tournaments are strong and rigid, and their compositions are intrinsically recoverable for n>=4. Strongness also follows immediately because each interval is strong and successive intervals meet. The n-2 triangle count is the minimum among ALL strong tournaments on n vertices (Moon, cited through this repository source), not the minimum among arbitrary tournaments; transitive tournaments have zero. No new proof of Moon's minimum or product theorem is asserted.

Closest mechanism: the multiplicative Moon chain, distinct from the Session268 two-parameter atom/product compiler. Hostile: these chain counts never have a factor 7. Corrected near miss: a modular target is not an integer target. Least-used coordinate: the exact 2-adic valuation of the difference between equal-size chain gadgets.

Concept board: physical strong chains; Hensel digit lifting; fixed-score fibres; partition entropy; local-versus-integral spectrum. The result joins all five without a global enumeration.

## 2. Equal-cost digit switches

For j>=1 choose two compositions A_j,B_j of a common total size c_j:

    j=1: A_1=(5), B_1=(2,3), c_1=5;
    j=2: A_2=(7), B_2=(3,4), c_2=7;
    j>=3: A_j=(2j), B_j=(j,j), c_j=2j.

Write a_j=H(IC(A_j)), b_j=H(IC(B_j)). Then

    (a_1,b_1)=(17,15), (a_2,b_2)=(65,45),
    a_j=1+2^(2j-1), b_j=(1+2^(j-1))^2 for j>=3.

Thus a_j and b_j are odd and

    v_2(b_j-a_j)=j.

For j=1,2 the differences are -2,-20. For j>=3 the difference is
2^j-2^(2j-2)=2^j(1-2^(j-2)), with odd parenthesis. Also both a_j,b_j equal 1 modulo 2^j. These exact valuations, rather than the absolute magnitudes of the counts, make the switches triangular.

## 3. Fixed-order, fixed-score dyadic controllability

**Theorem 1.** For every k>=3, there are exactly 2^(k-1) explicitly specified strong rigid tournaments, all of order

    n_k=k^2-k+7,

all with the same score sequence and the minimum possible n_k-2 cyclic triangles among strong n_k-vertex tournaments, whose Hamilton-path counts represent each odd residue modulo 2^k exactly once. A prescribed odd residue is decoded using k-1 exact integer switch decisions.

**Construction.** For epsilon=(epsilon_1,...,epsilon_(k-1)), concatenate A_j when epsilon_j=0 and B_j when epsilon_j=1; call the resulting chain T_epsilon. Its count is

    h(epsilon)=product_j a_j^(1-epsilon_j)b_j^epsilon_j.

Its order is 1+5+7+sum_(j=3)^(k-1)2j=k^2-k+7, including k=3 (the sum is empty). The inherited chain facts give strongness, rigidity, scores and triangles.

**Proof of residue bijection.** If epsilon and eta first differ at index j, their factors before j agree. Their j-th factor difference has valuation j. All subsequent factors in either word are 1 modulo 2^(j+1). All preceding factors are odd. Consequently

    v_2(h(epsilon)-h(eta))=j.

In particular the residues are distinct modulo 2^k. There are 2^(k-1) words and equally many odd residues, proving the assertion. Distinct residues also prove that the unmarked tournaments are pairwise nonisomorphic; no merely labelled count is substituted for an isomorphism count. There are exactly that many tournaments in the specified family, not necessarily in the entire fixed-score fibre. QED.

**Constructive inverse.** Start with h=product a_j. At stage j, the earlier stages ensure h=r modulo 2^j. If h differs from target r modulo 2^(j+1), replace h by h*b_j/a_j and select B_j; otherwise retain A_j. The current h is divisible by a_j because this factor has not yet been changed. The replacement changes exactly the next live bit, because

    h*b_j/a_j-h=(h/a_j)(b_j-a_j)

has valuation j. Later switches cannot disturb this bit. This proves the algorithm, with no tournament search and no oracle for H. Its k-1 decisions are arithmetic-operation counts, not unit-cost bit-complexity claims; operands have O(k^2) bits. Producing the adjacency matrix has its ordinary O(n_k^2) output cost.

The theorem is intentionally stated k>=3. Precision 1 is realized by C3, and precision 2 by the strong counts 3 and 5, at possibly different orders. The main fixed-order formula need not be extended to those boundaries.

## 4. The quadratic cost is the correct exponent in this architecture

**Theorem 2.** Suppose Moon interval chains on a common order n=N+1 have H values representing every odd residue modulo 2^k. Then

    N >= (3(log 2)^2/(2*pi^2)) (k-1)^2.

Here log is natural. Theorem 1 therefore has asymptotically optimal quadratic order within the entire Moon interval-chain architecture, up to an absolute multiplicative constant. It is not an optimality theorem for arbitrary tournaments.

**Proof.** The count of an interval chain depends only on the unordered multiset of its parts, an integer partition of N with parts at least 2. Thus the number of different counts is at most p(N), the number of unrestricted integer partitions. For any t>0, the formal partition product, evaluated at exp(-t), gives

    p(N) <= exp(tN) product_(i>=1)(1-exp(-ti))^(-1).

Absolute convergence and nonnegative terms justify the logarithmic expansion. Its logarithm, excluding tN, equals

    sum_(m>=1) 1/[m(exp(tm)-1)] <= (1/t) sum_(m>=1) 1/m^2
    = pi^2/(6t),

using exp(tm)-1>=tm. Set t=pi/sqrt(6N). Hence p(N)<=exp(pi sqrt(2N/3)). Since at least 2^(k-1) distinct counts are required, taking logarithms proves the bound. QED.

Moreover every chain with total part size N satisfies

    2^(N/2) < H < 2^N,

because 2^(p-1)<1+2^(p-1)<2^p and the number of parts is at most N/2. Thus quadratic vertex order also means quadratic bit-length in this family. Our residue representatives can be far larger than their targets. The construction does not hide a logarithmic-height exact-target algorithm.

## 5. An actual 2-adic isometry, with an integral obstruction

**Theorem 3.** For an infinite binary word epsilon, the scalar product

    Phi(epsilon)=product_(j>=1) (a_j if epsilon_j=0 else b_j)

converges in the 2-adic integers. Phi maps binary words bijectively onto the odd 2-adic integers, and if two words first differ at j, then

    v_2(Phi(epsilon)-Phi(eta))=j.

Fair independent bits push forward to the uniform (Haar) probability measure on the odd 2-adic integers.

**Proof.** The j-th factor is 1 modulo 2^j, so the products are Cauchy. The same first-different-bit proof used in Theorem 1 applies after passage to the limit. For an arbitrary odd 2-adic target, solve successive residues by the constructive inverse. The solutions are compatible: all factors with index >=k equal 1 modulo 2^k, so an additional factor cannot change the previous chosen bits. This gives a unique infinite word. Every odd residue cylinder of modulus 2^k has exactly one binary prefix of length k-1 as inverse image and therefore probability 2^(-(k-1)), proving the measure statement. QED.

The infinite object here is a LIMIT OF COUNTS, not a finite tournament with a 2-adic number of Hamilton paths. The source is a physical chain word, the map is its count modulo increasing powers of 2, the preserved predicate is the first differing digit, and the lost information is the actual integer target and full graph. Precision k and the full composition supply the finite-construction sidecar.

**Integral hostile, proved.** Every Moon chain count is coprime to 7. Indeed 2^a modulo 7 lies in {1,2,4}; none is -1. No factor 1+2^a is divisible by 7. In particular our construction approximates both 7 and 21 to arbitrarily high 2-adic precision but never realizes them, or ANY multiple of 7, exactly. The stronger architecture-specific hostile is 63: it IS attained by general tournaments at order 8 (inherited THM-1370, Corollaries; also `01-canon/theorems/THM-344-n8-h63-complete-omega.md`), but no Moon chain realizes it. Thus this constructor misses a genuine Chan–Pak target, not only the universal holes.

For example, at k=8 the inverse for target 63 gives the composition

    (2,3,7,6,4,4,5,5,12,14),
    n=63, H=12644053665875775 = 63 (mod 256).

The equality n=63 in this example is incidental; its Hamilton count is the displayed much larger integer. This is a particularly concrete local-to-global warning: a sparse family can have full 2-adic closure. It does not imply a new obstruction for arbitrary tournaments, which can have counts divisible by 7.

## 6. Reproduction and actual computational scope

Run `python3 reproduce.py > results.json`. All checks use explicit `require`/`raise`, not optimization-disabled assertions. `python3 -O reproduce.py > results-optimized.json` was also run, and `cmp results.json results-optimized.json` confirms byte-identical output. The standard-library program checks every odd target for k=3,...,15, exact switch valuations, literal subset-DP counts for both banks through j=7, 40 seeded small chain samples, and all four full-family k=3 witnesses. It independently checks strongness, scores and triangle counts on these literal graphs. Eighty target inversions at k=32,64,128,256 check exact arithmetic at larger precision; those are NOT raw graph DP checks. No enumeration of all tournaments, independent rigidity search, or audit of inherited Moon results is claimed.

For an actual adjacency witness use, for example,

    python3 reproduce.py --target 7 --bits 8 --graph

Replacing `--target 7` by `--target 63` gives the attained-spectrum hostile above. `results.json` includes target-63 examples at precisions 6,8,12,32,64.

It emits the composition and exact outneighbor-bit rows of a strong tournament whose H is 7 modulo 256, not H=7. The output explicitly records the exact integer H so that this distinction can be verified.
