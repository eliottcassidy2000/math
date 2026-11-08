# Independent audit of the order-five mobile module

**Verdict: PASS within the stated scope.** No completeness, density, prime-coverage, or unbounded-order census conclusion follows.

The portable `verify.py` imports only the Python standard library. It does not import producer code or inherited transfer code. Both ordinary and optimized Python runs execute explicit exception-based checks, not removable `assert` statements. `inputs/` freezes the audited producer report, code, and six banks; receipt hashes identify these inputs and the independent checker. The frozen producer scripts are provenance only: they are not the canonical verifier and their use of assertions, SymPy, or scratch-path imports is not covered by the optimized-verification claim. That claim applies only to the independent standard-library verifier.

## Literal checks and finite scope

The checker independently derives local transfers by full permutations split into two ordered paths. A separate subset dynamic program counts Hamilton paths, endpoint counts, and first-endpoint path-pair weights in literal adjacency-mask tournaments.

For masks [18,8,3,5,14], H=13, Q=104 and alpha=[18,14,18,15,26]. Exactly the subsets {4} and {0,1,2,3} have intermediate weights 13m with 0<m<7, giving m=2 and m=5. Thus they model binary 010 and 101 respectively. This is matching normalized two-bin sidecar data, not equality of full rooted states.

All r=4 untrimmed witnesses are checked literally: D=401,409,457, E=0, U=1118,d=13. Their gcd is 1. In fact every untrimmed fiber for r≥3 has gcd 1: its j=0 and j=1 scalars are odd and differ by 8. All r=7 trimmed witnesses are checked literally, including graph reachability from every vertex: D=23731,23779,23971,24739, E=d=13,U=35503; their gcd is 1. Full adjacency masks are retained in the receipts.

The independent checker exhausts every labeled tournament for n=1,...,6, derives all compatible (H,Q/H) counts, and compares complete sets of statistics signatures against the frozen banks. It independently validates every bank representative's H,Q,alpha and every qualifying root subset. The bank grouping is not isomorphism classification. Counts match the producer report. Consequently H=11 is excluded for binary-compatible modules only through order six; no stronger exclusion is established.

## Exact all-order algebra certificate

Use six states (0,0),(0,1),(1,0),(1,1),(2,0),(2,1), with the first coordinate recording empty/last-in/last-out for P and the second empty/good-first for R. Put

W = span{(0,1,0,-1,0,0), (0,0,0,1,0,-1)}.

This is exactly the subspace of vectors supported on indices 1,3,5 whose coordinates sum to zero (indices are zero-based).

The checker verifies S0 W⊆W and S1 W⊆W, each of the four triangle defects and both order-five defects annihilates W, and multiplying any such defect on the left by S0 or S1 maps its image into W. Therefore for every nonempty binary word w and every pair of these defects X,Y, X S_w Y=0. This is an inductive invariant-space proof for arbitrary word length, not a long finite spacer test. Empty spacers generally fail and must not be included.

In fact the certificate is universal. Categorize a block's local covers by (p,r), where p,r∈{0,1,2} record empty/good/bad endpoints. A nonempty block never has category (0,0). For two blocks or weighted blocks with equal H,Q,alpha, the eight category differences necessarily have the form

(0,1):x; (0,2):−x; (1,0):y; (2,0):−y;
(1,1):a; (1,2):b; (2,1):−x−a; (2,2):x−b.

The first two pairs follow from equality of H. Equality of alpha gives the seventh entry, and equality of Q gives the eighth. Applying the six-state update rule gives a transfer difference linear in x,y,a,b. The independent checker constructs its four coefficient matrices and verifies each annihilates W and enters W after either singleton. Thus every physical module satisfying Q/H=2^k and alpha/H=m has the same separation property relative to H times the corresponding k-bit word, provided 0≤m<2^k. The frozen `inputs/UNIVERSAL.md` states this general theorem. The checker also verifies its displayed four-parameter matrix against all four coefficient matrices independently reconstructed from category updates. Thus the displayed family and its invariant-space certificate are both bound to the receipt.

## Scalar formula and arithmetic boundary

For A=S1S0 the checker verifies A^4−6A^3+9A²−4A=0. With the 010 defect G, it verifies all sixteen initial values of ell A^k S1 G A^j v, j,k=0,...,3, equal 4^(j+1)−10+4·1_(j=0). Both sides obey that recurrence independently in j,k, so induction proves equality for every j,k≥0. Combining the baseline binary-word scalar gives the producer's full D formula for r=j+k+2. Distinctness follows directly from the increasing powers of four, with the j=0 boundary correction treated separately.

The trimmed strong construction is valid for 1≤j≤r−3. The deleted initial bit-0 singleton is a source and final bit-1 singleton a sink, so deleting them preserves D. The remaining first singleton is beaten by the root and remaining last singleton beats the root. All blocks lie on directed routes between these boundaries; this proves strongness, with the strong module itself independently checkable. The sidecar formulas follow by binary compatibility. The restriction on j ensures both boundary singletons really remain singletons.

For modules with counts H_i and at least one ordinary singleton between successive modules, expansion has only its baseline and single-defect terms:

D=(∏H_i)D_base + Σ_i (∏_(h≠i)H_h)δ_i.

The universal invariant-space certificate kills every term with two or more defects: consecutive chosen defects in an expansion have a nonempty intervening binary word, even if other module positions between them contribute their baseline words. The identity therefore holds for the entire transfer matrix before applying terminal vectors, not just D.

The general sidecar law is also exact. If block i has counts H_i,Q_i and permitted-first weight alpha_i, partition a root-started second path by its first nonempty block. This gives

U/d = 1 + Σ_i (alpha_i/H_i) ∏_(j>i)(Q_j/H_j), with d=∏H_i.

Matching each module's Q/H and alpha/H to its binary subword therefore preserves the original binary value plus one. This proof includes empty first paths and the entirely empty post-root path through the initial 1.

For arbitrary positive integer counts H_i, the greatest common divisor of the coefficients ∏_(h≠i)H_h is exactly (∏H_i)/lcm(H_i): at each prime its valuation is the sum of the valuations minus their maximum. The baseline coefficient ∏H_i is a multiple of that number. This proves the product/lcm divisor at all orders, with no bounded numerical inference. For t H3 modules and u H13 modules, t,u≥1, every D is therefore divisible by 3^(t−1)13^(u−1). More generally the common forced divisor is (∏H_i)/lcm(H_i). This is a lower bound on forced divisibility, not necessarily the exact gcd of a particular fiber. Repeated finite module types therefore retain an arithmetic obstruction. A single H13 module gives a gcd-one unbounded family of physical fibers but does not establish prime coverage or exponentially many prime-capable scalars.

## Reproduction

Run `python verify.py > CHECK.json` and `python -O verify.py > CHECK.optimized.json` from any directory, keeping `inputs/` alongside the checker. The census covers 33,867 labeled tournaments in total. Normal and optimized receipts should agree byte for byte.
