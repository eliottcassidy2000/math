# Every eight-label pair quotient is visible

**FINITE-EXACT QUOTIENT THEOREM, WITH ARBITRARY AMBIENT LIFT.** Let Q be any tournament on eight fixed quotient labels. Inflate each label to a fixed physical pair, each with either internal orientation. Comparing Q with its whole quotient converse while retaining every pair interior is always visible to the full named two-port response. For each fixed Q, all2^9 quotient/pair-converse states have distinct responses, also in every fixed ambient context.

This is not a theorem for every even quotient order, reconstruction of arbitrary sixteen-name tournaments, or even reconstruction between unrelated eight-label quotients. The order-eight completeness step below is an exhaustive exact certificate with a proved coverage argument. The ambient lift and reversal-cube consequence are all-order arguments.

## 1. Two sufficient tests with their inherited inputs

Let B_i={2i,2i+1}, i=0,...,7. Set T=Q[B_i] and U=Q-op[B_i], on the same names and preserving all child interiors. Each coarse bit a_i or c_i is expanded onto both physical names; x->y is fixed. Write Delta=F_T-F_U.

The established pair-inflation identity from `../even_pair_refinement/PROOF.md` gives

    [a_i a_j c_j]Delta = 4h_ij (mod 8),       i!=j,
    h_ij=sum_(S subset V-j, |S| even)(S_i(Q[S])+E_i(Q[S])) (mod 2).

The short-support version is inherited, not a new identity: equation(2) in the pinned [E1/P2 parity packet](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) specializes at order eight to

    h_ij=sum_(S subset V-j, |S|=2 or4) S_i(Q[S]) (mod 2). (1)

Thus a nonzero h entry yields an eight-actual-query finite-difference separator. The full physical response is queried; other cubics and quartics are not set to zero or discarded.

The second sufficient test is ordinary quotient endpoint imbalance. Put

    b_i=S_i(Q)+E_i(Q) (mod 2).

Then for every coarse mask c,

    Delta(0,c)=2 sum_i b_i c_i (mod 4).                   (2)

Here the quotient may have any positive order; the finite theorem will use eight.

To justify(2), a=0 makes x a universal source, so the response reduces to H(T+y). The classical odd-cycle expansion recalled in `../odd_pair_quotient/PROOF.md`, Section2, gives H=1+2(number of directed odd cycles) modulo four. Cycles not containing y have the same parity for T and U, since their scalar Hamilton counts agree exactly by whole converse and module run reversal. In a cycle through y, if any pair contributes one physical name, or two nonconsecutive physical names, exchanging those names is a free involution. Such terms vanish modulo two. The remaining cycles use whole consecutive pairs, and collapse to a quotient path P from s to t with contact c_t(1-c_s). Reversing P gives the opposite-quotient term, whose difference is c_t-c_s. Hence the half-contrast is the endpoint parity of all quotient simple paths. Inherited P2, coefficientwise E_i(S)=sum_(A subset S) S_i(A), identifies the all-simple-path start/end parity with full Hamilton end/start parity. This proves(2). A one-pair collapsed cycle is impossible under uniform attachment; singleton path terms cancel in the endpoint sum.

Consequently b_i=1 supplies the single actual query a=0,c=e_i. These tests are sufficient only; a zero test by itself is not a response collision.

## 2. Complete seven-name precomputation

For a seven-name tournament R, let A be its zero-diagonal0/1 adjacency matrix, regarded over F2, and d=A1. Define

    t(R)=d+A^3 1+diag(A^3) over F2.                      (3)

Equation(1) shows that column j of h(Q) is exactly t(Q-j), with its remaining physical labels retained. Indeed two-name paths starting at i contribute d_i. A three-edge directed walk can repeat a name only by returning to its start: loops and two-cycles are forbidden. Subtracting these closed three-walks leaves precisely the paths on four distinct names. This gives the last two terms in(3).

The checker encodes every seven-name tournament by its21 lower-to-higher edge bits, in lexicographic pair order. It evaluates(3) for all2^21 codes, retaining the complete seven-bit vector in one byte per code. Exactly50,688 vectors are zero. No isomorphism reduction or probabilistic sampling is used in this precomputation.

## 3. Exact extension coverage at eight names

Every eight-name tournament is uniquely specified by its deletion Q-0 and the seven bits for the directions0->j, j=1,...,7. Its28-bit code is therefore

    code(Q)=128 code(Q-0)+incidence_mask.

If all h_ij vanish, t(Q-0) must vanish. It is therefore necessary and sufficient to examine every one of the50,688 zero seven-name codes with all128 incidence masks, then test the seven remaining deletions Q-j. This checks6,488,064 candidate extensions and covers every possible h-zero tournament among all2^28 labelled quotients. A failed deletion is an exact nonzero h column, not a heuristic rejection.

The resulting complete labelled h-zero list has30,720 entries, of which26,880 are prime. A separate subset path-DP implementation, not the matrix formula, recomputes every shared-cubic column on all30,720 retained quotients and confirms it is zero. It also computes both full Hamilton endpoint parity vectors and proves, case by case over this complete list,

    S_i(Q)+E_i(Q)=1 for every i=0,...,7.                 (4)

Thus every quotient either has a nonzero h entry and is separated by the cubic test, or belongs to this complete residual and is separated by(2), at every anchor. This is the finite-exact order-eight theorem.

The retained list has exactly four relabelling orbits. Their minimum edge codes are600626,600653,602967,1472106, with orbit sizes1920,13440,13440,1920 and automorphism orders21,3,3,21. The middle two are prime; the outer two have a universal source or sink. Whole converse pairs the first with the last and the two middle orbits. These orbit claims are additional exact controls, not a substitute for the labelled coverage argument. The orbit checker constructs every relabelled code explicitly and compares full integer sets.

## 4. A physical prime residual control

The first prime residual is

    Q=(154,80,75,34,12,23,57,126).

Its shared-cubic matrix is zero, its earlier anchor-correction bits are all one, and its full start/end parity masks are155 and100. Their XOR is255, as(4) predicts. Its sixteen-name pair inflation has the actual physical query

    a=0, c={0,1}, x=16,y=17,

with source value5 and target value7 modulo eight. A literal18-name subset Hamilton DP checks both values independently of the quotient matrix test. This is a separator for the full observer, not a response collision or a relabelled comparison.

## 5. Fixed ambient lift and the two-level cube

Give both ports the core's original exterior incidences. The augmented core has18 names, hence even order. The companion `../mod8_transfer/PROOF.md` multiplies its Hamilton-count contrast by the same odd unit modulo eight in every fixed rooted ambient context. Every nonzero local contrast above remains nonzero. In particular the endpoint test's nonzero residue modulo four cannot disappear.

All coarse queries are exactly invariant under each pair-interior converse, by reversing its maximal path runs. They identify the root quotient bit even when all eight child bits are unknown. Once that bit is known, one isolated physical query per pair reads its interior orientation independently of the others, as in `../six_pair_reconstruction/PROOF.md`, Section4.

For a known fixed quotient Q, choose one nonzero shared-cubic coefficient when h is nonzero; its eight-query finite difference is a fixed root battery. If h is zero, one endpoint query suffices instead. Adding eight pair decoders gives at most16 adaptive queries. Including decoders for both quotient orientations gives at most24 fixed queries. These are nonoptimal upper bounds for the fixed-Q two-level family, not a universal battery for identifying unrelated quotients.

Hence all512 root/pair-converse states associated with each fixed eight-label quotient have distinct complete named responses in every fixed ambient context. The exterior is unchanged throughout. No claim for arbitrary even quotient orders follows from the finite census.

## 6. Reproduction

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The standard-library driver compiles three C++17 sources in a temporary directory, installs nothing, binds all source/data hashes and checks explicit failure exits. The complete seven-name vector table, complete labelled eight-name residual list, endpoint histogram, four orbit records and literal physical controls are retained losslessly. DATA.md gives byte and label conventions.
