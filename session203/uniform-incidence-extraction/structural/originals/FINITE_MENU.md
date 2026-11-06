# A finite target-dependent uniform-incidence menu cannot extract an affine count

Status: **PROVED CANDIDATE / INDEPENDENT AUDIT PENDING.** Joint argument with the author of `session203-complement-gadget/PROOF.md`. It uses that companion's physical strict pair-capacity lemma, also under audit. The eventual odd spectrum remains **OPEN**.

## 1. Constructors and quantitative statement

A uniform binary-incidence constructor consists of an arbitrary fixed core tournament A and, for each a in A, two bits e_a(0),e_a(1). Input b=(b_1,...,b_m) supplies a transitive chain x_1->...->x_m, with a->x_i exactly when e_a(b_i)=1. A core vertex is sensitive when its two incidence bits differ. Both polarities are allowed. Put

    u_m(b)=1+2^m+sum_(i=1)^m b_i 2^(m-i).

Fix q>=1 constructors and real alpha>0,beta. The chosen constructor may depend arbitrarily on the input b. The cores, their orders, internal arcs and incidence rules may be chosen afresh for every m, provided their number is at most q.

**Theorem.** If

    m-8 >= (q+1)(2q-1),                         (1)
    3^(m-8) > alpha*2^(m+1)+beta,                     (2)

then some b in {0,1}^m has

    H(T_j(b)) != alpha*u_m(b)+beta for EVERY j=1,...,q.

In particular no fixed finite menu can realize this affine target on all sufficiently large indices. The menu may select a different constructor at every target; this is stronger than excluding one uniform template on the entire cube. There is still no bound on the order of any individual core.

For the specific extractors 6u-1 and 6u+1, condition (2) holds for all m>=28. Thus it suffices to take

    m >= max(28, 8+(q+1)(2q-1)).                  (3)

Indeed 3^20=3486784401 > 6*2^29+1=3221225473, and increasing m preserves the strict inequality because the left side triples while the exponential part of the right side doubles. In particular any menu of at most three constructors fails on the block 2^28+1<=u<=2^29. The numerical bound is not claimed sharp; the companion's one-template bound m=13 is sharper when q=1.

## 2. Pin eight bits to force every two-sensitive constructor large

Restrict to the subcube of words with fixed prefix0011 and fixed suffix0011, leaving the middle t=m-8 bits free. The prefix and suffix are disjoint under (1).

Take any constructor with at least two sensitive core vertices a,c. Consider only the induced tournament on a,c and all m digit vertices. We construct at least 3^t Hamilton paths in this subtournament, for EVERY assignment to the free middle bits.

Use paths of the form

    L, a, M, c, R,

where each of L,M,R is an increasing subsequence of the digit chain. Choose two distinct prefix vertices as first(M) and first(R), with bits ensuring respectively

    a->first(M), c->first(R).

Each sensitive incidence prescribes one bit. Prefix0011 contains two copies of either bit, so distinct choices exist even if both prescriptions coincide.

Choose two distinct suffix vertices as last(L) and last(M), with bits ensuring

    last(L)->a, last(M)->c.

Again suffix0011 supplies distinct choices for either pair of prescriptions. Put every unused prefix vertex in L and every unused suffix vertex in R. Assign each free middle vertex independently to L,M or R, keeping each bin in increasing order.

The prescribed extrema remain exactly the four chosen vertices: L ends in its chosen suffix vertex, M begins in its chosen prefix vertex and ends in its chosen suffix vertex, and R begins in its chosen prefix vertex. Every other required edge is an increasing-chain edge. In particular all three bins are nonempty, and the internal apex arc a-c is not used.

Every one of the3^t assignments therefore gives a Hamilton path, and different assignments give distinct vertex words. Inserting the omitted core vertices by deterministic tournament path insertion preserves an injection. Consequently, throughout this SAME pinned subcube and for EVERY constructor having at least two sensitive vertices,

    H(T_j(b)) >=3^t.                                  (4)

The pinned bits do not depend on the core, its sensitive polarities, or the chosen constructor. Only the four selected extrema can differ between constructors. This uniformity is essential for a target-dependent menu.

Since u_m(b)<=2^(m+1), condition (2) implies that no such constructor can hit the target anywhere in the pinned subcube.

## 3. Strict pair curvature for the remaining constructors

For a constructor with exactly one sensitive core vertex, the companion's physical cut-capacity identity proves that H(T_j(b)) is a quadratic multilinear polynomial in the bits and every pair coefficient is strictly negative. More explicitly, removing the sensitive vertex leaves a fixed tournament U, and each pair coefficient is -c_(x_i,x_k), where c is a positive exposed-path capacity. Both sensitive polarities give the same negative sign. The physical bound is c_(x_i,x_k)>=2H(U-{x_i,x_k})>=2.

Fixing eight bits changes only the constant and linear terms and preserves every pair coefficient between the remaining t bits. Subtracting the affine target alpha*u_m+beta also changes only those terms. Thus the residual polynomial f_j has every pair coefficient strictly negative on the free subcube.

A zero-sensitive constructor has a constant count on all inputs. Its residual cannot vanish at two points having distinct u_m, because alpha>0.

## 4. A colored rectangle defeats a finite menu

Suppose every point of the pinned subcube is hit by the menu. By Section 2 choose a successful constructor with at most one sensitive vertex at every point, and color the point by its index, using at most q colors.

Set K=q(2q-1). By (1), choose two disjoint free coordinate blocks I,J of respective lengths 2q-1 and K, fixing any other free bits at zero. Order the coordinates in each block. For row r=0,...,2q-1 and column s=0,...,K, let b(r,s) have ones on the first r coordinates of I and the first s coordinates of J, besides the fixed prefix and suffix. This is a 2q-by-(K+1) grid.

In each column, let n_1,...,n_q be its color-class sizes, including zero sizes. Since their sum is2q,

    sum_j binom(n_j,2) >= q.

For example, this follows from sum n_j^2>=4q by Cauchy--Schwarz. Thus every column supplies at least q labels of the form (unordered row pair, common color). There are only

    q*binom(2q,2)=q^2(2q-1)=qK

possible labels. Across K+1 columns there are at least q(K+1)>qK occurrences. A label cannot occur twice in the same column, so two distinct columns share one. Their two rows give a rectangle with all four corners of one color j, hence all four residual values f_j are zero.

If j has no sensitive vertex, even two of these corners contradict its constant count and strictly increasing u_m along each grid coordinate.

Otherwise let I' be the nonempty set of coordinates added between the two chosen rows and J' the nonempty set added between the columns. These sets are disjoint. Since f_j is quadratic, its rectangle mixed difference is

    f_j(b_22)-f_j(b_21)-f_j(b_12)+f_j(b_11)
      = sum_(a in I', c in J') coefficient_(b_a b_c) <0.

All four terms vanish, a contradiction. This proves the theorem.

The argument does not infer a small success set from a large average. The two-sensitive part is a deterministic lower bound on a common subcube; the one-sensitive part is a finite colored-rectangle obstruction. Both features are needed to permit arbitrary target-dependent selection.

## 5. Quantitative consequence and scope

For fixed alpha,beta and sufficiently large m, any uniform-incidence menu covering every target in the full dyadic block must have size q satisfying

    (q+1)(2q-1) > m-8,
    equivalently q > (sqrt(8m-55)-1)/4.                            (5)

Hence menu size must grow at least on the order of the square root of the digit length. This is a necessary lower bound, not a sufficient construction theorem, and not a computational complexity lower bound for unrestricted tournament generation.

The result still permits positional incidence rules, changing chain-to-chain edges with the target, a different or nonlinear encoding of u, or an unbounded menu. It concerns extracting the chosen affine function of the given digit index, not general membership of a Hamilton-path count. In particular none of the individual target values in the obstructed block is proved absent from the unrestricted spectrum.

The specialized pumping packet `PROOF.md` additionally identifies the precise surviving coefficient H(B)H(S) of g for grafts retaining E(u) that avoid its explicit base-three pumps. That structural information is complementary to, but not needed for, this broader finite-menu obstruction.

## 6. Reproduction and inherited routes not excluded

Run `python3 verify_menu.py` and `python3 -O verify_menu.py` from this packet directory, and compare their JSON receipts. The verifier constructs 74,648 explicit distinct pinned paths across 504 word/polarity/arc cases, compares 120 literal graphs by independent subset Hamilton-path DP, checks all 1,200 relevant color-class compositions through q=5, and tests 800 reproducible colored-grid/negative-quadratic controls. These are bounded implementation checks; Sections 2--4 prove the all-order statement.

The independently written companion `session203-complement-gadget/anchor_controls.py` and its normal/optimized receipts separately check 120 literal graphs, 2,072 explicit paths and 46,148 color-class histograms through q=7. The base companion verifier checks the strict cut-capacity identity itself. The two verifier paths are intentionally different.

This theorem does NOT exclude the inherited twelve-form or forty-form occupancy conjectures. Those constructions use two freely chosen source indices, their correlated source coordinates, and a different physical interface; they do not extract the affine function of one prescribed uniform digit word. Nor does it forbid an unbounded target-dependent core menu, a positional encoder, or an actual graph rebuilding operation.
