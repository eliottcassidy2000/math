# Fixed-interface attachment hulls and correlation constraints

19 October 2026. Status: PROVED reductions and applications of classical correlation inequalities, root reviewed and approved. Exact normal and optimized controls pass. No new morphology exclusion, seven-thin theorem, or historical novelty is claimed.

This continues the exact [bounded-chain model](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md) and [suffix-state representation](../boundary_suffix_states_20261019/BOUNDARY_SUFFIX_STATE_REPRESENTATION.md). Its two gains are different: arbitrary passive attachments can be compressed while retaining a fixed active interface; and actual suffix coefficient arrays satisfy nonlinear correlation constraints beyond the earlier flow/capacity system. The lower correlation theorem is classical Fishburn, not a claimed new inequality; its integration with the exterior-state model is the point here.

## 1. Preserve the active interface before compressing

Keep the full critically oriented root core Q in a k-thin P=Q union A union B. On the lower side let T be the vertices of A incomparable with some core vertex, and let A0=A minus T be the passive part. Thus A0 is below all Q and no member of T lies below A0.

Define

    C={c in A0 : c is incomparable with at least one t in T},
    D=A0 minus C,
    R=T union C.

Every d in D lies below every t in T: it is comparable with all T, and cannot be above an active vertex. Further, no c in C can lie below d in D, because c<d<T would make c comparable with every member of T. Consequently every relation between D and R is d<r or incomparability. D is also below all Q. This is the one-sided gap-compression setting, without needing additional distinguished roots inside R.

For each t in T let t_Q(t) count its incomparable core neighbors and d_T(t) its incomparable active neighbors. Every member of C is counted at least once among passive incomparability incidences with T, so

    |C| <= sum_(t in T) [k-t_Q(t)-d_T(t)] <= (k-1)|T|,
    |R| <= k|T|.                                            (1)

The sharper incidence sum can be used whenever the interface is fixed.

## 2. One replacement preserves all deletion counts simultaneously

Fix a linear extension sigma of D. For each r in R let gamma(r) be the number of D letters following its final D predecessor. Then gamma(t)=0 for t in T, gamma is antitone on R, and gamma(r) is at most the number of D vertices incomparable with r.

Replace D by one chain of length max gamma (defined as 0 when R is empty), imposing exactly the same permitted terminal gaps against R. Also put the chain below all Q and B. Call the resulting lower exterior A'_sigma=R union replacement chain.

THEOREM 1. For every subset S of T, with the SAME family of replacements on the right,

    e(A minus S)=sum_(sigma in LE(D)) e(A'_sigma minus S).    (2)

In particular, all the H(S) coefficients needed by the suffix-state representation are preserved simultaneously. The statement even holds for subsets S that are not upper filters.

Proof. Condition on sigma in A minus S. Every r in R minus S must follow its last D predecessor and has no D successor. Thus, for every fixed extension word of R minus S, its interleavings with D are counted by the same terminal-gap formula as its interleavings with the replacement chain. Removing S changes neither sigma nor gamma for the remaining R labels. Sum over the R minus S words, then over sigma. The replacement is NOT recomputed separately for each S; gamma on the original R is fixed once.

The replacement poset stays k-thin together with the unchanged Q and B. Each r gains gamma(r) chain incomparabilities, no more than its old D incomparabilities. Every chain vertex misses a subset of the R vertices missed by the final original D vertex, at most k vertices. It is comparable with Q, B and its own chain. Every c in C has an incomparable neighbor in T, so gamma(c)<=k-1; gamma(T)=0. Hence

    |A'_sigma| <= k|T|+k-1.                                 (3)

The set T and all its internal relations and Q incidences are unchanged. Every added or retained passive vertex remains below Q, so the active interface is still EXACTLY T.

For seven-thin critical-root cores, |T|<=36 from the preceding note. Thus a lower replacement has at most 258 vertices. Compressing both passive sides gives at most

    |Q|+k(|T_A|+|T_B|)+2(k-1)
      <= 2k(k-1)^2+4k-2                                    (4)

vertices overall, at most 530 when k=7. These are generous universal bounds, not a claim of practical exhaustive enumerability.

## 3. Exact convex-hull statement, and what it does not improve

Normalize the fixed-interface coefficient array by h_A(S)=e(A minus S)/e(A). Dividing(2) by e(A) expresses it as a convex mixture of the bounded replacement arrays, with weights e(A'_sigma)/e(A). Every replacement is itself an actual allowed passive attachment to the same fixed T and Q.

Therefore the convex hull of all normalized H arrays from arbitrary finite passive attachments to this FIXED active interface equals the convex hull of the finite family of actual attachments having size at most k|T|+k-1. Equivalently, every homogeneous linear inequality in the H coefficients is valid for all such arbitrary passive attachments if and only if it is valid on that bounded family. The finiteness is up to relabeling passive vertices while fixing T and Q.

This does not assert that every convex mixture is one actual attachment. Individual realizability and factorization across the lower and upper sides remain separate issues.

Nor is this a stronger global linear core cone than the original <=26-vertex chain-atom hull. The first compression already implies that the convex hull over ALL actual exterior distributions equals the convex hull of its actual small chain atoms, since both inclusions hold. The present larger bound gains information by CONDITIONING ON T and retaining its H relations and product factors. Forgetting T again cannot improve a linear inequality that already uses the complete first atom hull. Its possible additional power is through fixed-interface reasoning, unions over discrete interfaces, and nonlinear realizability constraints before averaging.

## 4. Fishburn correlation inequality for two deletion branches

Let P0 be any finite poset of n>=2 vertices, and let a,b be distinct maximal vertices. Then

    (n-1)e(P0)e(P0 minus {a,b})
       >= n e(P0 minus {a})e(P0 minus {b}).                  (5)

In particular e(P0)e(P0 minus {a,b}) is strictly greater than the product of the two single-deletion counts. No bounded-degree assumption is needed. This is the two-maxima instance of Fishburn's inequality, credited to Fishburn (1984). The original attribution and normalized-ideal statement are recorded in [Chan–Pak, Multivariate correlation inequalities for P-partitions, Theorem 3.1](https://arxiv.org/pdf/2212.11954), checked at PDF page 6 on 19 October 2026. The following proof is included so that the precise hypotheses and factorials can be audited here.

Here is a finite-grid derivation from the classical FKG inequality, avoiding an unverified appeal to continuous FKG. The external input is [Fortuin–Kasteleyn–Ginibre, Proposition1, Communications in Mathematical Physics22 (1971),89–103](https://math.bme.hu/~balint/oktatas/perkolacio/percolation_papers/fortuin_kasteleyn_ginibre.pdf), primary statement checked19 October2026. It says increasing functions on a finite distributive lattice are positively correlated under a log-supermodular positive weight.

For a positive integer M, the weakly isotone maps P0 to {0,...,M} form a finite distributive lattice under coordinatewise minimum and maximum. Constant weight1 satisfies the required lattice inequality. Apply FKG to indicators that coordinates a and b exceed a specified upper threshold. Rescale by M and take M to infinity. Ordinary Riemann-sum convergence on this polytope gives positive correlation for those upper-coordinate events under the uniform measure on

    O(P0)={z in [0,1]^n : z_u<=z_v whenever u<v in P0}.

The polytope boundary and each fixed threshold hyperplane have Lebesgue measure zero, so these limits apply directly. Its volume is e(P0)/n!: the regions obtained by strictly ordering all coordinates are equal-volume simplices, and exactly the linear extensions are allowed.

Put E_a(epsilon)={z_a>1-epsilon}, and define E_b likewise. Since a is maximal, fixing z_a near1 leaves a section whose volume tends to e(P0 minus a)/(n-1)!. The excluded part relative to that limiting section lies in coordinate strips of total volume O(epsilon). Thus

    Prob(E_a)=epsilon*n e(P0 minus a)/e(P0)+O(epsilon^2).

The two maxima are incomparable. Fixing both coordinates near1 leaves a section tending uniformly to e(P0 minus {a,b})/(n-2)!, giving

    Prob(E_a intersect E_b)
       =epsilon^2*n(n-1)e(P0 minus {a,b})/e(P0)+O(epsilon^3).

Positive correlation, divided by epsilon^2 and followed by epsilon tending to zero, is exactly(5).

This does NOT say that the two “last vertex” indicators in a random linear extension are positively correlated: those events are mutually exclusive. The upper-corner events are continuous-coordinate events, and the factorial conversion is essential. Likewise, maximality cannot be omitted: for a two-element chain, using its two vertices would make(5) read1>=2, which is false.

## 5. Nonlinear constraints on the finite H interface

At a suffix state S, take two distinct available active maxima a,b in A minus S. They are also maximal in the full remainder. Put n=|A|-|S|. Equation(5) becomes

    (n-1)H(S)H(S+a+b) >= n H(S+a)H(S+b).                    (6)

Dropping the unknown size factor gives a useful homogeneous, size-free necessary condition

    H(S)H(S+a+b) >= H(S+a)H(S+b),                           (7)

with strict inequality for every actual finite remainder. For the bounded depth-r representation, use these rows only when |S|+2<=r unless more coefficient levels are explicitly retained.

This adds information beyond the earlier diamond IDENTITY for transition products. That identity simply said that two deletion orders reaching the same set have equal path weights. Equation(7) adds a genuine inequality between four different H values.

For example, the formal integer array for a two-element active antichain

    H(empty)=10, H({a})=H({b})=5, H({a,b})=2

has nonnegative stopping values(0,3,3,2). It satisfies branch denominator5 at the first level and denominator4 at the second, and it has the shared-set diamond identity. Nevertheless(7) would require20>=25. No passive attachment can realize it, for any remainder size. Thus the new constraint is not a consequence of nonnegative flow, shared-set ties and those degree-based branch bounds.

One interface making the parameters explicit has a five-antichain core {x,y,q,r,s} and maximum degree k=6. The incomparable active vertices a,b both miss s; a also misses q and b also misses r, while each lies below the other core vertices. Their core incomparability counts are2, giving branch bounds5 and then4. The impossible formal H array assigns only 4/10 of its mass to maximal gap 2 at s, even though its two terminal-active probabilities are both 1/2. The inequality detects this lost simultaneous-tail constraint.

Crucially, this very array is the sum of TWO ACTUAL arrays for that same fixed interface. Take a three-vertex passive chain, below all Q and below a, but incomparable with b. Its array is (5,4,1,1). Interchange a and b to obtain (5,1,4,1). Both ambient posets are six-thin and each attains equality in the exact-size lower bound with n=5. Their sum is (10,5,5,2). Thus the equally weighted mixture of their normalized arrays violates even the size-free lower inequality. This supplies a small exact demonstration that fixed-interface convexification loses genuine nonlinear realizability information. The special compression mixture in Theorem 1 comes from an actual attachment and therefore still obeys its inequalities; arbitrary mixtures need not.

## 6. Full filter-lattice form and an upper correlation bound

For a fixed original exterior size N=|A|, define

    V(S)=H(S)/(N-|S|)!.

The exact-size lower diamond (6) says V(S)V(S+a+b)>=V(S+a)V(S+b). These local inequalities imply full log-supermodularity on the upper-filter lattice:

    V(S intersect T)V(S union T)>=V(S)V(T).                  (8)

Here S and T are upper filters contained in the active interface; all four counts refer to the SAME original attachment. In the depth-r representation use only |S union T|<=r.

For completeness, take logarithms f=log V. The local diamond says that an addable element's increment does not decrease when another incomparable available element has already been added. If X is contained in Y and a is addable at both, build Y from X through legal additions. Every added element is incomparable with a: an element above a was already in X, and an element below a could not belong to the filter Y without a. Chaining the diamonds shows that the a-increment at Y is at least the a-increment at X. Now add the elements of S minus T in a legal order, once starting from S intersect T and once starting from T. Sum these increment comparisons to obtain (8).

This is exactly Fishburn's normalized-ideal theorem with complementary lower ideals. The factorial depends on the common ambient size N; do not normalize each H entry with an independently compressed exterior size.

There is also a complementary external theorem. By duality, [Chan–Pak, Correlation inequalities for linear extensions, Theorem 1.1](https://arxiv.org/pdf/2211.16637) gives

    H(S)H(S+a+b)<=2 H(S+a)H(S+b).                           (9)

Their primary paper, dated September 2024 and published in Advances in Mathematics 458 (2024), article 109954, states this for two minimal elements and n>2. The dual yields maxima; the n=2 antichain case is an immediate equality. This upper bound uses their combinatorial-atlas theorem, not the elementary FKG argument above.

Consequently each actual available-maxima diamond ratio lies in [n/(n-1),2], or in the size-free interval [1,2]. The formal array (10,2,2,1) passes the same nonnegative stopping and branch-denominator 5/4 tests, and passes the lower bound for every n>=2. It fails the upper bound because 10>8. These two hostile controls show independent missing directions in the earlier flow/capacity relaxation. Neither quadratic bound is asserted to survive arbitrary convexification.

## 7. More than two simultaneous maximal gaps

For any set W of m distinct maximal vertices of an n-vertex P0, repeated finite-lattice FKG applied to the upper-coordinate events gives

    (n)_m e(P0)^(m-1) e(P0 minus W)
        >= n^m product_(a in W)e(P0 minus a),               (10)

where (n)_m=n(n-1)...(n-m+1). The same upper-corner calculation supplies the leading joint coefficient. The m=1 case is equality; m=2 is (5).

Equality occurs when P0 is a disjoint union of chains and W contains their maximal vertices: multinomial extension counts verify every factor exactly. This includes the antichain and the earlier passive-chain-plus-independent-active-vertex example. No classification of equality cases is claimed.

If I_A(q) is an antichain of s vertices, it is also an upper filter, so its members are maximal in A. Combining (10) with the exact maximal-gap mass formula gives

    Prob[alpha(q)=s]
      >= s! * n^s/(n)_s * product_(a in I_A(q)) Prob[a is last]. (11)

This is an additional joint-mass restriction on an actual profile distribution. It is not obtained by assuming independence among last-position events.

## 8. Exact controls and remaining scope

Run from this directory:

    python verify_attachment_correlation.py
    python -O verify_attachment_correlation.py

Both runs have identical JSON receipts and 49,932 explicit guards. The selected compression controls include two D orders with different retained profiles, a stripped universal prefix, a deletion that is not an upper filter, an empty active interface and an empty passive part. Each fixed-order interleaving count is independently checked against direct extension enumeration for every S before the totals and normalized mixtures are compared.

The correlation controls cover all 407 naturally labeled posets of orders 1 through 5: 4,485 diamonds, 5,458 multiway rows and 34,278 full filter-lattice rows. This is an implementation guard, not the mathematical proof. The receipts also record the two actual six-thin constituents of the forbidden convex mixture and the upper-bound hostile array. Script and receipt hashes are frozen in this packet's manifest.

The fixed-interface H hull is finite, but its bounded generators may still be numerous and depend on the interface. The flow, capacity, exact support-max, rank-one and correlation constraints together are necessary; no proof here makes them a sufficient realizability description. A new balance or Gold exclusion still requires applying them to a genuinely surviving interface or morphology and covering the stated universe exactly.
