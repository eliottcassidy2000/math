# Growing five-vertex carriers: unbounded clock valuation and an exact tail inverse

2026-11-04. **PROVED:** normalized physical carriers with a five-vertex strong component have arbitrarily high dyadic coefficient valuation, including the required odd coefficient 63. A fixed-front/arbitrary-transitive-tail subclass nevertheless has an exact all-length and all-clock inverse at any prescribed constant. It does not repair the missing +4 branches. **INDEPENDENTLY VERIFIED:** the six positive SCC7/8 witnesses in the companion closure packet. No global spectrum, infinite branch closure, or native guard theorem is asserted.

## 1. Interface and retained graph

Use the inherited seven-coordinate signature `(h,p,F,G,J,e,s)` from S182 `induction/PROOF.md`. A carrier `(B,W)` has actual named rows and actual port subset W. The strong attachment consists of nonempty TT_t→B and root r, with W→r and r→TT_t and B\W. Its count is

    H_t = c + 2^t F,       c=h−J.

The inherited terminal-SCC port condition is the precise strongness guard. These are direct adjacency constructions, not programs of native Hamilton-complement moves.

Write T(u) for the transitive tournament whose vertices are the binary positions of positive u, marked exactly at its set bits. Its last vertex is marked. Its signature has `h=1`, `F=u`, `s=u mod 2`, and J=j(u), where the exact physical scan starts g=j=e=0 and at each low-to-high bit b performs

    j <- j+b(g+e),       g <- 2g+b,       e <- b.

For an arbitrary fixed marked front Q with signature `(h,p,F,G,J,e,s)`, the ordered sum Q→T(u), retaining both port sets, has

    F_new = F+p u,
    c_new = h−J−h j(u)−e(u mod 2)−G u.                 (1)

This follows from the inherited ordered-sum law, not a new substitution principle. The count is exact and the joint correction is explicitly transported. Because the tail's last vertex is marked, every such full attachment is strong. The front may be strong or not for this statement. With no tail, use Q itself and check its terminal guard separately.

## 2. Actual normalized five-vertex carriers with unbounded valuation

Take the strong five-vertex tournament

    rows = (2,12,17,21,3),       W={3}.

Rows are outgoing-neighbor bitmasks in physical name order 0,...,4. Direct endpoint subset DP gives

    sigma(Q,W) = (11,100,12,27,0,1,4).

For every integer m≥0 put

    u_m = (63·2^(20m+8)−3)/25.

This is a positive integer: 2^20=1 modulo25 and 63·2^8=3 modulo25. Attach the canonical transitive tail T(u_m). Equation (1) gives

    F_m = 12+100u_m = 63·2^(20m+10),
    c_m = 11−11j(u_m)−(u_m mod2)−27u_m.                (2)

Thus v2(F_m)=20m+10 is unbounded, while the odd part of F is exactly63. The first component is already marked, so there are no initial unmarked singleton components to remove. These are actual terminal-port carriers with growing tail order and actual growing port sets, not formal signatures or leading zero padding. The first mask is645, on ten tail vertices.

This is a concrete limitation of extrapolating the S184 singleton/C3 normalization or the SCC≤4 valuation bound to SCC5. A finite list of dyadic scales cannot be justified merely by removing leading unmarked singleton blocks once SCC5 is permitted. This does not prove that the fixed constants −35 or −29 occur at high scales: equation (2) keeps the defect, and c_m≤11−27u_m tends to minus infinity. No high valuation is mistaken for a successful return.

For source F=21, the required next coefficient has odd part63, so this family lives on exactly the coefficient ray relevant to the missing +4 siblings. Its constant fails the required values, which is the precise remaining obstruction in this family.

## 3. Exact inverse for one fixed front and arbitrary growing tail

Assume Q has a nonempty port set, hence G>0. For a prescribed constant C, equation (1) implies

    0 < u ≤ floor((h−J−C)/G).                         (3)

All omitted terms are nonnegative. Enumerate exactly these integers u and evaluate (1). This is a complete inverse for every positive-length transitive tail; also inspect the empty-tail case separately. If the bound is negative there is no positive tail.

For a fixed source coefficient F_source, a fixed-clock-shift digit return exists in this subclass precisely when one of these actual carriers has `c_new=3c_source+d` and oddpart(F_new)=oddpart(3F_source). The matching integer clock shift is determined uniquely from the powers of two. The candidate and source clocks must both be at least one. This test covers every fixed integer clock shift, not only a bounded shift search, because the constant itself bounds u.

The finite certificate applies this theorem to every marked strong front of orders1,3,4,5 and both target constants −35 and −29, requiring odd coefficient63. There are no hits. Every strong type through order5 and every actual port mask is retained in `bank.json`; `explore.py` generated these types by exhaustive labeled tournaments and canonicalization over every permutation. This is an all-length statement for **one marked front followed by an arbitrary transitive tail**, not for arbitrary words of SCCs≤5, not for unmarked fronts, and not for multiple large components.

As an auxiliary exact fixed-coefficient observation, arbitrary words over the complete SCC≤5 bank have finite search at every prescribed F. For any prefix P before a nonempty marked final component, F_total≥p_P≥2^|P|: each subset contributes at least one ordered two-path cover, and the final component's marked F is at least one. Therefore a carrier whose final SCC has order at most b has total order at most floor(log2 F)+b. Positivity in the ordered-sum recurrence justifies pruning when an accumulated F exceeds the target or prefix p exceeds it. This is an all-order inverse at fixed F, not an all-clock normalization. The exploratory fixed-coefficient searches are not needed for the theorem above.

## 4. Verification and independent escape controls

`verify.py` uses explicit exceptions and the inherited dependency-free direct subset DP copied with attribution in `inherited.py`. It checks the exact front signature, twelve members of (2), and the full literal suffix and strong attachment at the first member. It exhausts the constant-bounded inverse for all marked strong fronts through order5. Finally it independently recomputes all seven coordinates, strongness, and full attachment counts at carrier clocks1,2,3 for each of the companion closure packet's six physical SCC7/8 witnesses.

The companion witnesses demonstrate that the growing-tail exclusion is scoped: a different strong internal component really does escape. Their all-clock identities and address interpretation are proved in that packet. Here the literal rows, masks and resulting numbers are checked independently; no supplied signature is trusted.

The positive family (2) supplies unbounded physical clock coefficients but not regenerative digit closure. The next constructive requirement remains joint control of F, c, terminal ports and legal nonnegative clock, with actual return to a repeatable class. Strongness of direct constructions is not native-operation legality.
