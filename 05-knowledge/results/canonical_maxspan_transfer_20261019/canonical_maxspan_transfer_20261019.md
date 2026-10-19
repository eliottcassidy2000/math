# Canonical maximum-span roots and exact secondary-pair transfer

Date: 2026-10-19. Status: structural lemmas with a bounded exact impact assessment. The canonical-root screen has a narrower logical scope than existing any-root morphology exclusions. No new full morphology scan or full seven-thin theorem is claimed.

## 1. Majority order, criticality and the full root core

Assume a finite non-total poset P has no balanced pair. Orient each incomparable pair in its probability direction strictly above 2/3 and keep every comparability forward. The resulting tournament is transitive: a directed triangle would have cyclic comparison probabilities summing to more than two, while each linear extension satisfies at most two of those three indicators. Let M be this total order.

Choose an incomparable pair x before y of largest index span d in M. It is critical. Indeed, if z<x but z is not below y, then z is incomparable with y, giving a larger span than x,y. This proves D(x)⊆D(y); the successor argument is dual and gives U(y)⊆U(x).

Write

    L=D(y)\D(x), T=U(x)\U(y), Z=I(x)∩I(y),
    l=|L|, u=|T|, c=|Z|,
    Q={x,y}∪L∪T∪Z, n=|Q|=2+l+u+c.

This is the complete root core: every vertex incomparable with either root is in Q. Every exterior vertex is below both roots (call this set A) or above both (call it B). Hence A<B, no exterior vertex lies between x,y in M, and their projected span in M restricted to Q is still d.

## 2. Exact tie lemma, including exact common-neighbor counts

Let s,t be any incomparable pair in Q whose span in M restricted to Q is d. Projection cannot increase span, while global maximality bounds its full span by d. Therefore its full span equals d. In particular, no exterior letter was deleted from between s and t by the projection.

Every common incomparable neighbor of any global maximum-span pair lies strictly between its endpoints in M. A common neighbor outside their interval would form an incomparable pair with the farther endpoint of span greater than d. Applying this to s,t, every ambient common incomparable neighbor lies between them, where there are no exterior letters. Consequently

    I_P(s)∩I_P(t)=I_Q(s)∩I_Q(t).

Thus the common-neighbor count of a tied projected pair is **exactly** computable from Q, not merely bounded below by its core count. Oriented by M, the pair s,t is critical by the argument in Section 1. These facts hold for any selected maximum-span root; no canonical tie-breaking assumption is needed yet.

## 3. Canonical-root screen and its precise scope

Among all global maximum-span incomparable pairs, choose x,y to maximize c=|I(x)∩I(y)|. This rule is invariant under order duality. Such a root exists because the poset is finite.

For this canonically chosen root, every core incomparable pair s,t of projected span d must satisfy

    |I_Q(s)∩I_Q(t)|≤c.

Otherwise Section 2 gives a second global maximum-span pair with a larger ambient common-neighbor count, contradicting canonical selection. A projected order violating this condition can therefore be discarded when studying the canonical root.

This is a **canonical-root** exclusion, not an exclusion for every maximum-span root. A global proof may choose a canonical root initially and cover all remaining possibilities for that choice. Existing any-root exclusions remain valid for that root as well. The canonical screen must not silently be inserted into a proof asserting the stronger any-root statement.

### Why the naive core degree-sum tie-break adds nothing

One may further break ties by maximizing the ambient sum d_P(x)+d_P(y). But the proposed test using only a tied pair's core degree sum is automatically redundant. The root degrees are known exactly:

    d_P(x)+d_P(y)=2+l+u+2c=n+c.

For any core pair s,t, writing c_Q=|I_Q(s)∩I_Q(t)| gives

    d_Q(s)+d_Q(t)
      =|I_Q(s)∪I_Q(t)|+c_Q ≤ n+c_Q.

If c_Q=c, its core degree sum cannot exceed the roots' ambient degree sum n+c. Thus this test never excludes an additional order after the common-count screen. Additional certified exterior information might still prove a larger ambient degree sum; the redundancy statement concerns the core-only lower bound, not every possible use of a degree tie-break.

## 4. Any-root transfer to a previously excluded morphology

This separate mechanism does not require a canonical root. Let s precede t be a tied projected pair as in Section 2, and put

    c'=|I_Q(s)∩I_Q(t)|,
    l0=|D_Q(t)\D_Q(s)|,
    u0=|U_Q(s)\U_Q(t)|.

The ambient common count is exactly c'. We now bound the two ambient exclusive counts without guessing the exterior.

Let S consist of the original roots x,y and every core vertex of degree seven. A saturated vertex is comparable with every exterior vertex. No core vertex can lie below A or above B, since it is either a root or incomparable with a root. Therefore A<z<B for saturated z. Define

    F_A={v : some z in S satisfies z≤v},
    F_B={v : some z in S satisfies v≤z}.

Every A is below F_A, and every B is above F_B. Set

    α(v)=0 for v in F_A, and 7−d_Q(v) otherwise;
    β(v)=0 for v in F_B, and 7−d_Q(v) otherwise.

An exterior contribution to D_P(t)\D_P(s) can only be an A-vertex incomparable with s. Its number is at most α(s). Dually, an exterior contribution to U_P(s)\U_P(t) is at most β(t). Thus

    l0 ≤ l' ≤ l0+α(s),
    u0 ≤ u' ≤ u0+β(t).

Seven-thinness also gives l',u'≤6−c'. The previously proved [two-sided two-exclusive reduction](../mixed16_c0_exclusion_20261019/two_sided_two_exclusive_reduction_20261019.md) requires l',u'≥2 in a hypothetical counterexample. Therefore every actual secondary morphology belongs to the finite integer rectangle

    l' in [max(2,l0), min(6−c',l0+α(s))],
    u' in [max(2,u0), min(6−c',u0+β(t))],
    common count exactly c'.

If the rectangle is empty, the projected order is impossible. If every morphology in the rectangle has already been excluded by an independently established **any-root** theorem, the projected order is likewise impossible. This is a safe transfer of earlier results to another maximum-span pair, not an assumption that a relaxed exterior assignment is realizable.

In particular, c'≥5 makes the rectangle empty, since seven-thinness would force both exclusive counts at most one, contradicting the two-exclusive reduction. This conclusion already follows from a lower bound c_Q≥5, but Section 2 gives exactness here.

If s∈F_A then α(s)=0 and l' is exactly l0; if t∈F_B then β(t)=0 and u' is exactly u0. Thus both exclusive counts are exact under these two one-sided conditions. Requiring both endpoints in F_A∩F_B is sufficient but stronger than necessary. Saturated endpoints provide a useful special case.

## 5. Explicit dependency catalogue and avoidance of circularity

The bounded assessment below uses only these already completed any-root results, up to duality:

- At least two exclusive neighbors on each side: [packet 99 corollary](../mixed16_c0_exclusion_20261019/two_sided_two_exclusive_reduction_20261019.md).
- (2,2,c), c≤2: [packet 95](../mixed22_c02_exact_exclusion_20261019/mixed22_c02_exact_exclusion_20261019.md).
- (2,2,3): [packet 96](../mixed22_c3_exact_exclusion_20261019/mixed22_c3_exact_exclusion_20261019.md).
- (2,3,c), c≤1: [packet 107](../mixed23_c01_exact_exclusion_20261019/mixed23_c01_exact_exclusion_20261019.md).
- (2,3,2): [packet 111](../mixed23_c2_exact_exclusion_20261019/mixed23_c2_exact_exclusion_20261019.md).
- (2,4,0) and (3,3,0): [packet 121](../mixed24_mixed33_c0_exact_exclusion_20261019/mixed24_mixed33_c0_exact_exclusion_20261019.md).

Each of these results predates this transfer lemma and was established without it. The direction of dependence is these fixed results → transfer assessment, with no reverse edge. The ongoing (2,2,4) work is not assumed. Future use must freeze the available exclusion catalogue first; a theorem must not use its own unproved exclusion as a transfer premise. Canonical-only exclusions also cannot be substituted for any-root exclusions of the secondary pair without proving that it meets their canonical assumptions. There is one useful exact transfer case: under the (maximum span, maximum common count) canonical rule of Section 3, a tied pair with c'=c is itself canonical. Previously established canonical-root exclusions may therefore transfer to such same-c pairs, provided proof dependencies remain acyclic. A lower-c pair is not canonical under that rule. If additional tie-breakers such as degree sum are imposed, their conditions must also be checked separately. The current fixture implementation uses only the prior any-root catalogue, not this future canonical-to-canonical extension.

The baseline projected-order screens use the explicitly cited known six-thin theorem for spans≤3, as in the linked certificates: Marcin Peczarski, *The Gold Partition Conjecture for 6-Thin Posets*, Order 25 (2008), 91–103, https://doi.org/10.1007/s11083-008-9081-9. The tie and rectangle arguments above are otherwise elementary. Standard majority-order, saturation and counting techniques are not claimed as original; no historical priority claim is made.

## 6. Bounded exact impact, not population coverage

The fixtures contain the first 100 residual cores used in each earlier bounded probe, for four still-open larger branches. They were chosen before this tie screen and are ordered-prefix samples, not random samples. Each fixture contains the exact predecessor masks and expected number of orders retained after the prior degree, directed-defect, maximum-span and known-six-thin screens. The checker independently re-enumerates those baseline orders, verifies the count, and applies the new tests.

| Root morphology | Sample cores | Baseline orders | Canonical common-count exclusions | Remaining orders | Nonempty sampled cores fully removed |
|---|---:|---:|---:|---:|---:|
| (2,4,1) | 100 | 1,514 | 804 | 710 | 10 |
| (3,3,1) | 100 | 1,150 | 598 | 552 | 34 |
| (2,5,0) | 100 | 3,020 | 2,000 | 1,020 | 13 |
| (3,4,0) | 100 | 2,476 | 1,400 | 1,076 | 0 |

In total, the canonical screen removes 4,802 of 8,160 previously retained orders, leaving 3,358; 57 nonempty sampled cores lose all retained orders. These figures do not establish whole-core or whole-morphology exclusions outside the stated canonical sample setting.

The any-root rectangle-transfer test removes **zero** orders in these four samples with the fixed dependency catalogue above, including the sharper α/β bounds. It is mathematically valid but provides no additional observed reduction here. This negative result is recorded rather than implying that every new lemma improves these fixtures.

Run `python -O assess.py` from this folder. The standard-library checker verifies strict transitive masks, original root morphology, actual core degrees, criticality of every tied pair, exact common intersections, the degree-sum redundancy, all finite rectangles and catalogue membership. It records every canonical exclusion's order, tied pair and common-neighbor witness. Normal and optimized controls agree. `fixtures.json`, `impact.json`, `impact_normal.json` and `impact_details.json.gz` provide self-contained exact inputs and outputs; no numerical solver is used. No existing certificate engine or ongoing scan was modified.
