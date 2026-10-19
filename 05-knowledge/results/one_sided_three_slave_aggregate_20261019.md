# Aggregate balance with at most three slaves on one side

Date: 2026-10-19. Status: proved all-size structural inequality, with a direct exact-polynomial certificate. The four-slave extension has an explicit actual-poset obstruction. No claim of settling all seven-thin posets or of priority over the literature.

## 1. Theorem and exact scope

Let (x,y) be an ordered critical incomparable pair in a finite poset P. Assume it has no lower slaves and at most three upper slaves:

    D(x)=D(y),
    U(y) ⊆ U(x),
    T=U(x)\U(y),       m=|T|≤3.

Thus every a in T satisfies x<a and a∥y. Let

    Z=I(x)∩I(y),       c=|Z|.

The induced orders on T and Z are arbitrary. All other vertices may form an arbitrary exterior consistent with these root neighborhoods.

**Theorem.** For every such P,

    sum_{z in Z} [p(x<z)+p(z<y)] ≤ 4c/3.              (1)

The order-dual statement holds with at most three lower slaves and no upper slaves. The constant is sharp for m=3: the poset

    x<a1<a2<a3,      z<a3,      y isolated

has 24 extensions and aggregate numerator 32, so its sole summand is exactly 4/3.

**Maximum-span consequence.** In a hypothetical poset with no balanced pair, an incomparable pair of maximum index span in its >2/3 majority order cannot have no lower slaves and at most three upper slaves, or the dual configuration. For c>0 this follows from (1); a separate finite suffix certificate in Section 7 covers c=0.

For seven-thin posets, the new one-sided branches are (lower slaves, upper slaves)=(0,2), with 0≤c≤4, and (0,3), with 0≤c≤3, together with their duals. The zero/one-slave cases overlap earlier results.

## 2. Reduction from arbitrary exterior to bare cores

Let Q contain x,y,Z,T. Every exterior point lies below both roots or above both roots; call these sets A and B. We have A<B, no Q point lies below A, and no Q point lies above B. For any legal core word o,

    number of extensions restricting to o
      = e(A)e(B)L(pre(o))R(suf(o)),                  (2)

where pre(o) precedes the first root, suf(o) follows the last root, and L,R are their normalized extension counts. All A letters precede the first root and all B letters follow the last root; concatenation proves the factorization.

Define the nonnegative prefix and suffix increments

    U_u=L(u)-L(u without its last letter),
    V_v=R(v)-R(v without its first letter).

Appending or prepending the removed letter is a valid injection, because no core point is below A or above B. Each normalized feature is one plus the sum of increments along its prefix or suffix tree.

Expand the aggregate deficit 4cE0-3N0 in these increments. A coefficient fixes initial core word u and terminal core word v. If those endpoints are impossible, the coefficient is zero. Otherwise delete them to obtain a bare core R with k common neighbors and m'≤m upper slaves, still with no lower slaves. Every deleted common neighbor lies outside both roots and contributes one to its two comparison indicators. Therefore the coefficient is exactly

    (c-k)e(R) + [4k e(R)-3N_R].                     (3)

It suffices to establish the bare-core inequality for all m'=0,1,2,3. This family is closed under the required deletions. In particular, deleting upper slaves introduces no lower slave.

## 3. Two conditioned orders and cumulative thresholds

A bare core has only the following relations besides the internal orders of Z and T:

    x<a for every a in T,
    z<a for selected pairs z in Z, a in T,
    y isolated.

The opposite cross-direction a<z would force x<z and is impossible.

Partition its extensions by the induced common order w1,...,wc and the induced upper-slave order a1,...,am. Adding these two chains leaves exactly

    w1<...<wc,
    x<a1<...<am,
    w_(Lj)<aj when Lj>0,                            (4)

where 0≤L1≤...≤Lm≤c. Explicitly Lj is the cumulative maximum of the common indices originally below any ai with i≤j, or zero if there is none. This cumulative maximum matters: the newly imposed upper chain propagates the predecessors of earlier upper letters to later ones.

The conditional extensions are uniform extensions of (4), with y isolated. Different conditioned word-pairs are mixed with their actual extension counts, not with equal weights. Thus a homogeneous extension-count inequality for every threshold instance proves it for arbitrary bare common and upper-slave posets.

## 4. Delete the isolated root before counting

Fix a threshold instance with m upper slaves, and temporarily delete y. A remaining word consists of the common chain and the single special chain

    x,a1,...,am.

Write j_x,j_1,...,j_m for the numbers of common letters before those special letters. There are c+m+2 gaps in which to reinsert y, all allowed.

For this fixed word, summing the aggregate indicator over all insertions of y gives

    n(word)=c(3c+2m+5)/2
              -(c+m+1)j_x + sum_{h=1}^m j_h.       (5)

Indeed, the sum of the common-letter ranks before y is inserted is

    c(c+1)/2+(m+1)c-j_x-sum_h j_h.

There are c-j_x common letters after x. Summing the number before y over its insertion gaps now gives (5).

Let E_base be the number of these y-deleted words. Then

    E=(c+m+2)E_base,
    N=sum_word n(word),
    D=4cE-3N.

The contribution of one y-deleted word to D is

    d(word)=c(2m+1-c)/2
              +3(c+m+1)j_x-3 sum_h j_h.            (6)

These formulas eliminate the need to enumerate the possible positions of y separately.

## 5. A finite shuffle formula valid for every c

Set the nonnegative integer gap lengths

    v0=L1,
    vi=L_(i+1)-Li for 1≤i<m,
    vm=c-Lm.

For m=0, use the single length v0=c. The common-chain gap indices are partitioned into m+1 consecutive blocks. The first m blocks have v_i available gap indices; the last has v_m+1, including the gap after the last common letter.

Split the special chain x,a1,...,am into consecutive parts of lengths h0,...,hm, where their sum is m+1. Require that aj belongs to block index at least j. This enforces exactly j_h≥Lh. There are respectively 1,2,5,14 admissible cuts for m=0,1,2,3.

For a cut h, its exact extension weight is

    W_h = product_{i<m} binom(vi+hi-1,hi)
                       times binom(vm+hm,hm).       (7)

For a special letter that is kth within its block, its mean common rank is the block offset plus

    k(vi-1)/(hi+1),    i<m,
    k vm/(hm+1),      i=m.                          (8)

These formulas follow from weakly increasing gap indices. Adding 1,2,...,h to such indices converts them to a uniform h-subset, whose order-statistic means give (8). If vi=0 for i<m, the corresponding count is zero when hi>0 and one when hi=0. Formal means in zero-weight cases are finite and harmless; no occupied empty block is counted.

Since (5)-(6) are affine in the common ranks, replace each rank by its mean in a fixed cut. Then

    E=(c+m+2) sum_h W_h,
    N=sum_h W_h n(mean ranks),
    D=sum_h W_h d(mean ranks).                      (9)

This directly derives rational polynomials in the integer variables v_i. It is a finite symbolic identity, not interpolation or a bounded-size inference.

## 6. Nonnegative certificates for m=0,1,2,3

The standalone standard-library verifier expands (9) with exact rational arithmetic and records every coefficient.

### No upper slaves

The roots are exchangeable, so each summand in (1) is one. Algebraically D=cE≥0.

### One or two upper slaves

For m=1, the degree-four deficit has 13 positive monomial coefficients. For example, writing v0=p,v1=s,

    D=(p+s)(s+1)(4p²+3ps+10p+s²+5s+6)/4.

For m=2, the degree-five deficit has 47 positive monomial coefficients. All exact coefficients are in the accompanying JSON certificate. Because the gap variables are nonnegative, these polynomials are nonnegative.

### Three upper slaves: an integer-boundary certificate

Write (v0,v1,v2,v3)=(p,q,t,s). The degree-six deficit has some negative ordinary coefficients, so it is important not to claim positivity merely by inspection of all monomials.

Its boundary restriction p=s=0 is exactly

    B(q,t)=(t+1)(q+t-1)H(q,t)/24,                   (10)

where

    H(q,t)=6q³+6q²t+18q²+4qt²+14qt+12q
             +t³+5t²+6t.

Every coefficient of H is positive and H(0,0)=0. Since q,t are nonnegative integers, either q+t≥1, making every factor in (10) nonnegative, or q=t=0, making B=0. Hence B≥0 on the required integer domain.

The exact polynomial difference

    D(p,q,t,s)-B(q,t)

has 145 strictly positive monomial coefficients. Each such monomial contains p or s. The verifier derives D from the 14 admissible shuffle cuts, proves equality of its boundary with (10) coefficient by coefficient, and checks every coefficient of the difference. Thus D≥0 for all integer gap lengths.

The complete deficit contains 164 monomials: the 19 boundary monomials plus the 145 positive remainder monomials. The certificate is specifically an integer-domain certificate; no positivity claim is made for fractional gap lengths.

Equations (7)-(10), the conditional-word mixture, and the exterior lift (3) prove Theorem (1).

## 7. Maximum-span consequence

In a poset with no balanced pair, the >2/3 comparison tournament is transitive: a cyclic triple would have three probabilities summing to more than two, although at most two cyclic indicators can hold in any extension. Its order M extends P.

An incomparable pair x<_M y of maximum index span is critical. Every common neighbor z lies strictly between the roots in M; otherwise one of its spans with a root would be larger. Therefore each summand in (1) would exceed 4/3. When c>0 this contradicts (1).

### The c=0 endpoint: one-pair suffix certificates

When c=0, the core Q consists only of x,y and its m upper slaves. All lower exterior A lies below the entire core, since A<x<a for every upper slave a. Thus only the suffix factor in (2) varies:

    E0=sum_o R(suf(o)).

Let M be the actual majority order projected to Q. It is a linear extension of Q with x before y. No exterior point lies between x and y, so the root span in this projected order is exactly its full span. Projection can only decrease other core-pair spans. Hence every incomparable core pair has projected span at most the root span.

There are finitely many possible cores and projected orders, because m≤3. Naturally label the upper poset using one of its linear extensions. The exact counts are:

| Upper slaves m | Naturally labelled upper posets | Admissible projected orders |
|---|---:|---:|
| 0 | 1 | 1 |
| 1 | 1 | 2 |
| 2 | 2 | 6 |
| 3 | 7 | 48 |

For each of these 57 projected orders, the verifier identifies one oriented incomparable pair a<_M b with the following property. For every suffix-tree generator s, including the empty generator,

    3 N_s(a<b)-2 E_s ≤ 0,                           (11)

where E_s counts bare core extensions ending with s and N_s counts those also having a before b. A nonempty s is an ordered suffix of upper letters; every core word ending with it automatically has both roots earlier.

Expanding each R(suf(o)) into its nonnegative suffix increments V_s gives

    3 N_P(a<b)-2 e(P)
       = e(A)e(B) sum_s [3N_s(a<b)-2E_s] V_s ≤ 0.

Therefore p(a<b)≤2/3, contradicting its forward direction in the >2/3 majority order. This proves the c=0 exclusion.

The finite certificate is exact integer enumeration of the indicated tiny rooted family, not a search over all small posets. The receipt records each core's predecessor masks, each admissible projected order, its chosen oriented pair, and every generator's E_s, N_s, and margin. There is no linear programming, floating-point calculation, or exterior-size cutoff.

Together with the c>0 aggregate argument, this establishes the full stated one-sided maximum-span exclusion. It does not claim to handle other two-sided slave multiplicities.

## 8. Sharp limit and failure of baseline monotonicity

The four-upper-slave extension is false for this aggregate. Take

    x<a1<a2<a3<a4,     z<a4,     y isolated.

There are 35 extensions, with N(x<z)=28 and N(z<y)=20. Consequently

    p(x<z)+p(z<y)=48/35>4/3.

The pair (z,y) is nevertheless balanced, with probability 4/7. This is an obstruction to the aggregate extension, not to the 1/3–2/3 conjecture.

The more ambitious claim that adding common-to-upper relations decreases the disjoint baseline also fails already with two upper slaves. The same family at m=2 has aggregate 19/15, exceeding the disjoint baseline 5/4. Thus the proof above uses exact threshold certificates and does not assume attachment monotonicity.

The companion one_sided_slave_aggregate_obstruction_20261019.md gives the all-m obstruction formula

    4E-3N=-(m+4)(m-3)/2,

and its separate verifier checks the actual poset counts.

## 9. Reproduction and independent implementation checks

Run:

    python 05-knowledge/experiments/verify_one_sided_three_slave_aggregate_20261019.py

No third-party Python packages are required. The verifier:

- derives the four symbolic count/deficit polynomials directly from the 1,2,5,14 shuffle cuts;
- verifies the positive coefficients and the factored integer-boundary certificate;
- compares E,N,D with an independent chain-state recurrence for all 12,649 threshold instances with c≤20 and m≤3;
- checks 400 actual posets with arbitrary upper-slave/common orders and exterior attachments, through total order 18, using independent whole-poset subset DP;
- checks all 57 c=0 projected maximum-span orders with exact one-pair suffix certificates;
- verifies the exact m=4,c=1 obstruction.

The numerical checks are implementation audits. The direct symbolic identities and the reductions prove the arbitrary-size result. The full receipt and sparse polynomials are in one_sided_three_slave_aggregate_audit_20261019.json.

The publication verifiers retain explicit guards under `python -O` and emit the full JSON receipts to standard output. Run the obstruction verifier as `python 05-knowledge/experiments/verify_one_sided_slave_obstruction_20261019.py`.
