# An all-size aggregate inequality for a critical pair with one slave on each side

Date: 2026-10-19. Status: proved structural inequality with an exact symbolic certificate and independent combinatorial audit. No claim of solving the general 1/3–2/3 conjecture or all seven-thin posets; no priority claim.

## 1. Main theorem and its use

Write p(u<v) for the probability that u precedes v in a uniformly random linear extension of a finite poset P. Let (x,y) be an ordered critical incomparable pair:

    D(x) ⊆ D(y),       U(y) ⊆ U(x).

Suppose there is at most one lower exclusive neighbor b, meaning b∥x and b<y, and at most one upper exclusive neighbor a, meaning x<a and a∥y. Every other point incomparable with either root is either the other root or a common incomparable neighbor. Put

    Z = I(x) ∩ I(y),       c = |Z|.

**Theorem 1 (aggregate inequality).** For arbitrary c, arbitrary induced order on Z, and arbitrary exterior,

    Σ[t∈Z] (p(x<t)+p(t<y)) ≤ 4c/3.                  (1)

The constant is sharp: equality holds for the disjoint union of x<a, b<y, and any c-point poset Z. If one slave is absent, the stronger bound 7c/6 holds; if both are absent, the sum is exactly c.

**Theorem 2 (maximum-span exclusion).** In a hypothetical finite poset with no balanced pair, an incomparable pair of maximum index span in its >2/3-majority order cannot have at most one slave on each side. This includes c=0.

Theorem 1 supplies actual all-size coverage beyond the earlier finite c≤3 certificates. For a seven-thin poset, it eliminates the mixed-one-slave branches (d(x),d(y),c)=(3,3,1), …, (7,7,5). Together with the known at-most-one-total-slave screen, any surviving maximum-span root pair has at least two slaves on one side. In particular it has at most four common neighbors.

Theorem 1 alone does not say that every poset containing such a critical pair has a balanced comparison among these endpoints. The use of maximum span in Theorem 2 is essential to the conclusion stated here.

## 2. Exterior counts and nonnegative tree increments

Let Q consist of x,y,Z and the slaves that are present. Every exterior point is below both roots or above both roots; call these sets A and B. We have A<B. No point of Q is below a point of A or above a point of B. These assertions follow immediately from the specified incomparability neighborhoods and transitivity.

Fix a legal word o of Q. Its prefix before the first root contains only b and common neighbors; its suffix after the last root contains only a and common neighbors. The number of extensions of P restricting to o is exactly

    e(A)e(B) L(pre(o)) R(suf(o)),                    (2)

where L(u) counts extensions of A with the prescribed ordered prefix u, normalized by e(A), and R(v) counts extensions of the prescribed suffix v with B, normalized by e(B). Both empty features are 1. Concatenating these two interleavings with the intervening fixed core word proves (2) bijectively.

For every nonempty prefix u and suffix v, define

    U_u = L(u)−L(u with its last letter removed),
    V_v = R(v)−R(v with its first letter removed).

Appending or prepending that letter gives a valid injection, so U_u,V_v≥0. Consequently each feature is 1 plus the sum of increments along its prefix or suffix tree.

Expand the normalized aggregate deficit

    4c e(P)/(e(A)e(B))
       −3 Σ[t∈Z] (N_P(x<t)+N_P(t<y))/(e(A)e(B))

in these tree increments, with U_empty=V_empty=1. A coefficient of U_u V_v restricts the core word to begin with u and end with v. An impossible restriction contributes zero. Otherwise u is an ideal word, v is a disjoint filter word, and deleting both leaves an induced bare core R with both roots, k common neighbors, and zero or one slave on each side. Each deleted common neighbor contributes exactly one to its two comparison indicators. Therefore the coefficient is

    (c−k)e(R) + D(R),                               (3)

where

    D(R)=4k e(R)−3 Σ[t∈Z(R)] (N_R(x<t)+N_R(t<y)).

Thus it is enough to prove the aggregate inequality for every *bare* core, including the versions with missing slaves. The argument retains actual boundary constraints; it does not assume arbitrary weights on core words. For a target bound λc, the corresponding coefficient is (λ−1)(c−k)e(R) plus the reduced λ-deficit, so the stronger missing-slave bound also lifts.

This exterior-lifting identity and its nonnegative tree coordinates were independently developed and checked in the companion `root4_two_common_tree_balance_20261019.md`.

## 3. Conditioning only the common-neighbor order

Consider a bare core with both slaves. Its only relations beyond the internal order on Z are

    x<a,  b<y,  b<t for t in an upset of Z,
    t<a for t in a downset of Z,  and optional b<a.

Other directions would compare a root with a common neighbor or compare the roots.

Partition its extensions by the induced word w1,…,wc on Z. In each part, imposing that chain gives precisely a threshold core:

    w1<⋯<wc,
    x<a, b<y,
    b<w_(K+1) if K<c,
    w_L<a if L>0,
    optional b<a,                                  (4)

for 0≤K,L≤c. Here K is one less than the earliest index in b's common upset, and L is the latest index in a's common downset. When K<L, b<a is forced. Conditioning on the common word preserves uniformity among extensions of (4). The mixture weights are the extension counts of these conditional cores, not equal weights on common words.

This conditioning is used only after the exterior-lifting reduction. Conditioning the common order directly inside a general exterior would not make the remaining core-word weights uniform.

## 4. The threshold deficit as three shuffle polynomials

For a threshold core let E be its extension count, N the aggregate numerator, and

    D = 4cE−3N.

In a word let i_x,i_y be the numbers of common letters before x,y. Its aggregate indicator is c+i_y−i_x. Its contribution to D is therefore

    c+3(i_x−i_y).                                   (5)

### 4.1. Gap region K≥L

Put p=L, q=K−L, r=c−K, all nonnegative. Let D0 refer to b∥a and D1 to b<a. Then

    4D0 = A0 + q A1 + q² A2 + q³ A3,               (6)

where

    A0=(p+1)(r+1)[p³+r³+2pr(p+r)+3(p²+r²)+8pr],

    A1=p⁴+r⁴+4pr(p²+r²)+4p²r²
       +6(p³+r³)+17pr(p+r)+10(p²+r²)+22pr+5(p+r),

    A2=2(p³+r³)+3pr(p+r)+8(p²+r²)+6pr+6(p+r),

    A3=p²+r²+p+r.

Moreover

    D1=D0 + (q+1)(q+2)H/60,                         (7)

where

    H=2q³+[10(p+r)+14]q²
      +[5(p²+r²)+40pr+35(p+r)+24]q
      +15(p+1)(p+r)(r+1).

Every displayed coefficient is nonnegative, so both deficits are nonnegative. For c>0, D0 vanishes exactly when p=r=0; D1 is strictly positive.

### 4.2. Overlap region K<L

Put p=K, q=L−K−1, r=c−L. Here b<a is forced. Then

    4D=(p+1)(r+1)J,                                (8)

where

    J=p³+r³+4q(p²+r²)+2pr(p+r)+7(p²+r²)
      +7q²(p+r)+8pqr+25q(p+r)+16pr+20(p+r)
      +4q³+20q²+24q+8.

This deficit is strictly positive.

### 4.3. Direct derivation, rather than interpolation

Here is a finite shuffle formula that proves identities (6)–(8). Choose one of the six orders of b,x,y,a respecting b<y and x<a; when required, also impose b<a. For h1,h2≥0 with h1+h2≤4, let h3=4−h1−h2. Split this special-letter word into three consecutive parts of those lengths.

A segment consisting of n common letters and h ordered special letters, ending in a common letter, has

    binom(n+h−1,h)

interleavings. The mean number of common letters before its j-th special letter is j(n−1)/(h+1). A segment not required to end in a common letter has binom(n+h,h) interleavings and mean jn/(h+1). These follow by the symmetry of the h+1 weak-composition gaps. For n=0 the ending-common count is zero if h>0 and 1 if h=0; formal means in zero-weight cases are harmless.

In the gap region, the three common-block lengths are p,q,r. The first two segments end in a common letter; the third is unrestricted. Disallow a in the first segment. If b occurs in the third special part, force its entire special prefix through b to precede all r common letters. Write f for the number of free special letters remaining in that part; otherwise f=h3. The extension weight is exactly

    binom(p+h1−1,h1) binom(q+h2−1,h2) binom(r+f,f).  (9)

The forced third-part letters have common rank p+q. Free third-part ranks are p+q+jr/(f+1). Earlier ranks have the segment offsets and means above. Summing (9), multiplied by (5), over the finitely many orders and cuts yields (6) or (7).

In the overlap region, the common-block lengths are p,q+1,r. Disallow a in the first two special parts and require b in those parts. If b is in the second part, force its prefix through b before every common letter of that part, and let f be the number of remaining free second-part letters; otherwise f=h2. The weight is

    binom(p+h1−1,h1) binom(q+f,f) binom(r+h3,h3).   (10)

Forced second-part ranks are p, free ones p+jq/(f+1), and third-part ranks p+q+1+jr/(h3+1). Summing (10) times (5) gives (8).

This is a direct polynomial proof. The companion standard-library verifier performs these finite sums with exact rational multivariate polynomials, proves equality with the printed formulas coefficient by coefficient, and verifies nonnegativity. It does not fit polynomials to numerical data. An independent state DP checks all 20,832 threshold parameter choices through c=30, including c=0, both relation flags, and all zero-block endpoints.

## 5. Missing slaves

If b is absent, y is isolated in the bare core. After conditioning on the common word, the only attachments are x<a and w_L<a. Delete y temporarily and classify a word by i common letters before a, L≤i≤c, and j common letters before x, 0≤j≤i. Reinsert y in its c+3 possible gaps. Direct summation gives

    E=(c+3) Σ[i=L..c](i+1),
    N=(c/2) Σ[i=L..c](i+1)(3c+7−i),
    7cE−6N=cL(L+1)(c−L+1)≥0.                     (11)

Thus the aggregate is at most 7c/6. If a is absent, dualize (11), using L'=c−K. If both slaves are absent, swapping the roots gives p(x<t)+p(t<y)=1 for every common neighbor. These cases are closed under the deletions in Section 2, so the same stronger bounds hold with arbitrary exterior.

Equations (6)–(11), the common-order mixture, and the nonnegative tree coefficients prove Theorem 1. The equality baseline follows either from (6) with p=r=0 or from uniform shuffles of the two 2-chains with an arbitrary fixed common word.

## 6. Majority-order witness migration and the c=0 endpoint

If P has no balanced pair, orient every pair by its probability exceeding 2/3. A directed triangle is impossible: in any total order at most two of its cyclic comparison indicators hold, whereas three probabilities exceeding 2/3 would have sum exceeding 2. Hence this majority tournament is transitive and defines an order M extending P.

Choose an incomparable pair x<_M y of maximum index span. If t<x and t∥y, its span with y would be larger; if y<t and t∥x, its span with x would be larger. Therefore the pair is critical. Every common incomparable neighbor t lies strictly between x and y in M; otherwise one of its spans with a root would be larger. Thus each summand in (1) would exceed 4/3. For c>0, Theorem 1 contradicts this.

For c=0 with both slaves present, only {b,x,y,a} remains in the core. Let C,D≥1 be its normalized single-prefix and single-suffix features. Direct enumeration gives

    E0=2CD+C+D+ε,
    N0(b<x)=2CD+C,
    N0(y<a)=2CD+D,                                 (12)

where ε=1 for b<a and ε=2 for b∥a. Both displayed probabilities are at least 1/2. In a no-balanced poset they therefore orient b<_M x and y<_M a. All exterior points are below both roots or above both roots, so x and y are consecutive in M. Maximum span is thus 1 for *every* incomparable pair.

The incomparability graph is consequently a disjoint union of paths on consecutive M-ranks. Different components form an ordinal sum. Within a path component of m≥2 vertices, every pair at rank distance at least two is comparable, so each extension is uniquely a collection of disjoint adjacent swaps, equivalently a matching of that path. There are F_(m+1) matchings, and F_(m−1) contain a given endpoint edge, with F_1=F_2=1. The endpoint reversal probability is

    F_(m−1)/F_(m+1) ∈ [1/3,1/2].                   (13)

It is a balanced pair, contradiction. This closes c=0. If one or both slaves are absent, the established at-most-one-total-slave screen already gives a balanced pair.

## 7. Scope, literature, and remaining seven-thin cases

Peczarski's primary 2017 *Comments on the Gold Partition Conjecture* restates the known at-most-one-slave result and the one-sided-chain mechanism:

https://cdm.ucalgary.ca/article/download/62349/46785/178291

The existing six-thin theorem is Peczarski, *The 1/3–2/3 conjecture for 6-thin posets*, Order 25 (2008), 91–103, DOI https://doi.org/10.1007/s11083-008-9081-9. The complete original six-thin configuration table has not been recovered here. A focused source search did not establish whether the present all-size mixed-two-slave aggregate inequality appears elsewhere, so no novelty or priority claim is attached.

For seven-thin posets write ℓ=d(x)−c−1 and u=d(y)−c−1 for the lower and upper slave counts of the selected critical pair. Theorem 2 rules out ℓ,u≤1. Since ℓ,u≤6−c, this forces c≤4. The cases with at least two slaves on one side remain, as does the general conjecture. The two-sided asymmetry cannot be discarded merely by averaging over a fixed exterior word.

## Reproduction and artifacts

Run:

    python 05-knowledge/experiments/verify_two_slave_all_common_20261019.py

The standard-library verifier derives the three symbolic identities, matches the printed nonnegative formulas, checks 20,832 threshold cores by an independent recurrence, and checks 1,891 instances of (11). Its output is `two_slave_all_common_audit_20261019.json`.

The companion `two_slave_allc_reduction_independent_audit_20261019.md` independently checks the common-order conditioning, exterior-deletion lift, missing-slave identity, and shuffle derivation. The earlier finite-c packet is retained unchanged, including its stronger pointwise result for c=2 and its genuine pointwise obstruction at c=3.

Both verifiers retain explicit guards under `python -O` and write JSON receipts to standard output. Run the separate exterior check as `python 05-knowledge/experiments/verify_two_slave_allc_exteriors_independent_20261019.py`.
