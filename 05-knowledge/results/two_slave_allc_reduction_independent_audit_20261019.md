# Independent audit of the all-c chain reduction

Date: 2026-10-19. Scope: independent proof audit of conditioning, arbitrary-exterior lifting, missing-slave cases, and the exact shuffle-block polynomial construction. The complete all-c argument passes this audit; numerical checks are supplementary.

## 1. Conditioning on the common vertices is exact

Let a bare core consist of roots x,y, common incomparable neighbors Z, and optional slaves b,a. When present, b is below y and can only be below common vertices; a is above x and can only be above common vertices. The optional relation b<a is allowed. There are no other root-incident relations in the bare core.

Fix a linear extension W1,...,Wc of the induced common poset P[Z], and condition the full extension on inducing exactly that common order. This is equivalent to adding the chain relations W1<...<Wc. It is important to fix only the common order, not an order on the slaves as well.

Let

- K+1 be the earliest index i such that b<Wi; put K=c if there is no such i;
- L be the latest index i such that Wi<a; put L=0 if there is no such i.

Every original b-to-common relation becomes equivalent to the single threshold b<W_(K+1). Every common-to-a relation becomes equivalent to W_L<a. Relations with an out-of-range endpoint are omitted. Thus the conditioned core has exactly the form

    W1<...<Wc,
    b<W_(K+1), W_L<a,
    x<a, b<y,
    and possibly b<a.

If K<L, the threshold chain already forces b<a. The boolean recording an optional direct b<a is consequently redundant in that region. No added common-chain relation can introduce a root/common comparability: x has only its outgoing relation to the maximal peripheral point a, and y only its incoming relation from the minimal peripheral point b.

The conditioned extensions are uniform linear extensions of this chain-threshold poset. Different common orders generally have different numbers E_W of compatible full extensions. The unconditional probability is their mixture with weights E_W/sum E_W. Therefore proving a homogeneous count inequality separately for each chain-threshold core proves it for every bare core. Uniform averaging over common words would be unjustified, but is not needed.

## 2. The arbitrary-exterior lifting remains valid for all c

All points outside the roots' incomparability union split into A below both roots and B above both roots. For a fixed full core word o, the extension count factors as

    e(A)e(B)L(pre(o))R(suf(o)).

Prefix increments U_u=L(u)-L(u with its last letter removed) and suffix increments V_v=R(v)-R(v with its first letter removed) are nonnegative, by appending and prepending the removed letter respectively.

Expand the aggregate deficit

    Delta=4c E0-3 sum_t N_t

in the U,V coordinates. A coefficient fixes an initial core word u and terminal core word v. If no extension has those endpoints, the coefficient is zero. Otherwise u is an ideal word, v is a filter word, and they are disjoint. Delete their vertices to obtain a bare core R with k common vertices. A slave can disappear, but a new slave cannot appear.

Each deleted common vertex lies outside both roots and contributes exactly one to its two comparison indicators. The coefficient is consequently

    (c-k)e(R) + [4k e(R)-3 sum_{t in Z(R)} N_R(t)].

This identity is exact for every c. Thus a bare-core aggregate bound closed under deleting either slave and any common vertices lifts to arbitrary exterior without an additional analytic hypothesis. Missing-slave cores must be covered explicitly; Section 3 does this.

## 3. One missing slave: an exact stronger formula

Suppose b is absent. The chain-threshold core consists of W1<...<Wc, x<a, W_L<a, and an isolated y, where 0<=L<=c. Let E be its extension count and N its aggregate comparison numerator.

Delete y temporarily. Any remaining extension is uniquely indexed by

    L<=i<=c, 0<=j<=i,

where i is the number of common vertices before a and j the number before x. Inserting y into its c+3 possible gaps gives

    E=(c+3) sum_{i=L}^c (i+1).

For the fixed (i,j) word, the rank of Wr before y is inserted is

    r + 1[j<r] + 1[i<r].

Hence its contribution, summed over the c+3 insertions of y, to the aggregate numerator is

    c(3c+7)/2 + i-(c+2)j.

Summing j from 0 through i gives

    N=(c/2) sum_{i=L}^c (i+1)(3c+7-i).

The elementary sums yield the exact nonnegative identity

    7cE-6N = c L(L+1)(c-L+1) >= 0.

In particular N/E<=7c/6, which is stronger than 4c/3. Equivalently,

    N/E = 7c/6 - c L(L+1) / [3(c+3)(c+L+2)],
    E = (c+3)(c-L+1)(c+L+2)/2.

This includes c=0, where the comparison sum is empty.

If a is absent and b<W_(K+1), reverse the order, exchange the roots, and set L'=c-K. This gives

    E = (c+3)(K+1)(2c-K+2)/2,
    N/E = 7c/6
          - c(c-K)(c-K+1)/[3(c+3)(2c-K+2)].

If both slaves are absent, the roots are exchangeable. For every common t,

    p(x<t)+p(t<y)=1,

so the aggregate is c exactly. This remains true for any induced common poset, even before conditioning it to a chain.

## 4. Independent audit of the shuffle-block polynomial derivation

Write j_v for the number of common-chain elements before a special letter v. For every full extension,

    sum_t (1[x<t]+1[t<y]) = c+j_y-j_x.

Its aggregate deficit contribution is therefore c+3(j_x-j_y). This identity holds whether x or y appears first.

Fix the order pi of the four special letters b,x,y,a. It must have b before y and x before a, and must have b before a when the optional relation is present. There are at most six such orders. Assign each special letter a common-chain gap j_v in {0,...,c}, weakly increasing along pi. This is a bijection with the original chain-shuffles. The attachment conditions are precisely j_b<=K and j_a>=L.

The general counting fact used in the derivation is elementary. A weakly increasing sequence of h gap indices in {0,...,m-1} has binomial(m+h-1,h) possibilities. For m>=1 its kth coordinate has mean k(m-1)/(h+1). To see the mean, add k to the kth coordinate, obtaining a uniform h-subset of {1,...,m+h-1}; its kth order statistic has mean k(m+h)/(h+1), and subtracting k gives the stated expression. Equivalently use h+1 symmetric stars-and-bars gaps.

The polynomial convention at m=0 is safe: the binomial polynomial is zero whenever h>0, and is one when h=0. There is no occupied empty block. Formal mean expressions in such zero-weight terms are finite and multiplied by zero; they neither create nor delete an actual extension.

### 4.1. The nonoverlap region K>=L

Set p=L, q=K-L, r=c-K, so p,q,r>=0. Split the chain gaps into

    0,...,p-1;
    p,...,p+q-1;
    p+q,...,p+q+r.

Let h1,h2,h3 be the counts of consecutive special letters in these blocks, with h1+h2+h3=4. The condition on a rejects precisely the assignments putting a in the first block.

If b is in the first or second block, j_b<=K is automatic. If b is in the third block, b and every preceding special letter in that block are forced to the initial gap K=p+q. Let f be their number and let h=h3-f. The remaining h third-block letters may still equal K; they are not required to be strictly later. Their allowed weak gap sequence therefore has binomial(r+h,h) choices. The forced prefix has no multiplicity. This is exactly the correct count, without overcounting letters sharing gap K.

Thus each assignment contributes the product

    binomial(p+h1-1,h1)
    binomial(q+h2-1,h2)
    binomial(r+h,h).

The unforced means are given by the general counting fact after adding each block's starting gap. The forced means equal K. Multiplying this product by c+3(mean(j_x)-mean(j_y)) and summing derives the exact deficit polynomial. The same construction applies with or without b<a, by restricting the special orders accordingly.

### 4.2. The overlap region K<L

Set p=K, q=L-K-1, r=c-L. Split the gaps into

    0,...,p-1;
    p,...,p+q;
    p+q+1,...,p+q+1+r.

The condition on a requires a in the last block. The condition on b excludes b from the last block. If b is in the middle block, b and earlier letters in that block are forced to its initial gap K=p. The remaining middle letters can occupy any weak sequence in its q+1 gaps. Thus, if f letters are forced, the weight is

    binomial(p+h1-1,h1)
    binomial(q+h2-f,h2-f)
    binomial(r+h3,h3).

Again the forced letters are counted once, unforced letters may share the forced gap, and all block sizes may reach their allowed boundary values. The relation b<a is automatic in this region. The corresponding sum of weighted wordwise deficits is exact.

### 4.3. Symbolic audit result

I independently inspected and reran the direct shuffle-block derivation. Its three exact expanded deficit polynomials have respectively 43, 53, and 44 positive monomial coefficients:

- K>=L without b<a: denominator 4, 43 positive terms;
- K>=L with b<a: denominator 60, 53 positive terms;
- K<L: denominator 4, 44 positive terms.

All chamber variables p,q,r are nonnegative integers. Hence each deficit is nonnegative for every c, not merely sampled c. The symbolic derivation comes directly from the finite shuffle sums above, rather than from interpolation. Its extension polynomials have degree at most four, and its deficit polynomials degree at most five, as also follows immediately from the four special letters and the linear wordwise deficit.

The derivation's checks against independent chain-state dynamic programming pass on all 125 triples 0<=p,q,r<=4 in each chamber. These checks confirm the implementation but are not used to infer the all-c identity.

## 5. Independent actual-exterior implementation checks

The separate standard-library verifier verify_two_slave_allc_exteriors_independent_20261019.py uses its own transitive closure and whole-poset subset dynamic programming. It imports no code from the main proof verifier. It generated 400 actual posets with c=4,5,6,7, arbitrary common-poset density, and arbitrary lower/upper exterior relations, through total order 17. It verifies the exact root neighborhoods and the aggregate inequality. All checks pass.

It also independently confirms the missing-slave formulas of Section 3 for all 45 threshold instances with c<=8. Results are recorded in two_slave_allc_exteriors_independent_audit_20261019.json. These checks do not enlarge the theorem's scope; the proof already covers every finite c and arbitrary exterior.

## 6. Final audit conclusion

The common-order conditioning, the exact shuffle-block deficit construction, the missing-slave formulas, and the tree-increment exterior lifting are all valid. Together they establish, for every finite c, the universal aggregate bound

    sum_{t in Z} [p(x<t)+p(t<y)] <= 4c/3

for a critical pair with at most one lower slave and at most one upper slave, with arbitrary finite exterior.

When c>0, this eliminates such a pair as a maximum-span choice in a hypothetical >2/3 majority order with no balanced pair. The argument does not eliminate the c=0 branch by itself, and it does not prove the larger slave-multiplicity cases.
