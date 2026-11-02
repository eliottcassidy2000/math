# Independent audit: asymptotically sharp all-order stability

**Verdict: ACCEPT for the final n>=4 theorem.** This is a separate
continuation, not part of the frozen first Session138 checkpoint. The proved
conclusion is

    kappa_n/[1+n kappa_n^2/(2c0)] <= K_n <= kappa_n,
    c0=(1-exp(-4))/4, n>=4,

and consequently `K_n=8/n+O(n^(-2))`. Exact equality K_n=kappa_n at arbitrary
finite order remains unproved. No fairness or transitivity claim about all
near-extremizers is established below.

## 1. Pfaffian normalization and one-edge robustness

The key coefficient comparison is exact for each of the two signed
Pfaffian-square sums. Fix an old skew entry x=Kij. Every allowed principal
set containing i,j gives a Pfaffian `alpha_S x+beta_S`; its coefficient is,
up to sign, the Pfaffian after deleting i,j. The map deleting those endpoints
is a bijection onto the allowed principal sets excluding both endpoints.
For the plus margin these are the even old-vertex sets. For the minus margin
they are the even augmented sets containing the observer; removing two old
vertices preserves that condition. The normalization2^(1-r) is unchanged
within the original sum. Hence, in `f(x)=a x^2+b x+c`, the subsets excluding
both endpoints contribute exactly a to c, not a different power-of-two factor.

The sets containing just one endpoint contribute additional nonnegative
constants. Cauchy on the intercept/coefficient vectors gives

    4a(c-a)>=b^2, or c-b^2/(4a)>=a when a>0.

If a0, every coefficient Pfaffian vanishes and b0. In determinant coordinates
an independent normalization check is `a=(u_R+/-v_R)/4`; this remains correct
at an empty R if its formal coordinates u_R=v_R=1 are retained. One must not
mistakenly replace a by the whole pair-deleted margin.

Thus f>=a. Since endpoint interpolation L equals f+4q a, it follows directly
that `f>=L/(1+4q)` for both signs. Iterating independent edge rounding leaves
the fixed product of factors `(1+4q)^-1`. Nonempty binary terminal matrices
have both signed margins>=1, by the same odd-integer Pfaffian parity counts.
Therefore

    u-|v| >= product_e(1+4q_e)^-1 >= exp(-4V)

for every nonempty weighted tournament. No expectation/minimum interchange
is used incorrectly: each signed margin is bounded first, then their minimum.

Boundaries are essential and checked: order1 gives both margins1; a fair
order2 edge gives margins1/2 and1, with the second quadratic coefficient0.
The empty matrix has margins2 and0, so it is deliberately excluded from the
product lower bound. In particular the later gap argument cannot use an
empty pair-deletion remainder at n2.

## 2. Integrated gap, without a limiting argument

At n>=3, one-edge interpolation of Delta has correction
`q_e(u_R-|v_R|)`. The pair-deleted remainder is nonempty, and after previously
rounded edges are binary its variance is at most the total V_t of later
unrounded edges. Thus the correction is at least `q_e exp(-4V_t)` on every
conditioning branch. All q_e remain their original values until rounded.
The intervals `[V_t,V_t+q_e]` partition `[0,V]` in reverse order. Their
left-endpoint exponential rectangles dominate the integral, giving

    Delta>=E Delta(T)+(1-exp(-4V))/4.

This is stronger than the earlier `V exp(-4V)` draft bound. It is valid at
n3, but not at n2: Delta is identically0 there with possibly positive V.

## 3. Greedy weighted matching lemma

Each discarded edge is assigned exactly once to the selected maximum-weight
edge that removes it. In a group centered at uv of weight w, split assigned
weights by endpoint into sums A,B and total S. No assigned edge is incident
to both u and v. Its weight is at most w. Counting selected/assigned pairs
and pairs within each endpoint star gives

    W_group >= wS+(A^2+B^2-sum q_e^2)/2
            >= wS/2+S^2/4 >= S^2/4.

Possible wedges between the two stars were merely omitted, not subtracted.
Different groups have disjoint edge sets, so their counted unordered wedge
pairs do not repeat even if their vertices overlap. With at most n/2 groups,
Cauchy gives `W>=(V-M)^2/(2n)`. The zero-weight case is immediate and may be
handled without dividing by a zero number of groups. The lemma holds for
arbitrary nonnegative weights, not only variances bounded by1/4.

## 4. Matching-support transfer and the uniform bound

After rounding the complement of this fixed greedy matching, the surviving
variance is M on every branch. Arbitrary binary cross orientations are
allowed by the already-audited matching-support theorem. Its coefficient is
at least kappa_n: a_c decreases with the number of pairs; a perfect even
matching has b_c, and a smaller matching has coefficient a_c above the
full-dimensional b. Therefore

    Delta>=E Delta(T)+kappa_n M
         >=kappa_n(V-sqrt(2nW)).

At n>=4 the binary proper-cycle bound and nonnegative edge-rounding
corrections give Delta>=tau>=W. This is not asserted at n3 via tau.

Retaining the binary expectation sharpens the rate. Set B=E Delta(T).
The binary proper-cycle theorem gives B>=E tau(T)=tau(A)>=W for n>=4.
Consequently

    Delta>=kappa_n V+B-kappa_n sqrt(2nB)
         =kappa_n V+(sqrt(B)-kappa_n sqrt(n/2))^2-n kappa_n^2/2
         >=kappa_n V-n kappa_n^2/2.

No near-optimality assumption is used in this additive inequality.

For R=Delta/V, the nonlinear gap bound gives R>=c0 if V<=1. The maximum
of x/(1+n x^2/(2c0)) for x>=0 is sqrt(c0/(2n)), which is <=c0 for n>=4
because c0>1/6. This covers the small-variance case. For V>1, R>=c0/V,
so 1/V<=R/c0; R is positive by the nonlinear bound. Therefore

    R>=kappa_n-n kappa_n^2/(2V)
     >=kappa_n-[n kappa_n^2/(2c0)]R.

Rearranging gives the claimed positive rational lower bound uniformly over
all matrices with V>0. No extremizer existence, infimum attainment,
near-extremizer subsequence or limiting exchange is assumed.

Correction lineage is retained: the initial small-kappa domain condition
would also have admitted n2 despite its invalid positive-gap proof branch;
the final explicit n>=4 theorem repairs this. Dropping the baseline had
first yielded a subtractive lower bound, then a positive rational bound with
relative error O(n^-1/2). Retaining B and completing its square gives the
strictly stronger displayed denominator and relative error O(n^-1).
Since kappa_n=8/n+O(n^-2), the final error is O(n^-2), not merely O(n^-3/2).
Exact finite-order equality beyond the separately proved orders stays OPEN.

## 5. Independent exact execution controls

`check.py` is self-contained and uses Fraction Gaussian determinants,
Hamilton subset DP, and literal independent binary expectations. It imports
no producer kernel. All checks raise RuntimeError and remain active under -O.
It verifies:

*796 matrices: all ternary-edge matrices through order4, and8 fixed-seed
 six-fractional-edge rational examples plus one sharp fair-matching matrix at
 each order5..8;
*4566 one-edge quadratic checks, including exact pair-deletion normalization,
 completed-square inequality, affine/constant boundary, endpoint identity
 and multiplicative comparison;
*796 exact full-margin product bounds;
*1584 exact rational telescoping gap inequalities, using both ascending and
 descending variance orders and the stronger discrete corrections
 `sum q_e product_(later f)(1+4q_f)^-1`;
*765 exact baseline/triangle/wedge/matching/square-completion checks, including
 10 near-candidate checks (including the explicit sharp controls);
*384 arbitrary nonnegative rational weighted graphs,16 at each n2..25,
 with2467 greedy groups, checking the local and global squared inequalities;
*explicit r1,r2,empty and prohibited n2-gap controls.

The exponential integration and asymptotic conclusion are proved above;
no floating-point exponential comparisons are substituted for those proofs.
Normal and optimized receipts agree, with semantic digest
`5af0e882af37de26b39217b660b93e72b7692fe59cba189b3f8f3f97f133a7cf`.

## 6. Independent resolvent mechanism

`../../stability/RESOLVENT_ROBUSTNESS_CROSSCHECK.md` also passes proof review.
For real skew K, I+K is normal and its inverse norm is at most1. Both
`||(I+K)^-1 1||^2` and `||(I+K)^-T 1||^2` equal
`s=1^T(I+K)^-1 1`. Thus one skew-coordinate derivative satisfies
`|log D|'<=2` and `|log s|'<=1`, where D=det(I+K). Both signed margins have
logarithmic Lipschitz constant at most3. Nearest-binary coordinate distance
is2d with d<=2V, giving the independent bound `u-|v|>=exp(-12V)`.
The identical finite-interval integration gives
`Delta>=E Delta(T)+(1-exp(-12V))/12` for n>=3. Although weaker, this suffices
for the same asymptotic error rate through matching extraction. The result
therefore has a second mathematical robustness mechanism, not only a second
implementation of the Pfaffian argument.

## 7. Scoped near-extremal structure

For any fixed C>0 and a sequence with R<=C/n, the nonlinear bound implies
V>=c0 n/C once n is sufficiently large. Also
`W>=2V^2/n-V/4` and W<=RV give `V<=n/8+C/2`. Thus V=Theta_C(n), W=O_C(1),
and greedy matching loses only O_C(sqrt n) variance mass.

If R=(8+o(1))/n, the retained matching has r=Omega(n), since its weight is
Omega(n) and each variance is at most1/4. Its matching-support coefficient
is `(4/r)(1-o(1))`; combining with retained mass V-O(sqrt n) gives
`r>=(1/2-o(1))n`. This proves an almost-perfect matching of fractional mass
support, not that its edge weights are asymptotically fair or that binary
cross orientations are transitive. Those stronger assertions would need
additional arguments.

For the stronger near-candidate condition R<=kappa_n, preserving B before
square completion gives B<=kappa_n sqrt(2nB), hence

    B<=2n kappa_n^2=O(1/n),
    W<=B,  V-M<=sqrt(2nW)<=2n kappa_n=O(1).

Every nontransitive binary tournament of order>=4 has a proper directed
triangle, so its integer gap is at least1. Thus
`Pr(T nontransitive)<=B<=2n kappa_n^2`. This is a probability bound, not a
claim that every supported completion is transitive at finite n. The earlier
general R<=C/n estimates remain valid but are weaker for this narrower class.
