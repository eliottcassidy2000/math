# Asymptotically sharp all-order weighted stability

Status: PROVED, independently reviewed and accepted. Exact verification
receipts and reviewed-source pins accompany the continuation audit.

Let K_n be the infimum of Delta(A)/V(A) over weighted tournaments of order
n with V>0, where Delta=H-det(I+A)-|det A| and V=sum_e p_e(1-p_e).
Let r=floor(n/2) and define

    kappa_(2r)=(4/r)[1-(3/4)^r-(1/4)^r],
    kappa_(2r+1)=(4/r)[1-(3/4)^r].

Ordinal fair pairs (and a singleton for odd n) give K_n<=kappa_n.
For every n>=4, with c0=(1-exp(-4))/4, we prove

    kappa_n/[1+n kappa_n^2/(2c0)] <= K_n <= kappa_n.       (1)

In particular K_n=8/n+O(n^(-2)), and n K_n tends to 8. This does not
prove exact equality K_n=kappa_n at every finite order.

## 1. A multiplicative Pfaffian margin bound

Put f_+=u+v and f_-=u-v, where v=(-1)^n det A. Both are positive multiples
of sums of principal Pfaffian squares, respectively excluding or including
a fixed binary observer. Vary one old skew edge x=2p-1. Each relevant
Pfaffian containing both endpoints has the form x a_S+b_S, while those
excluding both endpoints include exactly the coefficient Pfaffians a_S.
Consequently f(x)=a x^2+b x+c, and, when a>0,

    c-b^2/(4a) >= a.

Indeed the squared norms of the coefficient Pfaffians contribute a to the
constant term through the subsets excluding both endpoints. Completing
the square in the remaining sum of squares leaves a nonnegative residual.
Additional subsets containing exactly one endpoint only increase the
constant term. If a=0, the sum-of-squares form implies b=0, and f is constant.

Write f=a[(x+g)^2+d], d>=1. The endpoint interpolation L satisfies

    L=f+a(1-x^2)=f+4qa,
    L/a=1+2gx+g^2+d >= 2-x^2=1+4q.

In particular f(x)>=a. Since L=f(x)+4qa, this gives f(x)>=L/(1+4q), for both signs. Round all edges independently, iterate
this inequality, and use f_+(T),f_-(T)>=1 for nonempty binary tournaments:

    u-|v| >= product_e (1+4q_e)^(-1) >= exp(-4V).          (2)

This is a bound for a nonempty tournament. Empty pair-deletion remainders
are not included in its use below.

## 2. A dimension-free nonlinear gap bound

For n>=3, the established one-edge margin interpolation gives a correction
q_e(u_R-|v_R|), where R deletes both endpoints. After previously rounded
edges are binary, the remaining variance of R is at most the total V_t
of all later unrounded edges. By (2), its correction is at least
q_e exp(-4V_t). Summing and retaining the terminal binary expectation gives

    Delta(A) >= E Delta(T) + (1-exp(-4V))/4.               (3)

To justify the last step, order the intervals [V_t,V_t+q_e] backwards.
They partition [0,V], and the decreasing exponential has integral at most
q_e exp(-4V_t) on each interval. No limiting partition assumption is used.

## 3. Weighted matching extraction

For nonnegative edge weights q, let W be the sum q_e q_f over all incident
unordered edge pairs. Greedily choose a maximum-weight edge, discard its
incident remaining edges, and repeat. Let M be the total selected weight.
For a selected edge uv of weight w, assign all edges discarded at that
step to it. Their weights split into two stars with sums A,B, total S.
Every assigned weight is at most w. Within this edge-group, wedge weight
is at least

    wS + (A^2+B^2-sum_assigned q_e^2)/2
        >= S^2/4 + wS/2 >= S^2/4.

The groups are edge-disjoint, hence these counted wedge pairs do not
repeat. There are at most n/2 groups. Cauchy-Schwarz gives

    W >= (V-M)^2/(2n),
    M >= V-sqrt(2nW).                                  (4)

The local weighted triangle inequality gives tau>=W. For n>=4 the binary
proper-cycle bound transfers under rounding and yields Delta>=tau>=W.

## 4. Retain the matching after rounding its complement

Round all edges outside the greedy matching first, discarding their
nonnegative corrections. On each conditioning branch, the matching-support
stability theorem applies, regardless of the binary cross orientations.
Its coefficient for the number of retained pairs is at least kappa_n.
The matching edge variances have not changed. Therefore

    Delta(A) >= E Delta(T)+kappa_n M
             >= kappa_n[V-sqrt(2nW)].                  (5)

The matching-support coefficient monotonicity follows directly from its
geometric-sum formulas: a_r decreases with r, and the full-even coefficient
b_r is below a_r. No assumption that the intermediate tournament has H=1
is being made.

## 5. Retaining the binary baseline improves the rate

Write B=E Delta(T). The binary proper-cycle theorem implies B>=tau>=W
for n>=4. Retain B in (5), rather than discarding it:

    Delta >= kappa_n V+B-kappa_n sqrt(2nW)
          >= kappa_n V+B-kappa_n sqrt(2nB)
          >= kappa_n V-n kappa_n^2/2.                  (6)

The last step completes the square in sqrt(B). This is an additive
finite-n inequality, with no assumption that A is close to extremal.

## 6. Uniform asymptotic conclusion

Fix A and write R=Delta/V and L=kappa_n/[1+n kappa_n^2/(2c0)].
If V<=1, (3) gives R>=c0. The maximum of x/(1+n x^2/(2c0)) over x>=0 is
sqrt(c0/(2n)), which is at most c0 since c0>1/6 and n>=4. Thus R>=L.
If V>1, (3) gives R>=c0/V and hence V>=c0/R. Equation (6) yields

    R >= kappa_n-n kappa_n^2/(2V)
      >= kappa_n-[n kappa_n^2/(2c0)]R.

Rearranging gives R>=L, proving (1) for every n>=4 without a threshold
condition. Since kappa_n=8/n+O(n^(-2)), the error term in (1) is O(n^(-2)).
This positive rational lower bound is not advertised as sharp at modest n.
Exact orders already proved and other bounds remain stronger in that range.
Orders one and two are excluded: r=0 makes the former formula undefined,
while the latter has Delta=0 despite possible V>0.

## 7. Near-candidate structure

If R<=kappa_n, retaining the baseline before completing the square gives

    B <= kappa_n sqrt(2nB),
    B <= 2n kappa_n^2,  W<=B,
    V-M <= sqrt(2nW) <= 2n kappa_n.

Thus in any such sequence, expected binary gap and triangle weight are
O(1/n), and the greedy matching loses only O(1) total variance. A binary
nontransitive tournament of order n>=4 has at least one proper directed
triangle, so Delta(T)>=1. Consequently

    Pr(T is nontransitive) <= B <= 2n kappa_n^2.

These are probabilistic and weighted-mass statements. They do not assert
that the original cross edges are all binary, that every pair is fair, or
that every supported completion is transitive at finite n. Exact global
finite-order equality K_n=kappa_n remains open beyond the proved orders.
