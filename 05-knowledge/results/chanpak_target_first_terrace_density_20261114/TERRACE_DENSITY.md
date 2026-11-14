# Exact coverage versus residue coverage

## 1. Seams and connected factors

For the normalized first triangle terrace, the gap word begins and ends
with a positive gap. The order-three cycle is included by the word
`(1)`. Use the interval packing from Session 337: select any subset of
gaps. A maximal selected run with a positive gap has weight equal to the
product of (2^r) for its positive gaps. A pure-zero selected run has
weight two if its length is even, zero otherwise. Sum the weights of all
subsets to get (H).

Fix the selection bits for the positive gaps (r_i). Between two
successive positives let (q_i) be the number of zero gaps, and write
(K_q(a,b)) for the total weight of subwords there when the bounding
positive selections are (a,b). A run touching a selected endpoint
has weight one on its zero sites; otherwise use the pure-zero rule.
All zero corridors are independent after the positive choices.

Write (A_q=K_q(0,0)), (B_q=K_q(0,1)=K_q(1,0)) and
(D_q=K_q(1,1)). With (Delta=1-z-z^2-z^3) their respective
generating functions are
[
 (1+z^2)/Delta,quad (1+z)/Delta,quad (1+z+z^2)/Delta .
]
Indeed a pure-zero word consists of even selected runs of weight two
separated by nonselected sites. With run series (W=2z^2/(1-z^2))
the word series is ((1+W)/(1-z(1+W))). Allowing one or both
selected endpoints gives ((1+zA)/(1-z)), then
((1+zB)/(1-z)): count the all-selected word once, or split at the
first unselected site from the new selected end.

When (q=0) the kernel is all ones, so it separates a product.
Cut between every consecutive pair of positive gaps with no zero
in between. The resulting components have (q_ige1) at all
internal links and each contributes one factor (H). Conversely
every finite nonempty product of connected component values is realized:
concatenate the component gap words with no intervening zeros.
Order is immaterial to (H), even if the tournament isomorphism class
changes. Thus the terrace spectrum is **exactly** the nonempty
multiplicative closure of the connected component values. These are
structural factors, not necessarily arithmetically irreducible.

This gives a target-first alternative to looking through orders:
factor the prescribed integer, test candidate connected values only
among its divisors, and use multiplicative dynamic programming. For
a connected scan, the lower bound below gives finite budgets in both
positive widths and zero corridor lengths, so candidate words may
be pruned directly at the target.

## 2. A convergent Dirichlet bound

Let (w_i=2^{r_i}) and choose each bit independently, with
(Pr(b_i=1)=w_i/(1+w_i)ge2/3). The connected count is
[
 H=prod_i(1+w_i) mathbb Eprod_i K_{q_i}(b_i,b_{i+1}).
]
Each factor inside the expectation is nonnegative and increasing
in the bits: (B-A=z(1-z)/Delta) and (D-B=z^2/Delta)
have nonnegative coefficients, using the nondecreasing
coefficients of (1/Delta).
Nonnegative increasing functions of independent bits satisfy
(mathbb Eprod f_igeprodmathbb E f_i). For two functions
prove this by induction, conditioning on the last bit and using
the one-bit covariance (p(1-p)(f_1-f_0)(g_1-g_0)ge0);
then iterate since a nonnegative product is again increasing.

Thus (Hgeprod_i(1+2^{r_i})prod_i E_{q_i}), where
(E_q=(A_q+4B_q+4D_q)/9) is the kernel expectation at probability
(2/3). Its series is ((9+8z+5z^2)/(9Delta)).
For (q=1,ldots,5), (9E_q=17,31,57,105,193).
Moreover (E_qge(9/5)E_{q-1}) for (qge1):
check the first three inequalities, and thereafter use the
tribonacci recurrence on the consecutive differences.

For (s=49/50) set
[
 R=sum_{rge1}(1+2^r)^{-s},qquad Q=sum_{qge1}E_q^{-s}.
]
Bounding the (rge9) tail by powers of (2^{-s}), and (qge6)
by powers of ((9/5)^{-s}) starting with the fifth value, gives
(R<793/1000, Q<1202/1000, RQ<954/1000).
The replay verifies these with integer power comparisons and rational
ceilings, not floating-point evidence. Even summing over all ordered
connected words, with repetitions counted, gives
(sum H^{-s}lesum_{kge1}R^kQ^{k-1}<infty).

Let (S) be the set of distinct connected counts. Every member is at
least 3; hence
(prod_{min S}(1-m^{-s})^{-1}<infty),
since (-log(1-x)le x/(1-3^{-s})) for (xle3^{-s}).
The product counts multisets, allowing multiple representations of
the same integer. The sum of (H^{-s}) over distinct terrace counts
is therefore finite. Each term for (Hle X) is at least (X^{-s}),
which proves the announced bound. The multiplicative collapse at
zero-length seams is essential: counting ordered arbitrary words
instead can count the same values exponentially many times.

## 3. Full dyadic residue coverage is compatible with the bound

Positive-only gaps realize (prod_i(1+2^{r_i})). To realize an
odd (u) modulo (2^k), use the factor 3 if needed to match modulo 4.
For each (r=2,ldots,k-1), multiply by (1+2^r) if the next bit
disagrees. Multiplying an odd number by (1+2^r) toggles bit (r)
without changing lower bits. If the list is empty use (1+2^k);
the (k=1) claim is immediate from any witness. The sum of the
positive gaps is (O(k^2)). If arbitrarily large witnesses are
desired, append a gap (rge k). This supplies an explicit integer
approximation at every dyadic precision, but the witnesses are not
bounded and need not ever equal (u). Thus a target-decreasing
tournament recursion cannot replace control of its absolute value
with even perfect dyadic lifting.

## 4. The 7/21/63 test

Single positives yield (1+2^r). With at least three positive
gaps in one component the bound is at least
(3^3(17/9)^2>96), and at least four certainly gives a larger
bound. For two positives, (qge4) gives at least (9E_4=105).
At (q=1) the count is (2(1+w)(1+v)-1); at (q=2,3) it is
(3+3(w+v)+4wv) and (5+6(w+v)+7wv) respectively, where
(w,v) are powers of two at least 2. None is 7 or 21:
only the first formula could yield 21, which would require a
product of two odd factors at least three to be 11.
For 63 the first would require such a product to be 32;
in the second either both (w,vge4) (too large) or
(63=9+11v) (impossible); in the third the minimum is 57
and its next value is 97. Products cannot yield 7,21,63 unless
one connected count is one of those numbers. None is a single
positive-gap value either.

For comparison take the ascending six-vertex spine, hubs (x	o y),
and word `100110`, whose upper and lower bits record arcs from
(x,y) respectively to the corresponding spine vertex. The literal
directed adjacency subset DP in the replay gives (H=63).
The third number in the multiplication-by-three chain is not supplied
even by the *whole* first triangle terrace.

## 5. Above the terrace

For strong (T) of order (nge5) sort the scores (s_i) descending.
Let (s_i^*) be the decreasing terrace scores, and (e_k) the
number of backward arcs across cut (k) minus one. Then
[
 e_k=sum_{ile k}(s_i^*-s_i)ge0 .
]
With (d=C_3(T)-(n-2)), summation by parts using the triangle
score formula gives
[
 2d=sum_{k=1}^{n-1} e_k
   igl[(s^*_k-s^*_{k+1})+(s_k-s_{k+1})igr].
]
The terrace gaps equal one at all interior cuts, zero at the ends.
Hence the sum of the interior (e_k) is at most (2d).
The largest two actual scores give (e_2ge2e_1), and the
smallest two give (e_{n-2}ge2e_{n-1}). For example
(s_1=n-2-e_1) and (s_2le s_1) imply the first inequality.
As the two interior cuts are distinct, the total (sum e_kle3d).

Reverse to forward all backward arcs crossing any cut with (e_k>0).
Their number is at most (sum_{k:e_k>0}(1+e_k)le6d).
The remaining reverse-arc intervals share no cuts. Group intervals
that meet at a common endpoint into maximal chains. Each occupied
interval is now exactly the backward-chain tournament of the terrace
classification (or a transitive pair if it has two vertices), and all
arcs between intervals go forward. Unoccupied vertices are singletons.
This proves the claimed bounded edit decomposition; it does **not**
bound the spectrum after even one such edit. A changed long arc can
probe joint endpoint data in many formerly independent components.
The next precise question is whether that probe destroys the
multiplicative collapse quickly enough to give an exact-target return.
