# A square-root inverse for the full binary gadget response

Session 197. This is an exact arithmetic algorithm for fixed coefficients. It does not select a physical source, prove cofinite coverage, or treat arbitrary integer states as physical.

## 1. Input and scope

The reviewed ROOTED.md indexes every ending-one, transitive-deletion gadget by L>=1. Its constants g(L), Q(L), Z(L) have initial values g(1)=Q(1)=Z(1)=1 and recurrences

    (g,Q,Z)(2l) = (g(l),Q(l),2Z(l)),
    (g,Q,Z)(2l-1) = (g(l)+l-(l mod2),Q(l)+Z(l),2Z(l))  (l>=2).

Here g is the Hamilton-path count, Q is the root-start two-cover response, and Z=2^bit_length(L-1). Their physical interpretation is inherited. The earlier SQRT_FIXED_PREFIX_INVERSE.md proves a block identity for g and an inverse for h g(L)+R L-e(L mod2). This extension retains both Q and Z.

**Theorem.** Fix integers

    A>=1, B>=A, 0<=C<=4A, D>=0, 0<=E<=A.

For a positive integer target T, all solutions L>=1 of

    F(L)=A g(L)+B L+C Q(L)+D Z(L)-E(L mod2)=T          (1)

can be found in O(sqrt(K) log(2+K)) exact arithmetic operations and O(sqrt(K)) integer words, where K=floor(T/B), with constant work when K<=1. The algorithm returns the complete solution set, including L=1 and both parities.

The inequalities are sufficient arithmetic hypotheses, not a claim that every five-tuple has a graph realization. For example (1) applies to a nonnegative linear readout of the physical one-seam output when its coefficients satisfy these bounds. An additional constant term can be moved into T, and a singleton-only correction can be tested separately. No unreviewed two-source graft theorem is required here.

## 2. Exact bit identities

For a width j>=0 and 0<=x<2^j, write revc_j(x) for the integer obtained by complementing its j binary digits and reversing them. Set revc_0(0)=0. With k=bit_length(L-1),

    Q(L)=1+revc_k(L-1),  Z(L)=2^k.                  (2)

The width k is the minimal width of L-1. Arbitrary leading padding is not allowed in (2). In terms of n=L-1, appending a binary1 leaves Q unchanged and appending a binary0 to a nonzero word adds the old 2^k. The leading-zero operation on the empty word is not used. These are exactly the recurrences above. Starting at n=0 proves (2). Splitting a fixed-width word into high and low parts gives

    revc_(p+q)(2^q x+b)=2^p revc_q(b)+revc_p(x).     (3)

The inherited block identity, valid for q>=1, x>=0 and 0<=b<2^q, is

    g(2^q x+b+1)=g(x+1)+u x+c(x mod2)+g(b+1)-1,     (4)
    u=revc_q(b), c=u mod2.

For completeness, if f(n)=g(n+1), then f(2n+1)=f(n) and f(2n)=f(n)+n+(n mod2), also at n=0. Append the q low digits of b. Their zero positions give coefficient u on x; only the first appended digit can contribute parity of x, with coefficient c. Setting x=0 fixes the constant. This proves (4).

Induction gives 1<=g(L)<=L. Odd L>=3 has g(L)>=3. Consequently F(L)>=BL and every solution satisfies L<=K. Moreover, for every L>=2,

    F(L)>=A+2B+C+2D=F(2).                           (5)

For even L, A g-E parity>=A. For odd L>=3 it is at least 3A-E>=2A. The other terms use L>=2, Q>=1 and Z>=2. Thus after testing L=1, targets below F(2) can be rejected immediately.

## 3. Dyadic-shell monotonicity

Partition L>=2 into shells with k=bit_length(L-1), so

    2^(k-1)+1 <= L <= min(2^k,K).

The shells k<=4 contain only finitely many L and are tested directly. For k>=5 choose

    q=ceil((k+3)/2), p=k-q, B0=2^q, H=2^p.

Then p>=1 and B0>=8H. Every L in the shell uniquely has L=B0 x+b+1, 0<=b<B0. Its exact x interval is

    max(0,ceil((2^(k-1)-b)/B0)) <= x
       <= floor((min(2^k,K)-1-b)/B0).                (6)

In particular 0<=x<H. Put u=revc_q(b), c=u mod2 and r(x)=revc_p(x). Equations (2)-(4) turn (1) into

    F_b(x)=A g(x+1)+(Au+BB0)x+Ac(x mod2)+C r(x)
       +A(g(b+1)-1)+B(b+1)+C(1+Hu)+D2^k-E((b+1) mod2). (7)

**Lemma.** Every nonempty branch (7) is strictly increasing on its integer interval (6).

Indeed g(x+2)-g(x+1)>=-x, the parity increment is at least -1, and r(x+1)-r(x)>=-(H-1). Since u>=c,

    F_b(x+1)-F_b(x)
      >= BB0-Ax+A(u-c)-C(H-1)
      >= A[B0-x-4(H-1)] > 0.                       (8)

The last inequality uses B>=A, C<=4A, B0>=8H and x<H. Neither g nor Q is assumed monotone. Widening the low-word block is precisely what absorbs the signed changes of the bit-reversal response.

## 4. Algorithm, completeness and cost

Test L=1. Apply (5). For every remaining shell, test the finitely many small-shell values, or do the following for k>=5:

1. Tabulate g(1),...,g(B0) by its recurrence, and r(x) for 0<=x<H.
2. For every b=0,...,B0-1, compute u by its q binary steps and the interval (6). Skip an empty interval.
3. Binary-search that interval for the first F_b(x)>=T. Test equality and, if it holds, return L=B0x+b+1.

Every candidate L<=K is in exactly one shell and exactly one branch. Inequality (8) proves the validity of binary search and permits at most one answer per branch. Equality testing discards all false candidates. The cutoff L<=K excludes every larger L. These observations prove exactness at all orders.

For a shell, B0 is at most a constant times 2^(k/2), H<=B0, and (7) uses O(1) arithmetic operations after tabulation. Computing the reversal tables by binary steps and all branch searches costs O(B0 log(2+K)) arithmetic operations. Summing the geometric shell sizes gives O(sqrt(K) log(2+K)); tables can be discarded between shells. Sorting the O(sqrt(K)) answers has the same bound. Storage is O(sqrt(K)) integer words.

Once (5) permits a non-singleton search, all five coefficients are O(T), and all search intermediates have O(log(2+T)) bits. The early tests also have to read the original input coefficients, which may be longer than T. A bit-cost bound therefore accounts for input reading and the cost of exact arithmetic on these lengths. The theorem is an arithmetic-operation bound, not polynomial time in the input bit length.

## 5. Reproducible controls and boundaries

The standalone verifier uses no third-party package. It separately counts g and Q on actual gadget adjacency matrices, tests both bit identities, verifies strict branch inequalities, and compares the complete inverse with exhaustive L scans across coefficient boundaries. These finite controls supplement, not replace, the proof above. The verifier explicitly rejects invalid coefficient hypotheses and tests zero targets, C=0, C=4A, E=0/A, dyadic boundaries, and L=1.

Three larger arithmetic controls have exactly the indicated single solution:

- (A,B,C,D,E)=(9,9,30,30,3), T=119881395: L=1234567.
- (39,54,156,216,6), T=5703260151: L=9876543.
- (100,300,400,1234,100), T=235297253052: L=100000007.

These are arithmetic examples. No graph provenance is asserted for these particular coefficient vectors.

Named dependencies: ROOTED.md SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e, for the physical gadget meaning; SQRT_FIXED_PREFIX_INVERSE.md SHA256 084b7c372524f41e99a4087179e0fd30baf35a0de58a0166e3bdcd84c5e407bf, for the earlier scalar inverse. The present proof restates all arithmetic identities it needs. The new content is the full-response shell decomposition, its monotonicity despite Q, and the complete inverse algorithm.
