# A square-root exact inverse for one final low-cycle gadget

Session 197. Status: all-order algorithm for a fixed prefix; no scalar-coverage or physical inverse-closure claim.

## 1. Inherited physical input and new question

The reviewed `ROOTED.md` one-seam law and `ARITHMETIC_INVERSE.md` use a rooted prefix with response `(h,e,R)`: all Hamilton paths, paths ending at the root, and ordered two-path covers whose second path starts at the root (the first path may be empty). Every such nonempty prefix satisfies

    h >= 1,  0 <= e <= h,  R >= h.

Every integer L>=1 indexes a physical gadget, and the output count is

    F(L) = h g(L) + R L - e (L mod 2).                 (1)

Here g(1)=1, g(2l)=g(l), and g(2l-1)=g(l)+l-(l mod2) for l>=2. The earlier exact inverse scans a congruence-restricted interval. This addendum provides a different, complete algorithm with O(sqrt(K) log K) arithmetic operations for K=floor(N/R). Neither algorithm is asserted to dominate the other for every prefix.

The new arithmetic theorem needs only the three displayed inequalities. A triple satisfying them need not be physically reachable. When the prefix is physical, every returned L gives its actual one-seam extension. It does not choose a suitable prefix for an arbitrary target.

## 2. A binary-block identity

The elementary bounds 1<=g(n)<=n follow by induction. For odd n=2l-1, if l is odd use g(l)<=l and subtract1; if l is even use g(l)=g(l/2)<=l/2<=l-1. In particular F(L)>=RL.

Let k>=1, B=2^k, 0<=b<B, and x>=0. Write b as exactly k binary digits, including leading zeros. Let q be the integer obtained by complementing those digits and reversing their order, and let c=q mod2. Then

    g(Bx+b+1) = g(x+1) + qx + c(x mod2) + g(b+1)-1.   (2)

To prove this, put f(alpha)=g(alpha+1). Appending a binary1 gives f(2alpha+1)=f(alpha). Appending a binary0 gives

    f(2alpha)=f(alpha)+alpha+(alpha mod2),             (3)

including alpha=0. Append the k digits of b to x, from most significant to least significant. A zero in position j contributes 2^(j-1)x to the coefficient of x. Only the first appended digit can contribute a parity term involving x: after that digit, the current parity is determined by the appended low digits, not by x. More precisely, a parity term involving x occurs exactly when the first appended digit is zero. Thus the total linear coefficient is q and the parity coefficient is c. The remaining constant is obtained by setting x=0, giving f(b)-f(0)=g(b+1)-1. This proves (2), also when b=0 or x=0.

## 3. Monotone low-word branches

For a positive target N, every solution has L<=K=floor(N/R). If K=0 return no solution. If K=1 test L=1 directly. For K>=2 choose

    k = ceil(bit_length(K-1)/2),  B=2^k.

Then K<=B^2, B<2sqrt(K), and B<=K. Uniquely write each L in [1,K] as L=Bx+b+1 with 0<=b<B and

    0 <= x <= floor((K-1-b)/B) <= B-1.

Since B is even, the parity of L is the parity of b+1. By (2), its count is

    F_b(x) = h g(x+1) + (hq+RB)x + hc(x mod2)
             + h(g(b+1)-1) + R(b+1) - e((b+1) mod2). (4)

**Lemma.** Each nonempty branch F_b is strictly increasing on its allowed integer interval.

Indeed g(x+2)-g(x+1)>=1-(x+1)=-x. The parity increment in (4) is at least -hc. Therefore

    F_b(x+1)-F_b(x) >= RB-hx+h(q-c)
                       >= h(B-x+q-c) > 0.          (5)

Here R>=h, q>=c, and a consecutive pair in the interval has x<=B-2. This proof does not presume that g itself is monotone; it is not.

## 4. Exact algorithm and complexity

Tabulate g(1),...,g(B) using the recurrence. For each b=0,...,B-1, binary-search the allowed x interval for the first F_b(x)>=N. Test equality, and if it holds return L=Bx+b+1. Sort the resulting list if desired.

Every admissible L occurs in exactly one branch. The strict inequality (5) makes binary search exact and gives at most one solution per branch. Equality testing removes every false candidate. Consequently the algorithm returns **all** L>=1 satisfying (1), not just a witness when one happens to be found. The L=1 boundary is included, and no parity class is dropped.

The table costs O(B) arithmetic operations and storage. Computing q by k bit steps and searching a branch costs O(log K) arithmetic operations; there are B<2sqrt(K) branches. Sorting at most B results has the same O(B log B) comparison bound. Hence the total is O(sqrt(K) log K) arithmetic operations and O(sqrt(K)) integer words of storage. The input integers and all intermediate values have O(log N) bits when K>=1. Thus a bit-complexity bound may be obtained by multiplying by the cost of arithmetic on O(log N)-bit integers; this is **not** polynomial time in the input bit length. The exact-operation bound is the claim here, not a unit-cost claim for arbitrarily large integers.

The algorithm can be intersected with the earlier congruence filter, but no further complexity gain is needed for this theorem. In particular, the method improves a general linear scan in K without resolving the separate problem of selecting a physically reachable prefix.

## 5. Independent finite controls supplied here

`verify_sqrt_fixed_prefix.py` is self-contained apart from checking the two named proof hashes. It implements the new block search, separately counts g on literal one-apex graphs, compares against exhaustive final-L scans, checks strict branch monotonicity, and reconstructs small physical one-seam outputs. All checks use explicit exceptions, so normal and optimized Python runs perform the same tests.

The three larger exact inverse controls are:

- Prefix(3,1,3), target486878921: only L=100000007.
- Prefix(5,1,7), target1383336817: only L=123456789.
- Prefix(13,2,18), target29960873285: only L=987654321.

The prefix graphs are explicitly provided and directly counted. The larger output counts follow from the inherited all-order physical transfer and the exact inverse, not from a claim of large adjacency subset-DP enumeration. These examples test the arithmetic method; they do not establish cofinite coverage.

Dependencies: `ROOTED.md`, SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e; `ARITHMETIC_INVERSE.md`, SHA256 4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7. Their physical one-seam law and gadget indexing are inherited. Identity(2), the monotone branches, and the square-root search are the subject of this separate addendum.
