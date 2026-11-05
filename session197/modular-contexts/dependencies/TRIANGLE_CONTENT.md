# Repeated triangle tests generate universal content; cyclic composition can lose it

Session197. This separate consequence of the reviewed central-factorial operator and content theorems has two parts. First, a one-letter physical wrapper bank computes the exact universal intact-module divisor, and a small modification uses only strong outputs. Second, an explicit physical family shows that the more restrictive cyclic divisor is not preserved under composition, even though the universal divisor is. No new prime/7p coverage or infinitude is asserted.

Let T be a nonempty n-vertex tournament. Write H(T) for its directed Hamilton-path count and B_T(x) for its ordered-bin polynomial. Let g(T), G(T), C(T) be, respectively, the ordinary bin-value gcd, the gcd of all whole-module context counts, and the gcd of H(C3[A,B,T]) over nonempty A,B. The independently reviewed content theorems give

    g=gcd_k gamma_k u_k,
    G=gcd_k k! gamma_k u_k,
    C=gcd_k (k!)^2 gamma_k u_k,                         (1)

where d=ceil(n/2), B_T=sum_(k<d) u_k b_k, and

    b_k=x product_(i=1)^k(x^2-i^2), gamma_k=(2k+1)!       for odd n,
    b_k=x^2 product_(i=1)^k(x^2-i^2), gamma_k=(k+1)(2k+1)! for even n.

These contents are odd and have the same prime support. Their gcd interpretation does not assert that the gcd itself is attained.

## 1. One repeated wrapper is a complete integer response bank

Define actual wrappers

    W(T)=C3[T,K1,K1],
    E(T)=C3[T,Tr2,K1].

The inherited full-profile operator for W is

    Df=x^2 f+x[f(x+1)-f(x-1)].

It preserves the parity of f. In the central basis its matrix is tridiagonal:

    D b_k=b_(k+1)+a_k b_k+beta_k b_(k-1),

where

    a_k=k^2+6k+3, beta_k=2k^2(2k+1)       for odd parity,
    a_k=k^2+6k+5, beta_k=2k(k+1)(2k+1)   for even parity.

At k=0 the last term is absent. These formulas, including the k=0 boundary, were proved in the named path-polynomial source. Evaluation at1 reads the coefficient of b_0, since b_0(1)=1 and b_k(1)=0 for k>0.

Put

    p_k=product_(i=1)^k beta_i=k! gamma_k,   p_0=1.      (2)

The product identity follows directly by separating the even and odd factors of(2k+1)!. Let M be the displayed tridiagonal matrix, with input basis elements as columns. Each term of(M^j)_(0,k) is a weighted path from height k to height0 using steps down, level and up. Such a path crosses every down edge i->i-1,1<=i<=k, at least once. Therefore

    (M^j)_(0,k) is an integer multiple of p_k;
    (M^j)_(0,k)=0 if j<k;
    (M^k)_(0,k)=p_k.                                  (3)

The last equality holds because the only k-step path from k to0 always moves down. No assumption about Hankel positivity for an arbitrary initial input is needed.

Consequently the d values

    H(T), H(WT), ..., H(W^(d-1)T)                      (4)

are an integer lower unit-triangular transform of the coordinates p_k u_k. In particular

    gcd_(j>=0) H(W^jT)
      =gcd_(0<=j<d) H(W^jT)=G(T).                     (5)

Later outputs are still integer combinations of p_k u_k by(3), so they cannot reduce the gcd further. These actual test graphs have orders n,n+2,...,n+2d-2. This is the same size bound as the earlier matching-context bank, but uses one fixed wrapper repeatedly. The bound is sufficient, not asserted optimal.

There is also an exact response-lattice statement. In any whole-module context, the maximal-run expansion has coefficients(r-1)! times integers on the ordered-cover coordinates f_r. The inherited Newton change of basis expresses each(r-1)!f_r as an integer combination of p_k u_k. Since the square transform in(4) has an integer inverse, every such context count is an integer linear combination of the d physical counts in(4). The coefficients may be negative. This is not a positive construction of an arbitrary target integer, and it is not a fixed finite state space independent of n.

The identity output j=0 matters. For T=K1, the universal gcd is1, whereas the gcd of all positive iterates j>=1 is3: W(K1)=C3 has universal content3, and the first two positive counts are3 and15. All subsequent counts are divisible by3. Excluding j=0 without another operation would change the theorem.

## 2. Restricting the test outputs to strong tournaments loses nothing

Both W(T) and E(T) are strong for every nonempty T, even when T is weak. The inherited scalar readouts are

    H(WT)=H(T)+B_T(2),
    H(ET)=3B_T(2)-H(T).

Thus

    3H(WT)-H(ET)=4H(T).                               (6)

All these actual Hamilton counts are odd. Define the finite strong bank

    {W^j(T):1<=j<=max(1,d-1)} union {E(T)}.             (7)

Its gcd divides4H(T) by(6), hence divides H(T) because it is odd. It also divides every nonidentity value in(4). Equation(5) therefore forces this gcd to divide G(T). Conversely G(T) divides every member of(7), because each is a whole-module context. Hence their gcd is exactly G(T).

It follows that the gcd over all strong whole-module outputs is the same as the gcd over all whole-module outputs, with no assumption that T itself is strong. Bank(7) has max(2,d) graphs and maximum order

    n+max(3,2d-2).

This is a sufficient bound, not a minimum-size assertion. At n=1 the two outputs have H=3 and5, recovering gcd1 without an identity output. Strongness and the odd-count argument are essential; the statement is not being transferred to arbitrary digraphs.

As another direct proof of equality for the infinite strong class, apply W and E to any whole-module output S. Both results are strong whole-module outputs of T, and(6) recovers4H(S). Their odd gcd must therefore divide every H(S).

## 3. Exact repeated-join hierarchy

Use the reviewed physical ten-vertex tournament

    V=[382,420,88,706,810,972,18,597,204,327]

with bin polynomial

    P(x)=x^10+64x^8+1000x^6+4020x^4+2772x^2,

and contents(g,G,C)=(9,27,81). Define J_m to be the ordered join of m copies of V, for every m>=1. In every ordered bin, paths in successive blocks concatenate uniquely, so

    B_(J_m)(x)=P(x)^m.                                 (8)

Pointwise powers multiply every prime valuation of the integer-value gcd, giving g(J_m)=9^m exactly.

The small evaluations are

    P(1)=7857, P(2)=156816, P(3)=1558521,
    P(4)=10412352, P(5)=52972425.

The first three have exactly four factors of3; P(4) and P(5) have exactly two, with

    P(4)/9=P(5)/9=2 (mod3).                            (9)

Write A_k=gamma_k u_k for the normalized central coordinates of P^m. The evaluation matrix of b_k/gamma_k at x=1,2,... is lower unit triangular: its entries are integers, it vanishes when x<=k, and equals1 at x=k+1. Thus each A_k is an integer combination of the first k+1 evaluations. The first three A_k have3-adic valuation at least4m, while A_3 has exactly2m because its P(4)^m term has coefficient one and its earlier terms have larger valuations.

Every A_k has valuation at least2m, since g(J_m)=3^(2m). For k>=3, k! is divisible by3. Combining these facts with(1), and using v_3(3!)=1, gives

    v_3(G(J_m))=2m+1,
    v_3(C(J_m))=2m+2.

The equality of prime supports in the content theorem excludes every other prime. Therefore the complete contents are

    (g(J_m),G(J_m),C(J_m))=(3^(2m),3^(2m+1),3^(2m+2)). (10)

This is an all-order physical family of strict gaps. The stronger contents are not multiplicative under ordered joins. For example,

    G(V->V)=243 != G(V)^2=729,
    C(V->V)=729 != C(V)^2=6561.

The ordered joins for m>1 are weak. Their strong wrapper companions are treated next; no strongness is silently inferred for J_m.

## 4. A strong wrapper lowers the cyclic content, and two cyclic stages escape it

Put S_m=W(J_m), a strong tournament of order10m+2, and Q_m=D(P^m). Every integer evaluation of Q_m is divisible by3^(2m). At x=4, equations(9) show

    Q_m(4)=16P(4)^m+4[P(5)^m-P(3)^m]

has exactly2m factors of3, because the two lowest terms sum to20*2^m after division by3^(2m), a nonzero residue modulo3. Hence

    v_3(g(S_m))=2m.

At x=1,2 the valuations of Q_m are at least4m. At x=3,

    Q_m(3)=9P(3)^m+3[P(4)^m-P(2)^m]

has exactly2m+1 factors of3. The same normalized central-evaluation argument now gives valuations at least4m for coordinates k=0,1 and exactly2m+1 for k=2. Every later coordinate has valuation at least2m, and its k! contributes at least one factor of3. Therefore

    v_3(G(S_m))=v_3(C(S_m))=2m+1.                     (11)

Only these3-adic assertions are made for S_m; no unproved exclusion of other primes is needed. Since C(J_m) has valuation2m+2, it does not divide C(S_m). Thus cyclic content is not monotone even under the actual triangle whole-module wrapper. This does not contradict the universal-content monotonicity, because the two universal valuations in(10) and(11) agree.

The failure is witnessed by an explicit second cyclic stage, not only by a gcd comparison. Let

    Z_m=C3[Tr2,Tr2,S_m],

which is strong of order10m+6. The reviewed cyclic kernel, or a direct three-bin run expansion, gives for any input bin profile Q

    H(C3[Tr2,Tr2,input])=4Q(3)+5Q(2)-11Q(1).

Substitute Q=D(P^m) and simplify:

    H(Z_m)=12P(4)^m+46P(3)^m-3P(2)^m-21P(1)^m.        (12)

The first term has valuation2m+1, strictly below all other terms. Thus H(Z_m) has exactly2m+1 factors of3 and escapes the original cyclic divisor C(J_m)=3^(2m+2) at every m>=1.

For m=1 the actual16-vertex graph has H=196004745, which is54 modulo81. This does not contradict the definition of C(V): Z_1 is a composition of two cyclic contexts, not a single context C3[A,B,V] with V as one whole block. Equation(5) also shows conceptually that the gcd over iterated cyclic contexts, including the identity, is G(T), not generally C(T). In fact the strong bank(7) already lies in this iterated-cyclic language. Thus the same gcd holds with the identity excluded and only strong outputs allowed, provided both W and E are permitted. This is distinct from the one-letter positive-iteration counterexample in Section1.

## 5. Reproduction, boundaries and dependencies

The standalone verifier checks the exact central lowering products, integer triangular response matrices, strong-bank gcd and reconstruction on actual small tournaments and signed formal profiles. It separately verifies the ordered-power contents and valuations, constructs the physical join/wrapper graphs, and directly counts the12- and16-vertex witnesses. Optimized Python and ordinary Python outputs must agree.

    python3 verify.py
    python3 -O verify.py

The file has no input dependencies. The all-order arguments above use the finite tridiagonal structure, integer triangularity and exact valuations, not extrapolation from the tested powers.

Named reviewed mathematical sources:

- Central-factorial D/E operators and physical wrappers: session197-odd-cycle-path-polynomials/PROOF.md, SHA256b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b.
- Universal-content theorem and normalized central evaluations: UNIVERSAL_MODULE_CONTENT.md in that directory, SHA2564f86f7ada5eff79e61bb2ef4b8139b76df8efedf6996a71dae327bc66b24d9cb.
- Exact cyclic content, physical V and its bin polynomial: CYCLIC_MODULE_CONTENT.md in that directory, SHA25615ec2425b33d2b16a0a7575592bda2f1a5f97b7b257656e25f008311c4272571.
- Cyclic Gram kernel used in(12): CYCLIC_GRAM_FORESTS.md in that directory, SHA2562fe61a1fcbd4133fbd9bddcdd8a337dfe521d10ca41514f69ecb25f358327c83.

The prior packets already prove the individual operator and content formulas. This packet adds the repeated-wrapper integer bank, the strong-only refinement and the all-order composition counterexample. No literature-priority claim is made.
