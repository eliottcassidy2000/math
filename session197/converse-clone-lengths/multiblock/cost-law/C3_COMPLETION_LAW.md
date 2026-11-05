# Exact clone-repair cost distribution for a cyclic three-block quotient

Separate Session197 corollary of the frozen multi-block completion theorem. The objects here are ordinary unmarked isomorphism classes of actual tournaments C3[Tr_a,Tr_b,Tr_c], with a,b,c>=1. All are strong. Let D be the fewest vertices that must be added, using arbitrary transitive cloning, to make the tournament self-converse. This is a symmetry-completion statistic, not a Hamilton-path-preserving transformation.

## 1. Exact bivariate generating function

The canonical quotient is C3. Its three anti-involutions exchange one pair of blocks and fix the third, so the frozen completion theorem gives

    D=min(|a-b|,|b-c|,|c-a|).

Only cyclic rotations identify the weight vectors. If the three weights are distinct, an unordered sorted triple produces exactly two ordinary classes, exchanged by converse. If at least two weights agree, it produces one class, and that class is self-converse.

Write sorted weights as x<=y<=z. At D=0 the two equal weights can be a=a with a distinct third value b, including the all-equal case. The third value and the repeated value determine the unmarked class uniquely, including when all three agree. Thus the order generating function for D=0 is

    S(z)=sum_(a,b>=1) z^(2a+b)=z^3/[(1-z)(1-z^2)].

For an integer d>=1, D>=d exactly when y-x>=d and z-y>=d. Put

    y=x+d+u,  z=x+2d+u+v,  x>=1, u,v>=0.

Each such sorted triple contributes two classes. Their total order is3x+3d+2u+v, so

    sum_n #{classes of order n with D>=d} z^n
      =2z^(3+3d)/[(1-z)(1-z^2)(1-z^3)].                            (1)

Subtracting the same expression with d+1 gives the exact-D series

    2z^(3+3d)/[(1-z)(1-z^2)],  d>=1.

Therefore the complete two-variable class series is

    F(z,t)=sum_[T] z^|T| t^D(T)
          =z^3(1+t z^3)/[(1-z)(1-z^2)(1-t z^3)].                  (2)

The exact coefficient counts are

    a_(n,0)=floor((n-1)/2),                 n>=3,
    a_(n,d)=2 floor((n-1-3d)/2),            d>=1 and n>=3+3d,

and zero otherwise. In particular the largest possible cost at order n is floor((n-3)/3). It is attained: choose x=1, the first gap equal to that floor, and put the remainder into the second gap. At n=3,4,5 the maximum is zero, as the formula says.

Setting t=1 agrees exactly with ordinary Burnside:

    U_n=(binom(n-1,2)+2[3 divides n])/3,  n>=3.                   (3)

The all-equal weights have a nontrivial automorphism group, but they are counted once here. These formulas do not replace ordinary class weights by labelled or automorphism weights.

## 2. An explicit limiting distribution

Choose an ordinary class uniformly among these U_n classes. For any fixed0<x<1/3, apply (1) with d=ceil(nx). The number of nonnegative solutions to3x_0+2u+v=N, with x_0>=0, is N^2/12+O(N): equivalently this is the coefficient of1/[(1-z)(1-z^2)(1-z^3)], whose pole at1 has order three and whose other poles have smaller order. Hence

    Pr(D/n >= x) -> (1-3x)^2,     0<=x<=1/3.                     (4)

At x=0 the probability is exactly one. At x=1/3 the probability is zero because D<=floor((n-3)/3). Beyond1/3 it is zero. Thus D/n converges to the distribution on[0,1/3] with density6(1-3x). The exact probability D=0 is floor((n-1)/2)/U_n, asymptotic to3/n, so there is no limiting atom at zero.

Because0<=D/n<=1/3, bounded convergence of the tail integrals gives

    E[D]/n -> integral_0^(1/3)(1-3x)^2 dx =1/9.                  (5)

The exact generating function of the unnormalized first moments is also obtained directly from (2):

    sum_n (sum_[T],|T|=n D(T)) z^n
      =2z^6/[(1-z)(1-z^2)(1-z^3)^2].                            (6)

Equation(4) is the full cost law, not an inference from (5). It quantifies a distinction: self-converse members become rare under ordinary-class sampling, while the optimal physical completion has a completely explicit nondegenerate linear-scale cost. None of these statements promises a prime or7p Hamilton count after completion.

## 3. Checks and source boundary

The standalone C++17 checker enumerates cyclic weight-orbits by literal rotation and compares their complete D inventory with (2) and (3). On small actual blow-up graphs it separately searches growing weight vectors in increasing added order and tests self-converseness by full degree-compatible anti-isomorphism backtracking, rather than assuming that equal block weights suffice. It verifies the exact maximum, tail series and first-moment coefficients.

    g++ -O3 -std=c++17 verify_c3_cost.cpp -o verify-c3-cost
    ./verify-c3-cost
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify_c3_cost.cpp -o verify-c3-cost-debug
    ./verify-c3-cost-debug

Both outputs must agree. The programs have no input-file dependencies. The all-order conclusion follows from the explicit sorted-triple bijection and the frozen minimum-completion theorem, not from the finite controls.

Dependency: PROOF.md in this packet, SHA256934f76db619e39187d26c83f8c944694ab47371971077d32f3cdf81fdc8a2be5. That source is frozen separately and is undergoing separate independent review; this corollary requires a distinct review binding before publication. No literature-priority claim is made.
