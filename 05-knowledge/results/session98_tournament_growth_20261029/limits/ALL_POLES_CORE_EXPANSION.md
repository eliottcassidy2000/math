# Exact meromorphic continuation from the inherited unit-tail core kernel

Status: PROVED by the argument below; finite exact controls are separately marked.
This strengthens the bounded-radius paragraph in `LIMITS_AND_NORMALIZATION.md`.
It is a direct analytic consequence of the classical Burnside formula and the
repository's already-proved THM-514 unit-tail decomposition. No literature
priority or faster counting claim is made. The all-poles extension is independently audited in this session.

## Theorem

Let a_n=A000568(n), including a_0=1, and b_n=n!a_n/2^(n(n-1)/2).
The generating function B(z)=Σ_(n>=0)b_n z^n, initially convergent for |z|<1,
has a single-valued meromorphic continuation to the entire complex plane.
Its poles are EXACTLY z=4^k, k=0,1,2,..., and the pole at 4^k has exact order
3k+1. In particular, the normalized sequence (b_n), not only (a_n), is
non-P-recursive.

### Exact core formula

A core μ is an unordered partition with every part odd and at least three;
the empty partition is included. Put

    m=|μ|, c=length(μ), a=m-c,
    z_μ=product_j j^(multiplicity_j) (multiplicity_j)!,
    E_μ=Σ_i(μ_i-1)/2 + Σ_(i<j)gcd(μ_i,μ_j),
    L_μ=binomial(m,2)-E_μ,
    C_μ=(m!/z_μ)2^(-L_μ).

Then the exact meromorphic expansion is

    B(z) = Σ_μ C_μ z^m / (1-2^(-a)z)^(m+1).                 (1)

For the empty core, m=c=a=E=L=0, z_μ=1, and the term is 1/(1-z).
Equation (1) converges normally on every compact subset avoiding its poles.

## Proof of the exact expansion

A full odd partition is uniquely μ union 1^f. Its pair-orbit exponent is
E_μ+cf+binomial(f,2), and its centralizer denominator is z_μ f!.
Writing n=m+f, its normalized Burnside contribution is therefore

    [n!/(z_μ f!)] 2^[-L_μ-(m-c)f]
      = C_μ binomial(n,m) (2^(-a))^(n-m).

Summing in n>=m uses the elementary identity

    Σ_(n>=m)binomial(n,m)t^(n-m)z^n=z^m/(1-tz)^(m+1).

This establishes (1) coefficientwise: only finitely many cores have m<=n.
For |z|<1, all original coefficients are nonnegative and b_n is bounded by
the independent support argument in `LIMITS_AND_NORMALIZATION.md`. Absolute
convergence then justifies rearrangement there. The normal convergence proof
below gives an alternative justification and supplies the continuation.

## Normal convergence on the whole punctured plane

For a nonempty core, c<=m/3, hence a>=2m/3. The permutation of core type μ
has no fixed points. It fixes no unordered pair, since an odd-order
permutation cannot swap a pair. Every pair orbit has size at least three.
Thus E_μ<=binomial(m,2)/3 and

    L_μ >= m(m-1)/3.                                      (2)

Fix a compact K and let R=max(1,sup_(z in K)|z|). For all sufficiently large
m, put x_m=R*2^(-2m/3)<=1/2. For z in K,

    |1-2^(-a)z|^(-m-1)
      <= (1-x_m)^(-m-1)
      <= exp(2(m+1)x_m).

The final expression tends to 1, so it is uniformly bounded for all such m.
At fixed m, Σ_(μ of size m)m!/z_μ is the number of fixed-point-free
permutations of m letters with all cycles odd. It is at most m!, because
these cycle-type classes are disjoint subsets of S_m. Combining with (2),
the sum of absolute core terms of size m is bounded by a fixed constant times

    m! R^m 2^[-m(m-1)/3].                                (3)

The ratio of consecutive expressions in (3) is
(m+1)R*2^(-2m/3), which tends to zero. Thus the tail converges normally.
The finitely many excluded core terms are rational and holomorphic on K
when K avoids their poles. This proves (1) defines a meromorphic function;
its only possible poles are the displayed denominator zeros. The set of
these poles has no finite accumulation point, so each neighborhood can be
handled after separating its finite rational pole-bearing terms.

## Exact location and order: no hidden cancellation

Since a=Σ_i(μ_i-1) is even, write a=2k. Necessarily

    m=2k+c <=3k.

For a fixed k there are only finitely many cores. More precisely,
μ_i=2d_i+1 identifies them with the ordinary partitions (d_i) of k.
Consequently the k-th pole bank contains exactly p(k) core types.
For k>=1, their sizes lie between 2k+1 and 3k.

The unique core attaining size m=3k is μ=3^k. Its term has denominator
(1-z/4^k)^(3k+1). Every other term with the same pole has smaller order.
All terms with different a are holomorphic in a neighborhood of 4^k by
normal convergence and the preceding finite separation. The leading term
therefore cannot cancel. For the top core,

    E_(3^k)=(3k²-k)/2,
    L_(3^k)=3k²-k,
    C_(3^k)=(3k)!/(3^k k!)*2^(-3k²+k).

In particular, including k=0 under the empty-core convention,

    lim_(z->4^k) (1-z/4^k)^(3k+1) B(z)
        = (3k)!/(3^k k!)*2^(3k²+k) >0.                  (4)

This proves the exact pole set and orders. A D-finite analytic germ can
have finite-plane singularities only among the finitely many zeros of the
leading polynomial of its differential equation. Infinitely many genuine
poles contradict this, so B is not D-finite and (b_n) is not P-recursive.

## First banks and exact coefficient polynomials

Let H_k be the finite sum in (1) over a=2k. Then

    H_0(z)=1/(1-z),
    H_1(z)=z³/[2(1-z/4)^4],
    H_2(z)=3z⁵/[32(1-z/16)^6]+5z⁶/[128(1-z/16)^7].

For arbitrary k define the polynomial

    P_k(n)=Σ_(μ:a=2k) C_μ 4^(km) binomial(n,m).

It has degree exactly 3k and is expressed in a positive binomial basis.
For k>=1 it vanishes at every nonnegative integer n<2k+1. Therefore the following
infinite-looking formula is a FINITE exact sum at every integer n>=0:

    b_n = Σ_(k>=0) P_k(n)4^(-kn),   P_0(n)=1.            (5)

For example P_1(n)=32 binomial(n,3), and
P_2(n)=98304 binomial(n,5)+655360 binomial(n,6).
This is the precise correspondence between symmetry-defect bands, their
finite core classifications, exponentially small corrections, and poles.
It restates the inherited core decomposition in normalized coordinates;
it does not make all core banks bounded-size or remove their gcd information.

After subtracting H_0,...,H_K, the remaining Taylor series has radius exactly
4^(K+1): it is holomorphic up to the next genuine pole. Consequently every
fixed finite bank truncation supplies a controlled asymptotic expansion.
It still does NOT permit exact rounding of raw a_n without an absolute
error bound after multiplication by 2^binomial(n,2)/n!.

## Inheritance, literature, and claim boundary

The closest proved mechanism is
`01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`
(source 2026-06-15): its exact integer core K[m,c] directly gives

    B(z)=Σ_(m,c) K[m,c]2^[-binomial(m,2)]
                     z^m/(1-2^[-(m-c)]z)^(m+1).

THM-514 already isolates the same finite cores and the soluble fixed-point
tail; the present addition is the normalized rational resummation, normal
convergence, and exact global pole proof. The gcd dependence of K is retained,
not newly eliminated. The adjacent inherited source
`01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md` gives the
divisor-profile add-a-part recurrence for the same core weights.

The classical external context remains E. M. Wright, *The number of
irreducible tournaments*, Glasgow Math. J. 11 (1970), 97–101,
https://doi.org/10.1017/S0017089500000914, especially its symmetry-correction
expansion on p.100. His use of irreducible means strongly connected, not
prime under modular substitution. The analytic statement here was proved
directly rather than attributed to that paper without evidence.

Targeted web searches on 2026-10-29 for tournament meromorphic generating
functions, normalized tournament poles, and A000568 poles did not locate an
exact prior statement of this pole theorem. This is a limited negative
search result, NOT evidence of originality. The formula is a natural direct
consequence of a classical enumerator and an inherited repository theorem.

## Finite exact controls

Run `python3 session98_growth/limits/verify_all_poles.py`.
The checker independently compares full odd-partition Burnside coefficients
with rational-core coefficient extraction through n=35; verifies the
support/loss bounds, core-bank counts p(k), maximal core uniqueness, and
leading pole amplitude through k=10; and checks explicit H_1,H_2 formulas.
The analytic normal convergence and infinite pole claim rest on the proof,
not on these finite tests.
