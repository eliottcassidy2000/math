# Arithmetic dichotomy for every binary block image of Thue–Morse

2026-10-19. **PROVED from the explicitly cited scalar p-adic Mahler input below.** Independently [reviewed](review/INDEPENDENT_REVIEW.md), including the cited-source hypotheses. The elementary positive-integer exclusion in [the elementary companion](THUE_MORSE_FINITE_SOURCE_EXCLUSION.md) is independent of this analytic strengthening. Neither theorem excludes all automatic carry words or solves the Mahler Z-number problem. No external-priority claim is made.

## 1. Statement and genuine safe instances

Let t=011010011001... be Thue–Morse. Let U,V be arbitrary nonempty binary words, let h(0)=U, h(1)=V, and put c=h(t). For a binary word c define its Mahler source by

    Phi(c)=-(1/3) sum_{n>=0} c_n (2/3)^n in Q_2.

Then exactly one of the following holds.

1. If UV=VU, the word c is periodic and Phi(c) is a rational number in the real interval [-1,0].
2. If UV!=VU, then Phi(c) is transcendental over Q. Its 2-adic irrationality exponent satisfies

       2 <= mu_2(Phi(c)) <= 2/log_2(4/3).

The same rational-versus-transcendental classification, and the same exponent bound in the transcendental case, hold after prepending any finite binary word or removing any finite initial segment. In the rational case the resulting binary word is eventually periodic and its rational source still lies in [-1,0].

When U,V are loops from the initial Mahler safe follower state to itself, each containing a reset edge, every concatenation gives a strictly safe path: infinitely many loop ends supply infinitely many resets. This is the exact reset-Buechi characterization in THM-4072, `01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md`, Section1, applied to a path beginning at state0. The theorem therefore closes the ordinary-integer entry gate on every such Thue–Morse-controlled two-loop bank. It does not say that arbitrary paths through those same loops fail integer entry.

For example U=000,V=100 gives

    Phi(c)=-(1/3) tau((2/3)^3),
    tau(z)=sum_{n>=0} t_n z^n.

Every suffix has real weighted carry tail at most18/19, since successive ones are separated by at least three places. These words stay in follower states0,1,2, yet their sources are transcendental. The rational periodic control (100)^infinity has source -9/19 and is equally safe. Thus a bounded safe follower and a uniform strict safety margin do not supply ordinary integer entry.

## 2. The precise external input and its hypotheses

**CITED:** Bugeaud and Yao, *Hankel determinants, Pade approximations, and irrationality exponents for p-adic numbers*, Theorem3.1 on manuscript page6 (PDF index5), the rational-argument remark on manuscript page13 (PDF index12), and the Thue–Morse verification on manuscript page15 (PDF index14):

https://irma.math.unistra.fr/~bugeaud/travaux/p-adic18.pdf

The scalar equation f(z)=A(z)/B(z)+(C(z)/D(z))f(z^d), with integer polynomials, gives transcendence at b=r p^w/s when infinitely many principal Hankel determinants are nonzero and B,C,D do not vanish at any b^(d^m). The remark permits nonzero coprime r,s, s>0, p not dividing rs. For Hankel index ratio limsup n_(i+1)/n_i=1 and max(r p^w,s)<p^(2w), its exponent bound is 2w/(2w-log_p max(r p^w,s)). Their Thue–Morse argument verifies the required ratio1. This is a scalar first-order theorem, not a general p-adic lifting theorem for automatic linear systems.

Here the complete substitution is

    p=2, d=2, w=L, r=1, s=3^L, b=(2/3)^L,
    A(z)=z, B(z)=1-z^2, C(z)=1-z, D(z)=1.

The integer L will be |U|+|V|. All conditions are visible:

- L>=2 and r,s are coprime, with s odd.
- b is nonzero and |b|_2=2^(-L)<1, so the bounded-coefficient series converges.
- Every rational b^(2^m) lies strictly between0 and1 in the real embedding. Hence neither1-b^(2^m) nor1-b^(2^(m+1)) vanishes.
- The coefficient sequence is Thue–Morse, exactly the sequence treated in the cited Hankel verification.
- max(2^L,3^L)=3^L<4^L, the inequality needed for the quantitative bound.

Thus tau((2/3)^L) is transcendental and has exponent at most

    2L/(2L-L log_2 3)=2/log_2(4/3).

The lower bound2 is the ordinary p-adic Dirichlet lower bound for an irrational p-adic number. No claim of equality with2 is made at the rational argument2^L/3^L.

## 3. The scalar equation is exact

The identities t_(2n)=t_n and t_(2n+1)=1-t_n give, as a formal power-series identity and in its convergence domains,

    tau(z)=tau(z^2)+z/(1-z^2)-z tau(z^2)
          =z/(1-z^2)+(1-z)tau(z^2).

For completeness t is not eventually periodic. If an eventual period p>=1 existed, then for every sufficiently large k with2^k>p we would have

    t_(2^k-p)=t_(2^k)=1.

But the k binary digits of2^k-p are the complement of those of p-1, padded to k positions, so

    t_(2^k-p)=(k-s_2(p-1)) mod2,

which alternates with k. This is a contradiction.

A rational generating function with coefficients in the finite set{0,1} satisfies an eventual finite linear recurrence. Its finite set of possible length-r state vectors, with deterministic recurrence, makes its coefficient word eventually periodic. Hence tau(z) is not rational. Equivalently, by the standard Hankel rationality criterion, its principal Hankel determinants are nonzero for arbitrarily large orders. For the quantitative ratio1 one uses the specific cited Thue–Morse Hankel verification, rather than assuming that arbitrary nonrational series have that property.

## 4. Unequal blocks become equal blocks without changing the word

Let sigma(0)=01, sigma(1)=10. Since t=sigma(t),

    h(t)=(h composed with sigma)(t).

Thus put

    E_0=UV, E_1=VU, L=|U|+|V|.

The same output word is a uniform block coding of t by the equal-length words E_0,E_1. For any finite binary word E define

    P_E(z)=sum_{j=0}^{|E|-1} E_j z^j.

Breaking the output into length-L blocks gives the exact identity

    F_c(z):=sum_{n>=0} c_n z^n
           =P_(UV)(z)/(1-z^L)
             +(P_(VU)(z)-P_(UV)(z)) tau(z^L).          (1)

This does not forget the order of the blocks: the entire polynomial difference is retained.

If UV!=VU, let Delta(z)=P_(VU)(z)-P_(UV)(z). This is a nonzero integer polynomial whose nonzero coefficients are1 or-1, so its leading coefficient is a unit. Consequently2/3 cannot be a root, by the rational-root theorem. Therefore Delta(2/3)!=0.

Evaluate(1) at2/3 in Q_2. It expresses F_c(2/3) as a rational number plus a nonzero rational multiple of the transcendental value tau((2/3)^L). Thus F_c(2/3), and hence Phi(c), are transcendental. A nonconstant rational affine transformation preserves a p-adic irrationality exponent: transforming a reduced rational approximant changes its numerator/denominator height by multiplicative constants in both directions. This transfers the displayed exponent bound.

The specialization guard is essential. Replacing the signed polynomial difference by merely the pair of block lengths would erase the exact coefficient that must remain nonzero.

## 5. The commuting boundary and finite contexts

The elementary word-commutation lemma says nonempty U,V satisfy UV=VU exactly when U=W^a and V=W^b for a common nonempty word W and positive integers a,b. One proof cancels the shorter common prefix and inducts on |U|+|V|. Therefore h(t)=W^infinity.

For any eventual binary period the formal geometric-series sum is rational and is the same rational when evaluated in R and Q_2. Its real value F_c(2/3) belongs to[0,3]; hence Phi(c) belongs to[-1,0]. This proves the rational branch, including all-zero and all-one controls.

If P is a finite prefix of length a, then

    F_(Pc)(z)=P_P(z)+z^a F_c(z).

If d is obtained by deleting the first a letters of c, then

    F_d(z)=z^(-a)(F_c(z)-P_(c[0:a]))(z).

At z=2/3 both are invertible rational affine transformations. They preserve transcendence and the irrationality exponent. They also preserve eventual periodicity, proving all claimed finite-context variants.

## 6. What is new here, and what remains open

The inherited arithmetic mechanism is the source-retaining carry isometry in THM-2228 and the repetition-height estimate in THM-4654. The elementary companion closes every fixed non-erasing binary Thue–Morse block template with a terminating source-specific mismatch receipt. This extension gives its exact arithmetic boundary: commuting blocks yield rational formal sources, and all noncommuting blocks yield transcendental formal sources, with a uniform approximation exponent independent of their lengths.

The connection to the recent Collatz image-guard descent is block coding plus an arithmetic guard on an actual retained source. The mechanisms are different: the Collatz result decodes disjoint terminal images and invokes a matched-block height theorem; here balanced substitution supplies an overlong periodic shadow, and the scalar extension keeps a nonvanishing signed block polynomial. No Collatz integer source is transported into a Mahler integer source.

The original deterministic task on arbitrary reachable safe states remains open. Safe paths need not be generated by Thue–Morse or by a first-order scalar Mahler equation. In particular, the result does not establish blanket nonautomaticity of positive Mahler integer codes. The general scalar/vector distinction must remain explicit.
