# Unique scalar Frobenius lifts and an exact enumeration

**Status:** PROVED. Addendum to Theorem3.1 of `SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md`. The previously frozen proof packet is unchanged.

Fix a prime p and m>=1. Put `N=p^(m-1)` and let

`phi(A)(X)=A(X^p)`.

The established existence theorem says that `f in F_p[X,X^-1]` admits an integer Laurent-polynomial lift A with

`A^p=phi(A) mod p^m`                                    (1)

if and only if f=0 or its support lies in one residue class modulo N.

## 1. Uniqueness modulo p^m

**Theorem.** For each admissible f, the residue class of A modulo p^m satisfying (1) is unique, among all Laurent-polynomial lifts, without a degree restriction.

**Proof.** Suppose A and B have the same reduction f and both satisfy (1). They are congruent modulo p. If `A=B mod p^e` with `1<=e<m`, integral binomial lifting gives

`A^p=B^p mod p^(e+1)`.

Combining this with (1) gives

`phi(A-B)=0 mod p^(e+1)`.

Spatial substitution is injective coefficientwise modulo any integer: it just places the coefficient of X^j at exponent p*j. Thus `A=B mod p^(e+1)`. Induction reaches e=m. ∎

Consequently the existence construction is a complete parametrization. If

`f=X^d g(X^N)`

and B is any integer lift of g, then the unique compatible lift is

`tau_m(f)=X^d B(X)^N mod p^m`.                          (2)

The formula is independent of the choice of B by integral lifting, and independent of shifting d by a multiple of N with the compensating change in g. Define `tau_m(0)=0`.

The uniqueness is at fixed precision. Different integer representatives of the same coefficient class remain possible. Nor does a nonmonomial f admit lifts at every precision: its exact terminal precision is established in the preceding theorem.

## 2. Multiplication survives; addition generally does not

The admissible reductions form a multiplicative monoid: two support cosets multiply into their sum coset. Uniqueness and the product identity imply

`tau_m(fg)=tau_m(f)tau_m(g) mod p^m`.                    (3)

When f is admissible at precision m+1, reduction of `tau_(m+1)(f)` modulo p^m equals `tau_m(f)`. These are compatible multiplicative lifts on successively smaller domains.

Addition does not obey the analogous identity, even when f, g, and f+g all lie in the domain. For p=2,m=2,

`tau_2(1)=1, tau_2(X^2)=X^2`,

but

`tau_2(1+X^2)=(1+X)^2=1+2X+X^2 mod4`.

The intermediate carry is essential. Thus the parametrization is not a ring homomorphism, and the union of support-coset sectors must not be silently treated as an additive linear space.

## 3. Exact count in every polynomial degree box

For D>=0 let `M_(p,m,D)` be the number of coefficient polynomials in `(Z/p^m Z)[X]` of degree at most D satisfying (1). Define

`k_r=#{j:0<=j<=D, j=r modN}` for `0<=r<N`.

Then

`M_(p,m,D)=1+sum_(r=0)^(N-1)(p^(k_r)-1)`.               (4)

Writing `D+1=q*N+t`, `0<=t<N`, gives the closed form

`M_(p,m,D)=1+t*(p^(q+1)-1)+(N-t)*(p^q-1)`.              (5)

**Proof.** A nonzero admissible reduction belongs to exactly one support-residue sector. The sector r contains `p^(k_r)-1` nonzero polynomials. Add the zero polynomial. Each has exactly one compatible lift by Section1, and the construction (2) preserves its minimum and maximum exponent, so that lift remains within the same degree box. ∎

This counts actual coefficient polynomials modulo p^m, rather than only their reductions or chosen representatives. In particular, when D<N,

`M_(p,m,D)=1+(D+1)(p-1)`.

At these degrees every nonzero compatible lift is monomial. At m=1, N=1 and formula(4) becomes `p^(D+1)`, as freshman's dream requires.

## 4. Existing independent controls

The producer's complete-lift enumeration already records both the number of compatible lifts and the number of admissible reductions. They are equal in every declared box, as the uniqueness theorem predicts. For example:

- p=2,m=3,D=4: exactly 7 compatible coefficient polynomials among 32,768 lifts;
- p=3,m=2,D=3: exactly 13 among 6,561;
- p=5,m=2,D=2: exactly 13 among 15,625;
- p=2,m=4,D=2: exactly 4 among 4,096.

The independent reviewer also enumerated the first three boxes and p=2,m=4,D=3, whose exact count is 5. The checker accompanying this addendum verifies every stored count against (4) and (5), independently enumerates the support-sector counts for bounded parameters, and checks the additive hostile. A fresh dense-convolution path then enumerates all 120,490 coefficient lifts in the four independent-review boxes, retaining the multiplicity of every reduction fibre. Every admissible fibre has size exactly one, and every inadmissible fibre is empty. Normal and optimized runs agree. Those finite controls support, rather than replace, the universal proof.
