# Complete scalar p-grid depth and the intrinsic Frobenius-lift obstruction

**Status:** PROVED, scoped to scalar integer Laurent-polynomial rules. Exact companion controls are in `verify_template_depth.py`. This is a continuation of `WRONSKIAN_FROBENIUS_TRANSFER.md`, whose power-lifting argument is repeated where needed. No priority claim is made over the classical theory of additive cellular automata or p-derivations.

## 1. The complete one-variable late-time classification

Let p be prime, `m>=1`, and `A in Z[X,X^-1]`. For an integer `n>=m-1`, reduce the kernel

`K_(n,m)=A(X)^(p^n) mod p^m`.

Write `f=A mod p`. The following trichotomy is exhaustive.

1. If f=0, then `K_(n,m)=0`.
2. If `f=aX^d` is a nonzero monomial, then `K_(n,m)` is a nonzero monomial at exponent `d*p^n`.
3. Suppose f has at least two terms. Put

`S=support(f),  d=min S,`

`s=min_(e in S, e!=d) v_p(e-d)`.

Then the exact p-adic granularity of the support differences of `K_(n,m)` is

`min_(u,v in support K, u!=v) v_p(u-v) = n-m+1+s`.       (1.1)

In particular, after a monomial shift and removal of the maximal spatial p-dilation already present in f, every nonmonomial scalar rule has the full template depth `m-1`. A drop in this normalized depth is impossible. This is a classification of the support's p-grid, not a formula for every nonzero coefficient or every possible odd/common spacing factor.

### Proof

We use two elementary identities:

- if `U=V mod p^e`, e>=1, then `U^(p^n)=V^(p^n) mod p^(e+n)`;
- for `n>=m-1`, `B(X)^(p^n)=B(X^(p^(n-m+1)))^(p^(m-1)) mod p^m`.

The first follows by a binomial expansion and induction; the second follows by applying it repeatedly to `B(X)^p=B(X^p) mod p`.

Thus at our times K depends only on f, and we may choose the canonical coefficient representatives `0,...,p-1`. If f=0, compare A to zero modulo p. If f is a monomial, compare A to its single monomial representative; its coefficient is a unit, proving the first two cases.

In the third case write that representative as

`A_0=X^d B(X^(p^s))`,

where B has a nonzero constant term `a_0`, all its nonzero coefficients are units modulo p, and at least one positive supported exponent is not divisible by p. Set

`N=p^(m-1), L=n-m+1+s`.

The fixed-template identity gives

`K_(n,m)=X^(d*p^n) B(X^(p^L))^N mod p^m`.               (1.2)

All support differences are therefore divisible by `p^L`. It remains to exclude `p^(L+1)`.

Let j be the smallest exponent of B not divisible by p, and write its coefficient as a_j. Every smaller positive supported exponent is divisible by p. A product of terms at such smaller exponents cannot sum to j modulo p. A term of exponent at least j can occur in a contribution to degree j only by taking exactly the j term once and constant terms in every other factor. Hence the following is an **exact integer coefficient identity**:

`[X^j] B(X)^N = N*a_j*a_0^(N-1)`.                       (1.3)

Its p-valuation is exactly m-1, so it is nonzero modulo `p^m`. The constant coefficient is also a unit. Equation (1.2) therefore has nonzero terms at exponents

`d*p^n` and `d*p^n+j*p^L`.

Their difference has valuation L because p does not divide j. This proves (1.1). ∎

### Boundaries and examples

- For `A=1+X`, s=0, so the full depth is m-1.
- For `A=1+X^(p^s)`, the visible grid is coarser by s powers of p, but after removing that original spatial dilation its depth is again m-1.
- For `A=X^7+X^11` and p=2, the support difference is 4, so s=2. Modulo `2^m`, late-time support differences have exact valuation `n-m+3`.
- Coefficients divisible by p outside the mod-p support do not alter the late-time statement: they are erased by the integral lifting congruence. They can matter before n=m-1.
- A seed S has not been included in (1.1). Multiplication by S, finite-torus aliasing, and nonunit content can change the evolved state's support. Apply the kernel statement before making a state-support claim.
- The proof requires one scalar commutative Laurent polynomial. The explicit noncommuting 2x2 lift in the first note prevents an arbitrary matrix-valued extension.

## 2. A lift-independent first Frobenius obstruction

For an integer Laurent polynomial A define its integral carry polynomial

`C_p(A)=(A(X)^p-A(X^p))/p`.

This is integral by the binomial theorem and Fermat's congruence on coefficients. Let `f=A mod p`. The polynomial `C_p(A) mod p` depends on the chosen integral lift A. Its class in the **additive quotient**

`F_p[X,X^-1] / F_p[X^p,X^-p]`                            (2.1)

does not. Denote that class by `delta_perp(f)`. The quotient is one of additive F_p-vector spaces; the denominator is a subring, not an ideal, and (2.1) is not claimed to be a quotient ring.

### Theorem 2.1. Complete first obstruction and existence of a compatible lift

For a nonzero Laurent polynomial f, the following are equivalent:

1. `delta_perp(f)=0`.
2. The support exponents of f all lie in one residue class modulo p.
3. `f=X^d g(X^p)` for some integer d and some Laurent polynomial g over F_p.
4. There exists an integer Laurent polynomial lift A of f such that

`A(X)^p=A(X^p) mod p^2`.                                 (2.2)

For f=0, the obstruction also vanishes and the zero lift satisfies (2.2).

This is an existence theorem for a lift. Even when these conditions hold, a specific previously chosen lift can have nonzero C_p(A) modulo p.

### Proof: change of lift

If `A_1=A+pB`, then binomial expansion gives `A_1^p=A^p mod p^2`, so

`C_p(A_1)=C_p(A)-B(X^p) mod p`.                          (2.3)

This proves lift independence of the class in (2.1). It also proves the equivalence of conditions 1 and 4: if C_p(A) is supported on p-divisible exponents modulo p, write it as `b(X^p)`, choose any integer lift B of b, and use A+pB. Conversely a compatible lift has zero class.

### Proof: derivative characterization

Over the integers, differentiation of the defining identity gives

`C_p(A)'=A^(p-1)A'-X^(p-1)A'(X^p)`.

Modulo p, this becomes

`C_p(A)'=f^(p-1)f'-X^(p-1)(f')^p`.                      (2.4)

A Laurent polynomial over F_p has derivative zero iff all its exponents are divisible by p. Thus condition 1 is equivalent to the vanishing of (2.4).

If `f'=0`, then `f=g(X^p)` and condition 3 holds with d=0. Otherwise divide (2.4) in the rational-function field to obtain

`f^(p-1)=(Xf')^(p-1)`.

Every `(p-1)`st root of unity in F_p(X) is a nonzero constant in F_p, since `T^(p-1)-1` splits into the corresponding distinct linear factors. Therefore

`Xf'=lambda*f` for some `lambda in F_p^*`.

Comparing coefficients shows that every supported exponent is lambda modulo p. This proves condition 2. Conversely condition 2 gives `Xf'=d*f`; substituting in (2.4) makes it zero because `d^p=d`. Finally, conditions 2 and 3 are equivalent by factoring any one support exponent. ∎

### Why this explains template-depth loss

A normalized nonmonomial f with at least one exponent not divisible by p has a nonzero transverse obstruction. No change of integral coefficient lift can remove its first off-grid carry. The minimal-exponent coefficient in (1.3) is an all-precision persistence certificate for the same obstruction.

If all support differences are divisible by p, the first transverse obstruction vanishes, but this means there was already a spatial p-dilation after shifting the monomial. Remove it and repeat. Section 1 shows precisely how many such removals are possible and why no further normalized depth reduction exists.

### Positive and hostile controls

- At p=2, `f=1+X` has `C_2(1+X)=X`, which is off the p-grid. No lift can make its first Frobenius defect vanish modulo4.
- At p=2, `f=1+X^2` has `C_2(1+X^2)=X^2`. The class vanishes, and the corrected lift `A=1+2X+X^2` satisfies `A(X)^2=A(X^2) mod4`.
- At p=3, the monomial lift A=2 has `C_3(2)=2 mod3`, so this fixed lift is not compatible. Its support is one residue class, and the corrected lift A=8 does satisfy `8^3=8 mod9`. The fixed-lift defect and the intrinsic obstruction are distinct objects.
- The additive quotient in (2.1) cannot be silently treated as an ideal quotient: for example X times a p-grid monomial need not remain on the p-grid.

## 3. Complete finite-precision existence theorem

### Theorem 3.1. Frobenius-compatible integral lifts to every fixed precision

Fix `m>=1` and `f in F_p[X,X^-1]`. There exists an integer Laurent polynomial lift A of f satisfying

`A(X)^p=A(X^p) mod p^m`                                  (3.1)

if and only if f is zero or its support lies in one residue class modulo `p^(m-1)`. Equivalently, for nonzero f,

`f=X^d g(X^(p^(m-1)))`                                   (3.2)

for an integer d and a Laurent polynomial g over F_p. At m=1 the condition is vacuous, as it must be. At m=2 it recovers Theorem 2.1.

**Sufficiency, including an explicit lift.** Put `N=p^(m-1)`, choose any integer Laurent lift B of g, and set

`A=X^d B(X)^N`.                                          (3.3)

Modulo p, freshman's dream gives exactly (3.2), so A lifts f. Also

`A^p-A(X^p)=X^(d*p) * [(B(X)^p)^N-B(X^p)^N]`.

The expressions inside the Nth powers are congruent modulo p. Integral power lifting gives divisibility by `p^m`. For f=0 choose A=0.

**Necessity.** Let `phi(A)=A(X^p)`. From (3.1), applying phi and raising to pth powers inductively gives

`A^(p^r)=phi^r(A) mod p^m` for every r>=1.                (3.4)

Take r=m-1; for m=1 the needed condition is already vacuous. The right side is supported on exponents divisible by `p^(m-1)`. If f is nonzero and nonmonomial, Section 1 says the exact p-content of the left side's support differences at time `p^(m-1)` is

`s=min_(u,v in support(f),u!=v) v_p(u-v)`.

Thus s>=m-1, precisely the one-residue-class condition. Zero and monomial f already satisfy the conclusion. ∎

### Exact obstruction depth

For a fixed nonmonomial f whose support differences have p-content s, a compatible lift exists through precision `p^(s+1)`, but none exists at precision `p^(s+2)`. This obstruction cannot be removed by altering coefficients, allowing extra terms divisible by p, or changing the integer lift's support. The necessity proof quantifies over every Laurent-polynomial lift.

The lift in (3.3) may depend on the requested precision. This theorem does not assert that a single lift works at every precision. If f is a nonmonomial, its finite obstruction depth proves that no such universal lift exists. For monomial f, the question of a single *integer* coefficient lift at all precisions also includes the equation `a^p=a` over the integers; allowing p-adic Teichmuller coefficients is a different coefficient ring.

### Degree-preserving construction and finite checks

For an ordinary polynomial f with support in `[d,D]`, choose g with minimum exponent zero in (3.2). Then (3.3) stays in the same exponent interval `[d,D]`; intermediate coefficients may appear. This permits an independent exhaustive check: enumerate all coefficient lifts modulo p^m up to degree D, record exactly which mod-p polynomials admit a compatible lift, and compare with the one-residue-class criterion. The companion does this in declared complete small universes, in addition to the general proof.

## 4. Connection contract and remaining scope

Source: the mod-p support of a scalar local rule. Target: its late-time p^m kernel granularity and the class of its first integral Frobenius defect. Maps: integral lifting, spatial substitution, differentiation of C_p, and the minimum non-p-divisible exponent. Preserved: the exact p-grid obstruction and the coefficient valuation m-1 witnessing it. Discarded: finer residue values, seed overlaps, finite-torus identifications, and noncommuting structure. Required sidecars: integral-lift choice for a fixed C_p, monomial shift, spatial p-dilation, coefficient precision, and the fact that (2.1) is additive.

The complete one-variable scalar p-grid classification and the finite-precision lift-existence iff are now proved. It does not classify finer prime-to-p spacing, all nonzero coefficient patterns, nonlinear automata, matrix automata, or arbitrary carry systems. Those remain separate questions. A potential controlled-commutator extension would need an explicit commutator valuation law; the scalar proof does not provide one automatically.
