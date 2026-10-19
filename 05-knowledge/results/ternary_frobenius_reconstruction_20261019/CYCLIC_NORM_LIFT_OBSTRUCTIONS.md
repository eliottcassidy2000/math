# Cyclic norm lift obstructions and the first Frobenius carry

Status: proved support-preserving theorem; separately finite-exact unrestricted cyclic lift census. 2026-10-19.

Use C, N, W_n, U_n, Phi from `TERNARY_FROBENIUS_RECONSTRUCTION.md`. A lift here must satisfy BOTH cyclicity Cf=f and the norm equation Nf=0. This restriction is essential: for the ambient polynomial lattice without cyclicity, p!=3 makes every mod-p norm-kernel element lift by the projector I-N/3.

## 1. Supported odd solutions have exponent-one torsion at every precision

PROVED. For every prime p>=5, every m>=2, and every odd n>=1, every f in

    W_n(Z/p^m Z) intersect U_n(Z/p^m Z)

has all coefficients divisible by p^(m-1). More precisely multiplication by p^(m-1) induces an isomorphism of F_p-vector spaces

    W_n(F_p) intersect U_n(F_p)
       -> W_n(Z/p^m Z) intersect U_n(Z/p^m Z).

Thus the right side is killed by p and has dimension floor(log_p n), but its coefficient-reduction map to W_n(F_p) is zero. No nonzero supported odd solution modulo p has a support-preserving lift to any higher precision.

### Additive-polynomial lemma

Over R=Z/p^m Z, m>=2, a polynomial a(T) with a(U+V)=a(U)+a(V) has the form

    a(T)=c T + p^(m-1) sum_{r>=1} c_r T^(p^r).

Indeed each nonzero coefficient of T^e must annihilate every middle binomial coefficient binom(e,j). If e is not a p-power and e>1, Lucas's digit criterion supplies a binomial coefficient not divisible by p, so the coefficient of T^e is zero in R. If e=p^r, every middle binomial coefficient is divisible by p, and

    binom(p^r,p^(r-1))
      = p binom(p^r-1,p^(r-1)-1)

has p-adic valuation exactly one: the final binomial is nonzero modulo p by Lucas. Hence the coefficient of T^(p^r) must be divisible by p^(m-1). Conversely these displayed polynomials are additive modulo p^m.

The same argument works with polynomial coefficients in a second variable, coefficient by coefficient.

### Proof of the theorem

The proof of Section 2 of the main note uses only that 2 and 3 are units until it invokes the additive-polynomial classification. Those units remain available in R. It therefore gives f=Phi(q), with q normalized on the axes; the symmetric odd sector is zero, and q is alternating. Its odd-in-first-variable part must be additive.

By the preceding lemma, writing Q_a=X^(n-a)Y^a-X^aY^(n-a), one has

    q = c Q_1 + p^(m-1) sum_{a=p^r<n, r>=1} c_a Q_a.

Terms with no available positive complementary exponent are simply absent. For n=1 the axis normalization gives q=0 immediately. For n>=3, direct expansion gives

    N_2 Q_1 = F_n := X^n-Y^n-(X-Y)^n,
    N_2 Q_a = F_n modulo p for a=p^r.

Consequently

    N_2 q = (c+p^(m-1) sum c_a) F_n modulo p^m.

If n is not a p-power, F_n has a coefficient prime to p. If n is a p-power, all its coefficients are divisible by p and at least one has valuation exactly one, by the same binomial calculation as above. If n is not a p-power, N_2q=0 forces c+p^(m-1) sum c_a=0 modulo p^m. If n is a p-power, it forces that scalar to be divisible by p^(m-1). In either case c is divisible by p^(m-1). Hence every coefficient of q, and therefore f, is divisible by p^(m-1).

Writing f=p^(m-1) g, the equations on f are exactly Cg=g and Ng=0 modulo p. Support is preserved coefficientwise. Conversely multiplying any supported mod-p solution by p^(m-1) gives a well-defined supported mod-p^m solution. This proves the claimed isomorphism.

### Integral Smith-spectrum corollary

Let M_n be N on the integral cyclic, two-variable-supported degree-n lattice, with n odd, and take its cokernel in the full integral degree-n polynomial lattice. For each p>=5, its p-primary torsion is exactly

    (Z/p Z)^(floor(log_p n)).

Indeed the rational kernel is zero, the mod-p kernel dimension is the classification in the main note, and the mod-p^2 theorem excludes every Smith invariant divisible by p^2. Thus, after discarding the 2- and 3-primary parts, the torsion order is the prime-to-six part of lcm(1,...,n), while its exponent is the product of primes 5<=p<=n. This is an all-degree structural corollary, not an observed Smith pattern. No assertion about the 2- or 3-primary parts is made.

As an additional finite control, `probe_supported_smith.py` computes the complete integer Smith diagonals for every odd n<=31. They match the proved p>=5 spectrum, including the two independent 5-torsion classes appearing at n=25.

## 2. A finite cyclic lift criterion without a support restriction

For f_p=Phi(X^(p-1)Y-XY^(p-1)), the integral identity is

    N f_p = F_p = X^p-Y^p-(X-Y)^p.

Every coefficient of F_p is divisible by p. Let M_p be the integer norm matrix whose columns are distinct cyclic monomial orbit sums of total degree p. A cyclic mod-p^2 lift of f_p has the form f_p+p g, with g cyclic modulo p, and it exists exactly when

    M_p g = -F_p/p modulo p

is solvable. Therefore it exists if and only if adjoining -F_p/p as a column does not increase the rank modulo p. This is an iff for all cyclic lifts, including corrections with genuinely three-variable support.

FINITE-EXACT: for every prime 5<=p<=61 (16 primes), adjoining the divided norm increases the rank by one. Hence f_p has no cyclic norm-zero lift modulo p^2 at any of those primes.

The production NumPy matrix and independent standard-library literal-substitution/sparse-column implementation agree. At p=59, the full mod-p solution space in degree 59 has dimension two, unlike the one-dimensional space at each other tested prime; nevertheless the specified f_59 still fails the lift test. Its restriction kernel is zero in that degree. The additional class is not another supported solution and is not invisible on the coordinate boundary.

OPEN: does this specified f_p fail to lift in the full cyclic lattice for every prime p>=5? The support-preserving theorem above does not settle this question; corrections with full three-variable support require their own obstruction.

## 3. Literal carry/logarithm bridge and its limit

Put C_p(A)=(A(T)^p-A(T^p))/p for an integral polynomial A. Then

    F_p(1+T,1)/p = ((1+T)^p-1-T^p)/p = C_p(1+T).

Modulo p this equals -sum_{k=1}^{p-1}(-T)^k/k, the classical finite logarithm at -T. The underlying finite logarithm is standard; see Elbaz-Vincent and Gangl, *On Poly(ana)logs I*, Compositio Mathematica 130 (2002), https://doi.org/10.1023/A:1013757217319, preprint https://arxiv.org/abs/math/0008089.

At T=1, its value is 2(2^(p-1)-1)/p. Vanishing of that ONE value is the Wieferich condition. It is only an evaluation of the divided norm. It does not make the divided norm the zero polynomial and does not imply solvability in the cyclic norm matrix. The missing coordinate is the full polynomial/cyclic compatibility class.

This is an exact substitution map between polynomial carriers. It is not an assertion that the ternary norm evolution is a cellular automaton, or that Wieferich primes decide a cellular-automaton conjecture.

## 4. Reproduction

    python3 -B probe_lift.py
    python3 -B verify_lift_independent.py
    python3 -B probe_supported_smith.py

Both scripts write their complete rank rows beside themselves, to `lift_census.json` and `lift_independent.json`. The independent verifier imports only the standalone standard-library polynomial routines from `verify_periods_independent.py`, not the production norm or matrix generator.
