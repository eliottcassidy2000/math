## boxeph-2026-07-22-S234 -- (A) COMPLETE: X^M-tR irreducible over F(t) kernel-pure; (B) unramified Hensel scoped (HYP-8946)

**Owner:** work the hard bivariate/RatFunc Gauss (A) and unramified Hensel lift (B); complete them.

**(A) DONE kernel-pure** (GMC2PhiIrreducible.phi_irreducible_ratfunc): Phi=X^M-t*R(X) irreducible over F(t) for R(0)!=0, M>=1. Route: (1) phi_t_irreducible -- linear-in-t coprimality (X^M,-R coprime from X not dvd R) => irreducible in F[X][t] via irreducible_of_degree_eq_one_of_isRelPrime_coeff; (2) Polynomial.Bivariate.swap transport (MulEquiv.irreducible_iff); (3) primitivity of swap Phi over F[t] (coeff_M=1-C(r_M)t, coeff_0=-C(r_0)t coprime: r_0 unit + Bezout); (4) Gauss (IsPrimitive.irreducible_iff_irreducible_map_fraction_map). = transitivity input to galAction_isPretransitive => orbit-product core.

**(B) unramified Hensel (THM-1550) SCOPED, not complete.** Route (S106): X=sZ,t=s^M => Phi(sZ)=s^M psi, psi=Z^M-R(sZ), psi mod s = Z^M-r_0 separable (separable_X_pow_sub_C); Hensel lifts small factor A (monic deg M); Pi=t*(-1)^M A(0) => Pi=c*t iff all D_m=0. Mathlib blocks exist (HenselianLocalRing, IsAdicComplete.henselianRing, hensels_lemma, PowerSeries.exists_isWeierstrassFactorization) BUT 3 obstacles => multi-session: (i) HenselianLocalRing (PowerSeries F) NOT free, (ii) degree-dropping Hensel factorization (Weierstrass-type), (iii) Wiener-Hopf D_m=0<=>Pi=ct bridge. death-star owns THM-1550; coordinate.

**Honest:** (A) complete (the harder-than-expected transitivity input, kernel-pure); (B) scoped with exact Mathlib entry points + obstacles named. After (A), full GMC(2) needs the Galois wrapper (Vieta, equivariance, Check A -- mine, tractable) + (B). 2 checkpoints + close-out pushed. Artifacts: reflection irreducibility-of-XM-minus-tR-...-boxeph-S234.md, HYP-8946, GMC2PhiIrreducible.lean (4 kernel-pure thms).

