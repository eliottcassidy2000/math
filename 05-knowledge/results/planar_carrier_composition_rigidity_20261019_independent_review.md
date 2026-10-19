# Independent review: polynomial realizations of finite universal-flow compositions

2026-10-19. **ACCEPTED ANALYTICALLY.** The arguments below independently validate both proofs of the finite-composition obstruction. Scope is finite compositions of the specified universal Hamiltonian flows in the original source plane, not arbitrary Keller maps.

## 1. Source completion and the external automorphism theorem

Each carrier S=c+p²DR has original source expression S-c=t⁵q⁴R(tq,xt²q), q=1+x²t. Its Hamiltonian derivation raises source t-order by at least four, its exponential exists in K[x][[t]], preserves the ordinary source Poisson bracket, and fixes t=0 pointwise. Finite compositions retain these properties. A polynomial realization therefore has determinant one and restricts to (x,0) on the line t=0.

[Janusz Gwoździewicz, Injectivity on one line](https://arxiv.org/html/alg-geom/9305008v1), Theorem 1.1 and directly applicable Lemma 2.3, were checked in the primary author text. They imply that this polynomial realization is an automorphism over every algebraically closed characteristic-zero field. This uses a proved special theorem, not the full Jacobian conjecture. For arbitrary characteristic-zero K, extend scalars to its algebraic closure; the identity ultimately proved there descends as a polynomial identity, without needing a separate descent of an inverse.

An automorphism fixing the line pointwise carries its principal ideal (t) to itself, so T=at for a nonzero scalar a. The determinant equation implies X=a^(-1)x+f(t). Restriction to the line gives a=1 and f(0)=0. Thus every hypothetical realization is a shear (x+f(t),t).

## 2. Valid boundary-evaluation route

Use the specific common ring L=K[s][[tau]][tau^(-1)], with x=s/tau, t=tau, p=s²+tau and y=s(s²+tau). The source substitution sigma:L to K[x]((t)), s maps to xt and tau maps to t, is well-defined, continuous and injective: its monomial exponent map (i,j) to (i,i+j) is injective, and the uniform lower bound on j makes every target coefficient a finite sum. Polynomial source expressions lie in L, and sigma intertwines their Hamiltonian derivations and convergent exponentials.

The boundary evaluation epsilon:L to K((s)), tau maps to -s², is also well-defined and continuous. Coefficients are polynomials in s, so a tail of tau-order N has image of s-order at least 2N. This statement would fail on an unrestricted coefficient field K(s) with unbounded negative s-orders; no map on all of K(s)((tau)) is asserted.

In L, S-c=tau(s²+tau)^4R(s²+tau,s(s²+tau)). Consequently both delta(s)=tau S_tau and delta(tau)=-tau S_s are divisible by (s²+tau)^3. For every element of L, the chain-rule formula for delta and the ring homomorphism epsilon yield epsilon(delta f)=0. Since delta raises tau-order, its exponential converges and epsilon(E f)=epsilon(f). Every finite composition, without any common invariant assumption, retains this equality.

For a polynomial realization, injectivity of sigma identifies its source-polynomial rewrites with the same common-ring series. The shear obtained above therefore satisfies

    0=epsilon(C(x)-x)=epsilon(f(tau))=f(-s²).

The substitution K[t] to K((s)), t maps to -s², is injective. Hence f=0 and the polynomial realization is the identity. No inference from a t-adic equality to original q-adic divisibility is made; q is a unit in the source t-adic completion. Boundary evaluation in the restricted common ring is the indispensable justification.

## 3. Independent completed-cusp intersection route

Let A=K[p,y], D=p³-y² and let Â be its unlocalized D-adic completion. The identity y²=p³-D gives the unique normal form

    sum over n>=0 of D^n(a_n(p)+y b_n(p)), a_n,b_n in K[p].

This follows already at every quotient A/(D^N), since A is free with basis 1,y over K[p,D]. Substitution into K[s][[tau]] is continuous because D=tau(s²+tau)². A nonzero element of D-order n has tau-leading coefficient

    s^(4n)(a_n(s²)+s³ b_n(s²)).

The parenthesis is nonzero: its even and odd monomial exponents cannot cancel. This also proves injectivity of the completed substitution. If the image lies in K[[tau]], its order-zero coefficient must be scalar. The normal form then permits subtracting that scalar. Any nonzero remainder has D-order n>=1, and its tau-leading coefficient is a nonzero polynomial divisible by s^(4n), impossible for a scalar coefficient. Thus Â intersects K[[tau]] precisely in K.

The universal derivative formulas give delta(A) contained in DA, delta(D) contained in D²A, and delta(x) in A. Therefore the exponential preserves Â, and E(x)-x belongs to Â because delta^j(x) belongs to D^(j-1)A for j>=1. Compatibility with the common-ring realization follows by continuity and matching polynomial partial sums. Induction shows that every finite composition preserves the affine set x+Â. For its hypothetical shear realization, injectivity of sigma yields f(tau) in Â. The intersection result makes f constant; f(0)=0 makes it zero. This proof also requires neither commutativity nor a common Hamiltonian.

## 4. Independent finite controls and scope

A separate symbolic calculation passed 95 checks: seven carriers tested triple boundary divisibility of both vector-field components and zero evaluation of derivatives on six common-ring inputs; twenty-five canonical normal-form terms tested the completed-cusp leading coefficient for D-orders zero through four. These finite controls support, but do not establish, the universal statements proved above.

The identity conclusion permits nontrivial words of flows that cancel exactly. It excludes nonidentity polynomial realizations of finite compositions of these universal carriers. It does not exclude rational realizations, arbitrary other Hamiltonians, unrestricted infinite compositions, or transformations after a different choice of source chart. It supplies no chart-entry theorem and no proof of planar JC.

## 5. Final proof-text verification

The producer proof `planar_carrier_composition_rigidity_20261019.md`, SHA-256 `d96a53fea47af0453f66cbb13a9f39cc21a524d4952e499abfdb1693db18015d`, was read in full. Both the main boundary-evaluation proof and the separate completed-cusp proof are accepted. The text explicitly confines the evaluation to the polynomial-coefficient common ring, provides the needed continuity and injectivity, uses the classical fixed-line theorem only within its stated hypotheses, and keeps original-source polynomiality distinct from cusp-coordinate polynomiality. No mathematical correction was required. The independent checks described above were run separately from the producer's checker.

## 6. Final replay receipt

The final producer proof has SHA-256 `cf79b7cc010be7042849e1066da5fcfdcf9e72c800fd865df1fad76b29dece3f`; its mathematical content is unchanged, with status and reproduction references updated. The standalone checker was independently rerun both normally and with `python -O`. Both completed successfully, reported 470 checks, and matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `d92c484ee92c4705410474940fc16c23657b355f8f9551d69649264eb93620bd`. Analytical acceptance and finite reproduction are both confirmed.
