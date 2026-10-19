# Independent review: monomial-carrier scalar-time trichotomy

2026-10-19. **ACCEPTED ANALYTICALLY** after reading the complete producer proof. The classification applies to H_k=x^k t²(1+x²t), integer k>=0, and the specified original-source scalar-time exponential over any characteristic-zero field. It is not a classification of arbitrary carriers or mixed-invariant compositions.

## 1. Source existence and individual algebraicity

Every H_k is divisible by t², so its derivation raises t-order and all scalar-time source exponentials exist. Their coordinate images are nonzero. The invariant equation is a nonzero positive-degree polynomial in either coordinate when the other is adjoined; the special k=0 equation is also correctly supplied. Therefore algebraicity of one coordinate implies algebraicity of the other. Excluding simultaneous algebraicity consequently establishes individual transcendence, not merely a transcendental pair.

## 2. Scaling from fixed time to formal time

The dilation weights are correct: x has dilation weight one, t has weight -2, H_k has weight k-4, and its Hamiltonian derivation shifts weight by k-3. This proves the two displayed flow-scaling identities term by term.

For b=k-3 nonzero, the substitution tau=lambda z^b is injective on polynomials in tau and gives a finite algebraic field extension after adjoining z. Thus fixed-time algebraicity gives algebraicity over K(x,t,tau) by transitivity. The coefficient-array comparison is essential and valid: every fixed t coefficient of the flow is polynomial in time, while every fixed time coefficient is polynomial in x,t. Clearing denominators yields a finite polynomial relation compatible with both arrangements. Injectivity of the time substitution and equality coefficient by coefficient transfer that relation to the full formal-time flow. This is not an unsupported evaluation between unrelated infinite-series fields.

## 3. Exactness by finite-cover trace

The proof correctly freezes the generic initial coordinates as constants and enlarges that constant field to an algebraic closure C. The target variable v is nonconstant because its initial time derivative -a xE is nonzero. Choosing eta^a=E and w=eta/X produces the target function field C(v,w), with w^a=v²(1+v). Its degree over C(v) is |a|, by Eisenstein at v+1; for a<0 use w^(-1).

Algebraic formal-time coordinates lie in a finite extension L/C(tau). The target curve field M=C(v,w) lies inside L, has transcendence degree one, and L is finitely generated, so L/M is finite. Characteristic zero makes the extension separable. The time differential equation d tau=-w dv/(a eta E) is correct. Tracing it gives a primitive of w dv in M. The proof needs no descent of exactness to the smaller original constants, since the ensuing obstruction holds over C itself.

## 4. Complete rational exactness test

Character projection for the cyclic Kummer extension places a primitive in the eigenspace w C(v). This also holds for |a|=1, when that eigenspace is the whole rational field. Differentiating wR(v) gives precisely

    a v(1+v)R'+(2+3v)R=a v(1+v).

A rational solution cannot have poles away from 0,-1. At these points the leading-pole conditions are ar=2 and ar=1 for positive pole order r. Neither is possible for a>=3 or a<0, so R is polynomial. Higher degree can occur only in the resonance a=-1, degree three. The linear-coefficient comparison excludes all other cases. In the resonant case the constant coefficient is zero, the linear coefficient of R must be -1, and the quadratic coefficient on the left is consequently -2 rather than the required -1. This excludes the resonance as well.

Thus a=1,2 are the only nonzero possibilities for algebraic formal time. The explicit primitives for those two cases are correct. The constructive algebraic cases are proved separately, so no unsupported converse to the trace implication is used.

## 5. The exceptional dilation case k=3

Here the induced scalar-time equation is dU/dlambda=U²(1+U), with source scalar-time series U=u+lambda u²+O(u³). The proof's descent of possible algebraicity from K(u,E) to K(u) is valid: after t=u/x², x is transcendental over K((u)), and hence so is E=u²(1+u)/x. Comparing the coefficients of a polynomial relation in E leaves a nonzero polynomial relation for U over K(u).

Flow invariance gives the displayed differential identity, and direct partial fractions give dh=dr/r for h=1/U-1/u and r=(1+U)u/[U(1+u)]. If h is nonconstant, the algebraic curve field is finite over C(h). Extending d/dh and taking the norm of r gives a nonzero constant logarithmic derivative for a rational function of h, impossible at infinity. Therefore h is constant; then r is constant. Substituting U=u/(1+cu) into r forces c=0 by coefficient comparison, so U=u, contradicting the nonzero-time first jet. This closes precisely the case not covered by scaling and establishes individual coordinate transcendence there.

## 6. Constructive cases and quintic irreducibility

The cited k=4 and k=5 results have already been independently reviewed. At k=4 only the first coordinate is asserted nonpolynomial; the polynomial second coordinate is accurately retained.

For k=6, the coordinates z=x³t,E give the stated rational inverse and delta(z)=z⁶/(E-z²). The rational time coordinate A_E(z)=1/(3z³)-E/(5z⁵) has derivative one along the vector field. The source exponential extends to fractions, so its existing source branch Z satisfies the cleared quintic identity. No arbitrary algebraic branch is introduced.

The linear-in-E irreducibility argument is correct. The two coefficients are coprime because at each fifth-root point Z=zeta z the second coefficient is 5z⁷(zeta²-1-3lambda z³), nonzero for lambda!=0 and transcendental z. The coefficient of Z² is a unit in K(z)[E], providing the second primitiveness step. Hence Z has degree five over K(x,t).

The source images are rational expressions in Z,E, while Z=X³T, so their joint field has degree exactly five. A rational X would give a cubic relation for Z. For T, I independently recomputed the remainder of (E-V²)³-TV⁸ modulo the quintic; its V⁴ coefficient is exactly 3E, independently of T. This makes the lower-degree contradiction valid and proves separate nonrationality. The proof does not claim unsupported individual extension degrees.

## 7. Scope and status

The theorem therefore has exactly the stated three-way outcome: rational at k=4, algebraic nonrational at k=5,6, and individually transcendental elsewhere. The source branch, base fields, and generic-variable roles remain explicit throughout. No mathematical correction was required in the reviewed draft. Finite computations corroborate the interfaces; the scaling, trace, pole and degree arguments supply the universal quantifiers.

## 8. Final text and independent reproduction

The accepted final proof SHA-256 is `98ba600bf22b72f62f605952adfe8a2521be6b2a0c9970ff959dbf8288cad26a`. Its added restricted double-series formulation explicitly records injectivity of tau=lambda z^b for either sign of b, and its Laurent-field explanation justifies evaluation of the k=6 rational time primitive. Both refinements are sound.

The checker source was inspected and independently run normally and with `python -O`. Both runs completed successfully on 118 controls and matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `856900fb0bbb8996f693f8d3072308f4282baff95d3546b0ea5328c93fdb49a0`. Finite k interfaces are explicitly distinguished from the all-k proof. This full review is complete, with no unresolved mathematical blocker.
