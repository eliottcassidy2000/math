# Independent review: boundary-tangent Hamiltonian compositions

2026-10-19. **ACCEPTED ANALYTICALLY.** The broader-family theorem and both the induced boundary dynamics and polynomial-rigidity deduction have been independently checked. Packet110 is unchanged.

Let K have characteristic zero, q=1+x²t, and take finitely many Hamiltonians H_i=c_i+t²q A_i(xt,t), with A_i polynomial in its two arguments. Every finite composition of their scalar-time source-completed flows that is polynomial in the original source variables is the identity.

## 1. Two convergent realizations

In K[x][[t]], each nonconstant H_i-c_i has t-order at least two. Consequently delta_i(x) has t-order at least one, delta_i(t) has order at least two, and the derivation raises t-order by at least one. Its exponential is a continuous automorphism, fixes t=0 pointwise, and preserves the ordinary source Poisson bracket. Finite compositions retain these properties.

Set s=xt, tau=t, h=tau+s² and use B=K[s][[tau]][tau^(-1)]. Then H_i-c_i=tau h A_i(s,tau), delta_i(s)=tau(H_i)_tau has tau-order at least one, and delta_i(tau)=-tau(H_i)_s has order at least two. Hence delta_i raises every integer tau-order by at least one. Its exponential defines a continuous automorphism of B. The injective continuous substitution sigma:B to K[x]((t)), s maps to xt and tau maps to t, matches both convergent realizations and finite compositions.

## 2. Boundary dynamics, without a closed-ideal shortcut

The continuous homomorphism epsilon:B to K((s)), tau maps to -s², is well-defined because the coefficients are polynomials in s and tau exponents are bounded below. A tau-order N tail has s-order at least 2N. This is not an evaluation on all of K(s)((tau)).

For H=tau h A, direct differentiation gives

    epsilon(delta(s))=s⁴ A(s,-s²),
    epsilon(delta(tau))=-2s⁵ A(s,-s²).

Define dbar=s⁴ A(s,-s²) d/ds on K((s)). It raises s-order by at least three, unless it is zero. The two identities imply

    epsilon composed with delta = dbar composed with epsilon

first on Laurent polynomials and then on B, by continuity of both sides and density. Thus

    epsilon composed with exp(lambda delta)
      = exp(lambda dbar) composed with epsilon.

Both exponentials converge in their respective valuations, so the passage from powers to the full exponential is justified. For a finite composition C, the same intertwining holds with the corresponding finite composition Cbar of boundary exponentials.

Since epsilon(h)=0, one obtains epsilon(C(h))=0. Since epsilon(tau)=-s², one also obtains epsilon(C(tau))=Cbar(-s²), which has leading term -s² and is nonzero. Therefore epsilon(C(q))=0, using q=h/tau and the invertibility of C(tau) in B. This establishes setwise preservation of the q-curve for any polynomial realization.

This proof does not assume ker(epsilon)=hB or infer preservation by taking limits inside an unproved closed ideal. It directly intertwines the derivations and their convergent exponentials. In particular original t-adic convergence alone cannot establish q-boundary restriction, since q is a unit in that source completion.

## 3. Polynomial rigidity

A hypothetical polynomial realization F has determinant one and fixes t=0 pointwise. Gwoździewicz's proved fixed-line theorem, as checked in the packet110 review against the primary author text, makes F an automorphism after algebraic closure. Preservation of the line ideal then gives T=at and X=a^(-1)x+f(t); pointwise line fixing gives a=1 and f(0)=0.

Evaluate q(F) on x=-1/s,t=-s². Setwise q-preservation gives

    0=2s f(-s²)-s² f(-s²)²
     =s f(-s²)(2-s f(-s²)).

The second parenthesis has constant coefficient 2 and cannot vanish in characteristic zero. Hence f(-s²)=0 and f=0. Identity descends to K as a polynomial identity. No general Jacobian conjecture or descent of an inverse is needed.

## 4. Boundaries and independent controls

If H-c=t³q² A(xt,t), its common-ring form is tau h² A. The induced boundary derivation vanishes, so the q-curve is fixed pointwise, a stronger property than the main family needs. Dropping the h-factor permits H=t²: its flow is the nonidentity polynomial shear (x+2 lambda t,t). This is a valid sharp boundary control.

A separate symbolic calculation passed 70 checks over seven choices of A: the two induced coordinate formulas, divisibility of delta(h) by h, and the intertwining relation on seven Laurent-polynomial inputs per carrier. These checks corroborate the all-A proof rather than replacing it.

The polynomial-coefficient restriction A(s,tau) matters for the stated common-ring convergence. No extension to arbitrary A(x,t), arbitrary localizations, rational realizations, infinite compositions, or other source charts is established here.

## 5. Actual producer-text review

The producer proof `planar_boundary_preserving_hamiltonian_cone_20261019.md` was read in full. The main proof correctly uses epsilon(C(h))=0 directly after obtaining T=t, rather than requiring a separate division by the transformed tau-coordinate. Its additional exact criterion for pointwise boundary identity of a single nonzero-time generator is valid: if v(s)=s⁴A(s,-s²) is nonzero with initial order m>=4, the leading displacement of its scalar exponential is lambda times that initial term, while all higher terms have larger order. The kernel of polynomial evaluation tau=-s² is the principal ideal (tau+s²), by monic polynomial division. Thus h divides A exactly for the pointwise-fixed subfamily. The explicit boundary solution s(1-3 lambda s³)^(-1/3) for A=1 differentiates to s(lambda)^4 and has the correct initial value.

One minor notation correction was requested and applied: subfamily membership is stated for H-c rather than H when additive constants are allowed. No mathematical proof change was needed.

## 6. Final verification and reproduction

The final proof SHA-256 is `340a748863602d213c7185f68cbfd03ea7043633ea081d181751959d07f5fb62`. Its corrected mathematical text is accepted; final changes only add status and reproduction links. The checker was independently rerun normally and with `python -O`. Both runs completed successfully on 354 checks and matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `5e7282701599c501a9aabf95ab0364ddb43f1d7b0f01b7d8488d7bb8ec6ca494`. These producer-script checks are additional to the 70 independent symbolic controls described above. Review is complete; no unresolved mathematical blocker was found.
