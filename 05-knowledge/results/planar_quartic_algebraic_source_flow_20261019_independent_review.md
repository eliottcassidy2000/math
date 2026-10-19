# Independent review: quartic algebraic source flow

2026-10-19. **ACCEPTED ANALYTICALLY.** The full proof `planar_quartic_algebraic_source_flow_20261019.md` was independently read. The claimed degree-four joint extension and individual nonrationality hold for every nonzero scalar time over every characteristic-zero field.

## 1. The actual branch and invariant coordinates

For E=x⁵t²(1+x²t), the displayed source derivatives are correct and increase t-order by at least one. Thus the specified scalar-time exponential is a genuine continuous automorphism of K[x][[t]], not merely a prospective algebraic root. Setting u=x²t gives E=xu²(1+u), with the stated rational inverse x=E/[u²(1+u)], t=u⁵(1+u)²/E². This proves algebraic independence and equality K(x,t)=K(u,E).

Direct differentiation yields delta(u)=-xE. For P(z)=z³/3+z⁴/4, P'(u)=u²(1+u), so delta(P(u))=-E² and the second derivative is zero. Exponentiation therefore gives P(v)=P(u)-lambda E² for the already defined v=C_lambda(u). The argument correctly avoids invoking an implicit-function theorem at u=0, where P' vanishes. The quartic relation alone is not used to choose a branch.

## 2. The two Gauss steps are legitimate

Let F=K(u). The polynomial P(V)-P(u) has the simple prime factor V-u in F[V], because P'(u) is nonzero in F. Thus (P(u)-P(V))/lambda has odd valuation at that prime and is not a square in F(V). The polynomial lambda W²+P(V)-P(u) is consequently irreducible over F(V) as a quadratic in W.

Its W-leading coefficient lambda is a unit of F[V], so it is primitive and irreducible in F[V,W]. Viewed as a V-polynomial over F[W], its leading coefficient 1/4 is a unit, so a second Gauss application proves irreducibility over F(W). Renaming the independent variable W as E is legitimate because u,E are algebraically independent. The resulting irreducible polynomial has V-degree four. Hence the specific source branch v has degree exactly four over K(x,t). No algebraic-closure or analytic assumption is needed.

## 3. Both coordinates individually, and their joint field

Applying the exponential to polynomial identities gives X=E/[v²(1+v)] and T=v⁵(1+v)²/E² in the common fraction field. The denominators are nonzero because the original source exponential is an automorphism and u,1+u,E are nonzero. Thus both coordinates are algebraic. Since v=X²T, their joint field equals K(x,t,v), proving the exact degree-four assertion.

If X were rational, Xv³+Xv²-E would be a nonzero annihilating polynomial of degree at most three, impossible. For T, the displayed remainder of V⁵(1+V)² modulo V⁴+(4/3)V³-4C was independently recomputed and agrees exactly:

    (4C+16/81)V³+(8C/3)V²+(4C/9)V-16C/27.

Here C=P(u)-lambda E² is nonzero by algebraic independence. If T were rational, subtracting E²T would preserve the nonzero V² coefficient 8C/3 and produce another annihilating polynomial of degree below four. Thus T is nonrational separately. The proof claims exact degree four for the joint field, without making an unsupported assertion that each coordinate individually has degree four.

## 4. Scope and comparison

The contrast with R=x² (transcendental coordinates) and R=x⁴ (explicit rational nonpolynomial first coordinate) is accurate. The R=x⁵ result is algebraic nonrational, with its branch realized in the source t-adic completion. It neither contradicts the Newton-cone theorem nor turns a formal inverse into a rational plane inverse. It does not classify other exponents, general carriers, or different-invariant compositions. No mathematical correction was required.

## 5. Final verification and independent reproduction

The accepted final proof SHA-256 is `f250854f169ddd245ed877deb8aa5da233ee3e15f67ef1a9f5a07bd5640783d7`. The standalone checker source was inspected and independently run normally and with `python -O`. Both runs completed successfully on 49 controls and matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `7102b901df11db22b29f3e396d26f60a8b426da92fd2f958d38e1e21bd4f25c2`. The source-jet truncations are valid because the derivation raises t-order; they corroborate the algebraic interfaces but do not substitute for the universal irreducibility proof. Review is complete with no unresolved mathematical blocker.
