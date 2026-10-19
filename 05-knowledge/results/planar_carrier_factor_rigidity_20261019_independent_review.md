# Independent review: planar carrier-factor rigidity

2026-10-19. **ACCEPTED ANALYTICALLY within the stated scope.** This review covers the new finite-order gate and scalar-time polynomial obstruction in `planar_carrier_factor_rigidity_20261019.md`, including its strengthened proof entirely in the original source completion. It does not audit the repository or assert the planar Jacobian conjecture.

## 1. The finite-support argument is valid

A constant nonzero Jacobian gives an étale, dominant, quasi-finite map; neither properness nor a global inverse is used. Dominance ensures the pullback of every nonconstant irreducible is nonconstant. Distinct irreducible plane curves intersect in a finite set over the algebraic closure. A common curve in two pullbacks would be contracted into that finite set, contradicting quasi-finiteness.

Squarefreeness also follows directly from the chain rule. If an irreducible source polynomial g occurred repeatedly in h(F), then g would divide both partial derivatives of h(F). Invertibility of the Jacobian matrix over the polynomial ring would imply that g divides both h_x(F) and h_t(F). Thus its curve would map into the singular locus of the reduced irreducible target curve h=0, a finite set in characteristic zero. This again contradicts quasi-finiteness. This agrees with the forward implication of Theorem 5.1 in [Jędrzejewicz's author text](https://arxiv.org/html/1109.2113), checked independently.

Consequently the r pullbacks have r nonempty disjoint factor supports inside the r-element complete support of H. Each support must be a singleton; squarefreeness then makes its exponent one. Comparing the invariant factorization gives the multiplicities. This does not assert preservation of arbitrary irreducibles by arbitrary Keller maps.

The passage to an algebraic closure is legitimate. The two marked polynomials t and 1+x^m t remain irreducible there: the latter is a primitive linear polynomial in t over K[x], with relatively prime coefficients 1 and x^m. No geometric irreducibility assumption on the other original factors is needed because the entire factorization is taken after base change.

## 2. The stabilizer and boundary controls are valid

After the factor permutation is killed by an iterate G, one has T=at and 1+X^m at=d(1+x^m t). Evaluation at t=0 yields d=1, cancellation of t yields X^m=a^(-1)x^m, and unique factorization yields X=bx. The determinant-one identity is ab=1, whereas the marked-factor identity is ab^m=1. Hence b^(m-1)=1. This supplies the inverse only after finite order has been proved.

The m=1 torus example, single-component shear, and m>=3 root-of-unity examples all have the claimed properties. In particular the restriction m>=2 is indispensable.

## 3. Formal-time obstruction: an independent original-source proof

The divisor identity S-c=t^5(1+x^2t)^4 R(t(1+x^2t),xt^2(1+x^2t)) is exact. The displayed birational inverse ensures that nonzero R remains nonzero. Thus write, directly in the original source variables,

    S-c=B(x)t^n+O(t^(n+1)),   n>=5, B(x) nonzero.

For delta={-,S}, one has

    delta x=n B(x)t^(n-1)+O(t^n),
    delta t=-B'(x)t^n+O(t^(n+1)),
    delta(t^j K[x][[t]]) subset t^(j+n-1) K[x][[t]].

The last inclusion also follows for every integer j in the Laurent ring. Therefore the exponential is topologically convergent and

    exp(lambda delta)(x)-x=lambda n B(x)t^(n-1)+O(t^n).

The terms starting at delta^2 x have valuation at least 2(n-1)>=n. The additive exponential group law replaces lambda by ell lambda for the ell-th iterate; characteristic zero makes this leading coefficient nonzero for every positive ell and nonzero lambda. Hence every such time has infinite order. This argument is independent of non-local-nilpotence and requires no comparison between the full D-adic completion and a Laurent field.

If both source images were polynomials, their invariant and ordinary Jacobian identities in the source completion would be polynomial identities. Applying the finite-order gate to H=S-c and m=2 then contradicts the preceding calculation.

The logarithmic leading-jet control and the final note's explicit comparison are sound: the map K[s][[tau]][tau^(-1)] to K[x]((t)), s=xt and tau=t, is well-defined and injective. The bounded lower tau order makes each target coefficient a finite sum, and distinct monomials map to distinct exponent pairs. The actual exponential series have polynomial s coefficients and bounded lower tau order. This does not give a map on all of K(s)((tau)) or identify that field with the full D-adic completion. Ordinary determinant one belongs to the original (x,t) chart, not to the ordinary (s,tau) Jacobian.

## 4. Finite checks and limits

The companion script was independently rerun normally and with Python optimization; its normal and optimized outputs match each other and the staged JSON, and it reports 448 exact checks. Separately, 24 direct-source symbolic controls were executed over R=1,p,y,p+y,D,D^2,(p-y)D,1+p+y. They checked n>=5, the leading coefficient nB of delta x, and the lower bound 2(n-1) for delta^2 x. These controls corroborate the algebra but do not replace its universal proof.

No blocking flaw remains in the reviewed statements. Rational scalar-time realizations for general R, compositions without a common invariant, entry into the source-normal chart, and the general planar Jacobian conjecture remain outside this conclusion. No priority assessment was performed.

## Final topology wording verification

The final proof with SHA-256 `dd78484770d5642f02c9d4a5e66be3770bd788c14cbb7bb6afad7b29edb58d87` was checked after the completion clarifications. Both formulations are accepted. For the actual unlocalized B₂ subset K[x,t], the inclusion D^N B₂ subset t^(3N)K[x,t] proves the required continuity on completion; no arbitrary localization or injectivity of completed rings is invoked. Independently, the bounded-pole polynomial-coefficient series ring and its injective substitution σ correctly compare the specified exponential series. Finally t need not be a B₂ generator: the source image of E(D)/E(p)^2 equals the source exponential on t by multiplicativity and D=p²t. The denominator is nonzero since E(p)=t+O(t²). These clarifications preserve the original-source polynomiality scope and do not add an assumption of a globally defined inverse.
