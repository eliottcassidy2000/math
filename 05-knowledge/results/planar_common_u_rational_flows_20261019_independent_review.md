# Independent review: rational flows preserving x²t

2026-10-19. **ACCEPTED ANALYTICALLY** for the rational-flow classification, with the field and completion distinctions below. The full producer proof has been read and checked.

## 1. Exact flow and source convergence

Let u=x²t, f in K[u], f'(0)=0, and work in characteristic zero. The Hamiltonian derivation satisfies

    delta(u)=0, delta(x)=x² f'(u), delta(t)=-2xt f'(u).

Writing a=f'(u), constant along the flow, solves these equations exactly:

    X=x/(1-lambda x a), T=t(1-lambda x a)².

The inverse is obtained by replacing lambda by -lambda, with u invariant. Since a is divisible by u, the denominator is 1+O(t), a unit in K[x][[t]]. Also delta(x) is divisible by t and delta(t) by t², so the source scalar-time exponential converges and agrees with the displayed rational formulas. It fixes t=0 pointwise.

The ordinary source Jacobian is one. For example, in the rational coordinates (x,u), the Jacobian of (X,u) with respect to (x,u) is (1-lambda x a)^(-2). Changing back to source coordinates multiplies by x²/X², giving one. This is an identity of rational functions; it does not claim the map is polynomial or defined everywhere. The invariant identities X²T=u and 1+X²T=q are exact.

## 2. Finite commuting words

The derivations for f_i(u) commute, since each is f_i'(u) times the same invariant-preserving derivation. Their finite product has the above formula with lambda a replaced by F'(u), where F=sum lambda_i f_i. If F is constant the product is identity. Otherwise D=1-xF'(x²t) is nonconstant and coprime to x: modulo x it equals one. Therefore X=x/D is not a polynomial. T=tD² is still a polynomial, so the rational pair is not a polynomial map. This proves the exact polynomial-output criterion; every such word remains rational and birational.

The derivative guard is important. If f'(0) is nonzero, the t-order-zero coefficient of X is x/(1-lambda x f'(0)), which is not in K[x] for nonzero time. Thus the same rational formula need not define a flow of K[x][[t]] or fix the original line pointwise. The stated hypothesis avoids this failure.

## 3. q-divisible subfamily and boundary scope

For f(u)=c+u²(u+1)g(u), one has

    H-c=t²q x⁴g(u).

This subfamily satisfies the derivative guard. Its coefficient R=x⁴g(x²t), when nonzero, is homogeneous of weight -4 for w(x)=-1,w(t)=2, outside the Newton bound -2. The Hamiltonian has weight zero; the derivation need not raise that weight.

Along q=0, use x=-1/s and t=-s², so u=-1. Since f'(-1)=g(-1), the exact rational boundary action is

    s maps to s+lambda g(-1).

In particular g=1 yields translation by lambda. This is a rational action on K(s), not a globally regular automorphism of the entire affine q-curve, whose coordinate ring is K[s,s^(-1)]. A nonzero translation may hit its missing parameter value zero. Likewise a nonzero scalar translation does not define an order-raising continuous exponential on all of K((s)); expanding its action on -1/s produces unbounded negative s-orders.

These distinctions explain why the rational family does not contradict the preceding joint-boundary completion theorem. Generic rational preservation of q and source t-adic convergence do not give the stronger joint-weight convergence used there.

## 4. Comparison and limits

For g=1, R=x⁴ supplies an explicit rational nonpolynomial case outside the Newton cone, while the earlier R=x² invariant family has individually transcendental coordinate outputs at nonzero times. There is no conflict: the invariants and completion mechanisms differ. The new pair has a genuinely rational nonpolynomial first coordinate, so it is not a Keller polynomial counterexample. No claim about all Hamiltonians or arbitrary different-invariant compositions follows.

## 5. Actual producer text and pointwise-boundary example

The full proof `planar_common_u_rational_flows_20261019.md` was reviewed. The displayed iterated derivatives of t terminate at order three, consistent with the exact polynomial second coordinate, while the iterates of x form the convergent geometric series. The rational differential-form Jacobian computation correctly avoids treating (x,u) as global polynomial coordinates. The criterion f(-1)=f(0), combined with f'(0)=0, is exactly equivalent to f-c being divisible by Z²(Z+1).

The additional choice f(Z)=Z²(Z+1)² has derivative 2Z(Z+1)(2Z+1). Its denominator is identically one on both t=0 and q=0, so this rational map is actually defined and pointwise identity along both complete affine curves. Its first coordinate nevertheless remains nonpolynomial at nonzero scalar time by the same coprimality argument. This is a valid demonstration of the distinction between rational boundary fixing and global polynomial regularity.

No mathematical correction was required. The stated birational domains and completion guards are accurate, and the result does not provide a polynomial Keller counterexample.

## 6. Final proof and reproduction receipt

The accepted final proof SHA-256 is `9d51b831ca10295546daa094879b7d04d56f7fb9b3e025a6c61bc0032de4ec92`. The standalone checker source was inspected and independently run normally and with `python -O`. Both runs completed successfully on 169 checks and matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `8f23e31d79d570d1738b9cda74a849db7f753fb3ee5bd5d4b7ff65652738ef66`. The finite controls corroborate exact derivative, rational inverse, Jacobian, composition and boundary identities; the universal classification rests on the proof above. Review is complete with no unresolved mathematical blocker.
