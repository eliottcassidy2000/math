# Independent audit of the maximal depth-preserving carrier

**Status: ACCEPTED analytically; exact replay recorded below.**
This audits [the carrier theorem](planar_jc_long_20260906_depth_carrier.md)
for every characteristic-zero field and every `S in K[p,y]`.

The necessity direction is exact: preservation of `P_0=K[p,y]` forces
`p|S_p,S_y`, because `gcd(p,p^3-y^2)=1`. The coefficients of powers zero
and one in `p` then give `S=c+p^2R`. No higher-depth condition or source
form was used in this direction. Conversely the displayed generator
images put `delta(p),delta(y)` in `P_0` and `delta(x),delta(u)` in `P_1`;
Leibniz proves preservation of every whole module `P_d`.

The source intersection is also an iff. Once `delta(A) subset A`, the
choice of `H` cannot cancel a pole of `delta(u)`, so “some H” and “every H”
are equivalent. The remaining numerator is `p^3y(4+E)R`, coprime to the
cusp divisor apart from the displayed last factor. On the quotient cusp
ring the eigenvalues of `4+E` are the nonzero integers `n+4`. This proves
the claimed divisibility for arbitrary polynomials, including cancellation
between different monomials representing the same cusp weight.

The local-nilpotence proof has satisfiable assumptions and correct scope.
If `delta` were locally nilpotent and nonzero, its nonzero invariant
`p^2R=t^2(1+x^2t)^2R(p,y)` would force both factors into its factorially
closed kernel. It would then kill `t,x`, a contradiction. Constants give
the zero derivation and are correctly excluded. This is not a statement
about the rationality of an individual time.

The independent program uses the rational inverse chart and a generic
triple `R,R_p,R_y`, checking every coefficient of both generator identities.
This differs from the producer's literal source differentiation. It also
checks the cusp Euler identity and the two carrier boundaries. Its finite
universe does not replace the all-polynomial proof.

```
python3 -B 04-computation/planar_jc_long_20260906_depth_carrier_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_depth_carrier_audit.py
```

Normal and optimized streams must agree; final byte pins are in the
[session manifest](planar_jc_long_20260906_manifest.json).
