# Independent audit of the native two-event persistent-bias certificate

Status: **SOUND** after source/code inspection, 2026-10-31.
Auditor: session128 feedback lane, independent of the paper-audit author.

Inspected `../paper_audit/verify_persistent_bias.py` and its exact JSON output,
and independently read the source paper pp.15--20 through Library.

1. Native source normalization is the cosine side of Eq.(60), consistent with
   Eq.(69). Differentiating Eq.(69) gives a positive common factor times
   `sum A_j[4sin(x-3theta_j)-sin(2x-6theta_j)]`, x=3theta, with
   A_1=-4/varphi^2 and A_2=9/varphi^3 after removing positive kappa*rho.
   The code's sine and cosine coefficients have these signs.
2. Every interval operation rounds outward to a rational denominator 10^450.
   Negative reciprocals use the correct reversed endpoints. Multiplication
   takes all four corners; squaring handles zero-crossing explicitly.
3. The Machin formula uses rigorous alternating-series brackets. sqrt(5)
   uses a strict integer-square bracket. The trigonometric polynomials are
   evaluated by interval Horner arithmetic, with a valid Taylor--Lagrange
   bound on the next omitted degree; no floating test is a proof gate.
4. The strict certificate R^2>5S^2 has the right mathematical consequence.
   If g=R sin(x-delta)+S sin(2x-epsilon), every zero lies in one of the two
   arcs |sin(x-delta)|<=S/R. On each arc, the first derivative's leading term
   has magnitude >=R sqrt(1-S^2/R^2)>2S, so g is strictly monotone there.
   There are at most two zeros. Two disjoint opposite-sign brackets certify
   two simple zeros, hence exactly two.
5. The signs are negative before the first root, positive between the roots,
   and negative after the second. Thus Pi=(r_2-r_1)/pi-1, exactly the
   expression enclosed by the code. The x=3theta change covers the circle
   three times and preserves its normalized sign integral.
6. Independently proved |Pi|<=1/3 for the complete harmonic carrier freezes
   both canonical seeds because h_1>1/3 and h_2<-1/3. This proves constant
   nonzero Pi and b_t=Pi[1-(1-lambda)^t], uniformly for 0<lambda<=1 and every
   positive gain sequence. It is not merely an observed long finite orbit.

The author incorporated the requested explicit ordered-root gate
`0<r1.lo<r1.hi<r2.lo<r2.hi<2*pi.lo`; inspected at checker line 101.
The final checker SHA-256 is
`590a160231f299e17a0c0b5cc3cbb818076afea6e7ef057e6cc457c6a536dc11`;
the JSON SHA-256 is
`41d5c847c0e67e9c7df30002c47c5190007bbf799c921a2a78ec3d6792af98e5`.
The detailed PAPER_AUDIT.md proof was also inspected after completion.

The certificate proves a counterexample to unconditional long-run J2
symmetry enforcement. It does not classify all fixed event sets, establish
aperiodicity, or refute any separately imposed schedule/symmetry hypothesis.

The source calls audits required but specifies no hard acceptance policy.
An additional filter could reject the witness; it would not prove that the
stated feedback equations correct it or eventually attain an accepted state.
