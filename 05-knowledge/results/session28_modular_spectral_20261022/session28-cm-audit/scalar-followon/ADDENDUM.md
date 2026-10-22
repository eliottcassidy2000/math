# Scalar Gaussian ratio: the self-dual point is a strict local maximum

Status: proved by exact symbolic reduction and rational interval arithmetic. This corrects the earlier numerical suggestion that the ratio might attain its global minimum at t=1. It does not determine the quartic spectrum or establish Gaussian energy optimality.

## Setup and exact reduction

Use q=exp(−2πt), x=2πt, D=q d/dq and

F=E4^11−2640 Δ E4^8+1813680 Δ² E4^5−244992000 Δ³ E4²,

Δ=(E4³−E6²)/1728.

The audited scalar ratio is

r(t)=[D²F−45 DF/x]/(990 Δ⁴).

The standard Ramanujan identities are D E2=(E2²−E4)/12, D E4=(E2 E4−E6)/3, D E6=(E2 E6−E4²)/2. Consequently differentiation along x=2πt acts on these expressions as ∂x−D. At t=1, modular transformation gives E6(i)=0 and E2(i)=3/π=6/x. Let y=4π²E4(i). Exact substitution yields

r(1)=8192(304325−17489844/y),

d²r(e^u)/du² at u=0 = 1024(23708225 y²−1361323584 y−69959376)/y.

The first derivative vanishes exactly, also following from r(t)=r(1/t). The Python script checks the vanishing symbolically, not by a floating-point tolerance.

## Rational sign certificate

The attached script bounds π through Machin's identity with alternating arctangent remainders. It bounds q=exp(−2π) using positive exponential Taylor sums with a geometric upper bound on the remainder, then inversion. To enclose E4(i), it sums the first twelve terms of 1+240 Σ σ3(n)q^n. Since σ3(n)≤n⁴ and (j+1)⁴≤16^j, the omitted positive tail is at most

240·13⁴ q^13/(1−16q).

All arithmetic through the sign test uses exact rational numbers. It proves

57.471215393911309459 < y < 57.471215393911309460.

On this interval P(y)=23708225 y²−1361323584 y−69959376 is increasing, and even its upper endpoint is negative. The computed upper bound is approximately −21413.62539. Hence r(e^u) has strictly negative second derivative at u=0: t=1 is a strict local maximum.

Illustrative, non-certificate numerical values are r(1)≈11919.94921586 and curvature≈−381539.73692644. Values at 0.999 and1.001 are both smaller. The sign proof does not rely on these samples.

## Rigorous consequences and remaining obligations

The previously audited asymptotic limit is r(t)→2723200 as t→∞ and, by reciprocal symmetry, as t→0. This is larger than r(1). Continuity therefore ensures the global minimum is attained in a compact interval. Since1 is a strict local maximum, a global minimizer t₀≠1 exists, and1/t₀ is a distinct global minimizer. This proves at least two reciprocal global minimizing temperatures, not exactly two and not their locations.

For an angular eigenvalue λ, Hessian sign is governed by r(t)+λ. If −r(1)<λ<−min r, this branch is positive at the self-dual temperature and both tails but negative at the off-center minima. Such a split instability pattern is conditional: no assertion is made that a particular rank88 lattice has an eigenvalue in this interval.

A full global lower bound still needs a certified compact-interval argument plus a tail estimate. One viable route is to use exact modular-polynomial expressions in E2,E4,E6 and interval enclosures for their absolutely convergent q-series, retaining derivatives and outward remainders. A proof of exactly two minima would additionally need derivative sign control between all stationary points. No global positivity or uniqueness of minima is claimed here.

## Reproduction

Run python3 reduce_ratio.py. The script requires SymPy and mpmath. Symbolic identities and rational sign checks precede purely illustrative floating-point samples. Files from the original frozen audit are unchanged.
