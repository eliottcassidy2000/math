# Entire MRB deformation and a sensitivity test

Status: PROVED analytic statements; numerical values VERIFIED by two floating-point methods, not interval-certified. This is an independent derivation of a parameterized version of Crandall's known eta-derivative formula, not a novelty claim for the MRB expansion.

## Definitions

Put h_k=log(k)/k and H=log(3)/3 (the maximum over integer k>=1). Define

C(t)=sum_{k>=1} (-1)^k (exp(t h_k)-1),
A_m(t)=(2/pi)exp(-t log(m)/m)sin(pi t/m),
D_m(t)=A_m(t)-C(t).

The usual objects are C(1)=C_MRB and A_m(1)=A(m). The parameter t changes the kernel, so this is a hostile deformation test, not a claim that the CM structure itself varies with t.

## Entire-function theorem

For every complex t,

C(t)=t eta'(1)+sum_{j>=2} [(-1)^(j+1) eta^(j)(j)/j!] t^j.

This is entire. Let S2=sum_{k>=1} h_k^2=zeta''(2)<infinity. For every integer J>=1, truncation after degree J has the uniform disk bound

|remainder| <= S2 H^(J-1) |t|^(J+1) exp(|t|H)/(J+1)!.

Proof: isolate the linear term, whose alternating series converges because h_k decreases to zero for k>=3. The nonlinear remainder obeys |exp(t h)-1-th| <= |t|^2 h^2 exp(|t|H)/2. Thus its sum is absolutely and locally uniformly convergent. Interchange the j>=2 sums; the coefficient is sum_k (-1)^k h_k^j=(-1)^(j+1)eta^(j)(j). For the degree-J tail use h_k^(J+1)<=H^(J-1)h_k^2 and the exponential remainder inequality. All differentiations are consequently legitimate. This proof avoids falsely claiming absolute convergence of the original j=1 sum.

At zero,

D_8'(0)=1/4-eta'(1)=1/4-gamma log(2)+(log(2))^2/2.

In particular the defect is not identically zero in this natural family.

## Numerical stress test

The companion script compares alternating nsum against the shifted Abel-Plana integral on real t:

C(t)=integral_0^infinity Im[exp(t Log(1+i x)/(1+i x))]/sinh(pi x) dx.

The integral has removable endpoint value t/pi. For derivatives in t, replace the numerator by Im[h(1+i x)^j exp(t h(1+i x))]. For j=1 the endpoint is 1/pi; j>=2 has endpoint zero.

At 60-digit working precision:

- D_8(1)=0.0000000189180712788007352705453393844198311058750464.
- D_8'(1)=-0.0896013973074147415174576834930579928139735691707305.
- D_8(0.999999)=0.000000108519379884141030911440680152674522476416724.
- D_8(1.000001)=-0.0000000706834147306861474969947527711579135594301343.
- A nearby numerical zero is t=1.00000021113585159506605080273174576616853879056740.

Interpretation: the small value is locally sensitive. A change of t by one part in a million changes its sign and exceeds the original defect in magnitude. This does NOT refute a special identity at t=1, establish coincidence probabilities, or supply interval-certified existence/uniqueness of the nearby numerical zero.

## What a proposed modular formula must explain

An exact bridge should specify how the full analytic kernel and summation functional are recovered, rather than only the branch jump at -8 or a trigonometric normalization. A claim about an analytic family can be tested against the exact eta-derivative coefficients and the uniform tail bound above. A claim solely at t=1 cannot be ruled out by the deformation test alone.

## Primary literature

Richard E. Crandall, *Unified algorithms for polylogarithm, L-series, and zeta variants* (2012), section 7.5, pp.28-29: https://www.marvinrayburns.com/UniversalTOC25.pdf . Its eta-derivative representation at t=1 is explicitly prior art. The locally uniform proof and parameter stress test here are derived directly, and no claim of literature-wide originality is made.

## Independent certificate upgrade

The independent audit in `analytic-audit/AUDIT.md` uses a differently parameterized Euler transform, a different Cauchy radius, and directed interval arithmetic. It proves a unique real zero in [0.999999,1.000001] and encloses it between

1.000000211135851595066050802731745766168538790567401655840
and
1.000000211135851595066050802731745766168538790567401655842.

Thus the earlier numerical zero is now independently certified; the original floating-point script remains labelled by what it actually does. This does not upgrade an unexplained modular relationship into a theorem.
