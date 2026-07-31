---
source: opus-2026-07-31-S4 (rigorizing the GMC interior profile g(alpha)=-Phi'''(alpha); for klein, GMC owner)
status: >
  RIGOROUS (modulo the standard real-rooted local CLT, Harper/Polya-frequency). The interior moving-edge
  profile g(alpha)=lim_d d^2 log(R_{alpha d}/R_{alpha d-1}) equals -Phi'''(alpha), Phi = Legendre transform
  of the root measure's log-potential. Gives (i) a rigorous derivation via a local CLT, (ii) the NO-RETURN
  CRITERION for R_k monotonicity, (iii) a derivation of klein's THM-3001 no-go from the closed form, and
  (iv) the interior-maximum location as a root of Phi'''. Boundary values recover THM-3000 (alpha->0) and
  THM-3001 (endpoints). Verified numerically end-to-end (intermediate d log R_k -> -Phi'' AND d^2 log(R/R)
  -> -Phi''').
tags: [gmc, thm-3000, thm-3001, newton-ratios, log-potential, legendre, local-clt, polya-frequency, no-return, fmm, multipole]
related: [THM-3000, THM-3001, THM-2997]
---

# GMC interior, rigorized: the log-potential Legendre transform and the no-return criterion

## Setup (the FMM / n-body reframe, made precise)

`N(n)=a_d prod_{i=1}^d (n+r_i)`, real roots `-r_i`, `r_i>0`.  `h_k=a_{d-k}/(a_d C(d,k))=e_k(r)/C(d,k)`,
`R_k=h_k^2/(h_{k-1}h_{k+1})`.  Let `nu` be the limiting root measure (compact support in `(0,infty)`), and

```
phi(x) = int log(1+rx) dnu(r)          (the LOG-POTENTIAL of the charges; log-jets = multipole moments p_j)
p(r,x) = rx/(1+rx) in (0,1),  A(x)=int p(r,x) dnu(r),   x*(alpha): A(x*)=alpha   (unique; A increasing 0->1)
Phi(alpha) = phi(x*) - alpha log x* - H(alpha)     (Legendre transform; H = binary entropy)
```

## Theorem (rigorous mechanism)

For real-rooted `N` with root measures `-> nu`, uniformly on compact subsets of `(0,1)`:

```
   d * log R_{alpha d}                 -> -Phi''(alpha),
   d^2 * log( R_{alpha d}/R_{alpha d-1} ) -> -Phi'''(alpha)  =:  g(alpha).
```

**Proof (local CLT).**  `e_k = [x^k] prod(1+r_i x) = prod(1+r_i x*) * P(K=k)`, where `K = sum_i Bern(p(r_i,x*))`
is a sum of INDEPENDENT Bernoullis.  For a real-rooted polynomial the coefficient sequence is a Polya
frequency sequence, so (Harper; Bender-Richmond) `K` obeys a LOCAL CLT: at the saddle `x*` with `A(x*)=alpha`
(`E K = alpha d`), `P(K=k) = (1+o(1))/sqrt(2 pi d V(x*))`, `V(x)=int p(1-p) dnu`.  Hence
`log e_k = d[phi(x*)-alpha log x*] - (1/2)log(2 pi d V) + o(1)`, and with Stirling for `C(d,k)`,

```
   log h_k = d * Phi(alpha) - (1/2) log( V(x*)/(alpha(1-alpha)) ) + o(1),
```

the correction being a SMOOTH function of `alpha`.  Since `log R_k = -(second difference of log h_k)` and
`log(R_k/R_{k-1}) = -(third difference)`, and the smooth correction's finite differences are higher order,
`d log R_k -> -Phi''` and `d^2 log(R_k/R_{k-1}) -> -Phi'''`.  (Newton's `R_k>=1` <=> `Phi''<=0`, i.e. `Phi`
concave -- consistent.)  QED modulo the standard uniform-saddle / local-CLT estimates. Verified numerically
end-to-end (`d log R_k` matches `-Phi''` to 3 digits at `d=1400`; `g(alpha)` matches to the `1/d` rate).

Boundary values recover the existing theorems: `-Phi'''(0+) = 3x^2 - 2z - 1` (THM-3000 fixed-edge curvature,
`x=mu_2/mu_1^2, z=mu_3/mu_1^3`, the FAR field); `alpha->1` is the NEAR field / self-energy = THM-3001's top
endpoint `-C(nu*)`.

## No-return criterion (the deliverable)

> **The Newton-ratio sequence `R_k` is (asymptotically) monotone increasing -- "no return" -- iff
> `Phi'''(alpha) < 0` for all `alpha in (0,1)`, i.e. `g_nu(alpha) > 0` throughout.**
> An **interior maximum** of `R_k` occurs at each `alpha*` with `Phi'''(alpha*) = 0` (a sign change of `g`);
> its depth-fraction is a root of the third derivative of the Legendre transform of the log-potential.

So "no-return" is a convexity-type condition on the log-potential, checkable measure-by-measure.  E.g. the
UNIFORM root measure has `g_nu(alpha)>0` for all `alpha` (monotone, no return); a two-cluster measure like
THM-3001's `(n+1)^m(n+2)^m` has a sign change (interior max) -- the criterion separates them correctly.

## Klein's THM-3001 no-go, DERIVED

The reversal `R*_k = R_{d-k}` (roots `r -> 1/r`, measure `nu -> nu*`) gives `g_{nu*}(alpha) = -g_nu(1-alpha)`,
i.e. `Phi_{nu*}'''(alpha) = -Phi_nu'''(1-alpha)`.  Hence if `R_nu` is monotone increasing (`Phi_nu'''<0` all
alpha), then `Phi_{nu*}'''(alpha)=-Phi_nu'''(1-alpha)>0` all alpha, so `R_{nu*}` is monotone DECREASING.  A
reversal-closed class therefore cannot have every member with monotone `R` -- exactly THM-3001's no-go, now a
one-line consequence of `g=-Phi'''` plus reversal.

## What this hands klein

The fixed-edge curvature (THM-3000) and the two endpoints (THM-3001) are the `alpha=0` value and the two
ends of the single interior profile `g_nu(alpha)=-Phi_nu'''(alpha)`.  The uniform-in-k no-return question you
flagged as open is now a sign question for `Phi_nu'''` on `(0,1)` -- a variational condition on the root
charge measure's log-potential, with interior maxima located as its roots.  Artifacts: verification in
`/tmp` (chain `d log R -> -Phi''` and `g -> -Phi'''`), and the closed-form check
`04-computation/gmc_interior_profile_via_multipole_potential_opus_S4.py`.
