---
id: HYP-2978
title: LRC14 Ramanujan exact-period projector packets
status: STUB / proof-interface reservation; computation and theorem test pending
source: codex-2026-06-24-ramanujan-projector
related:
  - HYP-2977
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2970
  - HYP-2969
  - HYP-2963
  - HYP-2901
  - HYP-2908
  - THM-572
  - OPEN-Q-108
---

# HYP-2978: LRC14 Ramanujan Exact-Period Projector Packets

This file reserves the exact-period Ramanujan-sum route requested in the
2026-06-24 Ramanujan prompt. The proposed object is not a proof yet. It is a
packet interface between:

- HYP-2972's rational twist ladder,
- HYP-2973's danger-count moments,
- HYP-2974's Fourier-Toeplitz dual,
- HYP-2975's endpoint taut-current layer, and
- HYP-2977's spectral-shadow route.

For a denominator `q`,

```text
c_q(n) = sum_{1 <= a <= q, gcd(a,q)=1} exp(2*pi*i*a*n/q)
       = sum_{d | gcd(q,n)} d*mu(q/d)
       = mu(q/gcd(q,n))*phi(q)/phi(q/gcd(q,n)).
```

The working idea is to replace scalar "q is blocked" data by primitive-unit
exact-period packets. For LRC14, the layer `q=14` is especially sharp:

```text
gcd(14,n)=14 -> c_14(n)= 6
gcd(14,n)= 7 -> c_14(n)=-6
gcd(14,n)= 2 -> c_14(n)=-1
gcd(14,n)= 1 -> c_14(n)= 1
```

Thus `c_14(r+s)` is a primitive-unit trace of the AP/Goddyn-Wong zero-credit
taut relation `r+s == 0 mod 14`, while `c_14(r-r')` records exact-period
residue coincidence with parity/seven-adic labels still attached.

The theorem-facing hope is:

```text
Any LRC14 residual that is invisible to raw qdiv, Haar/Baire boundary support,
danger-count moments, Toeplitz PSD, and spectral shadows must be visible in a
Ramanujan exact-period packet; otherwise it is an AP/GW equality atom or it
routes to the HYP-2908/THM-572 K33/H=7 state lift.
```

The computation still needs to decide whether the useful packet is a primitive
phase witness profile, a shifted Carmichael autocorrelation profile, a
Ramanujan-subspace projection of `N_S(t)`, or a supercharacter table over the
unit-group quotient.

Assumption challenge: the tournament vertices should not be fixed too early.
Candidate vertices are speeds, residues, primitive denominator phases,
denominators, Ramanujan modes, endpoint owner pairs, taut-current atoms,
Ramanujan-subspace components, and proof obligations. This route will start
with proof obligations / exact-period modes because they preserve the LRC
predicate "primitive weak witness, strict open mass, or boundary-only equality"
and destroy raw runner ownership unless endpoint/C27/K33 labels are reattached.
