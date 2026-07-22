---
id: THM-2101
title: "An additive orbit-residue proof bypasses the small-root product in one-variable DvdK"
status: >
  RESERVED / PROVISIONAL PROOF UNDER AUDIT. Total constant-term vanishing
  makes the small-root residue sum equal to one; transitive Galois incidence
  is intended to conflict with the zero full-root Lagrange sum. After the
  MISTAKE-243 repair, Check A, additive orbit incidence, and the full-root
  Lagrange identity are kernel-checked and root-imported. The analytic-germ-
  to-splitting-field subset bridge and DvdK1 wrappers remain absent. Do not
  cite the bypass itself as proved until those layers pass audit.
source: codex-2026-07-22-GMC2-additive-orbit-residue
related:
  - THM-1550
  - THM-1765
  - THM-2022
  - THM-2067
  - HYP-8946
  - MISTAKE-243
formalization:
  - 04-computation/lean/TournamentH7/TournamentH7/GMC2LaurentShiftCheckA.lean
---

# THM-2101 -- reserved additive orbit-residue bypass

**RESERVED / UNPROVED PENDING FINAL AUDIT.** The intended proof replaces the
small-root product and its `t`-adic valuation by the residue weights
`alpha^(M-1)/Phi'(alpha)`. Its algebraic formal core now builds without
`sorryAx`; the exact analytic-to-algebraic bridge is not stated and must be
proved before promotion.

## Audited paper mechanism

For a genuine straddling Laurent polynomial, shift the least charge to zero:

```text
Lambda(z)=z^(-M)R(z),             Phi_t(z)=z^M-tR(z).
```

At a sufficiently small regular complex basepoint `t`, the roots inside the
unit circle satisfy the residue identity

```text
sum_(|alpha|<1) alpha^(M-1)/Phi_t'(alpha)
 = (2*pi*i)^(-1) integral dz/[z(1-t Lambda(z))]
 = sum_(m>=0) t^m CT(Lambda^m).
```

Thus total constant-term vanishing makes this inside-root sum equal to `1`.
The sum over every root is `0` because the numerator degree `M-1` lies below
`deg(Phi_t)-1`. If the inside subset can be transported into one abstract
splitting field, Galois conjugation makes every translate sum equal to `1`,
while uniform orbit incidence and the zero full-root sum give a contradiction.
This would bypass THM-1550's small-root product/Hensel valuation route.

## Missing interfaces and scope

The load-bearing missing step is local analytic subset to global algebra:
choose one small non-discriminant basepoint, identify its analytic root germs
with roots in an abstract splitting field, transport the inside subset, and
prove equivariance under the Galois action. The following wrappers are also
absent:

- split off a neutral charge `q=0` by the immediate `m=1` argument;
- retain injective charges and nonzero endpoint coefficients so
  `R(0)R_lead!=0`, `M>=1`, and `deg(Phi_t)>M`;
- use the exact shift `M=-min(support)`; a larger algebraic Check A shift
  inserts a zero root and breaks the irreducibility/simple-derivative route;
- normalize the derivative/full-root identity from a monic nodal polynomial
  to the generally nonmonic `Phi_t` and prove the root-cardinality hypotheses;
- import the repaired module from `TournamentH7.lean` and connect the result to
  `GMC2DvdKInterface.DvdK1` without depending circularly on THM-2067 itself.

The one-sided hostile control `Lambda=u^(-1)` has every positive constant term
zero, but its unique small-root weight and its full-root sum are both `1`.
Thus `deg(Phi_t)>M`, equivalently a genuinely positive charge after shifting,
is indispensable.
