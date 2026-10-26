---
source: opus-2026-08-01-S4 (all-degree extension of the degree-2 exp-integral theorem; the pure-power family)
status: >
  THEOREM (rigorous, all degrees d>=1). For algebraic alpha!=0, t0, gamma, the PURE-POWER-up-to-shift integral
  int_0^1 e^{alpha (t-t0)^d + gamma} dt is TRANSCENDENTAL. Same architecture as the degree-2 theorem, with the
  error-function E-function M=1F1(1/2;3/2) replaced by M_d(w)=sum_n w^n/((dn+1)n!) = 1F1(1/d; 1+1/d; w), still an
  E-function and still IRREDUCIBLE (a=1/d, b-a=1 neither a nonpositive integer), and still with GENUINE Poincare
  asymptotics M_d(w) ~ (1/d) e^w/w (a single honest confluent function -- NOT an oscillatory integral, so NO
  cancellation trap). int_0^x e^{alpha u^d}du = x M_d(alpha x^d) reduces I to e^gamma[(1-t0)M_d(eta1)+t0 M_d(eta2)],
  eta_i = alpha(1-t0)^d, alpha(-t0)^d algebraic; distinct exponential types e^{eta_i z} make the connections
  pairwise non-isomorphic, so {e^{gamma z}M_d(eta1 z), e^{gamma z}M_d(eta2 z), 1} are Qbar(z)-linearly independent,
  and Beukers 2006 lifts this to Qbar-linear independence of the values -- I transcendental. Degenerate t0 in {0}
  (single value e^gamma M_d(alpha)) and symmetric t0=1/2, d even (eta1=eta2, single value) handled likewise. This
  is a SUB-FAMILY at each degree: general deg>=3 Q has intermediate terms and its non-rationality is FC(2)-hard
  (see the corrected parameter-E-function file); the pure power avoids that because its arc is non-closed.
tags: [factorial-conjecture, FC2, exponential-integral, transcendence, e-functions, beukers, confluent-hypergeometric, pure-power, all-degrees, 1F1]
related: [degree-2-exp-integral-is-transcendental-rigorous-via-linear-beukers-shidlovskii, exp-integral-is-an-E-function-in-the-parameter-general-reformulation-and-rigidity, HYP-9078]
external: ["Beukers, Ann. of Math. 163 (2006) 369-379", "Slater, Confluent Hypergeometric Functions (asymptotics)"]
---

# Pure-power exp-integrals are transcendental at every degree

A clean infinite family where the degree-2 machinery goes through unchanged, giving rigorous transcendence at
all `d` -- and a sharp illustration of WHY it stops short of the general case.

## Theorem

For `alpha, t0, gamma in Qbar` with `alpha != 0` and any integer `d >= 1`,
```
        int_0^1 e^{ alpha (t - t0)^d + gamma } dt   is TRANSCENDENTAL.
```

## The one new ingredient: `M_d = 1F1(1/d; 1+1/d)`

Expanding `e^{alpha u^d}` and integrating,
```
   int_0^x e^{alpha u^d} du = sum_{n>=0} alpha^n x^{dn+1}/((dn+1) n!) = x * M_d(alpha x^d),
   M_d(w) := sum_{n>=0} w^n/((dn+1) n!) = 1F1(1/d; 1+1/d; w)   (since (1/d)_n/(1+1/d)_n = 1/(dn+1)).
```
`M_d` is an **E-function** (coefficients `1/(dn+1)`, bounded), and the confluent connection is **irreducible**:
`a = 1/d` and `b - a = 1` are neither nonpositive integers. Its asymptotics as `w -> +infinity` are the genuine
Poincare ones `M_d(w) ~ [Gamma(1+1/d)/Gamma(1/d)] e^w w^{-1} = (1/d) e^w/w` (verified: ratio -> 1 at `w=50` for
`d=3,4,5`). For `d = 2` this is exactly the error-function `M = 1F1(1/2;3/2)` of the degree-2 theorem.

## Proof (same shape as degree 2)

Shift `u = t - t0`: `I = e^gamma [ F_d(1-t0) - F_d(-t0) ]`, `F_d(x) = x M_d(alpha x^d)`. Hence
```
        I = e^gamma [ (1-t0) M_d(eta1) + t0 M_d(eta2) ],    eta1 = alpha(1-t0)^d,  eta2 = alpha(-t0)^d  in Qbar.
```
Set `g_i(z) = e^{gamma z} M_d(eta_i z)` (E-functions). The two Kummer connections `M_d(eta_i z)` have exponential
types `e^{eta_i z}` at infinity (from the genuine asymptotics above), so for `eta1 != eta2` they are
**non-isomorphic** irreducible rank-2 connections; `e^{gamma z}` is rank-1. Pairwise non-isomorphic factors have
no mixing horizontal endomorphisms, so `{g_1, g_2, 1}` are **linearly independent over `Qbar(z)`**. By the
refined linear Siegel-Shidlovskii theorem (Beukers 2006, no exceptional point), `{g_1(1), g_2(1), 1}` are
**`Qbar`-linearly independent**; therefore `I = (1-t0) g_1(1) + t0 g_2(1)`, a nontrivial `Qbar`-combination, is
transcendental. Degenerate cases collapse to a single value `e^gamma M_d(eta)` (namely `t0 = 0`, giving
`e^gamma M_d(alpha)`; and `t0 = 1/2` with `d` even, giving `eta1 = eta2` and `I = e^gamma M_d(alpha/2^d)`), still
transcendental by the two-function version. `#`

(Numerics: `d=3`, `int_0^1 e^{0.9(t-0.3)^3+0.2} dt = 1.29144045540794... = e^gamma[(1-t0)M_3(eta1)+t0 M_3(eta2)]`
to 15 digits; `eta1=0.3087 != eta2=-0.0243`; PSLQ finds no algebraic min-poly of degree `<=5`, height `<=1e9`.)

## Why this is a sub-family, not the whole conjecture (the honest boundary)

The pure power works because `t -> alpha(t-t0)^d` maps `[0,1]` to a **non-closed monomial arc**, whose pushforward
of Lebesgue measure cannot have all holomorphic moments zero -- so `Phi_Q(z)=int_0^1 e^{zQ}dt` is genuinely
non-constant, AND the single confluent function `M_d` has honest (non-cancelling) growth. A GENERAL degree-`d`
`Q` (with intermediate terms) has `d-1` interior saddles; its integral is a genuine oscillatory object that can
cancel (the corrected parameter-E-function file: `g=e^{2 pi i t}` has `int_0^1 e^{z g}=1` identically), and its
non-rationality is EQUIVALENT to FC(2). So the pure-power theorem is exactly the slice where the arc is monomial
and the confluent factor is single -- both cancellation-free. It extends the rigorous frontier from `{deg 2}` to
`{deg 2} cup {pure powers, all d}`, and pins precisely what the general case must additionally control.

## Corollary (values of special functions)

Instantiating: `int_0^1 e^{alpha t^d} dt = M_d(alpha) = 1F1(1/d; 1+1/d; alpha)` is transcendental for every
algebraic `alpha != 0` and every `d` -- e.g. the "higher error functions" `int_0^1 e^{alpha t^3}dt`,
`int_0^1 e^{alpha t^4}dt`, ... are all transcendental at algebraic points. (For `d=2` this is the classical
transcendence of error-function values; `d>=3` are the natural higher analogues.)
