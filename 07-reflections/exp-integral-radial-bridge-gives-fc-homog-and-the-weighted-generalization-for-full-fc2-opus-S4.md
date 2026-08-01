---
source: opus-2026-08-01-S4 (the exp-integral transcendence claim; reconstructing the bridge HYP-9078 left open)
status: >
  RESULT (reconstructs the bridge) + inhomogeneous clarification. The external claim "int_0^1 e^{Q(t)} dt is
  transcendental for every nonconstant Q in Qbar[t]" implies FC(2)-HOMOGENEOUS, via the RADIAL bridge (NOT the
  failed log-bridge): FC(2)-homog = the [0,1] Lebesgue polynomial moment problem (Liu-Sun 2020 Thm 2.6), whose
  moment EGF is EXACTLY int_0^1 e^{s phi(a)} da with Q=s phi a POLYNOMIAL in a; if FC-homog fails (phi!=0, all
  moments 0) then int_0^1 e^{s phi} da = 1 for all s, so at an algebraic s0!=0 the integral int_0^1 e^{s0 phi} da
  = 1 (rational) for a nonconstant algebraic exponent -- contradicting transcendence. This is the implication
  HYP-9078 flagged "not reconstructed"; it is valid because the radial reduction gives a POLYNOMIAL exponent,
  unlike the log-bridge's log-exponent. Since FC-homog is already Liu-Sun, the content is consistency + the
  RIGHT bridge. FULL FC(2): my saddle-descent gives the leak EGF int_0^1 e^{s phi_D} C da with weight
  C=exp(phi_{D-1}/(D phi_D)), so full FC(2) needs transcendence of int_0^1 e^{poly + RATIONAL} da -- a
  GENERALIZATION of the poly-only claim. Degree 1 of the claim is proved (LW, concurrent opus HYP-9078); degree
  >= 2 is Siegel-Shidlovskii E-functions + Beukers, matching the friend's route.
tags: [factorial-conjecture, FC2, exponential-integral, transcendence, e-functions, beukers, siegel-shidlovskii, lindemann-weierstrass, radial-bridge, saddle-descent, HYP-9078, liu-sun]
related: [HYP-9078, fc2-saddle-descent-reduces-to-the-top-form-period, FC2-no-pole-homogeneous-is-the-lebesgue-moment-problem, fc-jc2-tandem-the-isolated-coincidence-and-the-fold-no-fold-duality]
external: ["Liu-Sun 2020 Thm 2.6 (FC homogeneous = [0,1] PMP)", "Lindemann-Weierstrass", "Siegel-Shidlovskii E-functions", "openai/ten-proofs (Lean)"]
---

# The exp-integral radial bridge gives FC(2)-homogeneous, and the weighted generalization the full case needs

## 1. The claim, and the bridge HYP-9078 left open

External claim (owner-relayed): `int_0^1 e^{Q(t)} dt` is TRANSCENDENTAL for every nonconstant `Q in Qbar[t]`.
HYP-9078 (concurrent opus) proved DEGREE ONE (Lindemann-Weierstrass), showed the real case is trivial, showed
the naive LOG-bridge `x=-log u` fails (it makes the exponent a polynomial in `log u`, not in `u`), and
concluded the asserted implication "this implies FC(2)" was "not reconstructed."

**It is reconstructed here, for the HOMOGENEOUS case, by the RADIAL bridge.**

## 2. The radial bridge (valid; polynomial exponent)

For homogeneous `f` of degree `d`, `L(f^m) = (dm+1)! int_0^1 phi(a)^m da`, `phi = f|_{a+b=1}` (my radial
reduction; = Liu-Sun 2020 Thm 2.6, the `[0,1]` Lebesgue PMP). The moment EGF is

```
   sum_{m>=0} ( int_0^1 phi^m da ) s^m / m!  =  int_0^1 e^{s phi(a)} da,     Q(a) := s phi(a)  (POLYNOMIAL in a).
```

(Verified: `int_0^1 e^{s phi} da = sum s^m <phi^m>/m!` exactly.) Now suppose FC(2)-homog FAILS: some nonconstant
`phi` (degree `<= d`, algebraic coefficients WLOG) has `int_0^1 phi^m da = 0` for every `m >= 1`. Then

```
   int_0^1 e^{s phi(a)} da = 1   for ALL s.
```

Pick any algebraic `s0 != 0`. Then `Q = s0 phi` is a nonconstant polynomial with algebraic coefficients and
`int_0^1 e^{Q(a)} da = 1` -- a RATIONAL number. This CONTRADICTS the transcendence claim. Hence:

> **[`int_0^1 e^Q` transcendental for nonconstant algebraic `Q`]  =>  FC(2)-homogeneous.**

This is a genuine implication, and it is NOT the log-bridge: the exponent `s0 phi(a)` is a polynomial in the
integration variable `a`, exactly the claim's object. (The log-bridge fails precisely because it produces a
polynomial in `log u`; the radial reduction avoids `log` by passing through the simplex.) Since FC(2)-homog is
already Liu-Sun 2020, the value is (i) the correct bridge, and (ii) that the claim is at least as strong as a
proved theorem -- a consistency check that also tells us where the claim's strength must go for the open case.

## 3. Full FC(2): the WEIGHTED (poly + rational) exp-integral

Inhomogeneous `f = sum_{d<=D} f_d` does NOT reduce to a single `[0,1]` moment problem. My saddle-descent
(`fc2-saddle-descent-...`) gives `L(f^m) ~ (Dm)! P_m` with the top-form leak `P_m = int_0^1 phi_D^m C da`,
weight `C = exp( phi_{D-1}/(D phi_D) )`. Its EGF is

```
   sum_m P_m s^m / m!  =  int_0^1 e^{ s phi_D(a) + phi_{D-1}(a)/(D phi_D(a)) } da,
```

an exponential integral of `e^{poly + RATIONAL}`. So:

> **Full FC(2) needs transcendence of `int_0^1 e^{poly + rational} da`, a GENERALIZATION of the poly-only
> claim.** The poly-only claim closes the HOMOGENEOUS case; the inhomogeneous coupling adds the rational
> weight, exactly the `phi_{D-1}/phi_D` depression term.

This pinpoints the friend's implication: to reach FULL FC(2) the poly-only claim is not, by itself, the whole
story -- either the weighted generalization is proved, or the friend's E-function route reduces the
inhomogeneous integral to the poly-only claim by machinery I have not reconstructed. Both are consistent with
"E-functions + Beukers lifting + horizontal-endomorphism rigidity."

## 4. The transcendence engine (degree >= 2): E-functions / Beukers

Degree 1 is Lindemann-Weierstrass (HYP-9078). Degree `>= 2`: `int_0^1 e^{Q}` for `deg Q >= 2` is an
E-function value (the exponent's exponential is annihilated by a rank-1 ODE `y' = Q' y`; the definite integral
is a period of the associated inhomogeneous / rank-2 connection -- "rank-2 factorial"). Siegel-Shidlovskii and
Andre-Beukers give transcendence of E-function values at nonzero algebraic points unless a "trivial"
(horizontal / flat-section) relation forces algebraicity. The claimed "horizontal-endomorphism rigidity"
theorem is exactly the statement that no such trivial relation occurs for a nonconstant `Q` -- forcing
transcendence. If that holds at all degrees, the poly-only claim is true, and (sec 2) FC(2)-homog follows for
free; the open work is the inhomogeneous weight (sec 3).

## 5. Honest scope

- Reconstructed: the RADIAL bridge `[int e^Q transcendental] => FC(2)-homog` (fills HYP-9078's gap; distinct
  from the failed log-bridge). FC(2)-homog itself is Liu-Sun 2020 -- no novelty there, only the bridge.
- Clarified: FULL FC(2) via saddle-descent needs the `poly+rational` (weighted) exp-integral transcendence.
- Not established: the friend's poly-only `=> FULL FC(2)` (the E-function/Beukers route is not reconstructed
  here); and the transcendence claim itself at degree `>= 2` (E-function theorem, external).
- Cross-checks that stand: degree-1 transcendence (LW), the EGF identity, the golden AMM constant
  `C = 1 + 2 log(phi)/log 5` (independent: klein THM-3027 = opus two-ray).

Handed to the concurrent opus (HYP-9078 owner): the radial bridge is the missing link for the homogeneous case;
the inhomogeneous case is the weighted generalization; the E-function engine is the shared degree-`>=2` tool.
