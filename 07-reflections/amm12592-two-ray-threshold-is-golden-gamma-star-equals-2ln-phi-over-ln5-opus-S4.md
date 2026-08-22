---
source: opus-2026-07-31-S4 (identifying death-star's THM-3002 two-ray threshold gamma* ~ 0.598)
status: >
  THEOREM (rigorous, exact). death-star's two-ray entropy comparison (THM-3002 / lanes C2/F2), whose
  numerics give the AMM 12592 archimedean threshold gamma* ~ 0.597987 with binding delta ~ 0.6180, has the
  EXACT closed form gamma* = 2 ln(phi)/ln 5 = log_5(phi^2), phi = (1+sqrt5)/2 the golden ratio. The binding
  demand-fraction sits at the GOLDEN SECTION of the unit interval: x* = (3-sqrt5)/2 = 1/phi^2 (so 1-x* = 1/phi
  is exactly death-star's delta = 0.6180). Hence the conjectured minimal deadline constant is
  C* = 1 + log_5(phi^2) = log_5(5 phi^2) = log_5((15+5 sqrt5)/2) = 1.5979874356654401497... Proof: the
  min-max Lagrangian collapses to two expressions for gamma, both equal to 2 ln phi/ln 5 at x*=(3-sqrt5)/2
  via ln x*=-2ln phi, ln(1-x*)=-ln phi, ln(1+x*)=(1/2)ln5-ln phi. Verified to 50 digits.
tags: [amm12592, minimal-C, biased-coin, extractor, entropy, threshold, golden-ratio, two-ray, THM-3002, HYP-9061, death-star]
related: [HYP-9061, THM-3002, amm12592-C-eq-1-szego-rigidity]
---

# AMM 12592: the two-ray threshold is golden -- gamma* = 2 ln(phi)/ln 5 (opus-S4)

## The object (death-star, THM-3002 / lanes C2/F2)

death-star's construction-side analysis of AMM 12592 (minimal deadline constant `C`, `T(n) <= Cn+D`) reduces
the exactly-fair dyadic-epoch closure to a capacity/entropy comparison. The asymptotic feasibility threshold
`gamma* = C*-1` is the smallest `gamma` for which a **two-ray entropy comparison** holds for every demand
fraction `x in (0,1)`:

```
   max_y [ gamma(1+y) H( (x-y)/(gamma(1+y)) ) + (x-y) ln 2 ]  >=  H(x),      H = natural-log binary entropy.
```

The supply is the better of two "rays": a capacity channel `gamma(1+y)` carrying entropy at fill
`frac=(x-y)/(gamma(1+y))`, plus a free-bit channel `(x-y)ln2`. death-star's solver returns
`gamma* ~ 0.597987` (`C_arch ~ 1.597987`) with **binding `delta ~ 0.6180`**
(`amm12592_archimedean_threshold_asymptotic.out`), and certified archimedean lower bounds climbing to it from
below (`rho(2048) > 3890/2437 = 1.596225`, `amm12592_archimedean_lower_bound_large_m.out`). Both sides point
at one constant; this note names it exactly.

## Theorem

> **`gamma* = 2 ln(phi) / ln 5 = log_5(phi^2)`, where `phi = (1+sqrt5)/2` is the golden ratio.**
> The binding demand fraction is `x* = (3-sqrt5)/2 = 1/phi^2` (the golden section of `[0,1]`), so
> death-star's binding `delta = 1 - x* = 1/phi = 0.6180339887...`, and the `y`-optimum has
> `frac = 1 - 1/sqrt5`. Hence
> ```
>    C* = 1 + gamma* = 1 + log_5(phi^2) = log_5(5 phi^2) = log_5( (15+5 sqrt5)/2 )
>       = 1.59798743566544014974502650204862157804...
> ```

## Proof (min-max Lagrangian collapse)

At the threshold the binding `(x*,y*)` is an interior saddle: `y*` maximizes the bracket, `x*` minimizes the
slack, and the slack is `0`. Write `frac = (x-y)/(gamma(1+y))`, `H'(u)=ln((1-u)/u)`. The three stationarity
equations are

```
   partial_x :  H'(frac) = H'(x) - ln 2                         (A)
   partial_y :  gamma H(frac) - H'(frac)(1+x)/(1+y) - ln 2 = 0  (B)
   value     :  gamma(1+y) H(frac) + (x-y) ln 2 - H(x) = 0      (C)
```

**(A) collapses the fill.** `H'(frac)=H'(x)-ln2` says `(1-frac)/frac = (1-x)/(2x)`, i.e.

```
   frac = 2x/(1+x)         =>   1-frac = (1-x)/(1+x).                         (A')
```

**(B) collapses to a clean ratio.** Using `(1+x)/(1+y) = 1 + frac*gamma` (from the `frac` definition) and the
identity `H(frac) - frac*H'(frac) = -ln(1-frac)`, equation (B) becomes
`gamma[H(frac)-frac H'(frac)] = ln2 + H'(frac)`, i.e. `gamma * (-ln(1-frac)) = ln2 + H'(frac)`. With (A'):
`-ln(1-frac) = ln((1+x)/(1-x))` and `ln2 + H'(frac) = ln((1-x)/x)`, so

```
   gamma = ln((1-x)/x) / ln((1+x)/(1-x)).                                     (B'')
```

**(C) gives a second expression.** With `P(x) = -2x ln x -(1-x)ln(1-x)+(1+x)ln(1+x)` and `H(x)`,
substituting (A') and `(1+y)` reduces (C) to

```
   gamma = (1+x) H(x) / [ (1+x) P(x) - 2x H(x) ].                             (C'')
```

**The golden root.** Set `x* = (3-sqrt5)/2`. Then (all exact)

```
   1-x* = (sqrt5-1)/2 = 1/phi,     1+x* = (5-sqrt5)/2 = sqrt5/phi,     x* = 1/phi^2,
   ln x* = -2 ln phi,    ln(1-x*) = -ln phi,    ln(1+x*) = (1/2)ln 5 - ln phi.
```

Plug into (B''):  `(1-x*)/x* = phi`, `(1+x*)/(1-x*) = sqrt5`, so `gamma = ln(phi)/ln(sqrt5) = 2 ln phi/ln 5`.

Plug into (C''): the three logs give `H(x*) = (1+x*) ln phi` and
`P(x*) = 2x* ln phi + (1+x*)(ln5)/2`, whence

```
   (1+x*) P(x*) - 2x* H(x*) = (1+x*)^2 (ln5)/2,     (1+x*) H(x*) = (1+x*)^2 ln phi,
   gamma = (1+x*)^2 ln phi / [ (1+x*)^2 (ln5)/2 ] = 2 ln phi / ln 5.
```

Both branches agree at `x*`, so `x* = (3-sqrt5)/2` solves the threshold equation `(B'')=(C'')` and the common
value is `gamma* = 2 ln phi/ln 5`. Numerically the raw two-ray optimizer has slack `< 1e-50` at `(x*,gamma*)`,
and `gamma_closed - gamma_numeric < 1e-50` (mpmath, 50 dps). QED. (The only surd identity `sympy` won't
auto-close is `(1+sqrt5)^2(3-sqrt5)/8 = 1`, i.e. `ln x* = -2 ln phi`; it is elementary.)

## Why golden, why base 5

`x* = 1/phi^2` divides `[0,1]` in the golden proportion `x* : (1-x*) = 1/phi^2 : 1/phi = 1 : phi`. The two
structural ratios at the saddle are `(1-x)/x = phi` (a golden split of the demand) and `(1+x)/(1-x) = sqrt5`
(the "stretch" of the capacity ray), and `disc(phi) = 5` is the base. So `gamma*` is a **golden-mean-shift
entropy in base 5**: the same reason `log_2 phi` is the entropy of `11`-free binary strings, here the
capacity/free-bit split is Fibonacci-geometric and its critical exponent is `log_5 phi^2`. This dovetails
with death-star's finding that the classical scheme is dyadic and the obstruction is purely archimedean
(THM-2976/THM-3002): the archimedean wall sits at the golden section.

## Scope (honest)

This is the **exact threshold of death-star's two-ray comparison** -- a clean theorem about THM-3002's
criterion, upgrading `~0.597987` (6 digits) to a closed form with the binding `delta=0.6180` identified as
`1/phi`. It equals the conjectured `C*` under death-star's reading (the `H=1` ample criterion is the true
construction threshold) and is corroborated by their certified lower bounds `rho(m) -> 1.596...`. It is NOT by
itself the unconditional `C*`: that needs the construction to achieve `gamma*` for all `R` (death-star's live
"periodic-orbit at `gamma~3/5`" target) matched by a lower bound reaching `gamma*`. The current rigorous
unconditional general-class anchor is `C*>=1` together with nonattainment of every `n+o(n)` deadline
(THM-3342); endpoint nonattainment does not imply `C*>1` (MISTAKE-368).
death-star's `rho(2048)>1.5962` is a much stronger certified lower bound whose limit this note pins.

## Prediction handed to death-star

If the certified archimedean lower bounds `rho(m)` and the ample threshold both converge to `C*`, they
converge to `log_5(5 phi^2) = 1.5979874356654401497...`, NOT to any nearby rational (`8/5 = 1.6` is ruled
out; the answer is transcendental). The next lower-bound rung should overshoot `1.5962` toward `1.59799`, and
the binding `delta` should stay pinned at `1/phi = 0.61803398...` as `m` grows. Both are cheap falsifiers.
