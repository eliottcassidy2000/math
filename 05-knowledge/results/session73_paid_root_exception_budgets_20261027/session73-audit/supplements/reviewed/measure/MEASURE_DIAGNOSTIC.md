# Diagnostic only: the beta clock boundary changes the source measure

**PROVED elementary continuation of inherited beta weights; not a repair by
reweighting the native task.** The original native probability is exactly
Beta(1,2). No claim below licenses replacing it. This optional note is separate
from the frozen core and the scalar/spatial supplement.

Inherited sources: `collatz_adaptive_mixture_flow_20261005.md` P1, P4 and P5;
`collatz_critical_beta_weights_20261005.md` B2 explicitly has total defect 2
for Beta(2,2), so normalization cannot be guessed; and
`collatz_mixed_refuel_weights_20261005.md` preserves counter/clock distinctions.
The two latter documents were inspected at the same commit as the core source
pins. This is a targeted inheritance check, not a literature novelty claim.

For alpha>0 and beta>1, put

```
W_(alpha,beta)(n)=B(K+alpha,L+beta)/B(alpha,beta)
```

on rooted sources, with the inherited root/unrooted conventions. Let
`lambda_(alpha,beta)=W_(alpha,beta) 1_(3|n)`. Its exact total mass is

```
c=alpha/(beta-1),  lambda_bar=lambda_(alpha,beta)/c.        (D1)
```

Indeed the nonroot inverse ROOT ray carries mass
`integral r/(1-r) dBeta(alpha,beta)=c`; summing the killed incoming balance
identifies this with the leaf deficit. The inherited P4 total-mass bound is
finite when beta>1, so its occupation proof applies unchanged. The native choice (1,2) has c=1. More generally c=1 along beta=alpha+1,
but other pairs on that line define different weights and do not license
changing the original measure. In particular Beta(2,2) has c=2, as the
earlier source already records.

**Sharp first literal-clock mean boundary:** under this different normalized
source law, the expected first-hit literal standard clock is finite **iff
beta>2**. The expected shortcut clock has the same boundary, whereas the
expected odd-step count is finite throughout beta>1.

For sufficiency, write `n=S^j(b)` at a minimal sibling base, with
`a(n)=a(b)+2j`, `a(b) in {1,2}`. At fixed price r, inherited P1 says
`sum_b g_r(b)<=2+2r-r^2`. The literal local cost is `1+a(n)`, so

```
sum_(n>1) (1+a(n)) f_r(n)
 <= (2+2r-r^2) [3/(1-r)+2r/(1-r)^2]-3
  = 6/(1-r)^2+3/(1-r)-(1-r)-5.                          (D2)
```

The subtraction removes the unexecuted root self-loop. Integrating and using
occupation gives the explicit finite upper bound, for beta>2,

```
E_(lambda_bar) clock_literal <= (beta-1)/alpha *
 [6(alpha+beta-1)(alpha+beta-2)/((beta-1)(beta-2))
  +3(alpha+beta-1)/(beta-1)-beta/(alpha+beta)-5].           (D3)
```

For example it is 67/2 at (alpha,beta)=(1,3), a loose sufficient constant,
not an optimized expected time.

For necessity, the lambda-positive ROOT rays are `S^j(1)` with
`j=2 mod3`; their literal clocks are `2j+3`. Their weights are
`B(j+alpha,beta)/B(alpha,beta)`. For all sufficiently large j, integrate
only over `1-1/j<=r<=1-1/(2j)`: the integrand's power of r stays uniformly
positive there, so the weight is at least a positive constant times
`j^(-beta)`. The resulting expected-clock subseries dominates
`sum_(j=2 mod3) j^(1-beta)`, which diverges for beta<=2.

The mean odd-step lifetime is finite by the inherited occupation mass for
beta>1. Literal clock equals shortcut clock plus odd-step count, proving the
stated common boundary for the two nonaccelerated clocks.

The phase-aware forcing mean is still at most
`(sum_(n>1)W_(alpha,beta))/(19c)` for every beta>1. Hence the strip
`1<beta<=2` has finite odd-step/phase-energy accounts but infinite physical
clock. The original native measure is exactly on its upper boundary.
This diagnoses the two accounts; it does not pay the original clock by
quietly substituting a new source distribution.
