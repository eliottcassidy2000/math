---
id: HYP-3718
title: THE OBSERVER ESCAPES AT THE CONVERGENT, NOT THE MEDIANT -- n/Phi_6(n)=[0;n-1,n] has semiconvergents j/(j(n-1)+1), j=1..n: j=1 = 1/n (the MEDIANT = the LRC threshold), j=n = n/Phi_6(n) (the CONVERGENT = the covering-min); the killer lcm(n-1,n)=n(n-1) BLOCKS the mediant escape (n(n-1)/n=n-1 in Z, so the killer sits at 0 at t=1/n, M=0, not lonely), while at each semiconvergent t=j/(j(n-1)+1) the lonely gap is EXACTLY that value (climbing monotonically to the convergent); so a COVERING set (which must kill resonance n via the killer) forces the observer's lonely-time escape to the CONVERGENT n/Phi_6(n) > 1/n, NOT the mediant 1/n. This is the cusp/off-cusp split: the mediant 1/n = the tight extremal {1..n-1} (measure-0 cusp, M=1/n, proven small n); the convergent n/Phi_6(n) = the covering regime (positive measure, M>1/n, the open n>=7 / apex genus>=1 case)
status: VERIFIED (the semiconvergent escape ladder for {1..12,182}: mediant 1/14 BLOCKED M=0, convergent 14/183 the max escape M=14/183; M at semiconvergent j = j/(13j+1) exactly; killer blocks mediant for n=4,7,10,14; structured global scan finds nothing below 14/183). The killer-blocks-mediant + convergent-escape mechanism is rigorous; the GLOBAL covering-min remains the open piece.
source: klein-2026-06-29-S29
depends_on:
  - HYP-3717   # the three-gap gives the covering-min; CF [0;n-1,n]
  - HYP-3715   # the covering-min = n/Phi_6(n) at the densest-core + killer
related:
  - HYP-2945   # prior Farey-mediant non-commutation work
  - HYP-3551   # M=j/D; 14/183; densest core + minimal killer
  - HYP-3597   # cusp (measure 0) vs off-cusp (positive measure)
results:
  - 04-computation/observer_escape_convergent_not_mediant_klein.py
  - 05-knowledge/results/observer_escape_convergent_not_mediant_klein.out
---

# HYP-3718 — the observer escapes at the convergent, not the mediant

## The semiconvergent escape ladder
`n/Phi_6(n) = [0; n-1, n]`; its **semiconvergents** are `j/(j(n-1)+1)`, `j = 1..n`:
- `j = 1`: `1/n` -- the **MEDIANT** (the lonely-runner conjecture threshold);
- `j = n`: `n/Phi_6(n)` -- the **CONVERGENT** (the covering-min).
For the covering set `{1..12,182}` (`n=14`), the lonely gap `M(t) = min_s ||s t||` at each semiconvergent
`t = j/(13j+1)`:
```
j :  1     2     3    ...   13     14
t : 1/14  2/27  3/40  ... 13/170  14/183
M :  0   2/27  3/40  ... 13/170  14/183   (M = the semiconvergent value itself; climbs to the convergent)
```
- `j=1` (mediant `1/14`): `M = 0` -- **BLOCKED**: the killer `182` sits at `0` (`182*(1/14) = 13 in Z`), so
  the observer is NOT lonely at the mediant.
- `j=14` (convergent `14/183`): `M = 14/183` -- the **maximum**, the observer's escape.

## The killer always blocks the mediant (the realizability node)
The minimal killer is `lcm(n-1,n) = n(n-1)`, and `n(n-1)/n = n-1 in Z`, so **at the mediant time `t = 1/n`
the killer is exactly at the observer (`0`)** -- the mediant escape is BLOCKED for every `n` (verified
`n=4,7,10,14`). A covering set is REQUIRED to kill resonance `n` (THM-523), and the minimal killer that does
so is precisely the point that lands at `0` at `t=1/n`. So **the covering condition itself blocks the mediant
escape**, forcing the observer to the next realizable node -- the **convergent** `n/Phi_6(n) > 1/n`. The
realizability node is the convergent (`j=n`), not the mediant (`j=1`).

## In observer terms (the picture)
The static observer at `0` looks for a lonely time. The mediant `1/n` would give exactly `M = 1/n` (the
conjecture threshold) -- but a covering set parks its killer at `0` there, so the observer is crowded out.
Climbing the semiconvergent ladder `2/27, 3/40, ..., 14/183`, the lonely gap grows monotonically; the
observer's actual escape is the top rung, the **convergent** `14/183`, where `M = 14/183 > 1/14`. The
covering floor's strict positivity above `1/n` IS the gap between the blocked mediant and the realized
convergent.

## The cusp/off-cusp split and the n>=7 regime
- the **mediant `1/n`** = the tight extremal `{1..n-1}` (the measure-0 cusp, `M = 1/n` exactly) -- the
  classical tight case (proven small `n`);
- the **convergent `n/Phi_6(n)`** = the COVERING regime (positive measure, off-cusp, `M = n/Phi_6(n) > 1/n`)
  -- the open case, the `n>=7` / apex-genus`>=1` regime (the apex prime `7`, where `genus(X_0(2p))` jumps to
  `1` at `N=14`).
So for `n >= 7` the relevant (covering) escape is the convergent, strictly above the mediant -- the LRC's
"uniform looseness" (HYP-2566) is exactly the convergent-over-mediant gap `n/Phi_6(n) - 1/n = (n-1)/(n
Phi_6(n))`.

## Global optimality (the other next target)
Structured scan over densest-core families (n=14): `min M = 14/183` at `{1..12,182}`; nothing beats it --
supports `14/183` as the global covering-min (HYP-3551). Still OPEN as a theorem (exotic, non-densest-core
coverings), but the convergent-escape mechanism explains WHY the floor sits at the convergent, not the
mediant.

## Net
The covering condition blocks the mediant escape (the killer parks at `0` at `t=1/n`), so the observer's
lonely-time escape is the convergent `n/Phi_6(n)`, not the mediant `1/n` -- giving `M > 1/n` directly. The
convergent (covering, off-cusp, `n>=7`/apex-genus`>=1`) vs the mediant (tight, cusp, small `n`) is the
realizability node of the whole bridge.
