---
id: HYP-3734
title: The covering-min Farey-neighbor question -- PROVED a_1=n-1 (M in (1/n,1/(n-1))) and the REDUCTION Farey-neighbor-of-1/(n-1) <=> binding D ≡ 1 mod (n-1) <=> on the semiconvergent ladder (gcd forces remainder 1); VERIFIED D ≡ 1 mod (n-1) for all covering-mins n=7..14 (answer YES for all computed; general proof OPEN). THE SMALL-DEPTH SPREAD FAMILY: the achievable rungs form an UP-SET [k_min, inf), k_min(n)=2,2,4,4,3 (n=7..11), no low rung for n>=12 (construction rung n takes over). The OBSTRUCTION is arithmetic: rung k demands covering radius floor(kD/(k(n-1)+1)) at EVERY modulus D -- the radius-0 moduli are EXACTLY D<=n-1 (= the THM-523 resonances!), the radius-1 moduli D in (n-1, ~2(n-1)] are the EXTRA demand that over-constrains the n-1 speeds, so low rungs fail as n grows
status: PARTIALLY PROVED. a_1=n-1 PROVED (THM-523 floor < M <= construction < ceiling). The reduction Farey-neighbor <=> D≡1 mod (n-1) PROVED (CF/gcd). D≡1 mod(n-1) VERIFIED n=7..14 (general proof OPEN). Spread family up-set + k_min=2,2,4,4,3 VERIFIED (IP). Obstruction mechanism (radius floor(kD/(k(n-1)+1)) per modulus; radius-0 = resonances <=n-1) VERIFIED.
source: klein-2026-06-30-S38
depends_on:
  - HYP-3732   # the Farey-rung invariant
  - THM-523    # primitive covering-min > 1/n; kill resonances <=n
related:
  - HYP-3718   # the semiconvergent ladder (D=k(n-1)+1)
  - THM-580    # the descent (the multi-modulus / radius-0 structure)
results:
  - 04-computation/farey_rung_spread_family_klein.py
  - 05-knowledge/results/farey_rung_spread_family_klein.out
---

# HYP-3734 — the Farey-neighbor question + the small-depth spread family

## The clean open question, resolved to a reduction
**PROVED: `a_1 = n-1`.** The covering-min `M(n)` satisfies `1/n < M <= n/Phi_6(n) < 1/(n-1)`:
- lower `M > 1/n` is THM-523 (primitive covering-min strictly exceeds the floor);
- upper: the construction `n/Phi_6(n)` is a valid primitive covering, and `n/Phi_6(n) < 1/(n-1)` since
  `n(n-1) = n^2-n < n^2-n+1 = Phi_6(n)`.
So `1/M in (n-1, n)`, hence the first partial quotient `floor(1/M) = n-1` ALWAYS. `M = [0; n-1, ...]`.

**PROVED REDUCTION: `M` is a Farey neighbor of `1/(n-1)` <=> binding `D ≡ 1 mod (n-1)` <=> `M` is on the
semiconvergent ladder.** Write `M = j/D` in lowest terms. `a_1 = n-1` gives `D = (n-1)j + r_1`,
`0 < r_1 < j`. The CF has length 2 (`M=[0;n-1,k]`, a Farey neighbor of `1/(n-1)`) iff `r_1 | j`; but
`gcd(j, r_1) = gcd(j, D) = 1`, so `r_1 | j` forces `r_1 = 1`, i.e. `D = (n-1)j + 1`, i.e. `D ≡ 1 mod (n-1)`,
i.e. `M = j/(j(n-1)+1)` (rung `k=j`).

**VERIFIED: `D ≡ 1 mod (n-1)` for ALL covering-mins `n=7..14`** (binding `D = 13,15,33,37,31,133,157,183`,
each `≡ 1`). So the covering-min IS a Farey neighbor of `1/(n-1)` for every computed `n`. **OPEN:** prove
`D ≡ 1 mod (n-1)` in general (equivalently, the covering-min is never an off-ladder Farey fraction).

## The small-depth spread family
The achievable rungs (IP, small speeds `Smax=4n`) form an **UP-SET `[k_min(n), inf)`** -- if rung `k` is
achievable so is every `k' > k` (a larger rung = larger `M` = more room). So the covering-min `= k_min(n)`:
```
n          : 7  8  9  10 11 | 12
rungs 2..6 : YYYYY YYYYY ..YYY ..YYY .YYYY | .....   (Y=achievable)
k_min(n)   : 2     2     4     4    3      | >6 (construction rung n)
```
The covering-min is the bottom of the up-set; the construction `n/Phi_6` is rung `n` (the top of interest).

## The obstruction (why low rungs fail as n grows) -- arithmetic, via the radius layers
Rung `k` (`t = k/(k(n-1)+1)`) demands, at EVERY modulus `D`, a covering of `Z/D` of radius
`r(D) = floor(kD/(k(n-1)+1))`:
- **radius 0** iff `D <= n-1` -- needs a speed `≡ 0 mod D`, i.e. a MULTIPLE of `D`. These moduli are
  **EXACTLY `D = 2,...,n-1`** = the THM-523 resonances. *So the covering condition "kill all resonances
  `<= n`" IS precisely the radius-0 layer of the Farey-rung demands.*
- **radius 1** at `D in (n-1, ~2(n-1)]` -- needs a radius-1 covering of `Z/D` (a speed within 1 of 0 for
  every rotation). For `D in (n, 2n-2]` these are NOT forced by resonance-killing -- the EXTRA demand.
- radius `r` at `D ~ (r/k)(k(n-1)+1)`, up to radius `k` at the binding `D = k(n-1)+1`.
Lower `k` keeps the radius small across MORE of the modulus range (radius grows `~ kD/(kn)= D/n`, slope `1/n`),
so the radius-1 band must be met together with a tight binding -- over-constraining the `n-1` speeds. As `n`
grows the radius-1 band widens, and low rungs become infeasible: `k_min` rises (`2 -> 4 -> 3 -> ...`), and by
`n=12` no small rung survives and the construction (rung `n`) wins. The exact `k_min(n)` is genuinely
arithmetic (which radius-1 coverings of `Z/D`, `D in (n,2n-2]`, are simultaneously realizable by `n-1`
primitive speeds).

## Net
The clean open question is reduced and largely answered: `a_1 = n-1` is PROVED, the Farey-neighbor property is
EQUIVALENT to `D ≡ 1 mod (n-1)` (PROVED equivalence), and that holds for all covering-mins `n=7..14` (general
proof open). The small-depth spread family is the up-set of achievable rungs `[k_min, inf)` with
`k_min = 2,2,4,4,3` (`n=7..11`); low rungs fail as `n` grows because rung `k` imposes radius-`floor(kD/(k(n-1)
+1))` coverings at every modulus -- the radius-0 layer is exactly the THM-523 resonances, and the radius-1 band
`D in (n,2n-2]` is the extra demand that over-constrains the speeds. The arithmetic of which radius-1 coverings
are realizable picks `k_min`.
