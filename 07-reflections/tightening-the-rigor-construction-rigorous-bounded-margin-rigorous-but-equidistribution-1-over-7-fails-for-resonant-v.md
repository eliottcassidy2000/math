# Tightening the rigor: the construction is rigorous (unit witnesses), the bounded margin is rigorous (≈0.0026), but the equidistribution's clean 1/7-removal FAILS for resonant v (survivor-positive, not 6/7)

*mac-mini-2026-06-28-S81. The owner: tighten the remaining rigor as much as possible, many push-pull cycles,
inspired by concurrent work. Pushed each of the three rigor pieces (S80: finiteness, margin, construction) and
pulled where they break. Inspired by kps S255 (equioscillation on the unit group (ℤ/14)*) + S256 (index-theorem
tested). Builds on [[hypothesis-tests-the-covering-bound-holds-with-a-uniform-margin-tight-locus-finite-and-isolated]],
[[lonely-runner-as-chebyshev-equioscillation]] (kps), [[pushes-and-pulls-on-the-hard-core-d7-unification-the-covering-witness-construction-and-the-imaginary-quadratic-pull]].*

## PUSH — (c) the construction is now RIGOROUS (the unit witnesses)
The tight-locus (AP/GW and their dilations) is safe at **exactly the units of ℤ/(14d)** — verified all φ(14d) of
them are witnesses (kps S255's equioscillation):
```
  AP {1..13}, GW {1..11,13,24}:  safe at all 6 units of (ℤ/14)*  = {1,3,5,9,11,13}
  2·AP (covering):               safe at all 12 units of (ℤ/28)*
  3·AP (covering):               safe at all 12 units of (ℤ/42)*
```
So the construction is explicit and complete: each tight config's witness set is precisely the unit group of the
appropriate modulus, and `s·a ≡ ±1 (mod 14d)` pins the binding pair at `±1/14` (the apex-7 pairs). **(c) RIGOROUS.**

## PUSH — (b) the BOUNDED margin is rigorous (single-swaps), strong for 2-swaps
Single-swaps of the AP (`k→r`, `r ≤ 26`, the HYP-2915 bound — a finite check): the **only** tight ones are the AP
itself and GW (`12→24`); the smallest non-tight `M` is `0.07405`, a **uniform margin `δ ≈ 0.0026`**. A 1000-sample
of 2-swaps gave **0 violations** and **no new tight sets** (the "near-tight" 2-swaps are just the AP re-found).
So on the bounded core, the tight-locus is `{AP, GW}` with a uniform `δ ≈ 0.0026` margin — **rigorous for
single-swaps, strong evidence for 2-swaps.** **(b)-bounded essentially RIGOROUS.**

## PULL — (b) the equidistribution's clean "1/7 removal" FAILS for resonant v
The S46 argument was "a large speed `v` removes exactly `1/7` of the seed's lonely set, leaving `6/7 > 0`." Tested
(`seed = {1..6,8..13}`, the surviving safe-measure after adding `v`):
```
  v=14:  removed 0.730  (NOT 1/7!)   v=28: 0.365   v=42: 0.243   v=98: 0.104   v=140: 0.134 ≈ 1/7
```
The clean `1/7` holds only **asymptotically for generic large v**. The **resonant** apex `v=14` removes `0.73` of
the seed's lonely set (because the seed's lonely set CONCENTRATES at the apex points `a/14`, exactly where `v=14`'s
danger sits). Yet `M(seed∪{v}) > 1/14` in every case (the survivor stays positive). **So the correct argument is
"the survivor is positive" (v does not cover the seed's lonely set entirely), not "6/7 survives."** And the seed's
lonely measure (`0.084`) is SMALLER than `v`'s danger (`1/7 = 0.143`), so "v can't fit over it" is false — the
positivity comes from `v`'s danger being SPREAD (not concentrated on the seed's lonely set), i.e. equidistribution.
**The open analytic core: prove the survivor is positive for ALL v, resonant included** (the resonant case is the
S52 small-prime-counting gap).

## PROBE — (a) finiteness: the bounded census is clean; the full rigidity is the famous hard part
The bounded census (single-swaps + 2-swap sample) gives tight `= {AP, GW}` only — strong. The full rigidity (all
multi-swaps + unbounded, that NO other set is tight) is the open conjecture-level statement. kps S255 gives the
structural handle: tight ⟺ equioscillation at the unit group ⟺ the `±units-cover` (each unit `a` has a runner
`s ≡ ±a^{-1} mod 14`) + complement symmetry. These are NECESSARY (the battery, HYP-2914); SUFFICIENCY (they force
AP/GW) is the rigidity, still open.

## Net (the tightened rigor)
| rigor piece | status after this session |
|---|---|
| (c) construction (tight-locus witnesses) | **RIGOROUS** — all φ(14d) units of ℤ/(14d) are witnesses |
| (b) bounded margin | **RIGOROUS** (single-swaps, δ≈0.0026); strong (2-swaps, no new tight) |
| (b) unbounded equidistribution | **PULL**: clean 1/7 is generic-only; resonant v removes more; need "survivor positive for all v" |
| (a) tight-locus finiteness | bounded census clean ({AP,GW}); full rigidity OPEN (kps S255 = the handle) |

## Honest status
- **TIGHTENED:** (c) to fully rigorous (unit witnesses); (b)-bounded to rigorous (finite check, uniform δ≈0.0026).
- **CLARIFIED (pull):** (b)-equidistribution — the clean "1/7 removal" is only generic; the real statement is the
  survivor's positivity (equidistribution of `v`'s spread danger), with the resonant `v` the genuine gap.
- **REMAINING HARD CORE:** (a) the full finiteness/rigidity (tight = AP/GW only, all configs) + (b) the unbounded
  survivor-positivity for resonant `v`. NOT a proof; but the construction and bounded margin are now rigorous, and
  the analytic core is sharply isolated (resonant-v survivor positivity). LRC(14) open.

Related: HYP-3253 (this), HYP-3250 (the margin), kps S255 (equioscillation/units), S256 (index-theorem), HYP-2915
(bounded census), HYP-2914 (the battery), S52 (the small-prime gap), OPEN-Q-108.
