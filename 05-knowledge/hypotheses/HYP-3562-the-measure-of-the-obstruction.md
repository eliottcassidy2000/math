---
id: HYP-3562
title: The MEASURE OF THE OBSTRUCTION -- every "X exists" theorem in the project (Redei H>=1, SC tournaments exist, the LRC lonely point exists, the metagraph is connected) is the NONVANISHING of ONE equivariant obstruction class, the Euler/Lefschetz class of the complement involution R, and the project has been computing its MEASURE under different measures: COUNTING (SC = trace(R) = P_n(-1) = Lefschetz number, the metagraph), LEBESGUE (meas(lonely) = the floor, the LRC), and EULER (chi_meas of the danger-cover nerve, HYP-3242); the sigma-EVEN/sigma-ODD split of every object IS the Lebesgue(measure)/counting(index) split, so the Borsuk-Ulam/sigma-odd index is exactly what DETECTS the obstruction at the extremal where the Lebesgue measure VANISHES; a disproof = the class is EXACT (a coboundary, all measures zero), a proof = it is ESSENTIAL (some measure positive), and essentiality is a TOPOLOGICAL = SET-INDEPENDENT invariant (the Gamma_0(N) story)
status: SYNTHESIS + verified grounding (SC=trace(R)=Lefschetz>0 => SC tournaments exist, n=3..7; the LRC extremal {1,2,3} has Lebesgue floor 0 but lonely-count phi(4)=2; covering {2,3,4} floor 1/8). A unifying framework, not a new proof.
source: mac-mini-2026-06-29-S23
related:
  - THM-587   # SC = P_n(-1) = trace(R) (the Lefschetz number = the metagraph obstruction's counting measure)
  - HYP-3544  # klein equivariant homology (SC = chi = sum (-1)^k b_k; R-odd Betti = the obstruction)
  - HYP-3242  # chi_meas(nerve) = the cap (the Euler-calculus measure of the obstruction)
  - HYP-3538  # the cap M_odd = the R-odd obstruction (the LRC counting/index measure)
  - THM-581   # the Borsuk-Ulam witness (sigma-odd) vs the floor (sigma-even)
  - THM-589   # the moment method computes the measure (1st = mean count, 2nd = W(n) concentration)
  - HYP-3553  # the Gamma_0(N) set-independence = essentiality of the class (topological)
reflections:
  - lrc14-first-obstruction-cocycle-generation-codex-s259.md  # the sheaf-Cech syndrome side
---

# HYP-3562 -- the measure of the obstruction

## The framework
Every existence theorem in the project is the **nonvanishing of one equivariant obstruction class** -- the
Euler/Lefschetz class of the complement involution `R` -- and the project has been computing its MEASURE
all along, under three measures:
| object | the obstruction | its MEASURE | existence forced because... |
|---|---|---|---|
| **metagraph** | `SC = P_n(-1) = trace(R)` = the **Lefschetz number** | COUNTING (`SC` = #fixed classes) | `SC != 0` => `R` has fixed points => SC tournaments exist |
| **LRC** | the lonely set (the hole in the danger cover), `R`-symmetric | LEBESGUE (`meas(lonely)` = the floor) | floor `> 0` => a lonely time exists |
| **nerve** | the cover nerve | EULER (`chi_meas`, HYP-3242) | `chi_meas != 0` => a hole exists |
Verified: `SC = 2,2,8,12,88` (n=3..7), always `> 0` -- the metagraph obstruction is always essential, so
self-complementary tournaments always exist (existence WITHOUT construction, the Lefschetz fixed-point
theorem). This is the finite, exact rehearsal of the LRC lonely-point existence.

## The two measures, and why the sigma-odd index matters
The obstruction carries **two measures**, and the project's `sigma`-even/`sigma`-odd split IS the split
between them:
- `sigma`-EVEN = the **LEBESGUE/measure** side (the bulk floor, the SOS/Brouwer part, `R`-even,
  `(A000568+SC)/2`).
- `sigma`-ODD = the **COUNTING/Euler/index** side (the Borsuk-Ulam index, the units, `R`-odd, `#NS`).
VERIFIED on the LRC: the extremal `{1,2,3}` (LRC4) has **Lebesgue floor = 0** but lonely **count =
`phi(4)=2`** (the units mod 4) -- the obstruction is nonzero via the COUNTING measure though the Lebesgue
measure VANISHES; the covering `{2,3,4}` has Lebesgue floor `1/8 > 0`. So **the Borsuk-Ulam / `sigma`-odd
index is exactly what detects the obstruction at the extremal where the floor measure vanishes** -- the
reason the witness (`sigma`-odd) is needed and not just the floor (`sigma`-even), made precise. `SC =
(R-even) - (R-odd)` is the alternating sum (the Euler/Lefschetz number) of the two measures.

## Disproof = exact, proof = essential; essentiality is set-independent
A **disproof** = the obstruction class is **EXACT** (a coboundary): all measures vanish, the hole is
filled, the lonely set is empty. A **proof** = the class is **ESSENTIAL** (not a coboundary): some measure
is positive. The decisive point: **essentiality is a TOPOLOGICAL invariant** -- it depends on the
homotopy type of the quotient / the congruence structure, NOT on the specific speed set. That is exactly
the `Gamma_0(N)` set-independence (HYP-3553): the obstruction lives in the equivariant cohomology of the
covering quotient, indexed by `N`, so its nonvanishing (floor `> 0`) is a set-free statement. klein-S4's
finding (the per-set `CV(N_R)^2` is unbounded) is the symptom of measuring the obstruction set-by-set
instead of reading its essentiality off the topology.

## The moment method computes the measure
Obstruction theory says WHEN existence is forced (the class is essential); the MEASURE says HOW MUCH (the
floor). The MOMENT method (THM-589) computes the measure: the **first moment** `E[#lonely] = E[#fixed] =`
the Siegel/Burnside mass (existence if `> 0`); the **second moment** `= W(n)` (the metagraph variance,
`CV^2 ~ 2/n`) gives concentration. So obstruction theory (topological, "essential"), the moment method
(probabilistic, "the measure"), and the equivariant homology (klein HYP-3544, "the `R`-odd Betti") are
three computations of the one measure-valued obstruction class.

## What it buys
A single frame: name the obstruction (the equivariant Euler class of `R`), and (i) the floor's positivity
becomes its essentiality -- topological, SET-INDEPENDENT (the `Gamma_0(N)` cure); (ii) the `sigma`-odd
index is precisely the counting measure that survives when Lebesgue vanishes (the extremal), explaining
why the witness is needed; (iii) the moment method is the computation of the measure, with the metagraph
(`SC` = Lefschetz, `W(n)` = variance) the exact finite rehearsal. The disproof is the class becoming a
coboundary -- which the set-independent essentiality forbids.
