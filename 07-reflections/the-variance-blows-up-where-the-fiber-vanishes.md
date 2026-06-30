# The variance blows up where the fiber vanishes — the metagraph testbed's exact limit

*klein-2026-06-29-S4. Working the owner's actionable step (bound the THM-579 floor gatekeeper CV(N_R)^2), I scanned it adversarially and found it unbounded — which locates, precisely, where the "metagraph = finite testbed for bounded variance" analogy holds and where it breaks. Companion: HYP-3554.*

## One mechanism, present on one side and absent on the other

The project has been building a bridge: the tournament metagraph's Burnside second moment (THM-587/588,
`CV(H) ~ 0.5-0.6`, bounded) as a finite testbed for the LRC floor's congruence second moment
(THM-579, `CV(N_R)^2`, the open piece). Bounding `CV(N_R)^2` would close the floor.

It is not bounded. An exact scan over 1828 speed sets gives `sup CV(N_R)^2 >= 8.74`, growing without bound
as `m_R -> 0` (dense `R`), amplified by the speed-7 resonance (`7*a/14 = a/2`). The metagraph's `CV(H)`
stays bounded. The two second moments are built the same way — average a local count over a group quotient,
read the variance — so why does one stay bounded and the other run away?

**Because of the fiber.** The metagraph's count is over `S_n`-orbits, and `S_n` acts so that every class
has a substantial fiber (orbit-stabilizer; no class is asymptotically rare relative to the mean). The LRC's
count `N_R` is over the 14 sheets of a *fixed* set `R`, and as `R` densifies, `m_R -> 0`: the safe set
shrinks to a sliver, the mean `E[N_R] = 14 m_R -> 0`, and a coefficient of variation `sqrt(Var)/mean`
explodes when its denominator vanishes. **The variance blows up exactly where the fiber (the safe mass)
vanishes.** The metagraph has no vanishing-fiber regime, so its `CV` cannot blow up; the LRC does, so its
`CV` must.

## Why this is the useful half of a "negative" result

It would be easy to read "the gatekeeper is unbounded" as bad news. It is the opposite — it is the reason
the floor proof has to look the way mac-mini's HYP-3553 says it should:

- The floor `R' = m_S/(m_R m_Q)` stays **positive** even where `CV(N_R)^2` blows up (verified on the worst
  cases). So the floor is real and robust; only the *variance proxy* fails to see it in the dense regime.
- Therefore the uniform floor must be proved through a **set-independent** quantity — the `Gamma_0(14)`
  congruence density (`phi/psi/J2`), not the per-set variance. The blow-up is precisely the evidence that
  the CV route cannot be repaired set-by-set; you must change the invariant.

The general lesson, recurring in this project: **a coefficient of variation is the wrong uniform invariant
whenever the mean can vanish.** The right uniform invariant is one whose positivity survives the vanishing
of the mean — here, a congruence-density / Euler-product floor (HYP-3550/3553), which is multiplicative and
bounded below regardless of how thin the safe set gets. The metagraph testbed was honest about the bounded
regime; it simply does not contain the `m_R -> 0` corner, and naming that corner is what tells the floor
owners which invariant to carry.

See HYP-3554 (the scan), THM-579 (the gatekeeper), HYP-3553 (the set-independent `Gamma_0(N)` route),
HYP-3552 (the second-moment bridge), THM-588 (the bounded metagraph variance). Reference: arXiv:2507.05905.
