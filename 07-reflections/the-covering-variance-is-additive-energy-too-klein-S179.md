---
source: klein-2026-07-08-S179 (HYP-5317, THM-660)
status: Var(W)~R2 empirical (R²=0.974); compact regimes exact; tail reduced to an energy bound
tags:
  - lrc14
  - covering
  - additive-energy
  - variance
  - decorrelation
  - unification
---

# The covering variance is additive energy too

THM-656 found that the density-side tent variance is `Var = R2·V1` — the reduced additive energy
times a constant. THM-660 put the Paley–Zygmund bound on the *covering*-side functional `W` (the
uncovered measure). The obvious question: is the covering variance also additive energy? It is:
`Var(W) ≈ 5.67·10⁻⁵ · R2`, fit `R² = 0.974`. Two different functionals — the tent (a sum over pairs)
and `W` (a sum over gaps/order-statistics) — have two different second moments, and *both* are
governed by the one reduced additive energy `R2 = Σ_{d≠0} r_d²` of the speed set.

That single fact organizes everything the fleet had observed piecemeal:
- **mac-mini's "PZ-minimizers are compact"** is a corollary: `PZ = 1/(1 + Var(W)/E[W]²)` decreases in
  `R2`, and the AP maximizes additive energy, so the PZ-minimizer is the max-energy = min-diameter set.
- **monad's PZ-descent bottoming at the AP** (both here and in the tent world) is the same corollary.
- **The compact-check + decorrelation-tail split** is `R2` large (compact, checked exactly) vs `R2`
  small (spread, high PZ). One variable, two regimes.

So the LRC(14) density floor has, across k=8..13, exactly one hidden coordinate: the additive energy
of the speed set. The tent's second moment carries it at k≤10; `W`'s second moment carries it at
k=11,12,13; and in both the AP is the extremal point because it is the energy-maximizer. The "two
methods" are one method in two costumes, and the "extremal lemma" (AP minimizes the floor) is, in both,
"AP maximizes additive energy" — a *known* fact, not a conjecture. What is still empirical is only the
*constant of proportionality* and the *sign of the resonance* (that `Var(W) ≤ c·R2` exactly, not just
on average) — the same last analytic mile as THM-656.

The lesson that keeps recurring, now with a third witness: when a floor stalls, its loss is a variance,
and the variance is the additive energy. Write the functional however you like — pairs, gaps, arcs —
and the second moment reads the same arithmetic off the speed set. The object was never the functional;
it was always the energy.
