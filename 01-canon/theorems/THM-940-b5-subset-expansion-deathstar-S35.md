# THM-940 — The B5 subset expansion and the discrete deviation identity (death-star-2026-07-16-S35)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCB5SubsetExpansion.lean`,
axioms propext/Classical.choice/Quot.sound; verify the build report in the session
log). The DISCRETE half of codex-S28's named singular-series obligation
(LRCB5RelationBudget's honest note). Source: HYP-7172.

## Statement

For any family `v` and ruler `q`, with `N_T := jointFail v q T` (multipliers where
every runner of `T` misses the safe band):

1. `momentS_eq_sum_jointFail`: S_d = Σ_{|T|=d} N_T (choose ↔ powersetCard + swap);
2. `B5_eq_subset_sum`: B5 = Σ_{k≤5} (−1)^k Σ_{|T|=k} N_T — the quintic functional
   grouped by SUPPORT, the exact discrete analogue of the singular series;
3. `equilibrium_binomial`: Σ_{d≤5} (−1)^d C(13,d)/7^d = **2052/16807** — codex's
   equilibrium constant DERIVED (truncated binomial at failure rate 1/7);
4. `B5_eq_equilibrium_add_deviation` (**the identification**): with the deviation
   ledger `D_T := N_T − (q−1)/7^|T|` (exact ℚ),
       (B5 v q : ℚ) = (q−1)·2052/16807 + Σ_{k≤5} (−1)^k Σ_{|T|=k} D_T — EXACTLY;
5. `B5_pos_of_deviation_debt`: |signed deviation sum| < (q−1)·2052/16807 ⟹ 0 < B5 —
   codex's `relationModel_pos_of_debt_lt` shape with the masses now DEFINED discrete
   quantities, not hypothesized reals. Also `deviation_empty`: D_∅ = 0.

## Where the trap enters (THM-939)

For q beyond the family's magnitude scale, a mod-q coincidence pattern of support
≤ 5 topping out above the dense pair forces an integer relation of the same
support — which A1/A2 forbid. The remaining analytic obligation is now sharp and
concrete: bound the ℚ-numbers D_T for trapped subsets (per-subset equidistribution),
rather than an unspecified "analytic bridge."

## Referee

`subset_expansion_quadcore_referee_deathstar_S35.py`: direct B5 = subset form =
equilibrium+deviation on 40 random (v, q), exact rationals — PASS.
