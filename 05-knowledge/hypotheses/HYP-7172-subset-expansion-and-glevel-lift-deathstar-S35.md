# HYP-7172 — The B5 subset expansion under the trap + the GLevel lift (death-star-2026-07-16-S35)

**Status:** CLAIMED — in progress (S35). Verify the build report in the session log.

**Owner directive:** "discharge codex's singular-series identification under the trap
and run the GLevel lift."

**Claim A (module `LRCB5SubsetExpansion.lean`):** the discrete half of the
singular-series identification, exact and kernel-pure:
(1) `B5_eq_subset_sum` — B5 v q = Σ_{|T| ≤ 5} (−1)^{|T|} N_T with N_T = the joint
band-failure count of the subset T (moment→subset via choose = #powersetCard);
(2) the deviation ledger — N_T = (q−1)/7^{|T|} + D_T defines exact ℚ-deviations;
B5/(q−1) = 2052/16807 + Σ (−1)^{|T|} D_T/(q−1) EXACTLY (the equilibrium constant
DERIVED, matching codex's `equilibrium` — the truncated binomial identity
Σ_{d≤5}(−1)^d C(13,d)/7^d = 2052/16807);
(3) `B5_pos_of_deviation_debt` — positivity from |deviation sum| < (q−1)·equilibrium,
consuming codex's `relationModel_pos_of_debt_lt` SHAPE with the masses now DEFINED
from discrete counts, not hypothesized. The trap's role (THM-939): for q beyond the
family's magnitude scale, mod-q coincidence patterns of support ≤ 5 with top above
the dense pair force integer relations — which A1/A2 forbid. The remaining analytic
obligation sharpens to: bound the D_T of trapped subsets (equidistribution), now a
concrete ℚ-number per subset.

**Claim B (module `LRCGLevelLiftDichotomy.lean`):** insert the dense TRIPLE
{w(j), w(j+1), w(j+2)} as one ≤6-block GLevel (kps-S22 `cite_glevel_lonely`, mass fee
3L/7 + (3/7)Σ1/wᵢ + out·(LΣwᵢ + 9) < L − out) when THM-938's pair dodge fails;
the dense core shrinks from TripleDenseCore toward the quad/7-wall normal form.
