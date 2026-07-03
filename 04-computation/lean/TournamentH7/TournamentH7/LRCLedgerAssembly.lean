/-
  TournamentH7.LRCLedgerAssembly — THE LEDGER-POSITIVITY ASSEMBLY (klein-2026-07-02-S117,
  HYP-4022).  The fee-aware arithmetic that discharges `cite_hunter_lonely`'s `hledger`
  (kps-S23, LRCRealRegions) from the two analytic inputs:

    * the SINGLES upper bound   `Σ_w rlength(I ∩ teeth w) ≤ c·(L/7) + F`   (kps `teeth_mass`);
    * the PAIR-FLOOR lower bound `pairCredits ≥ (c−1)·(L/49) − E`           (mac-mini JointRateCore
      per-cell obligation / kps pair-event run/gap analysis).

  Then the Hunter ledger `0 < L − Σ singles + pairCredits` holds whenever the fee budget clears
  the path-Bonferroni credit `L·(48 − 6c)/49` (POSITIVE for c ≤ 7, HYP-4021 `ledger_coeff`):

    L − Σ singles + credits ≥ L − c(L/7) − F + (c−1)(L/49) − E = L·(48 − 6c)/49 − (F + E) > 0.

  This is pure arithmetic (`linarith`), framework-agnostic in `L, singlesSum, credits`; it is the
  exact shape `cite_hunter_lonely`'s `hledger` needs, reducing the c ≤ 7 crux to EXACTLY
  {singles-bound (done) + pair-floor (active)}.  Sorry-free.
-/
import Mathlib

namespace LonelyRunner.LRC14.Ledger

/-- The path-Bonferroni credit coefficient `1 − c/7 + (c−1)/49 = (48 − 6c)/49`. -/
theorem ledger_coeff (c : ℝ) : 1 - c / 7 + (c - 1) / 49 = (48 - 6 * c) / 49 := by ring

/-- **The ledger-positivity assembly.**  From the singles upper bound and the pair-floor lower
bound, with the fee budget `F + E` below the path-Bonferroni credit, the Hunter ledger is positive.
`c ≤ 7` enters through the credit being positive (`48 − 6c ≥ 6 > 0`). -/
theorem hledger_pos_of_bounds
    (L singlesSum credits F E : ℝ) (c : ℕ) (hL : 0 < L) (hc : c ≤ 7)
    (hsingles : singlesSum ≤ (c : ℝ) * (L / 7) + F)
    (hcredits : ((c : ℝ) - 1) * (L / 49) - E ≤ credits)
    (hfee : F + E < L * (48 - 6 * (c : ℝ)) / 49) :
    0 < L - singlesSum + credits := by
  have hcR : (c : ℝ) ≤ 7 := by exact_mod_cast hc
  -- L - singlesSum + credits ≥ L - (c·L/7 + F) + ((c-1)·L/49 - E) = L·(48-6c)/49 - (F+E)
  have hkey : L - ((c : ℝ) * (L / 7) + F) + (((c : ℝ) - 1) * (L / 49) - E)
      = L * (48 - 6 * (c : ℝ)) / 49 - (F + E) := by ring
  nlinarith [hsingles, hcredits, hfee, hkey]

/-- Clean form with the credit named, matching `ledger_coeff`: for `c ≤ 7` the credit
`L·(48−6c)/49` is strictly positive, so any fee budget strictly below it clears the ledger. -/
theorem credit_pos {L : ℝ} (hL : 0 < L) {c : ℕ} (hc : c ≤ 7) :
    0 < L * (48 - 6 * (c : ℝ)) / 49 := by
  have hcR : (c : ℝ) ≤ 7 := by exact_mod_cast hc
  have h1 : 0 < L * (48 - 6 * (c : ℝ)) := mul_pos hL (by linarith)
  exact div_pos h1 (by norm_num)

/-- **The general ledger, parametrized by the pair overlap `P`.**  This subsumes
`hledger_pos_of_bounds` (worst-case `P = L/49`, valid only for `c ≤ 7`) and, crucially,
extends to ANY block size `c` — including `c ≥ 8` — provided the consecutive-pair overlap `P`
clears the `c`-dependent threshold.  The Hunter ledger `0 < L − Σ singles + credits` holds
whenever the singles fee plus the pair-error is below `L + (c−1)P`. -/
theorem hledger_pos_of_pairbound
    (L singlesSum credits P F E : ℝ) (c : ℕ)
    (hsingles : singlesSum ≤ (c : ℝ) * (L / 7) + F)
    (hcredits : ((c : ℝ) - 1) * P - E ≤ credits)
    (hthresh : (c : ℝ) * (L / 7) + F < L + (((c : ℝ) - 1) * P - E)) :
    0 < L - singlesSum + credits := by
  linarith

/-- **The explicit pair-overlap threshold for block size `c`.**  Ignoring the fee terms, the
ledger closes iff the consecutive-pair overlap `P` exceeds `L·(c−7)/(7(c−1))`.  KEY CONSEQUENCE:
for `c ≤ 7` the threshold is `≤ 0` (any `P ≥ 0` works — the worst-case floor `L/49` suffices);
for `c ≥ 8` it is strictly positive and grows toward `L/7` (e.g. `c = 13` needs `P > L/14`), so
`c ≥ 8` is NOT closable by the worst-case pair-floor `L/49` (`L/49 < L/14`) — it requires the
NEAR-EQUAL pair overlap, which exceeds `L/49` exactly when the block is tight. -/
theorem pairbound_threshold (L c P : ℝ) (hc1 : 1 < c)
    (hP : L * (c - 7) < 7 * (c - 1) * P) :
    c * (L / 7) < L + (c - 1) * P := by
  nlinarith [hP]

#print axioms hledger_pos_of_bounds
#print axioms hledger_pos_of_pairbound

end LonelyRunner.LRC14.Ledger
