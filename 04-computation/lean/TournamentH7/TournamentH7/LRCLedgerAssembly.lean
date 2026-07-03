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

#print axioms hledger_pos_of_bounds

end LonelyRunner.LRC14.Ledger
