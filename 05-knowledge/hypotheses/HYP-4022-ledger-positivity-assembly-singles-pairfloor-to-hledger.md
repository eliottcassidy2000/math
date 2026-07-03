---
id: HYP-4022
title: THE LEDGER-POSITIVITY ASSEMBLY -- the fee-aware arithmetic that discharges cite_hunter_lonely's hledger (kps-S23) from {singles upper bound + pair-floor lower bound}, reducing the c<=7 crux to EXACTLY those two analytic inputs. hledger_pos_of_bounds (LRCLedgerAssembly.lean, sorry-free, foundational-axioms-only): given singlesSum <= c*(L/7)+F (kps teeth_mass) and credits >= (c-1)*(L/49)-E (the PAIR-FLOOR, mac-mini JointRateCore per-cell obligation) and c<=7 and F+E < L*(48-6c)/49, then 0 < L - singlesSum + credits (the Hunter ledger = cite_hunter_lonely's hledger). Proof: L - singlesSum + credits >= L - (c*L/7 + F) + ((c-1)*L/49 - E) = L*(48-6c)/49 - (F+E) > 0 (nlinarith via the ledger_coeff identity 1-c/7+(c-1)/49=(48-6c)/49). credit_pos: L*(48-6c)/49 > 0 for c<=7. THE STATE: LRC(14) is reduced (klein S115 far-cut dispatch + kps-S23 cite_hunter_lonely + this assembly) to citation + window census + TWO analytic things: (1) the PAIR-FLOOR (pairCredits >= (c-1)(L/49)-err for near-equal teeth -- the COMMENSURATE case = 1/49 EXACTLY is already proven in LRCCommensuration.lean; the general/drifting case is mac-mini/kps active work); (2) the c>=8 near-equal blocks (credit 48-6c <= 0, needs a triple-Bonferroni credit or the large-scale limit fees->0). This file finishes the ARITHMETIC assembly of the c<=7 route
status: VERIFIED (Lean, LRCLedgerAssembly.lean, sorry-free, foundational-axioms-only; #print axioms hledger_pos_of_bounds = [propext, Classical.choice, Quot.sound]; registered, corpus green). Pure arithmetic (nlinarith) over the abstract ledger quantity (L, singlesSum, credits, F, E), framework-agnostic; the exact shape cite_hunter_lonely's hledger needs. HONEST: this is the ARITHMETIC that turns {singles-bound + pair-floor} into hledger for c<=7 -- it FINISHES the ledger assembly but does NOT prove the pair-floor (mac-mini/kps active) nor handle c>=8 (open, needs triples/scale). It reduces the c<=7 crux to exactly the pair-floor.
source: klein-2026-07-02-S117
depends_on:
  - HYP-4021   # S116: the path-Hunter inequality + ledger_coeff (this is the fee-aware discharge)
  - HYP-3980   # kps-S23: cite_hunter_lonely (whose hledger this discharges)
related:
  - HYP-3874   # mac-mini JointRateCore: the pair-floor (credits >= (c-1)(L/49)-err)
  - HYP-4020   # S115: the far-cut-7 dispatch (the hge7 leg cite_hunter_lonely feeds)
external: LRCCommensuration.lean (volume(danger P cap danger Q) = 1/49 exactly, the commensurate pair-floor)
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCLedgerAssembly.lean
---

# HYP-4022 — the ledger-positivity assembly (singles + pair-floor → hledger)

## The lemma (LRCLedgerAssembly.lean, sorry-free)
`hledger_pos_of_bounds (L singlesSum credits F E) (c : ℕ) (hL) (hc : c ≤ 7)`:
- `singlesSum ≤ c·(L/7) + F` — the **singles** upper bound (kps `teeth_mass`);
- `credits ≥ (c−1)·(L/49) − E` — the **pair-floor** (mac-mini JointRateCore per-cell obligation);
- `F + E < L·(48 − 6c)/49` — the fee budget below the path-Bonferroni credit;
⟹ `0 < L − singlesSum + credits` — exactly `cite_hunter_lonely`'s `hledger`.

Proof: `L − singlesSum + credits ≥ L − (c·L/7 + F) + ((c−1)·L/49 − E) = L·(48−6c)/49 − (F+E) > 0`, via the
`ledger_coeff` identity `1 − c/7 + (c−1)/49 = (48−6c)/49` (HYP-4021) and `nlinarith`. `credit_pos`: the credit
`L·(48−6c)/49 > 0` for `c ≤ 7`.

## Where LRC(14) stands now (the chain)
`LRC14Statement` ⟸ (klein S115 far-cut dispatch) ⟸ {citation} + {window census} + {≤6-far leg} + {≥7-far
leg}; the ≥7-far leg ⟸ (kps-S23 `cite_hunter_lonely`) given `hledger`; and `hledger` ⟸ (this assembly) given
{singles-bound + pair-floor} for `c ≤ 7`. So the whole proof reduces to **two genuinely-remaining analytic
things**:

1. **The pair-floor** `pairCredits ≥ (c−1)·(L/49) − err` for near-equal teeth. The **commensurate** case is
   already proven — `LRCCommensuration.lean` has `volume(danger P ∩ danger Q) = 1/49` **exactly**. The
   general/drifting case is mac-mini's JointRateCore per-cell / kps's pair-event run-gap analysis (active).
2. **The c ≥ 8 near-equal blocks** — the credit `48 − 6c ≤ 0`, so pairwise Hunter caps at 7. Crossing needs a
   **triple-Bonferroni** credit or the large-scale limit (fees → 0 with a sharper argument). Open.

## Honest scope
This **finishes the arithmetic assembly** of the `c ≤ 7` route — the exact bridge from the two analytic
bounds to `cite_hunter_lonely`'s `hledger`, sorry-free and framework-agnostic. It does **not** prove the
pair-floor or handle `c ≥ 8`; it reduces the `c ≤ 7` crux to exactly the pair-floor. The two remaining things
are now pinned precisely and each owned (pair-floor: mac-mini/kps active; `c ≥ 8`: open frontier).
