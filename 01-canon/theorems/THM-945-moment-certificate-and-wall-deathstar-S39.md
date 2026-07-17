# THM-945 — The optimal moment certificate, the moment wall, and the exact capped identity (death-star-2026-07-17-S39)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCMomentCertificates.lean`,
four standard-trio audits; verify the build report in the session log). Source:
HYP-7185. Closes THM-944's tail gap 328× → 1.077× and proves the remainder is not
moment-closable.

## Statement

1. `capped_cert_pointwise` (decide, 7 points): the OPTIMAL c ≤ 6 certificate,
   integer form: `30(C(c,3)+C(c,5)) ≤ 27 − 27c + 28·C(c,2) + 33·C(c,4)`, tight at
   the capped LP's vertex {1, 3, 4, 6} (λ₁ = −9/10 — certificate signs are free;
   the summation consumer is an equality substitution).
2. `B5_capped_floor`: on `CoverageCapped(6)` strata (no multiplier has ≥ 7 runners
   failing — the union-bound side of the 7-wall):
       **30·B5 ≥ 3(q−1) − 3S₁ + 2S₂ − 3S₄.**
   Equidistributed shortfall 0.0094(q−1) vs equilibrium 0.1221(q−1): THM-944's
   trivial route was 328× off; this is 1.077×.
3. `moment_wall`: EVERY certificate (any signs) valid on the full range c ≤ 13 has
   equidistributed floor ≤ 2052/16807 − 342/2401 < 0 — the explicit LEGAL histogram
   (450, 702, 1248, 1)/2401 on c = (0, 1, 3, 13) matches S₁/S₂/S₄ exactly with
   ledger −342/2401. **Moments alone cannot close the race; the coverage cap is
   exactly what closure requires** — and the cap is what the fragmentation/killer
   arc (THM-883, killer_box_thirteenth) provides on far strata.
4. `B5_eq_live_sub_deepSix` (**the exact identity**): on the capped stratum,
       **B5 = liveCount − #{p : bandCount = 6}.**
   Positivity is literally "live multipliers outnumber depth-six multipliers" — the
   7-wall in its cleanest discrete form. (At cap 5, B5 = liveCount identically:
   the race on capped strata is a COUNTING question, not an inequality question.)

## Referee

`moment_cert_referee_deathstar_S39.py`: pointwise PASS (tight {1,3,4,6}); capped
floor vs direct B5 PASS (200 families); wall witness exact; gap 328× → 1.077×.

## Remaining after this theorem

Supply `CoverageCapped(6)` on explicit far strata from the killer/fragmentation arc
(the S30–S32 ladder bounds exactly these simultaneous-failure counts) — wiring the
two arcs closes the loop; then the capped race is live-vs-deep-six census.
