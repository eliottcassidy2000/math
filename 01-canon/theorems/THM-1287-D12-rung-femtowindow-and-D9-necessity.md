---
id: THM-1287
title: "[VERDICT PENDING THE RUN — placeholder filled at close] The D=12 rung at N=9699691 (binder 23, primorial 9699690, femto-window 5.3·10⁻¹⁵ at ~10⁷ speeds) and the D=9 necessity sweep (no in-window F_9(N) off the L_9 progression, odd N ≤ 3000)."
status: PENDING — filled from the stored out-files at session close.
source: death-star-2026-07-19-S59g (HYP-7925; owner: run the D=12 rung at N=9699691)
depends_on: [THM-1276, THM-1270, THM-1271]
scripts:
  - 04-computation/lrc_D12_rung_N9699691_deathstar_S59g.py -> 05-knowledge/results/lrc_D12_rung_N9699691_deathstar_S59g.out
  - 04-computation/lrc_D9_necessity_sweep_deathstar_S59g.py -> 05-knowledge/results/lrc_D9_necessity_sweep_deathstar_S59g.out
---

# THM-1287 — [pending]

Setup verified before the run: Q = 116396303, gcd(23, Q) = 1, N = 9699691 ≡ 16
≢ 1 (mod 23), N ≡ 1 (mod L_12 = 9699690 = 2·3·5·7·11·13·17·19); all nine lower
rungs D′ = 3..11 dead by the THM-1270 lemma (gcd(2D′−1, Q_{D′}) =
5, 7, 3, 11, 13, 15, 17, 19, 21 respectively — the D=11 composite binder 21
killed via 3 | Q_11, its factors in the primorial). Window width 5.31·10⁻¹⁵.

Results: [FILLED AT CLOSE]
