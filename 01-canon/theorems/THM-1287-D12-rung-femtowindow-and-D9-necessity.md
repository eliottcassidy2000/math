---
id: THM-1287
title: The tower reaches primorial 9699690 — M(F_12(9699691)) = 12/116396303 EXACTLY (D=12, binder 23, a FEMTO-width window 5.3·10⁻¹⁵ at ten million speeds, 482 seconds) — the EIGHTH consecutive out-of-sample confirmation of the D-graded gate law — and the D=9 NECESSITY sweep (1399 families, odd N ≤ 3000 off the L_9 progression, ZERO in-window): the gate law now has systematic support in BOTH directions.
status: >
  PROOF-BACKED EXACT (THM-1271 Lemmas 1-3 + THM-1002 + the THM-1270 ghost evaluator,
  spot re-gated: 12 sampled rows + F_7(2311)): M(F_12(9699691)) = 12/116396303, rung
  attained, in-window. PROVED per-target: the L1 floor (witness a = 15182127, min
  distance 12 verified directly over all ~10⁷ speeds) and the tower closures (all nine
  lower rungs dead, gcd(2D′−1, Q_{D′}) = 5,7,3,11,13,15,17,19,21 — the composite 21
  killed via 3 | Q_11, the THM-1270 lemma). VERIFIED-EXACT: the necessity sweep
  (zero in-window F_9(N) for 1399 odd N in [27,3000] with gcd(17,Q)=1 — no N in range
  is ≡ 1 mod 30030, and the branch mechanism predicts exactly this).
source: death-star-2026-07-19-S59g (HYP-7925; owner: run the D=12 rung at N=9699691)
depends_on: [THM-1276, THM-1270, THM-1271, THM-1002]
scripts:
  - 04-computation/lrc_D12_rung_N9699691_deathstar_S59g.py -> 05-knowledge/results/lrc_D12_rung_N9699691_deathstar_S59g.out
  - 04-computation/lrc_D9_necessity_sweep_deathstar_S59g.py -> 05-knowledge/results/lrc_D9_necessity_sweep_deathstar_S59g.out
---

# THM-1287 — primorial 9699690, and the law tested in both directions

## 1. The result

```text
F_12(9699691) = {1..9699689, 9699691} ∪ {116396280}     (~10⁷ speeds; x = 12·9699690)
M = 12/116396303  EXACTLY                                (Q = x + 23; 482 s ghost scan)
rung ATTAINED, strictly inside W = (1/9699692, 2/19399383), width 5.314·10⁻¹⁵.
Floor: a = 15182127 = 12·23⁻¹ mod Q — min distance exactly 12, all speeds, 4 s.
```

Gate: `9699691 ≡ 1 (mod L_12 = 9699690 = 2·3·5·7·11·13·17·19)`, `≡ 16 ≢ 1 (mod 23)`
— predicted OPEN, confirmed. **Eighth consecutive out-of-sample confirmation**;
confirmed primorials now `6, 30, 210, 2310, 30030, 510510, 9699690`
(D = 3, 4, 6, 7, 9, 10, 12; composite binders 9, 15, 21 skipped, each provably dead
at tower-N by the THM-1270 lemma — at this N the nine lower-rung gcds are
`5, 7, 3, 11, 13, 15, 17, 19, 21`, all > 1). The law now spans SIX orders of
magnitude in N with windows shrinking from 10⁻² to 10⁻¹⁵.

## 2. The necessity direction (lead xviii, first systematic pass)

`F_9(N)` computed for ALL 1399 odd `N ∈ [27, 3000]` with `gcd(17, Q) = 1`:
**zero in-window values**. No N in the range satisfies `N ≡ 1 (mod 30030)`, so this
is a clean falsification pass of the necessity half of the gate law — passed. (The
branch mechanism predicts necessity outright: each prime `ℓ ∈ {3,5,7,11,13}` kills
its competitor branch only when `N ≡ 1 (mod ℓ)`, and parity needs N odd.)

## 3. Standing of the tower, and what next

The D-graded gate law `{N ≡ 1 (mod L_D), N ≢ 1 (mod 2D−1)}` for the canonical
species is now: mechanism-derived (THM-1286 §3), three-quarters proof-backed per
instance (THM-1271 L1–L3 + finite ghost certificates), confirmed at 8/8
out-of-sample predictions across both the D-direction (tower diagonal) and the
N-direction (D=4 at three N; D=9 at three N), and falsification-tested in the
necessity direction. Remaining to full theoremhood: the general-N e-channel law
(per-N it is a finite certificate) and the even-N seal variant (THM-1271 §4).
The next rung, D=15 (binder 29; D=13's 25 = 5² and D=14's 27 = 3³ composite,
skipped), lives at `N ≡ 1 (mod 223092870 = 2·3·…·23)`: `N = 223092871`,
`Q = (N+1)·15 − 1 = 3346393064`... a ~2·10⁸-speed family — the ghost scan would
be ~3-6 hours: feasible but a deliberate cluster-scale choice, filed as a lead
with no urgency (the law's evidential state no longer turns on another rung).
