# HYP-2738 — Consec-max of measS7 is IRREDUCIBLY AGGREGATE, not a single monotone lemma (THREAD B)

**kind-pasteur-2026-06-21.** Status: **REFUTED** (the monotone-lemma reduction) **+ PARTIAL clean fact**.

## The question (THREAD B, from opus Thread B + HYP-2735)
The LRC(14) cover bound's worst case is the consecutive run; the open crux is
"consec maximizes `measS7(E)=p_0(E)` among k-sets."  HYP-2735 found additive energy
`E(A)=sum_s r_A(s)^2` is the dominant scalar (`corr +0.58`).  TARGET (LEMMA-1): is
`measS7` MONOTONE increasing in additive energy (so consec-max reduces to a clean
variational/compression lemma)?

## Verdict: NO clean monotone reduction.

### (a) Single-step compression monotonicity — REFUTED.
Of 12165 add-energy-INCREASING single-element moves (random k=8 sets), **3201 STRICTLY
DECREASE `measS7`**.  Worst: `[2,3,4,5,8,12,13,14] -> [2,3,4,5,7,12,13,14]` raises
add-energy by 8 but drops `measS7` by 0.111.  So add-energy is not even locally aligned.
Script `04-computation/lrc_q108_threadB_addenergy_monotone_kpswf4.py`.

### (b) Within-fixed-span monotonicity — REFUTED.
Grouping all k=8 sets of `{1..15}` by span, there are **6300 inversions** (a higher-energy
set with strictly lower `measS7`) across 7 of 8 spans.  Add-energy is NOT a monotone order
even after fixing span.  (Codex independently corroborated with a strict inversion: AE=180
row `p0=0.028` < AE=176 row `p0=0.069`, `lrc14_forced_zero_energy_tournament_codex`.)

### (c) Rank correlation — moderate, NOT monotone.
Kendall-tau-b(ish) `(measS7, add_energy) = +0.43` over all 3432 primitive k=8 sets;
lexicographic `(add_energy, -span)` gives `+0.41`.  Positive but far from a total order.

## The structural reason (the genuinely new content) — survival-basis corner decomposition
The functional that ACTUALLY closes the cap is L_y (THM-534 moment-LP dual), NOT measS7
or add-energy.  Write `N = #missed inner sectors` (sector 0 always hit), `p_t=P(N=t)`.
The L_y weight vector decomposes (opus routeB, exact) as
```
   L_y(E) = 1 - E[N] + sum_{a=1..5} w_a C_a(E),   C_a(E) := E[(N-a)_+]   (deep-miss corners).
```
EXHAUSTIVE per-corner extremality (`04-computation/lrc_q108_threadB_corner_monotone_kpswf4.py`,
`..._clean_certificate_kpswf4.py`):

| corner | k=8 | k=9 | k=10 |
|---|---|---|---|
| `C_1=E[(N-1)_+]` | 1139 exceed consec | 481 exceed | 102 exceed |
| `C_2=E[(N-2)_+]` | 5 exceed | **consec-MAX** | **consec-MAX** |
| `C_3=E[(N-3)_+]` | **consec-MAX** | **consec-MAX** | **consec-MAX** |
| `C_4=E[(N-4)_+]` | **consec-MAX** | **consec-MAX** | **consec-MAX** |
| `C_5=E[(N-5)_+]` | **consec-MAX** | **consec-MAX** | **consec-MAX** |
| `E[N]` | NOT consec-min (1 lower) | NOT min (27 lower) | NOT min (106 lower) |

**Clean fact (PARTIAL, useful):** consec UNIQUELY maximizes the DEEP corners
`C_a = E[(N-a)_+]` for `a>=3` (and `a>=2` at k=9,10) — consec produces the heaviest
deep-simultaneous-miss tail.  This is monotone and exhaustively verified.  Consec is also
PARETO-OPTIMAL on `(-E[N], C_3, C_4, C_5)` (0 weakly-dominating shapes).

**Why this does NOT close the cap (the irreducibility):** the L_y certificate puts
NONZERO weight `w_1=1` on `C_1` where consec is NOT max, and a NEGATIVE weight `w_3=-1/5`
at k=8.  The consec-max of L_y is a SIGNED net balance: at k=8 the worst C_1-exceeding
shape `[0,1,2,3,7,8,9,10]` has C_1 surplus +0.258, but the `-E[N]` term contributes -0.494,
giving total `L_y(shape)-L_y(consec) = -0.228 < 0`.  Consec still wins, but only after the
deficit is outweighed.

**Impossibility of a clean nonneg certificate (T2/T4, decisive):** any test function
`phi` with `phi(0)=1, phi(t)>=0` gives an UPPER bound `p_0 <= E[phi(N)]`.  The clean family
`phi_clean = p_0 + sum_{a>=3} C_a` (nonneg, shape-dependence only via consec-MAX corners) is
itself MAXIMIZED at consec — i.e. the bound is WEAKEST (loosest) at consec, so it CANNOT
certify consec is the `p_0`-argmax.  To certify consec-max you must SUBTRACT a consec-max
quantity, which forces the alternating (inclusion-exclusion / Bonferroni) sign pattern.
The mixed signs of L_y are therefore STRUCTURALLY FORCED, not an artifact.

## Honest status
- **REFUTED:** "measS7 monotone in additive energy" (locally, within-span, globally) — no.
- **REFUTED:** "consec-max reduces to a single monotone-functional lemma" — no clean
  nonneg certificate can certify it; the closing functional L_y is irreducibly signed.
- **PARTIAL PROVED-by-exhaustion (clean rung):** consec uniquely maximizes the deep-miss
  corners `C_a, a>=3` (k=8..10).  This is the largest piece that IS monotone; it alone
  does not close the cap but is a rigorous component for any future proof.
- The remaining content (consec maximizes the SIGNED L_y) is exactly the Delsarte-LP
  saturation statement (HYP-2726): consec = the LP-tight code.  That is genuinely
  aggregate (a complementary-slackness / extremal-code fact), consistent with mac-mini's
  "plain Bonferroni FAILS, optimal LP closes" (THM-534).

## Files
- `04-computation/lrc_q108_threadB_addenergy_monotone_kpswf4.py` (+ `.out`) — (a)(b)(c).
- `04-computation/lrc_q108_threadB_corner_monotone_kpswf4.py` (+ `.out`) — corner extremality.
- `04-computation/lrc_q108_threadB_clean_certificate_kpswf4.py` (+ `.out`) — impossibility.

-> HYP-2735, HYP-2726 (Delsarte LP saturation), THM-534, THM-559, opus routeB, OPEN-Q-108.
