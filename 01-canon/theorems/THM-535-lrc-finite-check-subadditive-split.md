---
id: THM-535
title: The subadditive cap-split of the LRC(14) seven-sector finite check — cap_k ≥ (k−6)/7 (PROVED), meas(S7(consec_k)) < (k−6)/7 for k≥9, reducing meas(S7(consec_k)) < cap_k to EXACTLY three genuinely-tight rational checks (k=8,9,10); plus the closed form Φ(c,k)=c/(k−1) for c≤1/2 (PROVED) and the cap_k minimizer family {1}∪{top cluster}
status: MIXED. PROVED: the subadditive cap lower bound cap_k ≥ (k−6)/7 (each of the 13−k speeds forbids EXACTLY measure 1/7, p disjoint arcs of half-width 1/(14p)); the single-arc orbit-confinement closed form Φ(c,k)=c/(k−1) for c≤1/2 (support is EXACTLY [0,c/(k−1)), elementary proof). VERIFIED (exact rational): meas(S7(consec_k)) < (k−6)/7 for k=9..13 (so closes via the proved cap bound); the pair-Bonferroni UB(1,6) < (k−6)/7 for k=11,12,13 (an analytic, true-cap-free closure of those rows); meas(S7) recomputed by independent direct sweep (no IE); cap_k minimizers exhaustive over {1..13}. NET: the per-k finite check meas(S7(consec_k)) < cap_k (k=8..13) is reduced to EXACTLY THREE genuinely-tight rational comparisons (k=8,9,10), the rest now a clean proved chain. Does NOT prove LRC(14) (still needs HYP-2603 "AP maximizes meas(S7)" + the upstream THM-527/HYP-2602 glue).
source: mac-mini-2026-06-18-S6 (Angle G: closed form / analytic finite check)
depends_on:
  - THM-532    # the seven-sector relation-height split (defines S7, cap_k, the finite check)
  - HYP-2603   # codex seven-sector net-cap reduction meas(S7)<=cap_k => HYP-2602
  - THM-530    # cap_k = min_{|P|=13-k} meas(G_P); the global-witness floor
related:
  - THM-533    # finer-cover certificate (the COMPLEMENTARY side: uniform corr≤C·W bound for general E). THM-535 instead handles the consecutive-cluster finite check directly; the proved cap_k≥(k−6)/7 bound here feeds BOTH (it is the RHS both must beat).
  - THM-534    # sector-moment LP-dual certificate
  - HYP-2602   # the 1/7-spread crux this finite check serves
  - THM-531    # AP-orbit scale/translation invariance (meas(S7) is scale-invariant)
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case, 13 speeds). Steinhaus three-gap theorem; Bonferroni inequalities; subadditivity of measure.
---

# THM-535 — The subadditive cap-split of the seven-sector finite check

THM-532/HYP-2603 reduce the LRC(14) crux to the **finite check**
`meas(S7(consec_k)) < cap_k` for `k=8..13`, where
`cap_k = min_{|P|=13−k} meas(G_P)`, `G_P={x:‖px‖≥1/14 ∀p∈P}`, and `S7(consec_k)` is the
seven-sector cover set of the consecutive cluster `E={0,1,…,k−1}`. Until now this was
**six per-`k` numerical comparisons**. Angle G makes it a graded, mostly-proved statement.

## A. The cap lower bound cap_k ≥ (k−6)/7 (PROVED, fully rigorous)

For a single speed `p≥1`, the forbidden set `{x:‖px‖<1/14}` is `p` arcs of half-width
`1/(14p)` centred at `j/p`, `j=0..p−1`. They are **disjoint** (centre spacing `1/p`,
half-width `1/(14p) < 1/(2p)`), so the forbidden measure is **exactly**
`p·(2/(14p)) = 1/7` for *every* `p≥1` (verified `p=1..29`, all `=1/7`). By subadditivity,
> `meas(G_P) ≥ 1 − Σ_{p∈P} 1/7 = 1 − |P|/7`,  hence  `cap_k ≥ 1 − (13−k)/7 = (k−6)/7`.

This bound is **TIGHT at `k=12` (`6/7`) and `k=13` (`1`)** (one/zero speeds, no overlap to lose).
Verified on 200 random `P⊆{1..13}`: zero violations.

## B. meas(S7(consec_k)) < (k−6)/7 for k≥9 — closes k=9..13 via [A]

Exact rational values (independent direct sweep + IE agree):

| k | meas(S7) | (k−6)/7 | slack |
|---|---|---|---|
| 8 | 481/1470 | 2/7 | **−61/1470 (FAILS — k=8 is special)** |
| 9 | 2447/5880 | 3/7 | 73/5880 |
| 10 | 8899/17640 | 4/7 | 1181/17640 |
| 11 | 3419/5880 | 5/7 | 781/5880 |
| 12 | 121103/194040 | 6/7 | 45217/194040 |
| 13 | 14573/21560 | 1 | 6987/21560 |

So for **`k=9..13`**: `meas(S7(consec_k)) < (k−6)/7 ≤ cap_k` (the cap bound is *proved*, [A]).
The binding row of this clean bound is `k=9` (slack `73/5880≈0.012`).

## C. Pair-Bonferroni closes k=11,12,13 WITHOUT the true cap (analytic)

`meas(S7) = 1 − meas(∪_{j=1}^6 A_j)`, `A_j={sector j empty}`. The Bonferroni pair bound
`meas(∪A_j) ≥ meas(A_a)+meas(A_b)−meas(A_a∩A_b)` gives, for any pair `(a,b)`,
`meas(S7) ≤ 1 − meas(A_a) − meas(A_b) + meas(A_a∩A_b)`. With `(a,b)=(1,6)`:
`k=11: 2881/4410`, `k=12: 27155/38808`, `k=13: 6995/9702` — each `< (k−6)/7`. (The pair-avoid
measures `meas(A_j), meas(A_a∩A_b)` are sector-confinement measures; the inequality itself is
rigorous, its terms exact rationals.)

## D. The residual = EXACTLY three genuinely-tight rational checks (k=8,9,10)

Combining [A]+[B]+[C], the only rows where `meas(S7)` is genuinely close to `cap_k` (and the
subadditive bound does not suffice) are **`k=8,9,10`**, three exact rational facts:
- `k=8: 481/1470 = 1924/5880 < 2243/5880 = cap_8` (margin `319/5880`; here `meas(S7) > 2/7`, so the *true* cap `P*={1,5,7,8,9}` is needed — k=8 is the unique anomalous minimizer).
- `k=9: 2447/5880 < 1979/4004` (margin `65669/840840`).
- `k=10: 8899/17640 < 55/91` (margin `22913/229320`).

## E. Closed forms (the building blocks)

- **Φ(c,k) := meas{x: frac(ix)∈[0,c) ∀i=0..k−1} = c/(k−1) for c ≤ 1/2 (PROVED).** The support
  is exactly `[0, c/(k−1))`: (i) `x<c/(k−1) ⟹ ix<(k−1)x<c≤1/2<1 ⟹ frac(ix)=ix<c`; (ii)
  `x∈[c/(k−1),c) ⟹ i*=⌈c/x⌉≤k−1`, `i*x∈[c,2c)⊆[c,1)`, so `frac(i*x)=i*x≥c` fails; (iii)
  `x≥c ⟹ frac(x)=x≥c` fails. Verified `k=3..29`, all `c≤1/2`. This gives the inclusion-exclusion
  `L=1,2,3` (i.e. `c=L/7≤3/7<1/2`) main terms `= L/(7(k−1))` in closed form.
- **cap_k minimizer family:** `P*={1}∪{top consecutive cluster}` for `k=9,10,11,12`
  (`{1,11,12,13}, {1,12,13}, {1,13}, {1}`), with `meas(G_{P*})` denominators
  `7·lcm(top cluster)`. **`k=8` is anomalous:** `P*={1,5,7,8,9}` (not a top cluster).
  `cap_13=1` (empty `P`).

## Honest status

This is an **engineering/structural** result on the *finite check*, not new analytic input to
LRC(14). It removes the "check six values" character: rows `k=9..13` follow from a **proved**
cap bound `cap_k≥(k−6)/7` plus the inequality `meas(S7)<(k−6)/7` (k=11,12,13 even have a
true-cap-free Bonferroni closure), leaving exactly **three** genuinely-tight rational facts
`k=8,9,10`. The full LRC(14) proof still needs HYP-2603's "consecutive (AP) maximizes
meas(S7)" (so that `consec_k` is the extremal cluster) and the upstream THM-527/HYP-2602 glue.
**LRC(14) NOT proved.** No elementary all-k closed form for `meas(S7(consec_k))` exists
(its `(1−meas(S7))·k` drifts `≈5.4→3.9`, a log-modulated three-gap quantity).

Files: `04-computation/lrc14_finitecheck_closedform_macmini_0618s7g.py`,
`05-knowledge/results/lrc14_finitecheck_closedform_macmini_0618s7g.out`.
