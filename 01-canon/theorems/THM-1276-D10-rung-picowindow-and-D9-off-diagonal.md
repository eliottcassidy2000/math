---
id: THM-1276
title: The tower reaches primorial 510510 — M(F_10(510511)) = 10/5105119 exactly (D=10, binder 19, a PICO-width window 1.9·10⁻¹², 19 seconds at half a million speeds), and the D=9 gate is confirmed OFF the primorial diagonal (F_9(60061) = 9/540557, F_9(90091) = 9/810827) — out-of-sample confirmations #5, #6, #7 of the D-graded gate law.
status: >
  PROOF-BACKED EXACT (THM-1271 Lemmas 1-3 + THM-1002 + the THM-1270 ghost evaluator,
  re-gated this run: 142/142 table rows + both big members): all three values. PROVED
  per-target: the L1 floors (closed-form witnesses, direct residue checks) and the
  tower closures (every lower rung dead via gcd(2D'−1, Q_{D'}) > 1, the THM-1270
  lemma). The gate law itself remains verified-at-7-hits + mechanism-derived.
source: death-star-2026-07-19-S59f (HYP-7915; owner: run the D=10 rung at N=510511)
depends_on: [THM-1270, THM-1271, THM-1257, THM-1002]
scripts:
  - 04-computation/lrc_D10_rung_N510511_deathstar_S59f.py -> 05-knowledge/results/lrc_D10_rung_N510511_deathstar_S59f.out
---

# THM-1276 — primorial 510510, and the gate off the diagonal

```text
F_9(60061)   = {1..60059, 60061, 540540}     M = 9/540557    (2 s;  width 1.4e-10)
F_9(90091)   = {1..90089, 90091, 810810}     M = 9/810827    (3 s;  width 6.2e-11)
F_10(510511) = {1..510509, 510511, 5105100}  M = 10/5105119  (19 s; width 1.9e-12)
```

All three rungs ATTAINED, all three strictly in-window, all three predicted by the
gate law `{N ≡ 1 (mod L_D), N ≢ 1 (mod 2D−1)}` before computation. The two D=9
values are the gate's first confirmations at MULTIPLE N per level beyond D=4
(60061 = 2·30030+1 ≡ 0 mod 17; 90091 = 3·30030+1 ≡ 8 mod 17): the law's N-degree
of freedom is real, not a diagonal artifact. The D=10 run makes SEVEN consecutive
out-of-sample confirmations (4/247, 6/1271, 7/16183, 9/270287, 9/540557, 9/810827,
10/5105119); confirmed primorials: 6, 30, 210, 2310, 30030, 510510.

Per-target the run verifies: the L1 closed-form witness (a = 31798 / 572349 /
268691; min distance = D over every speed, directly), the tower-closure lemma
(gcd(2D'−1, Q_{D'}) > 1 for every lower rung — all dead), and gcd(p, Q) = 1.
N = 510511 = 19·26869 and N = 60061 = 17·3533 are divisible by their own binders
(N ≡ 0 mod p) — harmless, as the lemmas require only N ≢ 1.

The ghost evaluator's scaling is now demonstrated linear-ish in N (1 s at 30k,
19 s at 510k speeds). The next rung, D=12 (binder 23; D=11's binder 21 = 3·7
composite, skipped by the THM-1270 lemma at tower-N), lives at
N ≡ 1 (mod 9699690 = 2·3·5·…·19), N ≢ 1 (mod 23): N = 9699691 ≡ 16 ≢ 1
(mod 23), gate predicted OPEN → M = 12/((9699692)·12−1) = 12/116396303; the
family has ~10⁷ speeds and the scan ~10⁸ candidates: an hour-scale deliberate
future run, not a casual one.

Scope: these are gap-membership results for the canonical species at enormous N —
structurally decisive for the cross-N law, untouched by and not touching the
N=12/13 walls except through the ghost-channel lens already handed to those
threads (THM-1271 §4, boxeph/opus handoffs).
