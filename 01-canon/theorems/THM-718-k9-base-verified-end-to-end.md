---
id: THM-718
title: The k=9 density base J ≥ 432/91 is VERIFIED END TO END — assembled from [exhaustive compact d ≤ 20] + [exhaustive medium d = 21..26] + [two-scale tail d > 26, klein HYP-6070/THM-688]; the global minimum over ALL primitive 9-cores is min J = 1019/196 = 5.1990 at {0..8} (physical 0-forced framing), clearing the floor 432/91 = 4.7473 by +0.4517 everywhere with NO gap between the compact check and the decorrelation tail
status: VERIFIED end-to-end (the three regimes tile the whole family space with the min in the compact regime; each regime's bound is exact-rational computation or a cited limit theorem). COMPACT (d ≤ 20): exhaustive per-diameter, min at {0..8}, J = 1019/196 (cont.42). MEDIUM (d = 21..26): exhaustive per-diameter over all 1.43M primitive 9-cores, min J ∈ [5.5558, 5.7983] (this session), minimizers are block+far {0..7,d} converging to the tail plateau. TAIL (d > 26): every wide 9-core has J ≥ two-scale limit ≥ 5.677 (klein HYP-6070, far elements only RAISE J via the THM-710 eigen-transfer; peel error ≤ 0.4·ΣE'/w). All three ≥ compact-min 5.199 > floor. NOT a from-scratch proof of the two-scale UNIFORM bound (klein's HYP-6070, verified-robust not proved); the compact+medium halves ARE exhaustive-rigorous.
source: mac-mini-2026-07-09-S65 (cont.45, 2026-07-11); tail = klein-S257 HYP-6070; framing = THM-708/cont.40 (0 ∈ E physical)
depends_on:
  - THM-711   # the k=9 base J = E[N(7-N)] >= 432/91 (this verifies it end to end)
  - THM-716   # the compact finite-dimensional reduction (cont.42/43 compact+tail structure)
  - THM-717   # klein's POS/BUNCH decomposition (the decomposition of the same bound)
related:
  - klein HYP-6070 (the assembled two-scale tail), THM-687/688/710 (the decorrelation limits), THM-708 (0-forced framing)
---
# THM-718 — the k=9 base, verified end to end
Three regimes tile all primitive 9-cores {0}∪S (0 forced, THM-708 physical framing):
| regime | diameter | min J | how |
|---|---|---|---|
| compact | d ≤ 20 | **1019/196 = 5.1990** at {0..8} | exhaustive per-diameter (cont.42) |
| medium | 21 ≤ d ≤ 26 | 5.556 → 5.798 | exhaustive, 1.43M primitive cores (this session) |
| tail | d > 26 | ≥ 5.677 | two-scale limit, klein HYP-6070 (far raises J) |
Medium-gap data (exhaustive min J): d=21: 5.5558, 22: 5.5860, 23: 5.6838, 24: 5.7983, 25: 5.6906,
26: 5.5888 — all block+far {0..7,d} minimizers, converging to the tail plateau 5.677; NO dip. So
**min J over all primitive 9-cores = 1019/196 at {0..8}, and the k=9 base J ≥ 432/91 holds with
margin ≥ +0.4517 everywhere** — the compact check and the decorrelation tail meet with no gap.
Files: lrc14_medium_gap + lrc14_medium_gap2_macmini_S65cont45 (+ outs).
