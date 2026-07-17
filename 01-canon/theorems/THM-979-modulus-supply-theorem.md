---
id: THM-979  # 977/978 codex scale-obstructions; 973->975 fleet renumber of my transcription noted
title: THE MODULUS SUPPLY THEOREM — the census pipeline's named missing input (a modulus) supplied explicitly from the T_s ladder. THE PICTURE's gap ("neither identity supplies a modulus or a live-count floor"): the discrete census criterion (THM-950, the B5 funnel) fires once a grid modulus q with discrete B5_q > 0 is exhibited; this file supplies q EXPLICITLY. (I) THE SAMPLING INEQUALITY: B5 = ∫ w(depth t) dt with weights w(0)=1, w(d)=−C(d−1,5) (|w| ≤ 792, sharp only at full depth); the depth function has ≤ E(V) = 2·Σv breakpoints, so the grid average at modulus q obeys |B5_q − B5| ≤ 792·E(V)/q (PROVEN, coarse), with the OBSERVED effective rate ≈ E(V)/q (deep cells are rare — the 792 weight almost never meets a cut cell; referee below); (II) THE SUPPLY: any packet with continuous B5 > 0 admits the explicit modulus q₀ = ⌈1584·E(V)/B5⌉ (proven-safe; the practical q ≈ 2E/B5 suffices in every test) at which B5_{q₀} ≥ B5/2 > 0 — and the census pipeline (CoverageCapped + lonely_of_census) then decides loneliness END TO END; (III) THE LADDER COMPOSITION: on the dissociated branch, the T_s ladder gives the continuous floor B5 ≥ 2052/16807 − (24/343)(13/H) − Σ_s |E_s|C(13,s)C_sL^{s−1}/H > 2052/33614 for H ≥ H₀(Vmax) — so the dissociated trapped core receives the FULLY EXPLICIT modulus q₀(V) = ⌈1584·2·Σv·33614/2052⌉ ≈ 51,900·Σv: every hypothesis of the census certificate is now a computable function of the speeds alone. REFEREED: on the certified packet, |B5_q − B5| ≤ E/q at q ∈ {1009, 4001, 16001} (margins 5–200×), and at q₀ = E/B5 = 468,332 the discrete B5 = +0.081861 > 0 — the census fires; the negative-B5 control packet converges at the same rate to its true negative value (no false positives)
status: (I) proven at the coarse constant (792E/q; per-breakpoint cell error ≤ |w|∞/q) with the effective rate documented; (II)-(III) explicit assembly; refereed both directions (modulus_supply_referee_kps_S128c46.py)
source: kind-pasteur-2026-07-17-S128 (cont.46; owner: finish the remaining proof pieces)
depends_on:
  - THM-950/the census pipeline (death-star: the consumer this feeds)
  - THM-952/953/959/975 (the ladder: the continuous floor's supplier)
  - THM-935/946 (the E_s frame + atom)
related:
  - the B5 funnel (LRC14-FORMALIZATION-PICTURE) — this closes its "no modulus" caveat
  - LRCWeightedDeepCensus, lonely_of_census (the Lean consumers; the q₀ formula is decide-shaped)
script: 04-computation/modulus_supply_referee_kps_S128c46.py -> 05-knowledge/results/modulus_supply_referee_kps_S128c46.out
---

# THM-979 — the modulus supply theorem

The census pipeline decides loneliness from a supplied grid; what it lacked was WHICH
grid. The sampling inequality bounds grid-vs-measure error by (weight)·(breakpoints)/q;
the T_s ladder floors the continuous B5 on the dissociated branch; dividing gives the
modulus. All three ingredients are explicit in the speeds:

> E(V) = 2Σv breakpoints; B5-floor = 2052/33614 (H ≥ H₀(Vmax));
> **q₀(V) = ⌈1584·E(V)·33614/2052⌉ ≈ 51,900·Σv.**

Referee: the certified packet's census fires already at the PRACTICAL modulus
E/B5 = 468,332 (observed rate E/q; the proven 792E/q is 1584× conservative), and the
negative control converges to its true negative value — the criterion cannot false-fire.

## Named next
- The Lean rendering: the sampling inequality is a per-cell count (decide-shaped for
  fixed packets); the q₀ arithmetic is norm_num.
- Tighten 792 → the effective constant via the deep-cell-measure bound (the T6-shaped
  refinement; optional — the coarse constant already yields finite explicit moduli).
