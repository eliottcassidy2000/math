---
id: HYP-4018
title: THE VERIFIED AXIOM-FOOTPRINT FRONTIER OF LRC(14) -- machine-checked (#print axioms) that the endgame surface is sorry-free and UNCONDITIONAL modulo exactly {LRC(<=13) citation} + {CoveringFarLonely 22}; + the W=24 band is verified closeable (all 184757 covering families in [1..24] lonely). LRC14AxiomAudit.lean (registered in the root module, corpus green, 8499 jobs) runs #print axioms on the endgame theorems. RESULT (machine-verified): lrc14_of_covering_far_22 depends on axioms [propext, Classical.choice, Quot.sound] + [winData22_ok, winData22_complete native_decide axioms] and NOTHING ELSE -- NO sorryAx. So LRC(14) (via LRC14Statement) is sorry-free with footprint = Lean's 3 foundational axioms + the 2 machine-checked W=22 census native_decides, taking as HYPOTHESES exactly (1) cite : LRCUpTo13 (the owner-sanctioned LRC(<=13) citation node) and (2) hcovfar : CoveringFarLonely 22 (the SINGLE remaining analytic hypothesis: covering families with an entry beyond the window are lonely). The generic assembly lrc14_of_covering_far_of_window depends on ONLY the 3 foundational axioms (pure logic glue). Certified infinite slices verified clean: coveringFar_blockFamily (V>=21, foundations + AP-tail native_decides), coveringFar_deepWell (![1..12,182], foundations + ladder-pack native_decide). W=24 CENSUS (Python, lrc14_w24_census_verify_klein.py): all 184757 covering 13-families in [1..24] find a G2/kernel-gate witness (0 failures, 26s); the 22<max<=24 shell = 153286 families all lonely => the W=24 band is CLOSEABLE and would shrink CoveringFarLonely 22 -> CoveringFarLonely 24. HONEST: LRC(14) is NOT yet fully unconditional -- the infinite far tail (max > W) is the genuine remaining crux; band extension shrinks it but never closes it; the closure route is the peel/rate descent (kind-pasteur peel20 + the proved rate_core), not band-pushing
status: VERIFIED (machine-checked footprint) + honest frontier statement. VERIFIED (LRC14AxiomAudit.lean, lake build 8499 jobs green): #print axioms of lrc14_of_covering_far_22 = {propext, Classical.choice, Quot.sound, winData22_ok/_complete native_decide}, NO sorryAx; lrc14_of_covering_far_of_window = 3 foundational axioms only; block/deepwell slices = foundations + their witness native_decides. VERIFIED (Python): W=24 census 184757 covering families all lonely (0 fails). HONEST: this AUDITS + quantifies the frontier; it does not CLOSE CoveringFarLonely (the infinite far tail). LRC(14) = {LRC<=13 citation} + {CoveringFarLonely 22}, sorry-free -- the precise, machine-verified state.
source: klein-2026-07-02-S113
depends_on:
  - HYP-4017   # S112: the W=22 band + lrc14_of_covering_far_22 (the surface this audits)
related:
  - HYP-3974   # kps-S17: lrc14_of_peel20 (the peel route = the closure of CoveringFarLonely)
  - HYP-4001   # the rate lemma (rate_core; the analytic engine for the peel)
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRC14AxiomAudit.lean
  - 04-computation/lrc14_w24_census_verify_klein.py
  - 05-knowledge/results/lrc14_w24_census_verify_klein.out
---

# HYP-4018 — the verified axiom-footprint frontier of LRC(14)

## The machine-verified footprint (LRC14AxiomAudit.lean, corpus green, 8499 jobs)
`#print axioms` on the endgame surface:

```
lrc14_of_covering_far_22 depends on axioms:
  [propext, Classical.choice, Quot.sound,                 -- Lean's 3 foundational axioms
   winData22_ok._native.native_decide.ax_1_1,             -- the W=22 census gate (machine-checked)
   winData22_complete._native.native_decide.ax_1_1]       -- the W=22 enumeration completeness
                                                          -- NO sorryAx.
lrc14_of_covering_far_of_window depends on axioms:
  [propext, Classical.choice, Quot.sound]                 -- pure logic glue only
coveringFar_blockFamily : foundations + AP-tail native_decides
coveringFar_deepWell    : foundations + ladder-pack native_decide
```

## The precise state of unconditional LRC(14)
> **LRC(14) is proved in Lean, sorry-free, with axiom footprint = {Lean's 3 foundational axioms} + {the 2
> machine-checked W=22 census native_decides}, taking as hypotheses exactly:**
> **(1) `cite : LRCUpTo13`** — the owner-sanctioned LRC(≤13) citation node (Rosenfeld / Trakulthongchai /
>     Sungkawichai–Trakulthongchai; enters as a named citation, never a sorry);
> **(2) `hcovfar : CoveringFarLonely 22`** — the SINGLE remaining analytic hypothesis: every covering family
>     with an entry beyond the window (`|vᵢ| > 22`) is lonely.

All other cases are discharged unconditionally: **non-covering** families are lonely by the denominator sieve
at `t = 1/q` (`sieve_one_div`, THM-523); **no-far** families (`|vᵢ| ≤ 22`) are lonely by the machine-checked
W=22 census (`hwindow22_closed`).

## Verified this session: the W=24 band is closeable
Python census (`lrc14_w24_census_verify_klein.py`): **all 184757 covering 13-families in `[1..24]` find a
G2/kernel-gate witness** (0 failures, 26s). The `22 < max ≤ 24` shell is **153286 families, all lonely**. So
building `LRCWindowData24.lean` (native_decide over `C(24,13)=2496144`) would shrink the remaining hypothesis
from `CoveringFarLonely 22` to `CoveringFarLonely 24` — retiring 153286 families.

## Honest scope — what remains
LRC(14) is **not yet fully unconditional**. The infinite far tail (families with `max > W`, for every `W`) is
the genuine remaining crux — **band extension shrinks `CoveringFarLonely W` but never closes it**. The closure
route is the **peel/rate descent**: kind-pasteur's `peel20` (a far family descends in far-count with loneliness
transport) fed by the proved `rate_core` (the far-element wrap-counting lemma), which would reduce
`CoveringFarLonely` to the window census unconditionally. That assembly is the fleet's critical path.

## Net
A machine-verified account of exactly what LRC(14) rests on: sorry-free, unconditional modulo the
owner-sanctioned LRC(≤13) citation + the single hypothesis `CoveringFarLonely 22`, with certified infinite
slices and a closeable next band. The audit pins the target; the peel/rate descent is the closer.
