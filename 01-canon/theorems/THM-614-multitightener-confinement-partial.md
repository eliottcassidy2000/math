---
id: THM-614
title: (CONVERGENCE NOTE — SUPERSEDED BY THM-612 Lemma D.) Independent re-derivation of the m=2, f=2 confinement structure — the extremity dichotomy on the U-loose region and the "tighteners bounded by the even part" compactness — which mac-mini reached first (S32 anti-correlation + S33 Lemma D switch obstruction) and carried further to a finite per-U check. Retained only as an independent-convergence record; the canonical statement is THM-612 (Lemma D). Contributes an alternative closed-form compactness bound w_i ≤ 12·u_max and independent search confirmation.
status: CONVERGENCE / SUPERSEDED. The results here (extremity lemma; single-tightener components; w_i ≤ 12 u_max for f=2; f ≥ max(2, 7 meas R)) are PROVED and CORRECT but were reached independently and are subsumed by THM-612 Lemma D (mac-mini S32/S33), which additionally gives the switch-point divisibility w_i|w_j and the finite per-U reduction. Full confinement remains CONJECTURE. Do NOT cite this as a separate result; cite THM-612 Lemma D.
source: opus-2026-07-03-S61
depends_on:
  - THM-612   # Lemma D (mac-mini S33) — the canonical, stronger statement; this converges with it
related:
  - THM-613   # the Lipschitz slope idea behind the w_i ≤ 12 u_max form
  - HYP-4066  # opus-S61: the honest convergence record + non-closure
results:
  - 04-computation/lrc14_confinement_setup_opus_S61.py
  - 05-knowledge/results/lrc14_confinement_setup_opus_S61.out
---

# THM-614 — convergence note (superseded by THM-612 Lemma D)

**This is not an independent theorem.** Working the owner's "prove multi-tightener confinement" prompt, I
independently re-derived the `m=2, f=2` structure that mac-mini had already reached (S32 anti-correlation)
and extended (S33 **Lemma D**, now in THM-612). Recorded here only as an independent-convergence check;
**the canonical statement is THM-612 Lemma D**, which is strictly stronger.

## What converged (proved, but not new)
For primitive tight `S=2U∪F`, `q*=28`, `f=2`, on `R={g_U(2t)>1/14}` ((+1/2)-invariant):
- **Extremity dichotomy** — at every `t∈R`, exactly one tightener is `≤1/14`, the other `≥3/7`
  (`= ` mac-mini's S32 anti-correlation `{‖w_2t‖<1/14}={‖w_1t‖>6/14}` on `R`). Verified on 3728 points.
- **Single-tightener components**, types swapped by `+1/2`.
- **Compactness** — the tighteners are bounded by the even part. My closed form: applying the Lipschitz
  slope bound (THM-613) to the component at `U`'s global max (length `≥(M(U)−1/14)/u_max`, and `≤1/(7w_i)`),
  with `M(U)≥1/12`, gives `w_1,w_2 ≤ 12·u_max`. mac-mini's Lemma D gives the sharper
  `w_i < 4N/(7 meas R_U)` (`N=2·#lonely-intervals(U)`) plus the **switch-point divisibility**
  `w_i | w_j` (which I did not obtain) reducing `f=2` to a finite per-`U` check.
- `f ≥ 2` and `f ≥ 7·meas(R)` (elementary union/reflection bounds).

## Independent confirmation (the one genuinely useful artifact)
Exact-`M`/`q*` search over 938 structured even-block + odd-tightener families (`e=10,11,12`): **0**
primitive tight with `q*>14` — an independent check of mac-mini's confinement search on this slice.

## Honest status
Full confinement (`primitive tight ⟹ q*=14`) is **not proved** (THM-612): the residual is mac-mini's
"bound `v_max(U)`, the even part itself", plus all `m≥3`. This note claims no closure — see THM-612
Lemma D for the live, stronger line, and HYP-4066 for the convergence record.
