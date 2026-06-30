---
id: HYP-3552
title: METAGRAPH -> LRC transfer — the metagraph's mult(1)=0 / mult(2)=1 ("no first-order invariant, a single quadratic mode = cyclicity") is the structural reason the LRC cap's binding content is purely the PAIRWISE SECOND MOMENT S_2 (level 2) with the first moment forced to vanish; the analytic engine for that single mod-14 second moment is a congruence-conditioned Siegel-transform second-moment formula (arXiv:2507.05905)
status: OPEN / transfer hypothesis (klein-2026-06-29-S3) — the metagraph side is proved (THM-588); the LRC side (S_1=0, S_2 = whole binding content) is to be verified on mac-mini's cap matrix M, and the Siegel second-moment route is a proposed evaluation
source: klein-2026-06-29-S3
depends_on:
  - THM-588   # metagraph: mult(1)=0, mult(2)=1, Fiedler = cyclicity, gap=4/d
  - THM-587   # signed cycle index; P_n(-1)=SC (antipodal Euler number)
related:
  - HYP-3538  # the LRC cap M = M_even (+) M_odd; S_2 = -tr(M_odd)/2 is the obstruction
  - HYP-3543  # one R, three spectra (metagraph = cap = witness)
  - HYP-3551  # the LRC is anti-Littlewood: positive floor on the simultaneous-approx product
  - lrc-walk-on-metagraph-proof-attempt-s519  # the runner walk IS a closed metagraph walk
references:
  - arXiv:2507.05905  # Han-Lee, moment formulas of Siegel transforms with congruence conditions (dim 2)
---

# HYP-3552 — The LRC obstruction is the level-2 (pairwise) second moment; evaluate it via congruence-conditioned Siegel second moments

## The structural claim (metagraph side: PROVED)

THM-588: the dominance-reversal metagraph has **no first-order invariant** (`mult(1)=0`) and a **single
quadratic invariant** (`mult(2)=1`, the cyclicity / 3-cycle count, the Fiedler mode). So the lowest
nontrivial content of the metagraph — the slowest mode, the spectral gap, the whole "near-equilibrium"
obstruction — is one PAIRWISE (degree-2) quantity. The linear term simply does not exist; the binding is
purely quadratic.

## The transfer (LRC side: to verify)

By mac-mini's "one R, three spectra" (HYP-3543), the LRC cap co-emptiness matrix `M` shares the metagraph's
`R`-even/`R`-odd structure. The metagraph's no-linear/single-quadratic law predicts, for the cap:

1. **`S_1 = 0`** (the first moment / linear functional of the covering deficit vanishes) — the analogue of
   `mult(1)=0`. mac-mini already finds the FLOOR is purely `R`-even (the `R`-odd/`SIGN` part of the lonely
   *measure* vanishes); HYP-3552 sharpens this to: the linear term is identically zero by the same
   transposition/antipodal mechanism that kills `mult(1)`.
2. **`S_2` is the entire binding content** — the pairwise co-emptiness `S_2 = -tr(M_odd)/2` (HYP-3538) is
   the cap's cyclicity-analogue: a single quadratic/pairwise mode. The cap closure should therefore reduce
   to bounding ONE second moment, not an infinite hierarchy.

## The analytic engine (proposed): congruence-conditioned Siegel second moments

The cap `S_2` is a **second moment** (a variance / pairwise correlation of the covering count), and the
LRC lives mod `14 = 2*7` (congruence structure). arXiv:2507.05905 (Han–Lee) gives **first and second
moment formulas for Siegel transforms with congruence conditions** in dimension 2 — exactly second moments
of lattice-point counts restricted by `mod q` conditions. Proposal: realize the LRC pairwise co-emptiness
`S_2` as such a congruence-conditioned second moment (the speed set as a lattice, the danger zones as the
counted region, `mod 14` as the congruence), and evaluate/bound it by the Siegel second-moment formula.
This couples:
- the **metagraph** (says: the obstruction is level 2 = a single pairwise second moment, `S_1=0`),
- the **Siegel transform** (computes that second moment with the `mod 14` congruences),
- the **anti-Littlewood floor** (HYP-3551: positive floor on the simultaneous-approx product) as the
  positivity the second-moment bound must deliver.

## Why this is worth testing

It would collapse the LRC cap from "bound `M_odd` somehow" to a single, named, classical-analytic object:
a congruence-conditioned second moment. The metagraph's role is **diagnostic** — it proves (in a clean,
fully-computable sibling system) that there is no linear obstruction and exactly one quadratic one, so the
proof effort belongs entirely at the second-moment level.

## Next steps

1. On mac-mini's cap matrix `M`: confirm the first-moment / linear functional vanishes (`S_1=0`) and that
   `S_2 = -tr(M_odd)/2` is the sole binding quantity (no independent higher-moment obstruction).
2. Write the LRC pairwise co-emptiness as a 2-dimensional Siegel transform second moment with `mod 14`
   congruences; compare its main term + error to the `S75e` cyclotomic-cosine SOS bound.
3. Use THM-588's exact `gap = 4/d` and Fiedler=cyclicity as the model for "the binding mode is the lowest
   quadratic"; check whether the LRC per-level `rho_j` (THM-580 descent) is likewise gap-controlled by its
   lowest quadratic level.
