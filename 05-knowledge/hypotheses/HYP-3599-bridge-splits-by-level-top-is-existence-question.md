---
id: HYP-3599
title: The bridge from the apex skeleton (THM-590) to the full per-level rho_j (THM-580) SPLITS by descent level -- deeper levels (j>=1) are genuine measure-decorrelations bounded away from 0 (verified rho_j >= 0.56, increasing with depth toward independence), but the TOP level rho_0 = m_S/(m_{O_0} m_{S'}) literally contains m_S and vanishes at the boundary (verified rho_0 = 0.34,0.23,0.19,0.17,0.16,...->0 as m_S->0); so THM-580 reduction (a) (rho_j>=c uniform) HOLDS for j>=1 but FAILS for j=0, where rho_0>0 <=> m_S>0 <=> LRC(S) is the original existence question; and the full measure-rho_j is NOT bounded by the apex gap (it binds at the cusp end |O7|->7, opposite the doublet skeleton |O7|=2)
status: VERIFIED (391 descent levels + targeted boundary/cusp sweeps). Corrects HYP-3600 consequence-3 (the bounded-product floor is false as a MEASURE statement) and localizes THM-580(a)'s failure to level 0.
source: klein-2026-06-29-S18
depends_on:
  - THM-580   # the descent + reduction (a) rho_j>=c
  - THM-590   # the apex skeleton
related:
  - HYP-3597   # measure vs existence; inf over infinite families
  - HYP-3600   # the finite families (consequence 3 corrected)
  - HYP-3575   # mac-mini: rho_j = the Z_7 Gram gap
  - HYP-3576   # mac-mini: Gamma_0(14) averaging
results:
  - 04-computation/bridge_skeleton_to_rho_klein.py
  - 05-knowledge/results/bridge_skeleton_to_rho_klein.out
  - 05-knowledge/results/bridge_verify_klein.out
  - 05-knowledge/results/bridge_levels_klein.out
---

# HYP-3599 — the bridge splits by level; the top level is the existence question

## Finding 1: the measure-rho_j is NOT bounded by the apex gap

Actual `rho_j` (THM-580) does not track `g(O_j mod 7)` (THM-590). min `rho_j` by gap value: `g=0 -> 0`;
`g=0.198(doublet) -> 0.44`; `g=0.308 -> 0.68`; `g=1 -> 0.19`; `g=2 -> 0.71`. The doublet (the skeleton's
binding core) does NOT bind `rho_j`. Conditioning on cusp-proximity instead: min `rho_0` falls monotonically
`0.91,0.65,0.53,0.33,0.14,0.10` as `|O_j mod 7| -> 7`. **The full `rho_j` binds at the cusp end, opposite
the doublet skeleton.** The apex cyclotomic gap governs the discrete content, NOT the measure `rho_j`.

## Finding 2: the obstruction is LEVEL 0, and level 0 = the existence question

min `rho_j` by level: `j=0: 0.05 (->0); j=1: 0.56; j=2: 0.77; j=3: 1.00; j=4: 1.07`. Deeper levels are
bounded away from 0 and increase toward independence (smaller cores, lonely-richer, more overlap). Only the
TOP vanishes, structurally: `rho_0 = meas(lonely S)/[meas(lonely O_0).meas(lonely S')]` -- the numerator is
`m_S`. Verified along boundary-approaching coverings: `rho_0 = 0.34,0.23,0.19,0.17,0.16,...->0` as
`m_S->0`. So `rho_0 -> 0` is the tautology that the top-level odd/even overlap is the original lonely set:
`rho_0 > 0 <=> m_S > 0 <=> LRC(S)`.

## The honest bridge

No single inequality `rho_j >= skeleton`. The bridge splits:
- **j >= 1: genuine bounded measure-decorrelations** (`rho_j >= 0.56`). THM-580(a) holds; the skeleton is
  satisfied with room to spare. The descent has real content here.
- **j = 0: the existence question.** `rho_0 = m_S/(...)` vanishes at the cusp; no `m_S`-independent measure
  bound can hold. The bridge is the MEASURE -> EXISTENCE passage (HYP-3597): need `rho_0 > 0` (overlap
  NONEMPTY), not `rho_0 >= c`. And `rho_0 > 0 <=> LRC(S)`. The descent RESTATES, does not reduce, the top.

So the descent peels off the easy deeper levels (genuinely bounded) and isolates the top-level overlap = the
original existence question. The apex skeleton (THM-590) is the right object for the top level, but only in
the EXISTENCE sense -- certify the overlap is nonempty, not bound its measure (which provably vanishes).

## Corrections

- **HYP-3600 consequence 3** (bounded product `meas(lonely S) >= (4cos^2(3pi/7))^d . cap^d`): FALSE as a
  measure statement. `rho_0` is not `>= 4cos^2(3pi/7)`; it `-> 0` at the cusp. Holds only for the discrete
  skeleton, not the measure.
- **THM-580 (a)** (`rho_j >= c` uniform): holds for `j >= 1` (verified `>= 0.56`); FAILS at `j = 0`
  (`rho_0 -> 0`). The failure is exactly the original `m_S` re-appearing at the top.

## Net

The bridge does not exist as a measure inequality; it splits. Deeper decorrelations are genuine and bounded
(the descent's real content); the top level is `m_S` itself -- the existence question, where the measure
vanishes and only the discrete cyclotomic skeleton (in the existence sense) can act. The descent isolates
the hard core to the top level; it does not dissolve it. Same measure/existence dichotomy as HYP-3597,
localized to the top of the descent.
