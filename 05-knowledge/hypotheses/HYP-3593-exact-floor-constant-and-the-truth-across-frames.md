---
id: HYP-3593
title: The exact floor constant inf R' = 114382/332563 = 0.34394 (binding R={1..13}\{7}, Q={1,2}; denominator 7^2*11*617) clears the clean set-independent bound 1/(2 zeta(2))=3/pi^2=0.304 by +0.040; and the meta-truth across all 16 'most fundamental' frames -- the LRC(14) floor is the POSITIVITY OF ONE CYCLOTOMIC OBSTRUCTION ATOM at the apex-7 cusp (the doublet = cusp form f_14 = 4cos^2(3pi/7)=0.198) surrounded by the bulk floor 3/pi^2; the frames are a monotone zoom-in localizing that single atom
status: VERIFIED (exact rational inf R' over a 942-covering scan; clears 3/pi^2 with margin) + META-SYNTHESIS (the convergence of all prior frames onto one positive cyclotomic atom). The true inf over ALL coverings and the literal f_14 bound remain open (floor owners).
source: klein-2026-06-29-S14
depends_on:
  - HYP-3571   # the proof sentence; inf R'=0.344 (here exact)
  - HYP-3581   # the apex atom 4cos^2(3pi/7)
  - HYP-3587   # genus = obstruction = cusp form f_14
related:
  - HYP-3550   # 3/pi^2 = 1/(2 zeta(2)) Euler-product floor (the bulk number)
  - HYP-3592   # blue = SC spine (where the atom sits)
  - polynomial-method / Gamma_0(14) (mac-mini HYP-3553)
results:
  - 04-computation/lrc14_exact_floor_constant_klein.py
  - 05-knowledge/results/lrc14_exact_floor_constant_klein.out
---

# HYP-3593 — the exact floor constant, and the one truth across all frames

## The bound (exact)

`R' = m_S/(m_R m_Q)`, `S = R u 14Q`. Over the covering scan (942 size-valid coverings):
> `inf R' = 114382/332563 = 0.343941`, at `R = {1..13}\{7}`, `Q = {1,2}` (denominator `= 7^2 · 11 · 617`).
It clears the clean set-independent floor `1/(2 zeta(2)) = 3/pi^2 = 0.303964` by `+0.040`. The exact value
is messy and `R`-specific (the `7^2` is the apex / `Gamma_0(14)` `d=7` cusp; `11, 617` are the speed-set
primes); the CLEAN, frame-independent bound is `3/pi^2`. The per-set value is a chart; the bound is the
truth.

Two clean numbers survive every reframing:
- **`3/pi^2 = 1/(2 zeta(2))`** — the bulk / Eisenstein / `zeta(2)`-density floor (the local-global density,
  HYP-3550).
- **`4 cos^2(3 pi/7) = 0.198`** — the apex-cusp / doublet / cusp-form obstruction atom (HYP-3581).
The floor is the statement that the second does not sink the first.

## The truth across the 16 frames

The LRC(14) floor was, in turn, declared most fundamental as: covering/Diophantine; the tournament
metagraph; the staircase triangle; the 2-adic descent; the CV variance gatekeeper; the `R`-eigenspace
(complement = antipodal); the Euler product / anti-Littlewood; `Gamma_0(N)` / finite Siegel transform;
relations-not-things; ESSENTIAL × BOUNDED; the transitive collapse; the finite cyclotomic min; the
`X_0(14)` cusps = Klein group; genus = local-global gap / cusp form; genus = odd boundary / even-graph;
blue = SC spine. Read in order they are a **monotone zoom-in on one obstruction**: continuum measure ->
variance over all sets (wrong/unbounded) -> per-level decorrelation -> finite min over `2^7` cores ->
a single doublet -> one cyclotomic number at the apex cusp.

> **THE TRUTH: the LRC(14) covering floor is the positivity of a single cyclotomic obstruction atom at the
> apex-7 cusp** (the doublet = the genus-1 cusp form `f_14` = the σ-odd `M_odd` = the blue odd-boundary on
> the SC spine = `4cos^2(3pi/7)`), surrounded by the bulk floor `3/pi^2` that every local chart computes.

Three fingerprints held across ALL frames: (1) a split into a computable bulk and a missing obstruction;
(2) the obstruction is ONE-dimensional at `N=14` (genus 1, one cusp form, one doublet); (3) it is cyclotomic
at the apex `7` (`Q(cos 2pi/7)`). The recurring error was mistaking each chart for the territory; the cure
was the persistence test (it killed `b_1^-=7`, `89=F_11`, the CV route; it kept `28=T(7)=C(8,2)`,
`P_n(-1)=SC`, the apex `Z_7`, the doublet). What remained, frame after frame, is one positive number.

## Closing statement of the arc

> LRC(14) holds because a single cyclotomic number at the apex-7 cusp is positive (`4cos^2(3pi/7) > 0`
> survives the bulk), and the bulk around it is `3/pi^2`. The remaining proof is the bound on that one atom
> (the leading apex-cusp coefficient of `f_14`, HYP-3587). Everything else is charts.
