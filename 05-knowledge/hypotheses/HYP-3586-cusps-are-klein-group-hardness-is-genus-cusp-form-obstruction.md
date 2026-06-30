---
id: HYP-3586
title: The 4 cusps of X_0(14) ARE the Klein four-group (Atkin-Lehner W(14)=(Z/2)^2 acts regularly) = the n=4 tournament classes (THM-584), with d=1 bulk/T, d=2 doubling/+, d=7 APEX-hard/-, d=14 Fricke-full/S; LRC(2p) hardness = genus(X_0(2p)) (=0,0,1,2,2 for p=3,5,7,11,13, JUMPING 0->1 at N=14: LRC(6) genus0 solved, LRC(14) genus1 first-hard); nu_2=0 iff apex=3mod4 iff Paley; and the floor's obstruction is CONJECTURALLY the genus-1 weight-2 cusp form f_14 (the rank-0 curve 14a) component of the Gamma_0(14) 2nd moment (Eisenstein bulk (+) cusp form), which genus-0 LRC(6) lacks
status: VERIFIED modular data (genus/cusps/Atkin-Lehner/nu_2 by standard formulas; cusps=Klein structural). The hardness=genus and cusp-form-obstruction are well-grounded CONJECTURES/reframes. Complements mac-mini S30 (HYP-3585, the Z_7-core gap landscape) with the modular-curve geometry.
source: klein-2026-06-29-S10
depends_on:
  - HYP-3580   # mac-mini: the proof lives at the cusps of X_0(14)
  - HYP-3581   # klein: the floor = finite cyclotomic min 4cos^2(3pi/7) at the apex cusp / doublet
related:
  - HYP-3585   # mac-mini S30: the 128-core Z_7 gap landscape (arithmetic side); this is the geometric side
  - THM-584    # the n=4 Klein four-group = the cusps
  - HYP-3547   # Mersenne ∩ Heegner ∩ 3mod4 = {3,7} (the arithmetic "why 14"); genus is the geometric "why"
  - HYP-3553   # metagraph = finite Siegel transform / Gamma_0(N) (the automorphic frame)
  - THM-578    # the doublet = the apex cusp (width 2)
results:
  - 04-computation/lrc14_X0_cusps_atkinlehner_klein_klein.py
  - 05-knowledge/results/lrc14_X0_cusps_atkinlehner_klein_klein.out
---

# HYP-3586 — the cusps are the Klein group; hardness is the genus; the obstruction is a cusp form

## Verified (standard modular formulas; script attached)

- **4 cusps of X_0(14) = the Klein four-group.** `omega(14)=2 => W(14)={1,W_2,W_7,W_14}=(Z/2)^2` acts
  REGULARLY on the 4 cusps `{d=1,2,7,14}` (widths `14,7,2,1`, summing to `psi(14)=24`). So the cusps are a
  `V_4`-torsor = the n=4 iso-class structure `{T,+,-,S}` (THM-584). Dictionary: `d=1`(width14)=BULK=`T`
  (identity); `d=2`(7)=DOUBLING=`+`=`W_2`; `d=7`(2)=**APEX-7 hard**=`-`=`W_7`; `d=14`(1)=FULL-DENSE/covering
  =`S`=Fricke `W_14`. The "two order-2 structures" (parity vs doubling) are `W_7` and `W_2`; the complement
  is Fricke `W_14`. The binding **doublet** (HYP-3581/3585, `4cos^2(3pi/7)`) is the width-2 apex cusp.

- **genus(X_0(2p)) = 0,0,1,2,2** for `p=3,5,7,11,13` — JUMPS `0->1` at `N=14`. Among LRC-relevant apices
  `{3,7}`: `genus(X_0(6))=0` = LRC(6), SOLVED; `genus(X_0(14))=1` = LRC(14), FIRST HARD. Reframe (conjecture):
  **LRC(2p) hardness = genus(X_0(2p))** — the geometric companion of HYP-3547's arithmetic "why 14".

- **`nu_2(X_0(2p)) = 0 iff apex p ≡ 3 mod 4`** (0 for p=3,7,11; 2 for p=5,13) ⟺ Paley exists ⟺ the
  Borsuk-Ulam pillar (THM-581). The order-2 elliptic-point count IS the 3-mod-4/Paley condition.

## Conjecture (the cusp-form obstruction)

`X_0(6)` is genus 0 (rational, only Eisenstein series) — LRC(6) closes with the bulk. `X_0(14)` is genus 1
= the rank-0 elliptic curve `14a`, carrying a nontrivial weight-2 cusp form `f_14`. CONJECTURE: the
`Gamma_0(14)` second moment controlling the floor decomposes Eisenstein(bulk) ⊕ cusp-form(`f_14`); the
piece the metagraph/`CV(H)`/transitive rehearsal cannot see (the cusp, not the bulk — klein-S4) is the
`f_14` component, concentrated at the `d=7` apex cusp. Rank 0 (`L(f_14,1)≠0`) makes it a fixed
non-degenerate obstruction. So "what we're missing" = a genus-1 cusp form at the apex cusp.

## Reframes generated (the grid; ★ = big shift)

- floor = spectral gap of `D*D` (danger relation composed with itself) at the apex cusp = `4cos^2(3pi/7)` ★
- cusps = Klein `V_4` = n=4 classes ★★ (loop closure to THM-584)
- floor 2nd moment = Eisenstein ⊕ cusp form `f_14`; obstruction = the genus-1 cusp form ★★
- 4 cusps = "places"; floor = adelic Euler product (HYP-3550); hardness = local-global obstruction
  (genus 0 = Hasse/easy, genus 1 = elliptic-curve/Sha-like/hard) ★
- LRC(2p) hardness = genus(X_0(2p)) ★★
- descent = flow `W_2`-ward terminating at the apex cusp `d=7` ★
- doublet = planar-difference-set MINIMUM (worst) vs Fano QR `{1,2,4}` = flat/optimal (mac-mini S27/S30) ★

## Abnormalities to track

1. genus jump `0->1` at 14 (the hardness invariant). 2. cusps = Klein (arc closure). 3. `nu_2=0 ⟺ Paley`.
4. `14a` rank 0 (non-degenerate obstruction). 5. **the "6" cluster — RUN THE PERSISTENCE TEST**: `phi(14)=6`
(witnesses), 6 curves in the `14a` isogeny class, the cuspidal-group order — likely DIFFERENT 6's (a TRAP;
do not conflate without checking persistence across `N`). 6. the doublet = the width-2 apex cusp.

## The shift

Treat LRC(14) not as analysis-with-side-structure but as the **arithmetic geometry of the rank-0 genus-1
curve X_0(14)**: cusps = Klein group, Atkin-Lehner = the two order-2 structures, genus = hardness, cusp
form `f_14` = obstruction. The metagraph, antipodal map, signed cycle index, descent, doublet are charts on
this one curve. The floor closes when `f_14` is bounded at the apex cusp.
