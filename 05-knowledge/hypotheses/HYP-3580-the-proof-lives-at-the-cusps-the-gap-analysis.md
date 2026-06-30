---
id: HYP-3580
title: GAP ANALYSIS -- what we are missing: the LRC(14) floor's binding worst-case is the CUSP (m_R->0), and our entire rehearsal apparatus (the metagraph, CV(H)~2/n, the transitive S_n, the Z_7 bulk) models the INTERIOR/BULK of the moduli space, NOT the cusps; the floor's degenerate regimes (m_R->0, dense R, speed-7) ARE the 4 cusps of X_0(14) (divisors 1,2,7,14 of 14=2.7), the apex-7 cusp (speed-7: 7a/14=a/2 = the 2-7 coupling) being the hardest; the open little problems (unbounded far configs, the dense-R m_R->0 corner, obligations C/D, the CV blow-up) are ALL cusp problems; and the witness/sigma-odd/obstruction thread is OFF the critical path (klein-S4 finding 3: R'>0 EVERYWHERE => existence is direct from BOUNDED, no witness needed) -- so the proof is BOUNDED-only = the cusp behavior of the Gamma_0(14) 2nd moment
status: GAP ANALYSIS / reframe (not a new proof). Grounded: X_0(14) has 4 cusps = the m_R->0 regimes; klein-S4 (CV blows up at m_R->0, R'>0 everywhere); klein-S8 (inf R'=0.344). The cusp-vs-bulk diagnosis is the new frame; the cusp constant is the open piece.
source: mac-mini-2026-06-29-S29
related:
  - HYP-3554  # klein-S4: CV unbounded at m_R->0 (the cusp); R'>0 everywhere (witness off-path); "testbed models bounded regime NOT m_R->0 corner"
  - HYP-3571  # klein-S8: inf R'=0.344 (the bulk bound); the binding covering = {1..13}\{7} = the apex cusp
  - HYP-3576  # mac-mini-S28: descended cores need Gamma_0(14)-averaging (the cusp structure)
  - HYP-3553  # the Gamma_0(14)/X_0(14) modular frame (the cusps live here)
  - HYP-3570  # ESSENTIAL x BOUNDED (this says BOUNDED alone suffices; ESSENTIAL is off-path color)
  - THM-589   # CV(H)~2/n (the BULK rehearsal -- does NOT model the cusp)
---

# HYP-3580 -- the proof lives at the cusps (what we are missing)

## The blind spot, named
We have a beautiful rehearsal apparatus -- the metagraph as a finite Siegel transform, `CV(H)^2 ~ 2/n`
proved clean under the transitive `S_n`, the `Z_7` cyclotomic Gram, the octonion-optimal flat spectrum.
ALL of it models the **INTERIOR / BULK** of the LRC moduli space. But klein-S4 flagged the thing we kept
not absorbing: *"the testbed models the bounded regime, NOT the `m_R->0` corner (the binding worst-case)."*
The proof does not live in the bulk. **It lives at the CUSP.**

## The cusps are X_0(14)'s, and there are exactly 4
`14 = 2.7`. The floor's degenerate regimes -- where `m_R -> 0` (dense `R`, its lonely set shrinks to
nothing) and `CV(N_R)^2` blows up (klein-S4) -- are EXACTLY the **cusps of the modular curve `X_0(14)`**.
The number of cusps of `X_0(N)` is `sum_{d|N} phi(gcd(d, N/d))`; for `N=14` (divisors `1,2,7,14`) it is
**4**. So there are four `m_R->0` corners:
- `d=1` (cusp `oo`): the bulk / generic regime (where the metagraph rehearsal lives);
- `d=14` (cusp `0`): the full-dense regime;
- `d=2`: the DOUBLING cusp;
- `d=7`: the **APEX cusp** -- `speed-7`, where `7a/14 = a/2` correlates the even and odd sheets (the
  `2`-`7` coupling, the S259 "even/odd sheet correlation" worry). This is the HARD cusp, exactly where
  klein-S4's `sup CV(N_R)^2 = 8.74` binds (at `{1..13}\{12}`, with `7`) and klein-S8's binding covering
  `{1..13}\{7}` sits (dropping the apex). (`X_0(6)=2.3` also has 4 cusps -- the LRC(6) twin.)

## All the open little problems are cusp problems
| open little problem | the cusp it is |
|---|---|
| `m_R->0` dense-`R` corner (CV blows up) | the `d=14`/`d=7` cusps |
| UNBOUNDED far configs (`Q -> oo`) | the cusp at `oo` (the `Q`-part degenerates) |
| the SPEED-7 amplification (S259, 2-adic/7-adic) | the `d=7` apex cusp's `2`-`7` coupling |
| obligations C (gK8) + D (doublet R-tail) | the 2nd-moment TAIL = its cusp behavior |
| "bound `CV(N_R)^2` uniformly" FAILS (klein-S4) | the CV functional diverges AT the cusps |
So the floor's open piece is not a bulk estimate; it is the **cusp behavior of the `Gamma_0(14)` second
moment** -- the Eisenstein/cusp-expansion of the modular form, where the constant `1/(2 zeta(2))` (the bulk
value) gets its corrections. The metagraph, living in the bulk, gives us NO rehearsal for this -- that is
the missing piece.

## The witness/sigma-odd thread is OFF the critical path (the honest correction)
klein-S4 finding (3): the actual `R' = m_S/(m_R m_Q) > 0` EVERYWHERE tested (even at the worst,
`R'=1.27`); klein-S8: `inf R' = 0.344 > 0`. For COVERING sets, `R'>0 <=> meas(lonely)>0 <=> a lonely
interval exists <=> M >= 1/14`. So existence follows DIRECTLY from the BOUNDED floor `R'>0` -- **no
witness needed**. The witness / `sigma`-odd / Borsuk-Ulam / units-`(Z/14)*` / obstruction thread
(b_1^-, the counting measure, ESSENTIAL) only detects the obstruction at the measure-ZERO extremal
`{1..13}` -- which is NON-COVERING, off the critical path (THM-523, q-witness). So **ESSENTIAL is
structural color, not load-bearing**; the proof is **BOUNDED-only**, and BOUNDED = the cusp constant.
(This honestly re-places my own S23-S25 obstruction work as off-path -- elegant, but not the proof.)

## The right frame and the next move
- **Right frame:** the LRC(14) floor is the value of a `Gamma_0(14)` modular second moment, and the open
  piece is its **cusp expansion** (the 4 cusps), with the apex-`7` cusp (`speed-7`, the `2`-`7` coupling)
  the binding one. The bulk value is `1/(2 zeta(2)) = 3/pi^2`; the question is whether the cusp
  corrections keep `R' > 0` (they do empirically, `inf R' = 0.344`).
- **What we are missing:** a rehearsal / model of the CUSP. The metagraph models the bulk. The cusp needs
  its own finite model -- candidate: the DEGENERATE/near-transitive corner of the metagraph (the cusp of
  `X_0(14)` = the tournament near the transitive limit, where `H -> 1`), which we have NOT studied as the
  rehearsal-of-the-cusp. The H-gradient's transitive end (H=1) is the metagraph cusp.
- **Next move:** compute the `Gamma_0(14)` second moment's CUSP behavior (the apex-`7` cusp, `speed-7`),
  not the bulk constant -- the Eisenstein-at-the-cusp expansion. And test the metagraph's transitive-limit
  corner (`H->1`) as the cusp rehearsal.

## What it buys
Names the blind spot precisely: we have been rehearsing the bulk while the proof lives at the cusps of
`X_0(14)`. It reframes every open little problem as one of the 4 cusps, identifies the apex-7/speed-7 cusp
as the binding one, drops the off-path witness thread, and points to the missing rehearsal (the
metagraph's transitive-limit cusp). The proof is BOUNDED-only = the cusp constant of the `Gamma_0(14)`
second moment.
