---
id: HYP-3550
title: The resonance structure of cut(s|prefix Q) (the cumulative surviving-resonance mass over denominators b<=Q) is CRACKED as a TOTIENT-WEIGHTED FAREY SUM sum_{b<=Q} phi(b) delta_b, factoring as a CONTINUED FRACTION (additive: prefix = Farey sequence F_Q = mediant tree; per-channel delta_b = three-gap convergents of the speed ratios) times an EULER PRODUCT (multiplicative: density = zeta(2) product prod_p(1-p^-2)=6/pi^2); the 2-adic descent (THM-580) is the PRIME-2 factor and the zeta(2)/Farey floor the ODD-PRIME factor of ONE prime Euler product; Euler-product positivity (each (1-p^-2)>0) is exactly the floor-positivity gatekeeper, the residual being quasi-independence (HYP-3129). NOT a continued root or power tower.
status: SYNTHESIS + verified structure (M=1/(smallest surviving b); Euler product -> 6/pi^2; Farey F_Q adds phi(Q) mediants/level; three-gap convergents). The EXACT factorization (vs lower bound) requires channel-independence = the open gatekeeper.
source: mac-mini-2026-06-29-S17
related:
  - HYP-3549   # the additive(Farey)/multiplicative(units) interface of the disproof boundary
  - HYP-3548   # the floor R'>0 gatekeeper (this gives its multiplicative structure)
  - THM-580    # 2-adic descent = the prime-2 factor of the resonance Euler product
  - THM-523    # M(S)=1/(smallest surviving b); covering = kill all b<=14
  - HYP-3129   # resonance-killing / quasi-independence = the channel-decoupling residual
external: zeta(2)=pi^2/6 Euler product; Steinhaus three-gap; Farey/Stern-Brocot; HYP-2856 zeta2 floor
results:
  - 04-computation/cut_prefix_resonance_structure_macmini_20260629.py
  - 05-knowledge/results/cut_prefix_resonance_structure_macmini_20260629.out
reflections:
  - the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md
  - the-obstruction-combining-duality-additive-mediant-vs-multiplicative-I.md
---

# HYP-3550 -- the resonance structure of cut(s|prefix Q), cracked

## The object
`cut(s | prefix Q)` := the cumulative surviving-resonance (lonely-neighborhood) mass of a speed set
`s` over witnesses `a/b` with denominator `b <= Q`. Resonance `b` is **killed** iff some speed `== 0
mod b`; `M(s) = 1/(smallest surviving b)` (verified: AP{1..13}->1/14, {1..11,13}->1/12,
covering{1..11,13,84}->1/15). The floor mass sits near the small-denominator Farey points, so
> `cut(s | Q) = sum_{b<=Q, b survives s} phi(b) * delta_b`  (a TOTIENT-WEIGHTED FAREY SUM),
`phi(b)` = the count of Farey points `a/b`, `delta_b` = the surviving half-width at each.

## The crack: which nested form?
It is a **CONTINUED FRACTION x EULER PRODUCT**, not a continued root or power tower. Two factors:

**MULTIPLICATIVE (Euler product).** The totient density `sum_{b<=Q} phi(b) ~ 3Q^2/pi^2` gives the
floor constant `3/pi^2 = 1/(2 zeta(2))`, and `1/zeta(2) = prod_p (1 - p^-2)` (verified partial product
-> `6/pi^2 = 0.60793`). So the resonance density is the **zeta(2) Euler product over prime resonances**;
`cut(s|prefix Q)` in multiplicative form is the **partial Euler product over `p <= Q`**. KEY: every
factor `(1 - p^-2) > 0`, so the product is POSITIVE -- **Euler-product positivity IS the
floor-positivity gatekeeper.**

**ADDITIVE (continued fraction).** The prefix `Q` is the **Farey sequence `F_Q`**: `|F_Q| = 1 + sum_{b<=Q}
phi(b)`, and `F_{Q-1} -> F_Q` inserts exactly `phi(Q)` new fractions, each the **mediant** of
Stern-Brocot neighbours (the continued-fraction tree). `cut(s|prefix Q)` accumulates one Farey LEVEL
per step = continued-fraction depth. The **per-channel width `delta_b`** is governed by the **three-gap
(Steinhaus) theorem**: the gaps of the orbit `{k*alpha}` are `||q_i alpha||` at the CONVERGENT
denominators `q_i` of the speed ratio's continued fraction (e.g. `pi=[3,7,15,1,292,...]`, convergents
`3/1,22/7,333/106,355/113`); the cut transition (maxgap crossing `1/7`) happens at a convergent, and
the prefix truncates the continued fraction at `q_i <= Q`.

## The unification: ONE prime Euler product (2-adic descent x zeta(2) floor)
The 2-adic descent (THM-580) `meas(lonely S) = prod_j rho_j * prod_j meas(lonely O_j)` is the
**prime-2 factor** of the resonance product (it strips the `b` by 2-adic valuation, level by level). The
`zeta(2)` Farey floor is the **odd-prime factor** `prod_{odd p}(1-p^-2)`. **Together they are ONE Euler
product over ALL primes** -- so the two separate floor tools (kps's descent + the zeta(2) density) are
the prime-2 and odd-prime halves of a single multiplicative resonance structure. The functional-equation
dual `zeta(-1)=-1/12` is the cap/Bernoulli side (the resonance-killing reflection), making `floor>=cap`
the two sides of one zeta across `s<->1-s`.

## What it buys (and the residual)
- It identifies the resonance structure as **multiplicative (Euler product)**, whose **positivity is the
  floor gatekeeper** -- reframing OPEN-Q-108 as "show the actual coupled `cut(s)` is bounded below by its
  Euler product," i.e. control the **channel quasi-independence** (HYP-3129's resonance-decorrelation).
- It is the **same additive(Farey)/multiplicative(Euler) interface** as the disproof boundary (HYP-3549:
  units vs Farey pins) -- the LRC lives on the additive-multiplicative seam at every level.
- HONEST: this gives the STRUCTURE and the lower-bound mechanism; the EXACT factorization holds only
  under channel-independence, which fails for coupled covering speeds -- that coupling is the open work.
