---
id: HYP-3824
title: (Task 2) THE Gamma_0(N)-LOCALIZED MOMENT LP recovers the lonely measure EXACTLY at the atom-aligned modulus Q=183=Phi_6(14), where the S91 GLOBAL LP gave min m_0=0 -- but inf meas AT r=1/14 is 0 (extremal AP sets are single-point-lonely), so the floor is the LINEAR SLOPE inf_S|L_S(r)| ~ 0.26(1-14r) -> 0, not a positive constant at the critical radius. (Task 1) INDEPENDENT rebuild of opus-S30's PSL(2,7) left-right Cayley complex CONFIRMS its dims (Z^1=176,B^1=160,H^1=16) but CORRECTS the expander claim: rank(delta^0)=160=168-8 => 8 CONNECTED COMPONENTS each of size 21=|Aut(Paley_7)|=3*7 (single left+right gens commute => component(g)=<a>g<b><=21) -- the single-generator complex is INHERENTLY DISCONNECTED, NOT an expander; string defects give soundness ~O(1/L)->0. Genuine LTC needs a GENERATING SET per side (LPS-Ramanujan).
status: MIXED -- both CONFIRMED computationally. Task2: the localization is EXACT (Gamma_0-aligned Q=183 recovers |L|), the single-point-lonely / linear-slope finding is exact (extremal AP sets), the 'inf meas AT r=1/14 = 0' is a clarifying correction to the floor target. Task1: the dims match opus; the 8-component disconnection is a definite CORRECTION/SHARPENING (rank(delta^0)=160 => 8 comps, verified two ways); soundness ~O(1/L) confirmed. The exact sqrt21 H^1 class identification is deferred to the narrow-class structure (opus/kps).
source: mac-mini-2026-07-01-S92
related:
  - HYP-3822   # S91 the facility-location game + the GLOBAL moment LP (min m_0=0) -- this LOCALIZES it
  - HYP-3817   # fixed-point instruments: covering/LOCAL moment not transform -- Task2 confirms (local works)
  - HYP-3571   # the arithmetic 2nd-moment floor (Gamma_0(N)); Task2 is its facility-location realization
  - HYP-3823   # opus-S30 PSL(2,7) complex + kps-S25 (this rebuilds/corrects: 8 components, not an expander)
results:
  - 04-computation/localized_moment_lp_gamma0_macmini_20260701.py
  - 04-computation/psl27_cocycle_soundness_macmini_20260701.py
  - 05-knowledge/results/localized_moment_lp_gamma0_macmini_20260701.out
  - 05-knowledge/results/localized_moment_lp_subcritical_macmini_20260701.out
  - 05-knowledge/results/psl27_cocycle_soundness_macmini_20260701.out
  - 05-knowledge/results/psl27_components_macmini_20260701.out
---

# HYP-3824 -- localizing the moment LP, and the PSL(2,7) cocycle soundness/connectivity test

Two owner directives: (1) localize the S91 moment LP with the `Gamma_0(N)` congruence constraint and test
whether `min m_0 > 0`; (2) encode the certificate as an explicit cocycle on the `PSL(2,7)` complex and test
soundness.

## Task 2 -- the Gamma_0(N)-localized moment LP
S91 (HYP-3822): the GLOBAL LP `min m_0` s.t. `(sum m_k=1, sum k m_k=A, sum k(k-1)m_k=Phi)` gave `0` (a
same-moment distribution with zero loneliness exists). **Localize:** partition `[0,1)` into `Q` arcs
`I_c=[c/Q,(c+1)/Q)`; impose the PER-ARC coverage mass `A_c = int_{I_c} C dt` (the LOCAL 1st moment = the
arithmetic data). Within an arc, `min|{C=0} cap I_c| = max(0, |I_c|-A_c)` (local union bound), so
```
    min m_0(Q) = sum_c max(0, 1/Q - A_c)      (+ 2nd-moment sharpening  ... + Phi_c/Cmax_c)
```
a VALID lower bound on `|L|`, positive iff some arc is UNDER-COVERED, and `-> |L|` as `Q -> inf`.

**Result (computed, `S={1..12,182}`, `r=1/14`, `|L|=0.02390`):**

| `Q` | aligned | union bound | potential bound |
|---|---|---|---|
| 1 (global) | -- | 0.00000 | 0.00000 |
| 61 | `Gamma_0` | 0.00904 | 0.02191 |
| **183=Phi_6** | `Gamma_0` (atom) | 0.02312 | **0.02390 = |L| exactly** |
| 183000 | -- | 0.02390 | 0.02390 |

So the GLOBAL LP `= 0` but the **atom-aligned localization `Q=183=Phi_6(14)` recovers `|L|` EXACTLY** -- the
localized moment LP is the right instrument, and the binding modulus `Phi_6(14)=183` is its natural scale.

**BUT the floor over sets is subtle -- `inf meas` at `r=1/14` is `0`.** The extremal AP sets `{1..13}`,
`{2..26}` achieve `M=1/14` with a **single-point (measure-zero) lonely set** (`|L|=0`); the localized LP
correctly returns `0` for them. So there is NO positive first-moment floor at the exact critical radius. The
floor is the **LINEAR SLOPE** as `r` drops below `1/14` (computed, the AP is the minimizer):

| `r/(1/14)` | 1.00 | 0.99 | 0.97 | 0.95 | 0.90 |
|---|---|---|---|---|---|
| `inf_S |L|` (AP) | 0.000 | 0.0026 | 0.0078 | 0.0129 | 0.0315 |

`inf_S |L_S(r)| ~ 0.26 * (1 - 14 r) -> 0` linearly as `r -> 1/14^-` (the localized `Q=183` bound tracks it
exactly). **This CLARIFIES the `inf meas` target:** it is not a positive constant at `r=1/14`; it is the
linear collapse rate to the atom. (Reconcile with kps `inf meas >= 1/36`: that must be a sub-critical radius
`r ~ 0.906/14` or a different normalization -- flagged for kps.)

## Task 1 -- the PSL(2,7) cocycle: soundness AND connectivity
Independently rebuilt opus-S30's complex: `G=PSL(2,7)` (order 168), LEFT gen `a` (order 7, heptagon), RIGHT
gen `b` (order 3, Eisenstein); edges `g-ag`, `g-gb`; squares `{g,ag,gb,agb}`. **CONFIRMED opus's dims:**
`V=168, E=336, F=168` (each edge in exactly 2 squares = a surface code), `dim Z^1=176, B^1=160, H^1=16`.

**CORRECTION / SHARPENING (the finding):** `rank(delta^0)=160=168-8` => the complex has **8 connected
components**, each of size **`21 = |Borel| = |Aut(Paley_7)| = 3*7`**. Reason: with a SINGLE left generator
and SINGLE right generator, left and right multiplication COMMUTE, so the component of `g` is exactly
`<a> g <b>` (`<= 7*3 = 21` elements) -- the single-generator left-right complex is **inherently
disconnected**, no matter the generators. A disconnected graph has spectral gap `0`, so it is **NOT an
expander** (opus-S30's reported 2nd eigenvalue must be per-component). 

**Soundness:** even per-component, a **string** (a path of edges) violates only its `O(1)` endpoint squares
but is `Theta(L)`-far from `Z^1`, so soundness `~ (viol/|F|)/(dist/|E|) ~ O(1/L) -> 0` (surface-code, **not
`O(1)`-sound**). So the single-generator apex "LTC" is **degenerate** -- disconnected AND weakly sound --
exactly as the abelian tiling cube was degenerate (opus-S28). **Genuine expansion + soundness needs a
GENERATING SET per side (LPS-Ramanujan) + tensor local codes**, confirming and sharpening opus's own caveat.
The pretty part: each component IS the Frobenius-`21` = `Aut(Paley_7)` group -- the `sqrt21`/apex symmetry
appears as the connected component itself. (The ad-hoc QR cochain I tried is NOT a cocycle (94/168 squares
violated) -- the exact `sqrt21` `H^1` class needs the narrow-class projection, deferred to opus/kps.)

## The synthesis across both tasks
Both tasks say the SAME thing that HYP-3817/HYP-3822 have been saying: the working instrument is **LOCAL**,
not global. Task 2: the GLOBAL moment LP fails (`min m_0=0`) but the LOCAL (per-arc, `Gamma_0`-aligned) one is
exact. Task 1: the GLOBAL/abelian complex is degenerate; a genuine LOCAL code needs a real expander
(generating sets). The floor and the certificate are both **local objects** -- a per-residue moment and a
per-square check on an expander -- and the naive global versions of each provably fail.

## Honest scope
Task 2: the localization exactness (`Q=183` recovers `|L|`) and the single-point / linear-slope structure are
exact computations; the "`inf meas` at `r=1/14` `= 0`" is a real clarification of the floor target. Task 1:
the dims match opus; the **8-component disconnection is a definite correction** (verified two ways:
`rank(delta^0)=160` and the double-coset count); soundness `~O(1/L)` confirmed. The exact `sqrt21` `H^1` class
and the LPS-Ramanujan upgrade are the deferred next builds (opus/kps thread).
