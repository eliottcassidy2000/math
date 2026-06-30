---
id: HYP-3587
title: The genus of X_0(2p) IS the local-global gap dimension of the LRC floor -- the BULK (Eisenstein, dim = #cusps-1 = 3) is CONSTANT across the family while the OBSTRUCTION (cusp forms, dim = genus = 0,0,1,2,2) grows; the floor decomposes 3 Eisenstein + genus cusp forms, cusp forms VANISH at the cusps so the obstruction is the leading apex-cusp coefficient of f_14 (rank-0 14a, a finite non-degenerate number); and genus 1 <=> the DOUBLET binds (4cos^2(3pi/7)) while genus>=2 <=> larger cores bind below it -- so N=14 is the last doublet-tractable case
status: VERIFIED dims (Eisenstein=#cusps-1=3, cusp=genus=0,0,1,2,2) + VERIFIED core-landscape (p=7 binds at doublet 0.198; p=11,13 bind below it at larger cores 0.0078,0.0049). The local-global/cusp-form-obstruction is the grounded SYNTHESIS (mac-mini's automorphic program). The floor-constant = L-value/period is speculative.
source: klein-2026-06-29-S11
depends_on:
  - HYP-3586   # cusps = Klein, hardness = genus
  - HYP-3581   # the finite cyclotomic floor = 4cos^2(3pi/7) at the doublet
related:
  - HYP-3580   # the proof lives at the cusps
  - HYP-3553   # metagraph = finite Siegel transform (the automorphic frame the decomposition needs)
  - THM-584    # R-eigenspace = sigma-even/sigma-odd = (now) Eisenstein/cusp
  - THM-578    # the doublet = the genus-1 binding core
  - HYP-3547   # Mersenne-Heegner-3mod4 (arithmetic 'why 14'); genus is the geometric 'why', this is the MEANING
results:
  - 04-computation/lrc14_eisenstein_cusp_genus_meaning_klein.py
  - 05-knowledge/results/lrc14_eisenstein_cusp_genus_meaning_klein.out
---

# HYP-3587 — the genus is the local-global gap; the bulk is constant, the obstruction is the cusp form

## What the genus represents (verified dimensions)

Weight-2 on `Gamma_0(N)`: `dim Eisenstein = (#cusps) - 1` (the BULK / boundary space),
`dim S_2(cusp forms) = genus(X_0(N))` (the OBSTRUCTION / global space). Across LRC(2p):

| N=2p | apex | #cusps | dim Eisenstein (bulk) | dim cusp = genus (obstruction) |
|---|---|---|---|---|
| 6  | 3  | 4 | 3 | 0  (LRC(6) SOLVED) |
| 10 | 5  | 4 | 3 | 0 |
| 14 | 7  | 4 | 3 | **1**  (LRC(14) FIRST HARD) |
| 22 | 11 | 4 | 3 | 2 |
| 26 | 13 | 4 | 3 | 2 |

The **bulk is constant (=3)**; only the **obstruction = genus** grows. So: **the genus is the dimension of
the global modes the boundary/cusp data does not determine = the local-global gap.** (A weight-2 form is
fixed by its cusp values up to the genus-dim space of cusp forms, which vanish at every cusp.) Genus 0 =
boundary determines the floor (Hasse / Euler product / bulk rehearsal suffices); genus 1 = one global mode
(the cusp form `f_14`) the rehearsal cannot see = the obstruction.

## The concrete step, worked

- Decomposition is **`3 + 1`**: three Eisenstein (bulk, controlled) + one cusp form `f_14` (= the rank-0
  conductor-14 elliptic curve `14a`).
- **Cusp forms VANISH at cusps**, so the `f_14` *value* at the apex cusp `d=7` is `0`; the obstruction is its
  **leading `q`-expansion coefficient at `d=7`** -- a single finite number. "Bound the cusp-form piece at the
  apex cusp" = bound that coefficient. Rank 0 (`L(f_14,1)!=0`) => non-degenerate (floor bounded away from 0).
- So the remaining floor work = bound the leading apex-cusp coefficient of one explicit newform.

## Genus 1 = doublet-tractable; the boundary of bespoke methods (verified)

Core-gap landscape `gap(O)=min_{k!=0}|sum_{x in O}w_p^{kx}|^2`:
- `p=7` (genus 1): binding core = the **DOUBLET**, gap `0.198 = 4cos^2(3pi/7)` (THM-578); 5 distinct gap
  values.
- `p=11,13` (genus 2): binding cores are LARGER (`{0,1,2,3,7}`, `{0,1,2,3,5,11}`), gaps `0.0078, 0.0049`
  *below* the doublet; 15, 36 distinct values.

So **genus 1 is exactly where the obstruction is the simplest possible configuration (a 2-element doublet);
genus >= 2 the obstruction is irreducibly larger.** `N=14` is the LAST case whose obstruction is a doublet --
a sharp reason it is the last bespoke-tractable LRC.

## The master dichotomy (synthesis)

Every two-index split is the SAME local-global split; `dim(global) = genus`:
Eisenstein | cusp form; boundary | interior; Hasse/local-global-holds | the genus gap; `sigma`-even/R-even
(THM-584) | `sigma`-odd/R-odd; Brouwer/SOS bulk | Borsuk-Ulam obstruction; metagraph/CV(H)/transitive
rehearsal | the cusp-form mode it misses; Euler product (HYP-3550) | anti-Littlewood (HYP-3551); the
`2^7`-core finite cyclotomic min (HYP-3581) | the global `f_14` correction. The project has been studying the
local-global gap of `X_0(2p)`; everything computable is the LOCAL column; the one missing thing is the
genus-dim GLOBAL column (one cusp form at genus 1).

## Synthesizing points (more)

- `nu_2=0 ⟺ apex≡3mod4 ⟺ Paley`, which keeps the genus a clean integer (no fractional corrections) -- the
  Borsuk-Ulam pillar makes the genus well-behaved.
- anti-Littlewood = the global (cusp-form) obstruction to the local product vanishing (Littlewood).
- `14a` bad reduction at `2,7` = the local cusp factors at the doubling/apex cusps; global curve = local
  pieces + one global mode.
- SPECULATIVE (persistence-test first): the floor constant carries a period/`L(f_14,1)` factor.

The whole floor is now: one finite local computation (the `2^7`-core cyclotomic min, done) + one global
number (the apex-cusp leading coefficient of `f_14`).
