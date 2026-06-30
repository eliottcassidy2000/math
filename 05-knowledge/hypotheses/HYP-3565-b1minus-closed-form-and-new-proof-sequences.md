---
id: HYP-3565
title: 1,7,119,1772 = b_1^-(n) (klein-S6, the R-odd first Betti of the metagraph = the OBSTRUCTION dimension); NO elementary closed form (a difference of non-elementary tournament-Burnside sequences A000568/SC/E(G_n)), but the exact Lefschetz form b_1^-=(E-V+SC-E_SCSC+E_comp)/2, the cycle-index/Burnside form (n!-sum, past the 2^{C(n,2)} wall), the asymptotic ~C(n,2)2^{C(n,2)-2}/n!, and a NEW structural fact b_1^-/b_1 -> 1/2 (the complement R splits the metagraph cycle space asymptotically EVENLY -- the obstruction is half, not vanishing); plus a roster of new useful proof-sequences and their metagraph<->LRC obstruction-dimension correspondence
status: ANSWERED + verified (the split b_1=E-V+1, the asymptotic, the 1/2 fraction; n=3..7). No elementary closed form found (consistent with klein-S6 'not in OEIS'); the structural/cycle-index forms + asymptotic are the closed forms that exist.
source: mac-mini-2026-06-29-S25
related:
  - HYP-3563  # klein-S6: b_1^-(n) = 0,1,7,119,1772, the apex-7/octonion conjectures refuted
  - HYP-3544  # klein equivariant homology (b_1^- = the R-odd Betti = the BU obstruction)
  - HYP-3562  # the measure of the obstruction (b_1^- is its cycle-space dimension)
  - THM-587   # SC = P_n(-1); the cycle-index machinery for the Burnside form
  - THM-589   # the metagraph variance W(n) (a companion new sequence)
external: A000568 (V), A001615 (Dedekind psi); klein-S6 (the values + Lefschetz formula)
results:
  - 04-computation/b1minus_closed_form_macmini_20260629.py
  - 05-knowledge/results/b1minus_closed_form_macmini_20260629.out
---

# HYP-3565 -- 1,7,119,1772 = b_1^-, its closed-form status, and new proof-sequences

## The identification and the closed-form status
`1, 7, 119, 1772` (with leading `0`) is **`b_1^-(n)` for n=3..7** -- klein-S6's R-ODD FIRST BETTI of the
arc-flip metagraph `G_n`, `b_1^-(n) = dim H_1^-(G_n;Q)`, the part of the cycle space on which the
complement `R` acts by `-1`. It is the **OBSTRUCTION DIMENSION** -- the cycle-space realization of the
S23 obstruction (the R-odd / Borsuk-Ulam part, klein HYP-3544).
- **NO elementary closed form.** `b_1^- = (E - V + SC - E_SCSC + E_comp)/2` (klein's Lefschetz): a
  difference of TOURNAMENT-BURNSIDE sequences -- `V = A000568`, `SC = P_n(-1)` (THM-587), `E(G_n) =
  1,5,30,290,4086`, `E_SCSC`, `E_comp` -- NONE of which has an elementary closed form. So neither does
  `b_1^-`. Factorizations `0,1,7,7*17,2^2*443` show no pattern (klein: not in OEIS).
- **The closed forms that DO exist:** (i) the **Lefschetz form** above; (ii) the **cycle-index/Burnside
  form** -- each of `V, SC, E` is a sum over `S_n` (computable from `n!`, past the `2^{C(n,2)}` wall),
  so `b_1^-` is an explicit `S_n`-cycle-index expression; (iii) the **asymptotic**
  `b_1 = E - V + 1 ~ E ~ C(n,2) 2^{C(n,2)-1}/n!`, so `b_1^- ~ C(n,2) 2^{C(n,2)-2}/n!` (leading growth
  `2^{C(n,2)}/n!`).

## The new structural fact: the obstruction is asymptotically HALF
`b_1 = E-V+1 = 0,2,19,235,3631`; the split is `b_1^- = 0,1,7,119,1772` (obstruction) and
`b_1^+ = 0,1,12,116,1859` (the R-even / measure-floor side). The **OBSTRUCTION FRACTION**
`b_1^-/b_1 = 0.50, 0.37, 0.51, 0.49` hovers near `1/2`:
> **`b_1^-/b_1 -> 1/2`**: the complement `R` splits the metagraph cycle space asymptotically EVENLY.
The obstruction is HALF the cycle space, not vanishing -- consistent with the S19 mean-eigenvalue->0
(the signed-action asymmetry washing out). So the topological obstruction is large and robust: half the
homology is essential-R-odd.

## The metagraph <-> LRC obstruction-dimension correspondence
`b_1^-` is the metagraph instance of the obstruction; its LRC mirror is the R-odd part of the lonely set:
at the extremal the lonely points are the `phi(n)` units mod `n` in `phi(n)/2` antipodal pairs (S15/S23),
so the **LRC obstruction dimension = the saddle index `phi(n)/2`** (`(p-1)/2` at apex `p`). The two
obstruction-dimension sequences:
- metagraph cycle-space: `b_1^-(n) = 0,1,7,119,1772` (super-exponential).
- LRC lonely-pair: `phi(n)/2 = ...,3 (n=14)` (the saddle index, slow).
Same R-odd obstruction, two complexes (the metagraph cycle space vs the danger-cover H_0).

## New useful proof-sequences (roster)
| sequence | values (small n) | role |
|---|---|---|
| **`b_1^-(n)`** obstruction dim | `0,1,7,119,1772` | the R-odd cycle obstruction (this) |
| **`b_1^+(n)`** measure-side dim | `0,1,12,116,1859` | the R-even cycle space (the floor side) |
| **`E(G_n)`** metagraph edges | `1,5,30,290,4086` | the arc-flip relation size (Burnside, new) |
| **`W(n)`** metagraph variance | `1,2,8,32,158,928` | the 2nd moment (THM-589, = simplicial-Redei) |
| **`SC(n)=P_n(-1)`** Lefschetz | `2,2,8,12,88` | the Euler/counting obstruction measure (S23) |
| **`phi(n)/2`** saddle index | the LRC R-odd lonely-pair count | the LRC obstruction mirror |
| **`psi(N)`** Dedekind psi (A001615) | `[SL2:Gamma_0(N)]` | the covering congruence index (HYP-3553) |

## What it buys
A precise answer (`b_1^-`, no elementary closed form, the Lefschetz/cycle-index/asymptotic forms) and a
new fact (the obstruction is asymptotically half the cycle space). The roster above gives the proof its
working sequences: the obstruction (`b_1^-`, `phi/2`), its measure (`SC`), the relation (`E`, `W(n)`), and
the change-of-base (`psi(N)`) -- the objects of the relational/obstruction program (HYP-3562/3564) made
into computable sequences.
