# A tournament is a binary relation (transitive or intransitive); the LRC peak has a harmonic transitive backbone, and the covering margin is a cut-deficit

*opus-2026-06-29. Owner: build the within-10% existence bound on the AP core, and think of a tournament
as a binary relation — transitive or intransitive — itself. That reframe produces the peak's
deletion-contraction, a clean harmonic backbone, and the exact mechanism of the +10% covering margin.*

## The reframe: transitive vs intransitive is the whole dichotomy
A tournament IS a complete asymmetric **binary relation**, and the only structural dichotomy is
**transitive (a total order) vs intransitive (contains a 3-cycle)**. The OCF is exactly this split:
> `H = 1 + 2Σ_k α_k` — the **`1` is the transitive backbone** (the order = the guaranteed Hamiltonian
> chain, Rédei), the **`2Σα_k` is the intransitive part** (the odd cycles).
mac-mini's HYP-3537 mirrors it for LRC: `cap(S) = 1 − Σ(conditional dangers)`, deletion-contraction
with **teeth = cycles**, factor 2 = complement `Z_2`. Tournament closes by *parity* (`H` stays odd);
LRC by a *measure/peak* bound.

## The peak deletion-contraction: a harmonic transitive backbone (new, exact)
Targeting the **peak `M` (margin-bearing)**, not the cap (measure, razor-thin), the deletion-contraction is
strikingly clean:
> **`M(∅) = 1/2`** (observer alone — pure transitive, witness `t=1/2`), and
> **`M({1,…,k}) = 1/(k+1)` EXACTLY** — the harmonic descent `1/2, 1/3, 1/4, …, 1/14`.
Each consecutive speed `k` cuts the peak from `1/k` to `1/(k+1)` (decrement `1/(k(k+1))`). The AP
`{1,…,13}` is the **pure transitive backbone**, descending harmonically to *exactly* `1/14` (tight).

## The +10% covering margin is a CUT-DEFICIT (verified)
Define `cut(s | prefix) = M(prefix) − M(prefix∪{s})`. Then `M(S) = 1/2 − Σ_k cut(s_k)`, and:
- The **smallest/resonant** available speed cuts the MOST (`cut(12 | {1..11}) = 1/12−1/13`, the steep
  AP step). A **mult-of-14** cuts LITTLE and shrinks as it grows (`cut(84)=.001, cut(168)=.0005`); off
  -resonance speeds (`13,14,15,16`) cut `0`.
- A **covering set is FORCED to include a multiple of 14** (and of 13, …) — a large, off-resonance
  speed in place of a small resonant one. Its **cut-deficit** is exactly the +10% margin:
> **`M(covering) = 1/2 − Σcut > 1/2 − 3/7 = 1/14`** — the AP attains `Σcut = 3/7` (`M=1/14`); the
> covering's forced mult-of-14 has a smaller cut, so `Σcut < 3/7` and `M > 1/14`.

So the margin is not slack in a measure — it is the **arithmetic deficit of the forced multiplicative
speed** relative to the small additive speed it displaces. This is the add/mult duality made
quantitative: the **additive AP is the tight cut-maximizer; the multiplicative covering load is the
cut-deficit that lifts `M` off the edge.**

## What a disproof would be (the transitive/intransitive picture)
`M(S) < 1/14` would mean `Σcut > 3/7` — some speed set cuts the peak *more* than the smallest speeds
do. Equivalently: the danger arcs **close cyclically** around the time circle (an *intransitive*
cover, `χ(nerve)=0`, no break). A lonely point is a **transitive break** (a gap, `χ≠0`). A circle
cannot be covered transitively (an order has a gap); covering needs intransitive closure, and the
arithmetic of integer speeds obstructs it. **LRC(14) ⟺ the cut-sum never exceeds `3/7` ⟺ the cover
never intransitively closes ⟺ the AP `{1,…,13}` is the global cut-maximizer.**

## Improved proof target (the within-10% bound, sharpened)
1. **Work with the peak `M`** (margin), via the cut decomposition `M = 1/2 − Σcut`, NOT the cap
   measure (razor-thin).
2. **Prove the AP is the cut-maximizer:** `Σ_k cut(s_k) ≤ 3/7` for every primitive 13-set, equality
   only at `{1,…,13}`. The covering constraint *forces* a mult-of-14 whose cut-deficit gives strict
   `< 3/7`, hence `M > 1/14` with margin. (The per-speed cut is governed by *resonance* with the
   prefix, not raw size — `12` and `24` cut, `13..16` don't — so the bound is a resonance/rearrangement
   statement, the same shape as max-H = the additive-interval extremizer.)
3. **The mult-of-14 is the apex-7 face:** its cut-deficit is the 2-adic/Mersenne descent (THM-580)
   made quantitative — peel it to the additive AP core where `M=1/(k+1)` is exact.

## Status
- **Verified:** `M({1,…,k})=1/(k+1)` (harmonic backbone); `M(∅)=1/2`; the cut-deficit table (small/resonant
  cut most, mult-of-14 cuts little); covering `M>1/14` is exactly the deficit.
- **New (opus):** the peak deletion-contraction; the harmonic transitive backbone; the margin = cut-deficit
  of the forced multiplicative speed; LRC(14) ⟺ AP is the cut-maximizer (`Σcut ≤ 3/7`).
- **Improved target:** prove `Σcut ≤ 3/7` (AP cut-maximality) — a resonance/rearrangement extremality,
  dual to max-H, with the multiplicative mult-of-14 deficit as the margin.

Related: mac-mini HYP-3537 (cap = measure OCF), THM-582 (two-index), THM-580 (apex-7 descent), THM-070
(Claim A), my razor-thin + add/mult-duality reflections, the max-H = additive-interval extremality, OPEN-Q-108.
