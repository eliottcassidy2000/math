---
id: HYP-2347
name: lrc-divisor-tower-and-56-base
status: PARTIAL (THM-421 proves the single peeling step; the tower + window-fit is open)
date: 2026-06-06
session: claudebox-2026-06-06-S710
depends_on:
  - THM-421   # divisor-clock peeling (proved here)
  - THM-420   # prime-n shell companion
  - HYP-2346  # the mod-7 fiber window-fit (S643)
  - THM-381   # LRC <-> tournament reachability (arc-width), source of |Arc(n)|
  - THM-305   # v_2(A000568(n)) = (n-1)/2, the QR connection
provisional_id: true
---

# HYP-2347: The divisor-clock tower, the prime/composite dichotomy, and the 56-cell base

Extends the "poke" fiber-bundle suggestion (14 = 2.7, project onto divisors) past S643's single
mod-7 step into a full lattice picture, and ties the suggestion's "the half-turn leak misses only 56
cells" to the repo's recurring 56 = A000568(6).

## H1 (tower reduction). 
For composite `n` with prime factorization `n = p_1^{a_1}...p_k^{a_k}`, iterating THM-421 down the
**divisor lattice** of `n` reduces `C'(n)` to a tower of residual loneliness problems, one per
fiber `mod p_i`, the deepest being the radical/squarefree core. For `n = 14`: `14 -> 7` (peel mult-of-7,
mod-2 fiber) is depth 1, residual a `<= 4`-runner sub-config (verified dist {1:348,2:675,3:391,4:86}).
*Open:* prove the tower terminates in a uniformly dodgeable core.

## H2 (window-fit, the shared open core).
At the best clock the residual sub-config `{w_i = v_i/(n/p)}` is lonely within the perturbation window
`|s| < s_0`. Dichotomy (same as HYP-2346): **V' small** => the window is wide enough to host the
sub-config's own lonely time; **V' large** => the dominant runner forces Criterion B'. Proving this for
all composite `n` + proving the THM-420 transversal residual for prime `n` would close `C'` on the
whole frontier `n = 15,19,21,22,...`.

## H3 (the prime/composite partition).
`C'(n)` splits by primality of `n`: composite `n` via the divisor-clock tower (THM-421), prime `n` via
the `2n-1` multiplier shell (THM-420). The hard `n` are exactly those where BOTH handles are awkward:
`n` composite *and* `2n-1` a prime power that ramifies. `n = 14` (`2n-1 = 27 = 3^3`) is the first such;
the next squarefree-`n` double-trouble cases are worth enumerating (e.g. `n` with `2n-1 = q^e`).

## H4 (the 56-cell base = the combinatorial type-space of LRC(7)).
The poke suggestion referred to "the half-turn leak that misses only 56 cells." VERIFIED:
`A000568(6) = 56` = number of tournaments on 6 vertices = the **type-space of the LRC(7) base**
(7 runners = 6 movers; each combinatorial type = an order pattern of the floor-vector `floor(v t)`,
i.e. a tournament on the 6 movers via THM-381). The fiber bundle `LRC(14) -> LRC(7)` has this
56-cell base; the half-turn (mod-2) fiber is the **parity coordinate over each of the 56 base cells**.
The "leak misses 56 cells" reads as: the cover leaves a section over the *entire* base (S643 verified
the cover always leaks => loose), so the leak rides all 56 base types, not a proper subset. The
genuinely-hard residual is the arc-confined sub-family `|Arc(7)| = 6` of the 56 (THM-381 / S514).
QR fingerprint (THM-305): `v_2(56) = 3 = (7-1)/2 = |QR_7|` — the base count carries the 2-adic seam
that controls even-vs-odd `n` throughout the project.

## Why this matters
THM-421 makes H4 more than an analogy: the divisor projection is literally the bundle map, the mod-p
fiber is literally the residue coordinate, and the base type-space is literally A000568(n/p ... ). The
"everything is one object" line (coloring-partition-unification) gains a concrete LRC instance: the
divisor lattice of `n` IS the cover poset of the loneliness problem, glued by CRT.

## Tests to run next
- H1: implement the recursive tower for `n` with `>=2` distinct prime factors (e.g. `n=30=2.3.5`);
  measure residual sizes at each level.
- H3: enumerate `n <= 60` with `2n-1` a prime power (the THM-420 ramified set) and cross with composite
  `n` (the double-trouble frontier).
- H4: for the LRC(7) base, count the floor-vector cells exactly and confirm the 6 arc-confined of 56
  are precisely the residual-hard types under the half-turn cover.
