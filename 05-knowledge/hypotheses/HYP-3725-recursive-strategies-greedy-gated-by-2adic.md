---
id: HYP-3725
title: RECURSIVE PROOF-STRATEGY LEDGER for the open obligations -- finding: the densest-core GREEDY = the single-skip covering-min = n/Phi_6(n) for nearly all n, with SPORADIC small exceptions at n=4 (skip-2: 2/7<4/13) and n=8 (skip-6: 4/29<8/57); n=14=2*7 (and n=16) WORK. So the recursive-greedy covering-min is a strong heuristic but NOT universal (sporadic exceptions, NOT cleanly the powers of 2 -- 16 works); n=14 is in the good regime. The existence obligation is a recursive FAREY/semiconvergent climb (q_j=q_{j-1}+(n-1), the Euclidean/Stern-Brocot recursion from the blocked mediant to the convergent); the descent (THM-580) is the per-level recursion
status: VERIFIED (single-skip covering-min vs n/Phi_6(n), n=4..16: equal EXCEPT n=4,8 where skip-(n-2) wins; n=14, n=16 both work). HONEST CORRECTION: the exceptions are NOT all powers of 2 (n=16 works) -- they are sporadic small cases (n=4,8). The greedy is a strong-but-not-universal strategy; the Farey-climb existence recursion is clean; the strategies are STRATEGIES, not completed proofs.
source: klein-2026-06-29-S32
depends_on:
  - HYP-3724   # Phi_6 = Sylvester map; the Egyptian tower
  - HYP-3718   # convergent-not-mediant; the semiconvergent ladder
related:
  - HYP-3715   # M = covering radius; densest-core + killer
  - HYP-3551   # 14/183 the tightest (grid-supported global min for n=14)
  - THM-580    # the 2-adic descent (the per-level recursion)
results:
  - 05-knowledge/results/recursive_strategies_klein.out
  - 05-knowledge/results/greedy_fails_at_powers_of_2_klein.out
---

# HYP-3725 — recursive strategies for the open obligations; the greedy is 2-adically gated

## The open obligations (recap)
(O1) GLOBAL covering-min: no covering beats `n/Phi_6(n)` (`14/183` at `n=14`).  (O2) the 1D-three-gap `<->`
2D-Kershner metric bridge.  (O3) the conditional reduction `rho_j = Z_7`-core gap.  (O4) the `f_14` cusp-form
bound.  (O5) the top-level existence `rho_0 > 0`.

## Strategy for (O1): greedy / Sylvester -- strong but NOT universal (sporadic exceptions)
The densest-core greedy (skip `n-1`, killer `lcm(n-1,n)`) gives `n/Phi_6(n)`. Is it the (single-skip)
covering-min? Verified `n=4..16`:
```
n : 4    5  6  7   8    9 10 11 12 13 14 15 16
greedy=min & =n/Phi_6 ?  NO  Y  Y  Y  NO  Y  Y  Y  Y  Y  Y  Y  Y
```
**FAILS only at `n=4` (skip-2 gives `2/7 < 4/13 = greedy`) and `n=8` (skip-6 gives `4/29 < 8/57`).** For
every other `n` in range -- INCLUDING `n=14 = 2.7` AND `n=16 = 2^4` -- the densest-core greedy IS the
single-skip covering-min `= n/Phi_6(n)`. HONEST CORRECTION: the exceptions are NOT the powers of 2 (`n=16`
works); they are sporadic small cases (`n=4, 8`), where `skip-(n-2)` wins because `lcm(n-2,n) = (n-2)n/2` is
a small killer (`n-2, n` share the factor 2) -- but this advantage does NOT persist (`n=6,10,12,14` even,
greedy still wins). So the recursive-greedy is a **strong heuristic, not a universal recursion**; `n=14` is
in the good regime (greedy = `14/183` among single-skips; HYP-3551's grid search supports it as the full
global min). The lesson: the global-min (O1) is genuinely `n`-dependent and not closed by naive greedy --
the sporadic exceptions are a real obstacle.

## Strategy for (O5) existence: the recursive FAREY/semiconvergent climb
The observer's escape (HYP-3718) is reached by a **linear recursion**: the semiconvergents `p_j/q_j =
j/(j(n-1)+1)` satisfy `q_j = q_{j-1} + (n-1)`, `p_j = p_{j-1} + 1` -- the Stern-Brocot / Euclidean recursion,
each rung the mediant of the previous and the convergent. The climb starts at the **blocked mediant**
`1/n` (`j=1`, killer at `0`) and ascends monotonically (`M_j = p_j/q_j` increasing) to the **convergent**
`n/Phi_6(n)` (`j=n`), the escape. So existence (`rho_0 > 0`) is a recursive ascent up the semiconvergent
ladder -- the Euclidean algorithm IS the existence proof, bottoming at the convergent.

## Strategy for (O2)/(O3): the descent as the per-level recursion
The 2-adic descent (THM-580) `S -> S' = E/2` IS the recursion: it reduces a covering to a smaller one,
bottoming at the apex (the odd core). (O3) the apex projection recurses down the descent; (O2) the covering
radius is computed by the binary peel `g = ceil(x/2)` (HYP-3723). The descent's depth `d <= log2(max speed)`
bounds the recursion. The cusp obstruction (HYP-3599) localizes to the TOP level, where (O5)'s Farey climb
takes over (existence, not measure).

## The honest assessment
- The recursive-greedy (O1) WORKS for nearly all `n` (single-skip min `= n/Phi_6(n)`), with sporadic small
  exceptions `n=4, 8` (NOT the powers of 2 -- `n=16` works); `n=14` is in the good regime. A strong heuristic,
  not a universal recursion or completed proof; the sporadic exceptions show the global-min is genuinely
  `n`-dependent.
- The Farey-climb (O5) is a clean recursive existence structure (the Euclidean algorithm).
- The descent (O2/O3) is the per-level recursion; the cusp localizes the open core to the top, handled by
  the Farey climb.
- (O4) the `f_14` cusp-form bound has its own (Hecke) recursion -- not addressed here.

## Net
Recursive thinking applied to the obligations yields: a greedy/Sylvester strategy for the covering-min
(single-skip min `= n/Phi_6(n)` for nearly all `n`, sporadic exceptions `n=4,8`, `n=14` OK), a clean
recursive Farey/Euclidean climb for the existence escape (blocked mediant `1/n` -> convergent `n/Phi_6(n)`,
`q_j = q_{j-1}+(n-1)`), and the 2-adic descent (THM-580) as the per-level recursion (cusp localizes to the
top, handed to the Farey climb). The Farey-climb and descent recursions are clean; the greedy is strong but
has sporadic exceptions -- the honest signal that the global covering-min (O1) is `n`-dependent and not
closed by naive recursion.
