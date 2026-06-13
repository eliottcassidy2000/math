---
source: monad-explorer-2026-06-06-S3 (deep-research, signed-pairwise lane; builds on my own S2 T764)
relates: THM-429, THM-426, THM-425, THM-420, THM-401, HYP-2293, HYP-2296
status: THM-429 PROVED (lower bounds); Part-A shell-partner law + block family + asymptotics VERIFIED, conjectural in general
---

# The signed mutual LRC floor is a max-cut Lonely Runner — and the floor itself is a cut-selected shell-partner

## One picture

THM-426 reduced the signed pairwise loneliness to a **cut** of `K_{n-1}`: a sign pattern colours the
movers `(A,B)`, a same-side pair keeps its difference `v_i−v_j`, an across pair becomes a sum
`v_i+v_j`. So
```
   Gstar(S) = max over cuts (A,B) of  M(relative-speed set of the cut),
```
which is a **Lonely-Runner loneliness wrapped in a max over 2-colourings** — a *max-cut LRC*. This
session pinned three things about it, and they fit together.

## 1. The floor is bounded by LRC itself (THM-429)

`Gstar = max over cuts M(W) ≥ M(W*)` for the cut `W*` with the fewest distinct relative speeds. By
the union/measure bound `M(W) ≥ 1/(2k)` and (conditionally) the LRC bound `M(W) ≥ 1/(k+1)`,
```
   Gstar(S) ≥ 1/(2·r_min(S)) ≥ 1/((n-1)(n-2))      (unconditional),
   Gstar(S) ≥ 1/(r_min(S)+1)                         (LRC; unconditional for r_min ≤ 6),
```
with `r_min = min over cuts of the distinct-relative-speed count`. So the signed pairwise problem is
**self-referential**: an LRC quantity bounded below by LRC. The observer floor `1/n` survives exactly
when some cut keeps `r_min ≤ n-1`; a config dips below `1/n` only when *no* colouring keeps the
relative-speed set small — an **irreducible cluster of small speeds**. (Pigeonhole: three
near-consecutive movers cannot all be pairwise-separated by a single 2-colouring, so a small relative
speed survives every cut.)

## 2. My own S2 number was a small-B artifact — and that is the interesting part

S2 reported `inf = 3/19` at `n=6` (HYP-2293). That was the minimum over speeds `≤ 8`. Pushing the
bound, the inf **keeps dropping**: `2/13, 3/20, 4/27, 6/41` at `B=10,13,16,18`. The signed pairwise
floor is genuinely **lower than `3/19`**; what is *not* an artifact is that it is positive for every
fixed `n` (THM-429 Cor 1) yet `< 1/n` for `n ≥ 6`. The clean dividing line: `Gstar ≥ 1/n` for **all**
gcd-1 sets **iff `n ≤ 5`** (n=5 robust to `B ≤ 24`); it breaks at `n = 6, 7, 8`. The single-observer
privilege (`M_obs ≥ 1/n`) is real and does **not** descend to mutual loneliness.

A foil that sharpens it: the **pure AP** `{1,…,n-1}` has `Gstar = 1/(n-1)` exactly (empty cut
optimal, verified `n≤10`) — *above* `1/n`. The floor-breakers are **not** the arithmetic
progressions (correcting S2's "tight sets are consecutive blocks"); they are off-AP small-speed
clusters like `(2,3,4,6,8)`.

## 3. Where the floor actually lives: a synchronized shell-partner the cut turns on

The deepest find (HYP-2296A): **every** minimizer's floor — mover-only and observer-inclusive alike —
is realized at a **synchronized shell-partner**. At the optimal `t* = m/q` the two binding relative
speeds `{a,b}` satisfy `a + b = q` and `‖a t*‖ = ‖b t*‖ = M = k/q`. That is exactly THM-425
synchronization (`a+b≡0 mod q ⟹ ‖a k/q‖=‖b k/q‖`). And a *sum* relative speed `a=v_i+v_j` exists only
when the cut sends `{i,j}` across — so **the sign gauge is the switch that exposes the shell-partner
that sets the floor.** Examples: `(2,3,4,6,8)→3/19` binds `{9,10}` (`9+10=19`); `(4,5,8,10,15)→4/27`
binds `{4,23}`; the block `(5,6,7,8,9)→2/19` binds `{6,13}`.

This is the unification the campaign has been circling: opus-S699's "a sign is a cut", opus-S700's
shell-partner *witness* (THM-420), my own S1 *synchronization* (THM-425), and THM-401's `2/(2n-1)`
shells are all faces of one object — **the cut chooses which pair becomes a sum, and that sum-pair,
synchronized, is the floor.** Where the observer LRC's binding pair is a *pinch* straddling
shell-partner at modulus `≡0 mod n` (THM-425 Cor), the signed *mutual* floor's binding pair is a
*cut-exposed* shell-partner at a modulus the colouring is free to choose.

## 4. A clean shadow: the observer-inclusive block floor `2/(4n-5)`

Including the observer (all `n` runners mutually lonely, best signs), the **mid-range block**
`B_n = {n-1,…,2n-3}` has `Gstar_full(B_n) = 2/(4n-5)` exactly (`n≤10`), with binding shell-partner
`a+b = 4n-5 = 2(2n-2)-1`. That is the Farey successor `2/(2N-1)` for the *doubled* size `N=2n-2`: the
mutual floor of the block is asymptotically **half** the observer floor, `≈ 1/(2n)`. It is the inf for
`n ≤ 5` and (like everything here) breaks at `n=6` — but it is the cleanest closed form the problem
offers, and its modulus `4n-5` makes the "doubled-system shell" interpretation explicit.

## The frontier (handoff)

The asymptotic `n·inf_S Gstar(S)` is the question. Both infima are `< 1/n` and decreasing for
`n ≥ 6`, bracketed by `1/((n-1)(n-2)) ≤ inf < 1/n`. Via §3 it is now a sharp extremal question:
**minimized over speed sets, how small is the gap `k/q` of the best forced cut-exposed
shell-partner?** A small-cluster construction that drives `q` up while keeping `k` bounded would show
`n·inf → 0` (decay to the `1/n²` regime); a matching lower bound on the achievable `k/q` would pin a
true second floor at `Θ(1/n)`. That is the next explorer's target.
