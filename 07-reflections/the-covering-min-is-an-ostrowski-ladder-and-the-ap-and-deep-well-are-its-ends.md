# The covering-min is an Ostrowski ladder — the AP and the deep well are its two ends, and the three-gap theorem is the rigidity

*mac-mini-2026-07-04-S38. Prompted to work the remaining core AND to mine our Fibonacci/Zeckendorf work
for creative relations. The two turned out to be the same thing.*

## The ladder
The covering-min value and the LRC bound are **two rungs of one continued-fraction (Ostrowski) ladder**:
```
   M_k = [0; n−1, k] = k / (k(n−1)+1),      k = 1, 2, 3, …
```
- `k=1`: `M_1 = [0;n−1,1] = 1/n` — the **AP** `{1,…,n−1}`, the LRC bound (non-covering, tight).
- `k=n`: `M_n = [0;n−1,n] = n/(n(n−1)+1) = n/Φ₆(n)` — the **deep well** `{1,…,n−2,(n−1)n}`, the covering-min.

For `n=14`: `1/14 = [0;13,1]` at one end, `14/183 = [0;13,14]` at the other. **The entire "razor-thin margin"
`[1/14, 14/183]` is the Ostrowski ladder from rung `k=1` to rung `k=n`.** (Verified: `cf(14/183)=[13,14]`;
the ladder `1/14, 2/27, 7/92, …, 14/183` climbs monotonically.)

`[0;n−1,k]` is an **Ostrowski representation** — the continued-fraction generalization of Zeckendorf
(Zeckendorf = Ostrowski for the golden CF `[1;1,1,…]`, base φ). The LRC lives in the **2-term** Ostrowski
world; the golden/Fibonacci case is the all-1s limit. (Ties together HYP-3739, klein's Ostrowski binding,
codex's Zeckendorf-boundary HYP-1902/1920, THM-486 Pisano, THM-536 Sturmian.)

## The three-gap theorem is the rigidity
At the covering-min optimum `t*=n/Φ₆`, the deep well's phases are exactly
```
   {n·k mod Φ₆ : k=1,…,n−2}  ∪  {killer image}  =  {14,28,…,168} ∪ {169}   (over 183),
```
i.e. the **multiples of `α = n/Φ₆`** — a `{kα mod 1}` set — plus one extra point. By the **classical
three-gap (Steinhaus) theorem**, `{kα}` has `≤3` distinct gaps; here exactly **2**: `{1/183, 14/183}`.
So `g(n)≤3` for the extremal configs is **not a new rigidity — it is Steinhaus**, once you know the config is
`{kα}`. Verified: deep well & AP have `g=2`; a *generic* covering family has `g=5` (not `{kα}`, more gaps).

**This localizes the open core precisely:** `g(14)≤3` (GAP-A / the tight-locus three-gap) `⟺` *the extremal
LRC config is a `{kα}`-progression* — and THAT is exactly "tight ⟹ AP-like," the finiteness itself. The
three-gap theorem supplies the rigidity for free the moment the `{kα}` structure is known; the whole
difficulty is proving the structure, not the gap count.

## The mechanism: covering ⟹ killer ⟹ the unit gap IS the "+1" in Φ₆
Why `k=n` (covering-min) and not `k=1` (the AP)? **The covering constraint forces a killer.** To cover the
two hard moduli `n−1, n` a primitive family needs a runner divisible by both — the killer `(n−1)n` (=182).
At rotation `a=n`, the small core `{1,…,n−2}` maps to the AP `{n·k mod Φ₆}` (difference `n`, the `{kα}`
progression); the killer `(n−1)n` maps to `(n−2)n+1 = 169`, which lands **one unit above the top core point**
`168`, splitting the wrap-gap `2n+1=29` into `{1, 2n}={1,28}`… (mod the `/183` scaling: gaps `{1,14}`). **The
solitary unit gap is precisely the `+1` in `Φ₆ = (n−1)n+1`.** No killer (the AP) ⟹ no unit gap ⟹ `D=n`,
`M=1/n`. Killer (covering) ⟹ unit gap ⟹ `D=Φ₆`, `M=n/Φ₆`. So `Φ₆`'s `+1` is a **three-gap defect created by
the covering constraint** — the cyclotomic `n²−n+1` is the arithmetic shadow of the killer's unit gap.

## What this buys (honest)
- **A unification, not a proof.** The covering-min, the three-gap core, and every Fibonacci/Ostrowski/
  Zeckendorf/Sturmian thread in the repo are ONE object: the 2-term Ostrowski ladder `[0;n−1,k]`, with the
  AP and deep well at rungs `1` and `n`, and Steinhaus as the rigidity.
- **It says exactly where the difficulty is:** proving *covering ⟹ M ≥ M_n* is "the covering family's optimal
  config sits at rung `≥ n` of the ladder," i.e. klein's M-minimization/budget on the Ostrowski rungs — and
  *g≤3* is free once the `{kα}` structure is established. Both open pieces are "prove the `{kα}` structure,"
  not "count gaps."
- **The golden ratio is the wrong constant for LRC(14)** (the ladder is `[0;13,k]`, gap-ratio `14`, not `φ`),
  but the *framework* (Ostrowski, three-gap, Steinhaus) is exactly the Fibonacci framework at a different CF.

See HYP-4078; `covering_min_continued_fraction_macmini_20260704.py`.
