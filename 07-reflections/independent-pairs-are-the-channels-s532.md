---
source: oracle-2026-06-01-S532
status: result (n=4 claim verified) + the multi-channel framework (independent pairs = floor(n/2))
tags: [LRC, independent-pairs, matching, orthogonal-roots, channels, CRT, parity-law, n14, multi-channel]
---

# Independent pairs are the channels: the multi-channel metric is the matching number floor(n/2)

**Prompt (user):** for the multi-channel generalization, the metric is probably the
amount and state of INDEPENDENT PAIRS. For the 4-tournament (2 independent pairs)
the iso class is determined by flipping just 2 arcs (the matching arcs) with the
rest fixed.

This is right, and it unifies three things that were floating separately: the
covering-character channels (S531), opus-S524's CRT classes, and the source-sink /
orthogonal-root geometry.

## The user's n=4 claim — verified

An **independent pair** = two disjoint arcs = **orthogonal roots** in `A_{n-1}` (the
90-degree "independent" arcs of the polygon-simplex trinity). A maximal set of them
is a **perfect matching**, size `floor(n/2)`. For `n=4` that is `2` pairs.

Take the matching `M = {(0,1),(2,3)}`; the other four arcs (the `K_{2,2}` cross
frame) are fixed. Computed (`independent_pairs_channels_s532.py`): over the 16 frame
orientations, the number of distinct iso classes reached by the 4 settings of the 2
matching arcs is distributed `{1:4, 2:4, 4:8}` — and **8 of the 16 frames give a
bijection: the 2 matching bits hit all 4 iso classes of `A000568(4)`.** So with a
suitable fixed frame the iso class of a 4-tournament IS exactly the state of its 2
independent pairs — two bits, one per pair. Confirmed.

## Independent pairs = floor(n/2) = the channels

```
 n           3   4   5   6   7   8   ...  14
 indep pairs 1   2   2   3   3   4   ...   7     = floor(n/2)  (the matching number)
```

The punchline: **`floor(14/2) = 7` independent pairs = opus-S524's 7 CRT classes**
(the mod-7 split of 14 = `6 pairs {i, i+7} + the singleton {7}`). The CRT channels
were never about 7 being prime — they are the 7 independent pairs of a perfect
matching of the 14 vertices. The singleton `{7}` is the special pair: the speed
`n/2`, the runner whose orbit is a 2-cycle (the diameter / antipodal direction,
S530). So:

> the multi-channel structure of LRC(n) has **`floor(n/2)` channels = the independent
> pairs of the observer-source tournament**, one of which (speed `n/2`, even `n`) is
> the degenerate diameter pair.

## How the parity law (n=4) is the one-channel degeneration

S531: the n=4 inside debt vanishes iff `a+b+c` is odd. Decompose over the 2 pairs —
pair the observer `0` with one runner, the other two runners together:

```
a + b + c = a + (b+c) = parity(pair_1) XOR parity(pair_2).
```

So "sum odd" is exactly the **XOR of the two independent-pair parities** — a single
combined channel. n=4 has 2 pairs but they fuse into one Z/2 obstruction, which is
why one bit settles half the cases. This is the user's "amount AND state": amount =
2 pairs, state = their joint parity.

## Why n>=6 is genuinely multi-channel (verified)

For `n=6` (3 pairs) the covering character is supported on `k` with `3 ∤ k`, spread
over two residues mod 3 (and the full weight is a character mod 6). Probe: of 120
primitive 5-runner sets, **every single one admits a full-support resonance** — i.e.
the inside debt is **always active**; there is no `a+b+c`-style single congruence
that switches it off. The vanishing condition is now a **joint state of the 3
independent pairs**, not one bit. That is exactly why n=4 was the last clean case
(S531) restated in the user's language: at n=4 the pairs fuse to one channel; from
n=6 on they stay independent and you must control their joint configuration.

## The multi-channel generalization (the program)

LRC(n) = the `floor(n/2)` independent pairs never reach a joint state that covers
the time-circle (no lonely gap). Concretely:
- **amount** = `floor(n/2)` channels (the matching);
- **state** = each pair's contribution to the resonance / to the blocking-arc cover;
- n=4: the 2 states XOR to one bit (parity law, half the cases proved);
- n=14: 7 channels (opus's CRT); the "musical chairs" handoff is the rotation of
  which of the 7 pairs is the last blocker — i.e. the **joint state of 7 independent
  pairs is what must leave a gap.**

So the generalization of the parity law is: **find the joint-state condition on the
`floor(n/2)` independent pairs under which the inside debt vanishes (or is bounded).**
For n=6 it is a condition mod 6 over 3 pairs; for n=14 it is opus's 7-way correlation
read as 7 pairs. The next concrete step: compute the n=6 inside debt as a function of
the 3-pair state and find which states make it `>= 0` — the 3-channel analogue of
`a+b+c odd`.

## Anchor
`04-computation/independent_pairs_channels_s532.py` (+ `.out`): n=4 bijection (8/16
frames), `indep pairs = floor(n/2)` table (= opus CRT at n=14), n=6 always-active
inside debt. Builds on S531 (parity law), S530 (diameter pair), opus-S524 (CRT),
the polygon-simplex trinity (orthogonal roots).
