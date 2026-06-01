---
source: oracle-2026-06-01-S533
status: result — divisibility verified; clean frame-factorization holds only n<=4; breakdown = the coupling
tags: [tournaments, A000568, independent-pairs, frame, 2-adic, coupling, multi-channel, LRC]
---

# The almost-fixed frame: it works exactly at n<=4, and its failure is the LRC coupling

**Prompt (user):** consider `2^x` tournament iso classes as `x` independent-pair
flips on a fixed frame; understand how an ALMOST fixed frame can satisfy this.

The honest answer has a clean positive half and an illuminating negative half, and
the negative half is the same wall as the LRC difficulty.

## The algebra is real: `2^{floor(n/2)}` divides `A000568(n)`

Verified for n=3..12: `2^{floor(n/2)} | A000568(n)`, with quotients
`1, 1, 3, 7, 57, 430, 11971, ...`. The `floor(n/2)` independent-pair flips
contribute a genuine factor of `2^{floor(n/2)}` to the tournament count. (The
2-adic valuation even *exceeds* `floor(n/2)` for even `n >= 8` — an extra pair-flip
2 at n=8,10,12 — so the pair-flip factor is a lower bound, not the whole 2-part.)

So "x pair flips" really is a `2^x` factor of the iso-class count, for every n.

## The geometry is clean only at n = 3, 4

The quotient `A000568(n)/2^{floor(n/2)}` is the number of "frames." It equals **1**
exactly at n=3 and n=4. There, a *single fixed frame* (the non-matching arcs) plus
the `floor(n/2)` pair-flips realizes **all** iso classes: A000568(3)=2=2^1,
A000568(4)=4=2^2. (n=4 verified S532: with a suitable fixed frame the 2 matching
bits index all 4 classes.) This is the user's picture, and it is exactly true for
the triangle and the square.

## At n >= 5 the frame cannot be (almost) fixed — verified

For n=5, `12 = 3 x 4`, so a clean factorization would tile the 12 iso classes by
**3 disjoint** pair-flip blocks (3 frames, each contributing a free `2^2` block).
Computed (`almost_fixed_frame_s533.py` + partition search):

- 136 of 256 frames give a *full* block of 4 distinct iso classes (the pair-flip
  action is generically free — the local `2^x` structure is there);
- but **no 3 of those 136 blocks are pairwise disjoint** — there is NO partition
  of the 12 iso classes into disjoint pair-flip blocks. The blocks genuinely
  overlap; a minimal cover needs 5 frames (not 3), Hamming-spread 4 of 8 arcs.

So an exactly-fixed frame fails (quotient > 1), and an *almost*-fixed frame fails
too: the frames a cover needs are far apart, and no disjoint tiling exists. The
`2^{floor(n/2)}` factor survives in the *count* but not as a clean geometric
product of "frame x pair-flips."

## Why: the failure IS the inter-pair coupling = the LRC multi-channel wall

The reason the frame won't stay fixed is that flipping one independent pair changes
how the *other* pairs sit relative to the frame — the pairs are **coupled** through
the shared frame arcs. At n=3,4 there is at most one effective pair-degree (the two
n=4 pairs fuse to one Z/2, S532), so nothing to couple — the frame is fixed. From
n=5 the `floor(n/2)` pairs interact, and the frame must move to absorb it.

This is *precisely* the obstruction from the LRC side:
- the n=4 parity law works because the 2 pairs fuse to one channel (S531/S532);
- it fails at n>=6 because the `floor(n/2)` channels stay independent and couple;
- equivalently, the inside debt becomes a JOINT state of the pairs (S532).

> **"The frame cannot stay fixed for n >= 5" and "the LRC parity law fails for
> n >= 6" are the same fact: the independent pairs are coupled through the frame.**
> An almost-fixed frame can satisfy the picture only where the pairs do not couple —
> the triangle and the square. Past that, the residual frame motion is the coupling,
> and bounding it is bounding the LRC inside debt.

## What survives, as a usable handle

The `2^{floor(n/2)}` divisibility is a real invariant: it says the iso-class set
carries a free `(Z/2)^{floor(n/2)}` shadow (the pair-flips) even when no fixed frame
displays it. For LRC this is the right count of channels (S532). The program is not
"find the fixed frame" (it doesn't exist past n=4) but **measure the coupling** —
the deviation of the true iso/loneliness structure from the free `2^{floor(n/2)}`
product. That deviation, channel by channel, is the multi-channel generalization of
the parity law.

## Anchor
`04-computation/almost_fixed_frame_s533.py` (+ `.out`): `2^{floor(n/2)} | A000568`
(n=3..12); n=4 single-frame; n=5 136 full-block frames but NO disjoint 3-partition;
n=6 likewise. Builds on S532 (independent pairs), S531 (parity law), S529
(inside/outside).
