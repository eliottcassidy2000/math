# n=17 back to n=14: the skipped half-gate becomes predecessor-of-apex

codex-2026-06-02 S560.

The prompt was to return from the `n=17` attempt to `n=14` and ask whether any
concepts carried over in surprising ways.

They did.  The `n=17` concept that carried was not "prime" but "gate packet."
At `n=17`, the closest structured rows had the form:

```text
{r} union {17*q : 1 <= q <= 16, q != 8}.
```

That is a literal half-gate skip: `8=(17-1)/2`.

## The n=14 transfer

At `n=14`, the same scan over whole gate packets gives:

```text
n=14 scale 7  skip 6: gap/th = 5/924
n=14 scale 14 skip 6: gap/th = 5/1848
n=14 scale 28 skip 6: gap/th = 5/3696
```

The skipped label is not the literal half/apex `7`.  It is `6`, the predecessor
of the apex.  The old seven-ladder was already telling us this:

```text
(1) union {7*q : 1 <= q <= 13, q != 6}.
```

The S380 gate ladder repeats it:

```text
(1) union {14*q : 1 <= q <= 13, q != 6}.
```

The double-gate ladder repeats it again with scale `28`.

So the n=17 skip-half idea carries over with a composite correction:

```text
prime row:     skip the half-gate.
2*q row:       keep the apex q, skip its predecessor.
```

## Why this matters

The `q=7` speed in `n=14` is not an ordinary middle gate.  It is the apex:
the co-observer, the mod-7 singleton, the even-fold bridge, and now the natural
sum-multiple shield in the HYP-2061/THM-396 language.  Removing the apex is not
how the best near-cover behaves.  The near-cover keeps the shield and lets the
leak move to the predecessor gate.

That reframes the seven-ladder.  It is not just "many multiples of 7 plus a
breaker."  It is:

```text
retained apex shield + missing predecessor corridor + exported endpoint debt.
```

The scale chain `7 -> 14 -> 28` then halves the visible gap each time while
increasing boundary/debt:

```text
gap/th: 5/924 -> 5/1848 -> 5/3696
boundary count: 84 -> 168 -> 336
```

This matches the older debt-export story, but the n=17 session supplied a much
cleaner name for the missing packet.

## Assumption challenge

I did not force Tournament Analysis vertices to be runners.  I considered:
runners, gate multiples, skipped labels, fixed walls, endpoint leaves,
small-pinch pair-cells, and whole proof obligations.  The quotient used whole
gate-packet rows because the carried predicate is packet-level:

```text
how close does this gate repair come to an open cover?
```

That quotient preserves gap/debt ordering and skipped label, but destroys
endpoint ownership and small-pinch shield incidence.  It is useful for finding
the invariant, not for closing the proof by itself.

## Next proof move

The next n=14 attempt should combine this with HYP-2061:

1. Treat the retained apex `q=7` as a sum-multiple shield.
2. Treat the skipped predecessor `q=6` as the leak corridor.
3. Prove the scale chain exports endpoint debt faster than it closes the open
   gap.

If that can be made exact, the old seven-ladder near-counterexample stops being
a scary search artifact and becomes a certificate family: every time it shrinks
the visible gap, it pays by doubling labelled boundary debt.
