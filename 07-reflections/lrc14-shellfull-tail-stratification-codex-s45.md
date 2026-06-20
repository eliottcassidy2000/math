# LRC14 Shell-Full Tail Stratification - Codex S45

The useful correction is that HYP-2670's `p1/4` far-tail instinct was too
optimistic.  At B30 it looked clean; at B36 it fails by one row:

```text
E'=(0,1,2,4,8,14,26,34), w=38
Delta^+/p1 = 966562/3357319.
```

The failure is not chaotic.  It is low-fold and doubled-odd:

```text
extras = 2*(7,13,17),  w = 2*19,
fold_recip = 1/34.
```

That matters because it says the next proof target should be an exception
ledger, not another scalar guess.  HYP-2671 owns the dyadic new-speed `1/3`
block.  HYP-2672 now adds a second addressed exception type: doubled-odd tail
packets.

The revised ladder is:

```text
B13 finite pocket      -> 2/5 ledger
dyadic m=4 block       -> 1/3 ledger
intermediate 21..24    -> finite >1/4 ledger
doubled-odd tail       -> exception ledger
remaining post-dyadic  -> <3/10 decay target
```

A focused doubled-odd scan supports this split.  Among `2912` exact rows of
the form `{0,1,2,4,8,2a,2b,2c}`, odd `a<b<c<=29`, `w=2d`, the B36 row is the
unique one above `1/4`, and none exceed `3/10`.

So the geometric lesson is: the shell-full tail is not controlled by fold mass
alone and not by a clean `1/4` scalar.  It is controlled by addressed packet
types.  The doubled-odd packet feels like a tail analogue of the dyadic block:
a sparse, low-fold, arithmetically aligned phase packet that survives one more
scalar threshold.

LRC(14) is not proved, but the post-gate constant now has fewer nameless
places to hide.
