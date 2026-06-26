# LRC14 Pair-Good Decoy Generators

The decoys from the binding-pair switch carrier are not a new wild family.  They
are modular tooth events.

At a pair-crossing time `t=p/q`, runner `c` blocks exactly when

```text
14 * min(c*p mod q, q - (c*p mod q)) < q.
```

That is the generator.  A decoy is:

```text
good source pair + blocker residue tooth
```

The finite audit makes the qualitative point sharp.  Across the named rows
there are `6043` decoy times, but `5643` are core-only blockers and `5472` are
ordinary deep-tooth hits.  Large-tail rows look noisy because many pair
crossings hit the same AP-core tooth cover, not because they create thousands of
separate proof obligations.

The proof move is therefore to keep generator classes, not counts:

```text
source lane, source shell, blocker role, tooth bucket, zero-tooth flag
```

This meshes with HYP-3018's normal fan: decoy blocker teeth are the negative
side analogue of active bottleneck owners for positive safe bars.
