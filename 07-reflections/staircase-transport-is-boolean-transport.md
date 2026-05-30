# Staircase Transport Is Boolean Transport

Session: codex-2026-05-30

The useful surprise in THM-353 is that the staircase geometry does not need a
new transport theorem.

The only geometric work is to certify that the coordinates are finite:

```text
StTile n ~= {(hi,lo) : Fin n x Fin n // lo + 2 <= hi}.
```

After that, `StTiling n` is just the Boolean cube on those coordinates.  Every
mask-line is an xor line.  Every nonzero mask is fixed-point-free.  The entire
transport checksum falls out of THM-352.

What remains geometric is not the conservation law but the labels:

- single-tile masks mean local interval insertions;
- the complement mask means the global all-bit flip;
- all nonzero masks mean the complete Boolean line geometry;
- `goodCutBucket` means the interval-union height;
- the top `goodCutBucket` is strong connectivity.

That is a good separation.  The row checksum is bookkeeping; the staircase
model supplies semantics.  This is exactly the kind of boundary that should
make later computations less fragile.  If a quotient matrix fails the row
identity, the bug is not in the interpretation.  If it passes, the interesting
questions can move to where the off-diagonal mass goes.

The good-cut quotient also gets cleaner.  Instead of treating `goodCutCount`
as an unbounded natural number and remembering external bounds, Lean now uses

```text
goodCutBucket : StTiling n -> Fin (n+1).
```

The bottom bucket is all-down.  The top bucket is SC.  Between them lies the
interval gas from earlier sessions, but now its transport rows are certified
by the same finite-cube law as every other quotient.

References: THM-351, THM-352, THM-353,
`04-computation/lean/TournamentH7/TournamentH7/StaircaseBucketTransport.lean`.
