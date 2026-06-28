# LRC14 Two-Branch Wall-Signature Atlas Reflection

Date: 2026-06-28

HYP-3425 made the right move by turning relocation into a two-color interval
obstruction.  But positivity alone is still too soft for a proof.  The next
thing to name is the boundary of the survivor interval.

The concurrent HYP-3426 mirror audit is the immediate reduction layer: it
removes branch ambiguity and shows endpoint-owner supports stay tiny on the
expanded bank.  HYP-3427 should be read as the next certificate layer on top
of that, not as a competing HYP-3426 claim.

The HYP-3427 scout gives a useful finite certificate language.  In the tight
row `{1..11,13,84}`, the four survivor windows are exactly bounded by `E:84`
and odd branch walls `5` and `7`.  That is much closer to a lemma than the
statement "good measure is `1/105`": it names which walls prevent the bad core
from covering the even-safe set.

The broader audit is encouraging but also a warning.  There are only `27`
global signature types across `5524` windows in `67` rows, but the midpoint
binders still involve even, odd, `7`, and `14Q` roles.  A proof should retain
the wall word and binder roles, not collapse them to one measure or one energy
scalar.

The incoming S302 coordination note makes the same point from the additive
side: energy has to travel with a sheet sidecar.  HYP-3427 therefore ranks
`energy_sheet_sidecar_join` as a repair/join below the exact wall certificate,
not as a replacement for the survivor-window boundary word.

The strongest next route is a bounded wall-alphabet lemma:

```text
survivor window exists with legal left/right walls
```

where illegal or growing wall alphabets become named sidecar debts.  This fits
HYP-3424's transfer rule and HYP-3425's energy-sheet warning: scalar data may
rank packets, but the interval proof needs exact walls.
