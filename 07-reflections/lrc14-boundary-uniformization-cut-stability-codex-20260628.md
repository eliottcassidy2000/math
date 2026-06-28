# LRC14 Boundary-Uniformization Cut Stability

The new request could easily have scattered into named constants and beautiful
maps.  The useful discipline was to force every imported idea to answer the
same question:

```text
What is this quotient allowed to forget, and how does the lost coordinate
return before a theorem exit changes?
```

That keeps the Bring radical, Schwarz-Christoffel maps, BDH, Krasner, and the
Ramanujan-Soldner / Meissel-Mertens constants from becoming decoration.

The strongest synthesis is this:

```text
boundary geometry supplies the channels,
Menger cuts decide which channels are indispensable,
Krasner stability explains local collar invariance,
BDH/Mertens controls average residue-owner leakage,
recursive chirality catches mirror collapse.
```

This lands right on HYP-3405 and HYP-3406.  HYP-3405 has a local boundary
collar with a named unit-height leak.  HYP-3406 has an enlarged-bank residue
leak whose persistent form is owner-side.  HYP-3407 says these are two faces
of the same labelled packet theorem: if residue or height forgets the
load-bearing coordinate, an owner cut or local disk exit must resurrect it.

The most concrete next computation is not analytic yet.  Build the endpoint
owner-support graph for:

```text
petal 13->26
single swap 1->26
single swap 3->26
single swap 5->26
single swaps into 40 and 54
```

Then compute a real min-cut separating `unit-petal-named` from
`positive-Haar-open`.  If that min-cut is stable across the HYP-3406 bank
growth, it becomes much closer to a theorem than another scalar separator.

The paired local task is AP versus `13->27` in the HYP-3405 collar.  The
correct Krasner-style statement is not that p-adics prove the real LRC
intervals.  It is narrower: inside a labelled local disk, theorem exits do not
change; the first coordinate that leaves the disk is the unit-height lift.

The nicest surprise is that recursive chirality is no longer a separate
tournament hobby.  If owner support is the missing enlarged-bank coordinate,
then deleting endpoint owners and recording tail/tip child decks is exactly
the right next sidecar.  Chiral signature is how a quotient says, "I remember
which side of the boundary this owner lived on."

So the proof route I would hand to the next agent is:

```text
1. Build the owner-support Menger graph on HYP-3406 leaks.
2. Build the HYP-3405 local disk/height-exit collar table.
3. Add recursive chiral child decks to owner-support words.
4. Only then ask for a BDH/Mertens mean-square theorem over the enlarged bank.
```

The order matters.  Average analytic control should come after the exceptional
fibers are named, not before.
