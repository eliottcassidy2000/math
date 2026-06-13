# LRC as an A000568 Projection-Defect Problem

The user hypothesis was that LRC is ultimately analogous to a problem on
tournament isomorphism classes and A000568.  The tempting simple version is:

```text
turn a runner configuration into a tournament;
take its isomorphism class;
read loneliness from that class.
```

That simple version is false in the interesting way.  The better version is a
tower:

```text
runner movie
  -> labelled tournament movie
  -> observer-pointed tournament class
  -> unmarked tournament isomorphism class
  -> A000568 base
```

LRC does not live as a function on the bottom.  It lives as a section problem in
the fibers above the bottom.

During the checkpoint rebase, a concurrent S509 pair-cell operation-grid thread
landed and claimed HYP-1976.  That was not noise; it was the missing middle
layer.  The pair-cell grid says what many of the fiber coordinates are.  This
projection-defect reflection is therefore renumbered HYP-1977 and should be read
as the base/fiber counterpart to HYP-1976.

## What the Computation Says

The script `04-computation/lrc_a000568_iso_analogy_s509.py` uses the half-turn
phase tournament with the stationary observer included.  Time is cut by both
half-turn walls and LRC endpoint walls.  Each open cell is sampled exactly, then
projected to:

- the unmarked tournament isomorphism class;
- the observer-pointed tournament class;
- standard fingerprints: `H`, score sequence, `c3`, SCC count, observer score.

The key result is projection defect.  Safe and unsafe cells can land in the
same unmarked class and even in the same pointed class.

Examples:

```text
N4 sparse (1,3,4):
  safe t=41/96 and unsafe t=7/48 have the same unmarked class,
  H=5, score=(1,1,2,2), stationary_score=2.

N5 prime-ish (2,3,5,7):
  safe t=3/28 and unsafe t=67/560 have the same regular class,
  H=15, score=(2,2,2,2,2), stationary_score=2.

N7 mixed (1,4,6,9,10,15):
  safe t=38/105 and unsafe t=29/840 share an unmarked class,
  H=33 and the same score sequence.
```

Globally in the sampled cells:

```text
N=4: unmarked mixed 1, pointed mixed 2
N=5: unmarked mixed 3, pointed mixed 6
N=6: unmarked mixed 1, pointed mixed 4
N=7: unmarked mixed 5, pointed mixed 5
```

So A000568 is not wrong.  It is too low in the stack to hold the LRC witness bit
by itself.

## Twelve Routes Through the Maze

1. **Clock-word in G_N.**  A speed set traces a word in the isomorphism class
graph.  Pairwise half-turn wall crossings are arc flips upstairs and quotient
steps downstairs.

2. **Pointed A000568 lift.**  LRC has a stationary observer.  The first repair
is not A000568 but pointed tournament classes, where vertex `0` is fixed.

3. **Alpha threshold bundle.**  Half-turn tournaments see the `1/2` gap.  LRC
needs `1/N`, so the correct object is a bundle of threshold gauges over the
same class path.

4. **Endpoint-pressure sheaf.**  Endpoint protection cores and pressure DAGs
are not extra decorations; they are the fiber coordinates that decide whether
a base class admits a safe section.

5. **Burnside fiber mass.**  A000568 is Burnside orbit counting.  LRC should
ask for the distribution of forbidden endpoint fibers over those orbits, not
only the number of orbits.

6. **G_N wall-crossing corridors.**  Hard LRC rows are short corridors through
adjacent half-turn cells.  The relevant unit is a path segment in `G_N`, not a
snapshot.

7. **A000568 entropy barrier.**  The quotient is tiny compared with labelled
phase space.  If an obstruction survives projection to A000568, it is structural
rather than accidental.

8. **Even-graph cycle-space dual.**  The even-graph metagraph is the cycle-space
dual of the tournament quotient.  LRC endpoint covers are also cycle/cut
objects, so the dual may carry the missing fiber data.

9. **Hamiltonian-path H meter.**  `H=1` detects half-turn semicircle collapse,
but LRC safety can vary at fixed `H`.  H is a coordinate on the base path, not
the section.

10. **Projection-defect obstruction.**  The same class can contain safe and
unsafe cells.  The proof target becomes: show a fully unsafe fiber selection
cannot persist along the whole speed-set movie.

11. **Residue bucket transport.**  Endpoint residues are buckets over the class
path.  The quotient-balance laws already used in the metagraph should have an
LRC endpoint analogue.

12. **Staircase extremal geodesic.**  Initial segments visit very few classes
and have only boundary witnesses.  They look like a special geodesic in the
staircase quotient, not like generic LRC behavior.

## The Shape of the Conjectural Proof

A possible proof has this form:

```text
Assume an LRC counterexample.
Then its time movie gives a quotient walk in G_N.
Every fiber over every visited class must be endpoint-forbidden.
But the bucket balance / endpoint pressure / projection-defect structure
forces either a safe boundary section, a safe open-cell section, or a pressure
cycle that violates the terminal-core constraints.
```

That is exactly analogous to A000568 because A000568 is the count of base
states.  It is not merely analogous because LRC needs the fiber geometry over
those states.

## Next Experiments

Compute the pointed-class counts by Burnside, not brute force.  Then build a
colored A000568 count with one observer and one threshold band.

Run the same projection-defect audit over all primitive speed sets in small
boxes for `N<=7`: how often is a safe section already forced by a mixed bucket?

For n=14/n=18, where full canonical classes are too expensive, use fingerprints
as class surrogates: `H`, score sequence, `c3`, observer score, pressure SCCs,
and endpoint labels.  Search for a corridor where every surrogate bucket is
unsafe and pressure-acyclic; that would be the first genuinely counterexample
shaped class path.

Build the endpoint sheaf explicitly: base vertex = tournament isomorphism class,
fiber = endpoint-owner/protector labels observed in cells mapping to that class.
The LRC theorem then asks whether every speed-set section of this sheaf has a
safe point.
