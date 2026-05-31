---
source: codex-2026-05-31-S452
status: integration note
tags:
  - lonely-runner
  - tournaments
  - zeckendorf
  - polygonal-numbers
  - pairwise-distance
---

# LRC Runner-Distance Tournaments

The usual Lonely Runner move chooses one runner as stationary and watches a
star of distances from that runner to all others.  That is the right
one-frame problem, but it is a lossy coordinate system for the whole moving
configuration.  This session explored three lifts:

```text
star     = distances from one chosen runner to the others
cycle    = adjacent circular gaps, hence the two nearest neighbors
complete = every pairwise distance, one coordinate for each edge of K_n
```

The cycle lift is the cleanest new object, now recorded as THM-370.  Let the
`n` runner positions be ordered on the circle and let `g_i` be the gap from
position `i` to the next position.  A runner at vertex `i` is lonely exactly
when its two incident gaps are both at least `1/n`.

So, in a static configuration:

```text
safe gap = gap >= 1/n
lonely vertex = two adjacent safe gaps
no-lonely configuration = safe gaps form an independent set in C_n
```

The empty safe-gap mask is not feasible, because all gaps shorter than `1/n`
would sum to less than `1`.  Therefore no-lonely configurations are exactly
the nonempty independent sets of the circular gap graph, after imposing the
metric gap lengths.

This feels like the missing bridge between the Zeckendorf thread and the LRC
thread.  Zeckendorf is independence on a path at fugacity `x=1`; tournament
OCF is independence in an odd-cycle conflict graph at fugacity `x=2`; the
two-neighbor LRC obstruction is independence in the cycle `C_n`.  Cutting the
cycle at a chosen unsafe gap turns the obstruction into a path problem, so
Fibonacci/Zeckendorf normal forms are not imported poetically.  They appear as
the canonical accounting of safe-gap masks.

## What Changes Dynamically?

For integer speeds, the circular order is constant between collision times.
Inside one chamber, all gaps are affine functions of time.  The chamber changes
by adjacent swaps at collisions.

That makes the dynamic LRC problem look like:

```text
time chamber = circular order of runners
state word   = safe/unsafe gap mask in that chamber
bad interval = mask has no adjacent 1s
proof target = every legal chamber walk must hit adjacent 1s for some frame
```

This is very close to a tournament flip graph.  Adjacent runner swaps are the
geometric analogue of local arc flips, while safe-mask transitions are the
boundary crossing events.  The ordinary endpoint-cover proof machinery sees
when one stationary star has no uncovered time.  The cycle lift asks whether
the local chamber word can avoid adjacent safe gaps forever.

## Complete Pairwise Lift

The complete lift stores all distances

```text
d_ij(t) = ||(v_i - v_j)t||
```

for every pair of runners.  There are

```text
binom(n,2) = T_{n-1}
```

coordinates, exactly the triangular count of tournament arcs.  This is where
the polygonal-number inspiration belongs: pairwise LRC data naturally lives in
a triangular coordinate system, and the first nontrivial curvature statistic is
already triangular, namely directed 3-cycles of the circular distance
tournament.

Define a distance tournament by orienting `i -> j` when `j` lies within the
clockwise half-circle from `i`, omitting antipodal/collision ties.  This adds
score sequences, cut imbalance, and directed-triangle counts to the usual
nearest-distance statistics.

The S452 computation records the contrast:

```text
initial n=14 boundary: close-degree hist 0^14, score hist 6^14, 3-cycles 70
seven-ladder n=14:     close-degree hist 0^2 1^8 2^4, 3-cycles 110
S380 n=14:             close-degree hist 0 1^12 2,   3-cycles 106
initial n=16 boundary: close-degree hist 0^16, score hist 7^16, 3-cycles 112
```

The near-counterexamples are not just "one small gap."  They have a visible
global crowding profile in the pairwise graph, and the tournament scores detect
how that crowding is distributed around the circle.

## Reframing the Usual Routes

The usual scalar route asks for the maximum of the minimum distance from the
stationary runner.  The cycle route asks for the two incident gaps around that
runner.  This is a better local invariant: a runner can have a large nearest
distance in the stationary star while the global gap word is still paying debt
elsewhere.

The endpoint-cover route becomes a row/column shadow of the chamber-word route:
endpoint rows are the times when a gap hits exactly `1/n`, and speed columns
are ways of protecting those boundary crossings.  The new question is whether
the protected crossings can support a nonempty independent safe-gap mask in
every chamber.

The complete pairwise route is "all Galilean frames at once."  If runner `i`
is made stationary, the observed distances are the `i`-star inside the
complete pairwise matrix.  Thus a pairwise proof would not repeatedly choose a
stationary runner; it would prove a global tournament-distance invariant and
then read off a lonely star.

## Polygonal and Zeckendorf Hints

The polygonal number theorem says, morally, that high-dimensional integer
mass can be decomposed into a bounded number of structured polygonal pieces.
For LRC this suggests a search direction rather than a theorem:

```text
gap masks        -> cycle/path independent-set pieces
pairwise distances -> triangular K_n coordinates
triangle stress  -> directed 3-cycle curvature
higher stress    -> k-gon chamber patterns / polygonal decompositions
```

For `n=16` the cycle fugacity count is especially loud:

```text
I(C_16,2) = 2^16 + 1 = 65537.
```

This is not a proof, but it is a useful coordinate warning.  The pure dyadic
case is not just a p-adic tree; its two-neighbor obstruction has a Fermat-like
cycle count at the same time that its complete pairwise lift has
`T_15 = 120` triangular coordinates.  The dyadic and polygonal views meet at
the exact point where the repo already expects `n=16` to be most rigid.

## Practical Next Experiments

1. Build the chamber automaton for small speed sets: nodes are circular orders,
   edges are adjacent swaps, and intervals carry safe-gap masks.
2. For each known near-counterexample, record the safe-mask path and test
   whether masks are cuts of a path, a Fibonacci cube, or a cycle quotient.
3. Compute tournament-distance features over time: score histograms,
   close-degree histograms, directed 3-cycles, and cut imbalance.
4. Try a Farkas/Hall certificate in the cycle language: unsafe gaps must be an
   edge cover of `C_n`; show endpoint protectors cannot maintain such an edge
   cover over every chamber without exporting endpoint debt.

The session's slogan:

```text
Do not only ask how far runner 0 is from everybody.
Ask what circular gap word the whole pack is spelling.
```
