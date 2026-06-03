---
id: HYP-2121
status: OPEN proof-route correction; S586 gives exact tournament-count and extension-map evidence
source: codex-2026-06-03-S586
related: [HYP-2120, HYP-1977, HYP-1978, HYP-1981, HYP-2113, HYP-2109, HYP-1824, HYP-1825, THM-260, THM-381, THM-385]
---

# HYP-2121: the perspective gap is an observer-coupling defect, not an orphan-class defect

## Claim

The old tournament perspective curiosity should be read as a quotient-stack
warning for LRC.

Let

```text
U(n) = number of unlabelled tournaments on n vertices,
P(n) = sum over n-tournament classes of their vertex-orbit counts.
```

Then `P(n)` is the number of pointed/rooted tournament classes.  The small
coincidences

```text
P(3)=4=U(4),
P(4)=12=U(5)
```

come from a rank-shift count between an observer-coupled object at level `n`
and an observer-blind object at level `n+1`.  They are not a natural extension
bijection.  The first failure

```text
P(5)=48,
U(6)=56,
gap=8
```

is the first visible place where a rooted perspective does not carry enough
payload to determine the next unrooted structure.

For LRC, this says:

```text
observer-blind class
  -> observer-pointed/rooted class
  -> root plus incident threshold word
  -> endpoint-owner / gap-pressure sheaf
  -> proof-obligation automaton.
```

The safe-box predicate lives above the second or third arrow, not at the
observer-blind base.

## Evidence

S586 rechecked the repository occurrences and the exact counts.

The user's literal small examples occur in S25g:

```text
n=3: transitive 3 + cyclic 1 = 4 perspectives.
n=4: orbit-count contributions 4 + 4 + 2 + 2 = 12.
```

The exhaustive count table is:

```text
n  U(n)  P(n)  U(n+1)  U(n+1)-P(n)
1     1     1       1             0
2     1     2       2             0
3     2     4       4             0
4     4    12      12             0
5    12    48      56             8
6    56   296     456           160
```

Orbit-count distributions:

```text
n=3: {1:1, 3:1}
n=4: {2:2, 4:2}
n=5: {1:1, 3:4, 5:7}
n=6: {2:5, 4:10, 6:41}
```

The old extension-map language had a hidden flaw.  If one fixes an `n`-class,
adds a new vertex, and ranges over all `2^n` incident-edge words, then the
chosen root is not used.  S586 checked:

```text
n=3->4: rooted=4,  target=4,  orphans=0, root-sensitive parent classes=0
n=4->5: rooted=12, target=12, orphans=0, root-sensitive parent classes=0
n=5->6: rooted=48, target=56, orphans=0, root-sensitive parent classes=0
```

Thus the `gap=8` is not eight six-classes unreachable from five-classes.  All
six-classes are reached.  The defect is that a rooted five-perspective does not
remember the incident word/cocycle needed to distinguish the six-level
extension.

This corrects and sharpens S369/S370.  Their connection between the same `8`
and the eight LRC alpha stencils remains valuable, but the mechanism should be
phrased as projection failure of observer-coupled payload, not as orphan
six-tournament classes.

S509/HYP-1977 supplies the LRC negative evidence: unmarked A000568 classes, and
often even observer-pointed classes, can mix safe and unsafe cells.  Therefore
rooting is necessary but not sufficient.

HYP-2120 supplies the complementary exact positive slice: the LRC-relevant
source-perspective recursion is gapless because deleting a source is canonical.
This hypothesis is about the full rooted-perspective extension map: it explains
why the old full-count `gap=8` is not an orphan-class fact and why any LRC use
of full perspectives must retain the incident threshold fiber.

## Clocks That Matter

These clocks can be used as caches, heuristics, or coarse stratifications, but
not as exact safety predicates:

```text
raw unmarked A000568 class,
raw H / score sequence / c3 alone,
plain rooted class when the root does not alter the switch,
extension maps that range over all incident words and therefore ignore the root.
```

These clocks must remain in the proof state:

```text
observer-pointed class,
observer-source threshold arcs,
endpoint/gap labels,
endpoint-owner residues,
pressure SCCs and good-cut cores,
middle-state automata,
proof-obligation gates with Cprime/Phi terminals.
```

For the question "does this orbit hit the safe box even once?", the efficient
model should be:

```text
periodic or equidistributed time orbit
  -> wall-crossing word in a coarse quotient
  -> labelled fiber over each quotient state
  -> accepting condition: some observer-source threshold state with compatible
     endpoint/gap labels.
```

The unmarked class may tell us where the walk is, but the labelled fiber tells
us whether the safe box was actually hit.

## Proof Route

1. Treat observer-blind A000568 classes as the base graph only.
2. Lift each visited base state to observer-pointed threshold states.
3. Attach incident words as cocycles: for each observer, record the thresholded
   relation to every runner, not only the unmarked tournament class.
4. Attach endpoint-owner and gap-pressure labels.  This is where HYP-1977's
   mixed pointed buckets split.
5. Compile the orbit into a finite automaton: states are labelled quotient
   fibers; accepting states are exact observer-source/safe-box hits.
6. Use THM-381/THM-385 for exact observer-source semantics, HYP-2113 for the
   tournament speedup stack, and HYP-2109 for middle-state transitions.

## Assumption Challenge

Tournament vertices need not be runners.  S586 used quotient models as the
Tournament Analysis vertices:

```text
unmarked_A000568_class,
observer_pointed_class,
root_plus_incident_word,
observer_source_threshold,
endpoint_owner_sheaf,
proof_obligation_automaton.
```

This quotient preserves the predicate "existence of an LRC safe-box hit after
projection" only once endpoint/residue fibers are retained.  It destroys exact
phase duration and arithmetic if those labels are forgotten.

## Predictions

1. Any future projection-defect audit that keeps only unmarked or pointed
   tournament classes will continue to find mixed safe/unsafe buckets.
2. Adding threshold incident words should split many mixed pointed buckets.
3. The remaining mixed buckets should be separated by endpoint-owner labels and
   pressure SCCs.
4. The eight-stencil LRC residue should be expressible as a small family of
   observer-coupled cocycles over the five-to-six perspective gap.
5. A successful n=14 proof engine will not enumerate raw time.  It will search
   a labelled quotient automaton and discharge residual middle states with
   Phi/Cprime/endpoint gates.
