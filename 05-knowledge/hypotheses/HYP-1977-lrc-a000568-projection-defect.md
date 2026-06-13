---
id: HYP-1977
status: OPEN
source: codex-2026-06-01-S509b
related:
  - HYP-1932
  - HYP-1940
  - HYP-1951
  - HYP-1967
  - HYP-1968
  - HYP-1970
  - HYP-1971
  - HYP-1972
  - HYP-1973
  - HYP-1974
  - HYP-1975
  - HYP-1976
  - THM-372
  - THM-373
  - THM-374
---

# HYP-1977: LRC is a projection-defect problem over the A000568 quotient

## Statement

The Lonely Runner Conjecture is plausibly analogous to a problem on tournament
isomorphism classes counted by A000568, but not as a scalar property of a plain
unmarked tournament class.

The proposed structure is:

```text
runner time movie
  -> labelled half-turn / threshold tournament movie
  -> pointed tournament class movie
  -> unmarked tournament isomorphism class movie in G_N
  -> A000568 base count
```

LRC safety is a section or avoidance condition in the fiber over this A000568
base.  The obstruction is projection defect: the same unmarked, and often even
observer-pointed, tournament class can contain both safe and unsafe LRC time
cells.  Thus A000568 is the correct base quotient, while the proof-relevant
object is a sheaf of marked threshold, gap-length, and endpoint-protection data
over that base.

## Evidence

`04-computation/lrc_a000568_iso_analogy_s509.py` samples exact time cells for
small primitive speed families.  Each cell is cut by both half-turn walls and
LRC endpoint walls.  At the cell midpoint it builds the half-turn phase
tournament including the stationary observer, canonicalizes both the unmarked
class and the class with observer `0` fixed, and records the LRC safe bit.

Tournament Analysis declaration:

```text
pairwise observable: circular half-turn phase difference
switch/gauge: delta in (0,1/2), ties by increasing vertex label
tie Hamiltonian path: 0 -> 1 -> ... -> N-1
```

The probe finds mixed safe/unsafe fibers inside A000568 buckets:

```text
N=4: unmarked mixed classes 1, pointed mixed classes 2
N=5: unmarked mixed classes 3, pointed mixed classes 6
N=6: unmarked mixed classes 1, pointed mixed classes 4
N=7: unmarked mixed classes 5, pointed mixed classes 5
```

Concrete examples:

```text
N4 sparse (1,3,4):
  same unmarked class has safe t=41/96 and unsafe t=7/48
  both have H=5, score=(1,1,2,2), stationary_score=2

N5 prime-ish (2,3,5,7):
  same unmarked regular class has safe t=3/28 and unsafe t=67/560
  both have H=15, score=(2,2,2,2,2), stationary_score=2

N7 mixed (1,4,6,9,10,15):
  same unmarked class has safe t=38/105 and unsafe t=29/840
  both have H=33 and the same score sequence
```

The initial consecutive families have no safe open cells but do have safe
boundary events, matching the tight LRC behavior:

```text
N4 consecutive: 2 safe boundary events
N5 consecutive: 4 safe boundary events
N6 consecutive: 2 safe boundary events
N7 consecutive: 6 safe boundary events
```

The route-comparison tournament in the script treats twelve possible analogy
routes as vertices.  Its majority gauge over six attributes has
`H=65`, `c3=6`, and `scc_count=5`, with top routes:

```text
projection-defect obstruction
alpha=1/N gauge bundle
endpoint-pressure sheaf
G_N wall-crossing corridors
residue bucket transport
staircase extremal geodesic
```

## Interpretation

The computation refutes the strongest naive version:

```text
LRC safe/unsafe is a function of the unmarked A000568 class.
```

It also refutes the next naive version in many sampled cells:

```text
LRC safe/unsafe is a function of the observer-pointed half-turn class.
```

But this strengthens the more convoluted analogy.  LRC looks like a
projection-defect problem: the A000568 quotient is the right coarse base, and
LRC asks whether every speed-set movie admits a safe section in the richer
fiber above the quotient path.

## Predictions

1. A full proof route should use `G_N` wall-crossing corridors plus labelled
   endpoint fibers, not just scalar `H` or score sequence.
2. The useful A000568 analogue is likely a pointed/colored orbit count: one
   observer, one threshold `1/N`, endpoint signs, and possibly gap-length
   chambers over each unmarked tournament class.
3. Counterexample searches should look for quotient walks in `G_N` whose entire
   fiber over each visited class is endpoint-forbidden.  Mixed fibers are
   repair opportunities.
4. Tight initial segments should be treated as boundary sections, not open-cell
   sections.  Their A000568 class is a special boundary fiber of the staircase
   geodesic.
5. Hard n=14/n=18 rows should be projected to the same data shape by replacing
   full canonical classes with fingerprints: `H`, score sequence, `c3`, SCCs,
   observer score, endpoint labels, and pressure fibers.

## Sources

- `04-computation/lrc_a000568_iso_analogy_s509.py`
- `05-knowledge/results/lrc_a000568_iso_analogy_s509.out`
- `07-reflections/lrc-a000568-isoclass-analogy-s509.md`
