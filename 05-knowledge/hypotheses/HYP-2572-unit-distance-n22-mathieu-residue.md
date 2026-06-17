# HYP-2572 - The n=22 unit-distance frontier should split by the M22 -> M21 projective-plane residue

**Status:** OPEN proof program with exact finite residue scout.

**Source:** codex-2026-06-17-S4.  Responds to the user's prompt about the
unit-distance problem for 22 points, the Mathieu chain `M24 -> M23 -> M22`,
and the break where `M21` is `L_3(4) = PSL(3,4)`.

## Claim

The current small unit-distance frontier

```text
u(21) = 57,        60 <= u(22) <= 61
```

should be attacked as a one-point extension problem over a structured
21-point residue, not as another raw graph-count problem.

In any hypothetical 61-edge 22-point unit-distance graph, deleting a
minimum-degree point gives one of two ears:

```text
degree 4 ear -> a 57-edge 21-core,
degree 5 ear -> a 56-edge 21-core.
```

The Mathieu observation supplies a candidate side-channel for those ears.
The group `M22` is attached to the Steiner system `S(3,6,22)`.  Fixing one
point leaves a point stabilizer `M21 = PSL(3,4)` acting on the `21` points of
`PG(2,4)`, and the hexads through the fixed point become exactly the `21`
projective lines of size `5`.

Thus the deletion ear should be classified inside `PG(2,4)`:

```text
degree 5: line-hexad, near-line, or scattered 5-set
degree 4: punctured-line, 3-collinear, or 4-arc
```

The conjectural proof split is:

1. Coherent ears (`line_5`, `line_4`) should be small enough to attack with
   unit-circle chord caps and the existing Moser cap-endpoint ledger.
2. Scattered ears should be fed into an obstruction library via their secant
   profiles, in the same spirit as the graph-only 62-edge coimage candidates
   that are killed only after totally-unfaithful / geometry side channels are
   retained.

This is not a claim that arbitrary planar 22-point unit-distance drawings have
Mathieu symmetry.  The claim is that the Mathieu residue gives the right
proof-carrier vocabulary for the deletion/extension side-channel that raw
edge counts forget.

## Exact Evidence

Script:

```text
04-computation/unit_distance_n22_mathieu_residue_codex.py
```

Stored output:

```text
05-knowledge/results/unit_distance_n22_mathieu_residue_codex.out
```

The script builds `PG(2,4)` directly over `GF(4)` and verifies:

```text
points = 21
lines = 21
line_size = 5
point_line_degree = 5
pairs covered by a unique line = 210/210
```

It also checks the Witt-design arithmetic:

```text
|S(3,6,22) hexads| = C(22,3)/C(6,3) = 77
hexads through one fixed point = C(21,2)/C(5,2) = 21
```

Deleting the fixed point from those `21` hexads gives the `21` lines of
`PG(2,4)`.

Degree-4 ear-neighbor subsets in the 21-point residue:

```text
arc_4_no_three:     count=2520, secants={6:2520}, profile=((0,7),(1,8),(2,6))
line_4:             count=105,  secants={1:105},  profile=((0,4),(1,16),(4,1))
three_collinear_4:  count=3360, secants={4:3360}, profile=((0,6),(1,11),(2,3),(3,1))
unit-circle cap:    at most 3 unit chords among 4 neighbors of one extension point
```

Degree-5 ear-neighbor subsets:

```text
arc_5_no_three:       count=1008,  secants={10:1008},          profile=((0,6),(1,5),(2,10))
line_5:               count=21,    secants={1:21},             profile=((1,20),(5,1))
near_line_4_plus_1:   count=1680,  secants={5:1680},           profile=((0,3),(1,13),(2,4),(4,1))
three_collinear_5:    count=17640, secants={6:7560,8:10080},   profile=((0,5),(1,8),(2,7),(3,1))
unit-circle cap:      at most 4 unit chords among 5 neighbors of one extension point
```

The line cases are rare and coherent (`21` full lines, `105` punctured lines).
The scattered cases are numerous but carry many secants.  This suggests a
finite obstruction ledger: either the ear is a design line/punctured line and
must pass a sharp circle-cap / Moser-endpoint test, or it spends many secants
and should trigger an unfaithful-geometry obstruction.

## Tournament Analysis

This session explicitly did not use points or unit edges as the default
tournament vertices.

Vertices:

```text
line_5_hexad_ear
line_4_punctured_hexad_ear
near_line_4_plus_1
three_collinear_scattered_ear
arc_scattered_ear
raw_21_core_extension
Moser_27_quantum
graph_only_F_free_coimage
```

Pairwise observable:

```text
design coherence
ear relevance
unit-circle cap sharpness
geometry retention
proof burden
scalar-forgetfulness penalty
```

The majority tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

Ranking:

```text
line_5_hexad_ear
line_4_punctured_hexad_ear
near_line_4_plus_1
Moser_27_quantum
three_collinear_scattered_ear
arc_scattered_ear
raw_21_core_extension
graph_only_F_free_coimage
```

Interpretation: the proof-bearing vertices are not runners, points, or even
graphs.  They are side-channel proof carriers: coherent design ears, scattered
ear profiles, Moser endpoint ledgers, and graph-only coimages.

## Assumption Challenge

Default assumption challenged:

```text
The natural quotient vertices for n=22 unit-distance should be points,
edges, or graph isomorphism classes.
```

Alternate vertex sets considered:

```text
deleted-point ears
PG(2,4) lines
punctured lines
secant profiles
hexads
unit-circle caps
Moser endpoint ledgers
graph-only coimage obstructions
proof obligations
```

Predicate preserved by the chosen quotient:

```text
whether the degree-4/5 extension ear is coherent with the M21 residue line
structure or scattered across many projective secants.
```

Information destroyed:

```text
actual Euclidean coordinates, exact unit-distance embedding data, and the
specific identity of the 21-core graph.
```

Reason this destruction may be acceptable:

```text
existing n=22 work already shows raw graph coimages are too lossy; the new
quotient is meant to be a side-channel triage before the Euclidean embedding
and unfaithful-subgraph tests are reattached.
```

## Next Proof Targets

1. Attach the `PG(2,4)` ear type to the stored 21-core candidates from the
   `n=22` frontier scripts and measure which types occur under degree-4 and
   degree-5 extension attempts.
2. Prove or computationally certify a coherent-ear exclusion: a `line_5`
   hexad ear or `line_4` punctured-hexad ear cannot supply the missing
   61st edge unless it repairs the known Moser cap-endpoint defect.
3. Build a scattered-ear obstruction library indexed by secant profile
   (`4`, `5`, `6`, `8`, `10` secants) and compare it against the
   totally-unfaithful subgraph killers from the graph-only 62-edge coimage.
4. Test whether the `M21` line ledger can be aligned with the exact `57`-edge
   21-cores and the `P_2^- / P_2^+` Moser spine-ladder rows.

## Cross-Links

`HYP-2176`, `HYP-2188`, `HYP-2203`, `HYP-2224`, `HYP-2467`, `THM-408`,
`T840`.
