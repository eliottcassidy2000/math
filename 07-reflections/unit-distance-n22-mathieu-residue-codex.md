# Unit-Distance n=22 and the Mathieu Residue

codex-2026-06-17-S4

The useful thought in the prompt is the asymmetry in the Mathieu chain.
`M24`, `M23`, and `M22` arrive by stabilizing points in the big sporadic
object.  Then the next stabilizer is no longer just "one more sporadic layer":
`M21` is the classical group `L_3(4)`, acting on the projective plane
`PG(2,4)`.

That feels exactly like the 22-point unit-distance frontier.

The published small-`n` frontier says `u(21)=57` and `60 <= u(22) <= 61`.
So the unresolved question is whether one extra point can add enough edges to
hit `61`.  Any hypothetical 61-edge graph has a low-degree deletion point:

```text
degree 4 -> 57-edge 21-core
degree 5 -> 56-edge 21-core
```

So the problem is already asking for a 22-to-21 residue.  The Mathieu lens says
not to treat that residue as a featureless 21-point set.  Under `M22`, fixing a
point exposes `M21 = PSL(3,4)`, and the residual structure is `PG(2,4)`: 21
points, 21 lines, each line of size 5.

The scout made that concrete.  It built `PG(2,4)` over `GF(4)` and classified
possible neighbor sets of the deleted/extension point:

```text
degree 5 ears:
  line_5                 count 21
  near_line_4_plus_1     count 1680
  arc_5_no_three         count 1008
  three_collinear_5      count 17640

degree 4 ears:
  line_4                 count 105
  arc_4_no_three         count 2520
  three_collinear_4      count 3360
```

The coherent cases are tiny: a full projective line, or a punctured line.  The
scattered cases are plentiful, but they have visible secant debt.  That is the
new proof split.

For coherent ears, the Euclidean geometry is sharp.  If one new point is
adjacent to `d` old points, those `d` neighbors lie on a unit circle centered at
the new point.  Among `d=4` or `d=5` points on that circle, the script checks
the maximal internal unit-chord counts:

```text
d=4: at most 3 internal unit chords
d=5: at most 4 internal unit chords
```

So a line or punctured-line ear has to pass a strict circle-cap test before it
can help reach 61 edges.  This should be paired with the existing Moser
cap-endpoint ledger: the known 22-point 60-edge lane has a cap defect, and a
61-edge ear has very little slack to repair it.

For scattered ears, the finite-plane geometry gives a different handle.  A
5-arc has 10 secants; three-collinear 5-sets have 6 or 8.  A 4-arc has 6.  That
is not a unit-distance contradiction by itself, but it is a structured
obstruction index.  It tells the next script where to look for the same kind of
side-channel that killed the graph-only 62-edge coimage candidates:
unfaithful embedding constraints that raw edge counts cannot see.

The big caveat is important.  This is not saying a planar 22-point
unit-distance drawing secretly has `M22` symmetry.  It almost surely does not.
The point is more modest and more useful: the Mathieu residue is a vocabulary
for the deletion side-channel.  The proof may need to ask whether the
degree-4/5 ear behaves like a line, a punctured line, or a scattered projective
set before it returns to Euclidean embedding data.

That also matches the repo's recurring no-leak lesson.  The graph-only
quotient is too coarse; the proof lives in a retained side channel.  Here the
side channel is:

```text
deleted point -> ear size -> PG(2,4) line/secant type -> circle cap / Moser endpoint / unfaithful obstruction
```

The next useful computation is therefore not another broad anneal.  It is an
ear-type annotator for the existing 21-core and 22-point candidates.  For each
degree-4/5 extension attempt, label the neighbor set by its `PG(2,4)` type and
record which Euclidean obstruction kills it.  If all coherent ears die by
circle-cap/Moser endpoint constraints and all scattered ears die by a finite
secant-profile obstruction library, the Mathieu thought will have become a
real proof scaffold.
