# The five-vertex origin: H15 to H33, then a fork at every height

The regular tournament on five vertices with outgoing-row
bitmasks (6,12,24,17,3) generates the entire mixed-direction
phase grid. Choose q=4 as forward-leaf port and p=1 as
reverse-leaf port, with the operations defined in OPPOSITE_FANS.md.
Writing H(a,b) for a forward and b reverse leaves,
the exact values are

- H(0,0)=15.
- H(a,0)=16·2^a+1 for a≥1.
- H(0,b)=16·2^b+1 for b≥1.
- H(a,b)=16·2^(a+b)+3 for a,b≥1.

The proof is the same two-visit skeleton interpolation
in OPPOSITE_FANS.md, now with base table

```
15   33   65
33   67  131
65  131  259
```

It follows that at the origin each designated move has defect
+3, and either direction reaches H33. From a nonorigin
point on an axis, keeping the same direction has defect −1
and changing direction has defect +1. In the interior
both defects are −3. Thus moving a single leaf assignment
between the two ports places an exact two-count phase
shift at every total height at least two; it is not an
unbounded return from the interior. The H129 construction
in OPPOSITE_FANS.md is the point (3,0), with the forward
port relabelled to be last. In the order used there, it has
backward arcs 3→0, 7→0 and 7→1.

This gives physical initial H33 and a one-time branching
construction, but does not imply that H33 was absent from
all chronological one-parent histories or that arbitrary
prescribed phases can be reached.
