# Independent audit of zero-first, nonneutral repair

**PASS — full proof/source review, normal/optimized/frozen replay, and an
independent named-graph and observer reconstruction.** This audits the three
finite programs in [DIRECT.md](DIRECT.md). It accepts a genuine escape from
the previously excluded recipe, while retaining the exact failed-clock-lift
boundary. It does not assert an unbounded returning class.

## 1. Source, order and the change of operation

The source is the actual eleven-name prepared tournament, physical root8,
H1429 and F1563. I independently reconstructed the two C2 births at old
auxiliary9. The intermediate has `(H,F)=(2383,2605)` and the common thirteen-
name zero child has `(4287,4689)`. All three graphs are strong. The old names
and root8 remain unchanged; the permanent newborns are11 and12.

Each repair begins at that same literal thirteen-name graph. It changes
old arcs and sends H4287 to4289 or4291. In particular this is not an H-neutral
preparation before a nonzero C2 reset. The operation's order bookkeeping
therefore differs from ORDER_OBSTRUCTION's hypothesis: only two permanent
names have been spent when the repaired, nonzero-ready target is reached.
There is no contradiction with the accepted retained-order theorem.

## 2. All three physical programs and their joint interfaces

I read the complete source, its pinned imports, explicit permutations,
old-pair schedules, observers, C2 consumers and suffix controls. The6,7,9
pair changes are all implemented by the accepted helper/two-complement
compiler. The source reconstructs each ready-language target and identifies
it on the existing names; this identification is not a free relabeling of
the input. The native program implements the entire adjacency difference.

A separate verifier with no producer imports rebuilds every helper birth
and every directed guard, then performs each complement by complementing
all arcs and restoring precisely the consecutive guard arcs. It checks
all44 complements, all22 helper-only deletions, and that every completed
step changes exactly the declared old pair. Independent reachability checks
verify all66 helper/native states and every completed adjacency state are
strong. A fresh forward endpoint recurrence recomputes H,F at old root8,
and the root defect, throughout those states. It also checks the common
permanent-birth prefix and all immediate final children.

The verifier separately constructs the three canonical funnel/birth/converse
graphs and checks the full physical maps, including root8 and auxiliaries
12,12,11. In each target the auxiliary's sole exit is old root8. Direct
counts of the graph, its converse and both deletion graphs verify:

|Repair|H|F=G|S|L|K|Sigma|Tau|
|---|---:|---:|---:|---:|---:|---:|---:|
|A|4289|7653|4291|1253|6400|2|0|
|B|4291|6403|4292|1027|5376|1|0|
|C|4291|8371|4293|1715|6656|2|0|

The accepted sole-exit identity establishes `G=L+K=F` at the actual marks;
it is not inferred from H alone. In particular the changed F during repair
is visible, and no F-neutrality assertion is made.

Both following C2 births are strong in every case. The actual fifteen-name
children have `(H,F,defect)`

```
A: (12871,22959,10088),
B: (12875,19209,6334),
C: (12877,25113,12236).
```

Every old permanent name survives. Temporary helper order is14; the ready
targets have order13, and the next C2 children order15. Thus A really gives
a strong numerical(+2,+4) path from the old source, through a different
operation from the earlier H-neutral preparation recipe. The proof is
precise about this distinction.

## 3. Exact unbounded test and its failure

All flipped pairs lie outside old vertex0. The module expansion replacing
that vertex by TT_k is well typed: root r beats the module, and the module
beats every other outside vertex. A strong contracted graph lifts to a
strong expanded graph because an outside return path through the module
connects any two of its vertices, regardless of their internal order.
This proves the final target's strongness at every k>=1. No uniform native
guard program for every expanded clock is asserted here.

I independently recomputed all seven suffix coordinates from literal
endpoint tables. With the **incoming-to-root** port convention, they are

```
A: (11,6400,3821,5651,3364,11,11),
B: (9,5376,3197,3587,2112,9,9),
C: (13,6656,4179,6513,4080,13,13).
```

The inherited prefix assembly formula therefore gives exactly the three
H and F formulas displayed in DIRECT §4. The auxiliary still has only r
as exit. After removing it, every Hamilton path ending at r must traverse
the entire prefix first, so its deletion-ending count remains2,1,2.
Consequently the positive gap and zero joint deficit persist for every k.

The requested first-digit maps from the original clock have dyadic
coefficient2268, whereas these three targets have3821,3197,4179. Their exact
errors are1553,929,1911 times `(2^k-2)`. They vanish exactly at k1. The
independent verifier constructs the first inflated clock literally and
recomputes every joint observer, obtaining H11931,10685,12649 against the
requested8825,8827,8827. The all-clock failure follows from exact affine
formulas, not fitting these controls.

The draft independent checker initially used the opposite port convention
in this inflated construction. The literal k2 comparison rejected it; the
checker was repaired to the declared incoming-root convention before
freezing. No producer or mathematical statement required correction. This
makes the gauge distinction explicit in the final audit source.

## 4. Reproduction, pins and accepted boundary

From the session worktree:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_direct.py --output /tmp/nonzero_direct_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_direct.py --output /tmp/nonzero_direct_audit_optimized.json
python3 -B 04-computation/nonzero_preparation_20261103/surgery_direct_audit.py
```

Both producer modes reproduce the frozen82,808-byte receipt and191-byte
stdout byte for byte:600 gates,3 programs,44 complements,66 strong native
states. The independent verifier passes378 checks and reproduces
[DIRECT_AUDIT.out](DIRECT_AUDIT.out). Every DIRECT_MANIFEST length and hash
was checked; producer files were preserved.

|Artifact|Bytes|SHA256|
|---|---:|---|
|surgery_direct.py|7746|10c5e39d6439057818622607de83d97d666eddf187c3ff30099e051702538db9|
|direct.json|82808|316d0dc774d6399587b2269cd29403052e812f2ac98cc04d71cbbfab7453f2df|
|direct.out|191|732042888eb6750d76085122be7a17cbea4f187aedd6786685c6d303e815f29b|
|DIRECT.md, audited prepromotion version|6912|c8ed805ff6c87ce2856f0f4287d421f241a64498f4265f988b903b755302574c|
|surgery_direct_audit.py|7060|15b7f3fba5ea12e1ea0fc8e3895da1d69223e2bb36981f0e5256990e8654c254|
|DIRECT_AUDIT.out|127|54a9ae3cd730710fee01e42a29087ae067eec2ef2a6a099a5d64c3d7ced68559|

No correction remains. The accepted positive result is three actual finite
all-strong programs, with retained names/root and specified positive joint
output guards. Their immediate next C2 steps are paid. Return of those next
outputs, or a different repair valid for all clocks, remains unproved.
Status promotion and integration belong to the parent.
