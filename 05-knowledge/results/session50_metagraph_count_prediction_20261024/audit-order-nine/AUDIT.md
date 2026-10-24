# Independent audit: the corrected order-nine prediction

24 October 2026. **PASS.** The structural proof and independently reproduced
finite-exact computation give the ordinary tournament count 191,536 at order
nine. This is not a proof that unmarked metagraph adjacency or scalar prior
counts determine the answer. It needs the complete ordinary order-eight atlas
and its vertex actions. No mathematical error was found. A wording ambiguity,
“nontransitive” versus “non-vertex-transitive” for the regular-seven hostile,
was reported and fixed by the producer before this audit was frozen.

## Independent computational path

`audit.cpp` reads only the literal archived merged TSV, not the producer's
canonical atlas, actions, algorithms, regular matrix, refinement outputs, or
known target count. It uses a different canonicalization boundary: exhaustive
orders within degree buckets, rather than stable-refinement buckets. Every
order-eight row is checked to be a tournament; duplicate ordinary IDs must
agree, and all ordinary IDs and merged IDs are consecutive.

For each original-label parent it independently obtains the full ordinary
automorphism group by testing all degree-preserving vertex maps. Exhaustive
degree-bucket canonical keys prove all 6,880 classes distinct. Their labelled
orbit masses sum to 268,435,456 = 2^28. Distinctness and full labelled mass
certify completeness, rather than a trusted class-count assertion. The seven
automorphism-order multiplicities are 1:6556, 3:292, 5:10, 7:2, 9:16, 15:2,
21:2. There are 3,528 merged records and 176 self-converse records.

Actual original-label vertex permutations form all 256-subset orbits for each
parent. Child refinement starts from one cell and splits by **incoming**
neighbor counts, independently of the producer's outdegree/outneighbor
implementation. Every resulting stable cell is checked odd. The full pass
has 1,716,608 pointed states and upper trace 1,724,582. Its root-cell histogram
is 1:1712964, 3:3491, 5:48, 7:20, 9:85.

The regular sector is independently reconstructed from the 85 almost-regular
parents. Parent degrees force the root adjacency. All 765 deletions are
identified using order-eight degree-bucket canonicalization. There are 85
additional key-building calls in this independent implementation, hence its
reported 850 regular-sector canonical calls; this does not conflict with the
producer's 765 calls, which reuse already canonical input rows. No call
canonicalizes nine vertices. The root-changing matrix has trace 135, fifteen
components, and total error 9(85-15)=630. Rows coincide throughout each
component and their destination multiplicities sum to nine. The regular
vertex-orbit profiles, with class multiplicities, are:

- nine singleton orbits: seven classes;
- 1,1,1,3,3: two classes;
- 3,3,3: three classes;
- 9: three classes.

Only after the independent computation, `check.py` matches the producer's
parents by independent degree-canonical keys. All 6,880 per-parent attachment
state counts, upper traces, and five-bin root-cell histograms agree. All
supplied parent actions are distinct and preserve arcs; their cardinalities
match the independently exhaustive groups, so the supplied actions are
complete, not merely same-size candidate lists. Every entry of the independent
85-by-85 regular reroot matrix agrees after relabelling its states by parent
canonical keys.

## Exhaustive structural proof audit

The complete odd partitions of nine are 9; 7+1+1; 5+3+1; 5+1^4;
3+3+3; 3+3+1^3; 3+1^6; and 1^9. These exhaust the cases because each equitable
cell induces a regular tournament and therefore has odd size.

A nine-cell is precisely regularity, accounted for by the separate regular
pointed calculation. A five-cell has the unique vertex-transitive regular
five-tournament internally; all its other crossings are uniform by coprime
cell sizes. Its internal automorphisms therefore extend. With at most two
triple cells, nonuniform matching support is a forest and admits nonzero
compatible rotations on each connected component. Both categories have zero
error. Singletons contribute no error.

A seven-cell has uniform crossings to its two singleton cells. The quotient
is a tournament on three blocks, with one distinguished seven-block. There
are three transitive placements and one directed-cycle placement up to
isomorphism. The first three are exactly the ordered SCC padding cases; the
last is the unique strong substitution case. Their stable cells remain
separate. Each has the internal seven-cell error. Only the non-vertex-transitive
regular-seven class has positive error, namely 14. This category therefore
contributes 4*14=56, split as 42 padding plus 14 strong. There is no overlap
with regular nine or three-triangle errors.

For three genuine stable triple cells, each induced triangle permits only
rotations. Positive error requires a signed support cycle, hence all three
crossings are nonuniform. Their majority orientations form a three-vertex
tournament. A cyclic majority orientation makes every vertex have degree
four, so its actual stable partition is one nine-cell: it belongs to the
regular sector and must be excluded here. A transitive majority orientation
gives the three cell degrees 3,4,5, naming the cells intrinsically. Every
lower-to-higher crossing is a matching. The affine matching signs admit four
negative products. Each has full-rank phase-translation action, so one class
per sign pattern, with trivial automorphism group. The four distinct classes
therefore contribute 4*18=72. The twelve balanced classes have zero error.

These cases are disjoint and exhaustive. Thus total error is 630+56+72=758,
and the parent-only upper trace gives (1,724,582-758)/9=191,536. This is an exact
restricted correction calculation, not an unrestricted all-order recurrence.

## Independent hostile and completeness controls

`check.py` generates every labelled regular-seven tournament by deciding the
remaining neighbors of each vertex under its residual degree requirement.
The only feasibility pruning compares a remaining vertex's achieved degree
and still-undecided incident pairs against degree three. Its 2,640 outputs
are partitioned using all 7! maps into three orbits of sizes 720, 1,680 and
240, with automorphism orders 7, 3 and 21. Their errors are 0,14,0. Each
orbit passes exact orbit-stabilizer mass and full closure/disjointness checks.

The same independent script exhausts **all 14^3=2,744 equitable crossings** of
three fixed directed triangles: either uniform direction, or one of six
matchings or complementary matchings on every pair. Of these presentations,
434 are regular and excluded, including the two uniform cyclic quotients;
2,310 have actual stable partition 3+3+3. Exactly 648 presentations have
positive error. Reducing them by degree-ordered cell names and the 27 within-
cell rotations gives exactly four ordinary rigid classes, each error 18.
This checks the entire nonregular triple-cell exception, not only a selected
216-presentation family.

Normal and Python `-O` checks each pass 44,610 explicit requirements and produce
byte-identical JSON. Independent C++ `-O2` and `-O3 -DNDEBUG` full runs each
pass 19,195,918 explicit checks, with identical JSON and detailed parent and
regular-matrix tables. The compiler is not asked to trust `assert`; the
requirements remain live under optimized builds.

## Provenance, reproduction and boundary

`pins.json` records the frozen producer bytes used only in comparison and
proof review. Historical Markdown and executable source copies in `archive/`
are `.txt`, not live theorem or executable surfaces. The original TSV has
SHA256 61e4e25a0e83bbde3cdb756bd4a71febb4d0f2b8c5237143184e782738d02835.
`MANIFEST.json` is a closed, exact-path, SHA256 inventory of this audit packet.
Use `bash reproduce.sh` from this directory. A copied packet is replayed
independently before delivery, with no references to the original producer
workspace needed at runtime.

No Git or producer file was changed by the auditor. There is no authorization
blocker, target-order class census, full 2^36 tournament enumeration, known-a9
input, or unreported target canonicalization. The target value appears in
results and reporting only. The independent audit does not claim new priority
for affine signed holonomy or the already inherited 16-class family.
