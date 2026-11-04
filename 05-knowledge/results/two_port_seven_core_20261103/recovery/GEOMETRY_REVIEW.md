# Independent review of the seven-private-point extension theorem

Status: **ACCEPTED — PROVED with an exhaustive finite classification input**.
Independent recovery/audit lane,2026-11-03. Reviewed
[geometry/PROOF.md](../geometry/PROOF.md), its complete census and metric
verifier, and the parent's independent C++ core/list controls. No blocking
mathematical or implementation defect was found.

The conclusion is correctly quantified: every proper five-colouring of
the specified fixed Moser spindle extends across every placement of at
most seven distinct additional planar points, with all induced unit edges.
It strengthens this carrier's necessary private size from seven to eight.
It does not imply a global vertex lower bound for six-chromatic unit graphs.

## 1. Classification and complete universe

The reduction at a private vertex of degree at most two is valid. Its
deletion leaves a unit graph on at most six private vertices, which is
three-choosable by the earlier proved geometric classification. Every
private list has at least three colours by the old spindle's contact bound.
After colouring the remainder, the deleted vertex avoids at most two
privately used colours. Thus a new seven-vertex obstruction has minimum
private degree three; this reduction does not assume the resulting lists
are arbitrary realizable boundary patterns.

Both local geometric filters are necessary for the full induced unit
graph. Two circles give codegree at most two. A vertex's neighbours lie
on its unit circle, where adjacency is rotation by plus or minus60 degrees.
Each rotation orbit has six points, so induced neighbourhood components
are paths or complete six-cycles. This remains true if components lie
in different angular orbits. The filter is used only as a necessary
condition, never as sufficient metric realizability.

The Python Gray-code loop visits all2^21 masks bijectively, updating both
adjacency rows and degrees for the unique toggled edge. For every survivor,
all7! permutations are expanded. The exact equality of the survivor set
and the union of the six permutation orbits is checked, not inferred from
a count. The independent C++ engine enumerates masks in ordinary increasing
order and independently expands the six full permutation orbits before
checking set coverage. Their counts agree:236,926 minimum-degree-three
graphs,8,190 filtered survivors, orbit sizes420,2520,1260,840,630,2520.
All integer ranges are tiny compared with the implemented32/64-bit bounds.

## 2. The two geometric exclusions and rigidity

For mask122616, fixing the equilateral triangle forces the other common
unit neighbour of each side. Injectivity excludes reusing its third old
vertex. The three reflected caps have pairwise distance two. No fourth
point is one unit from all three: each pair of radius-one circles is
tangent, and the three tangent points are distinct. This eliminates the
whole realization class rather than a chosen coordinate candidate.

For mask380278, deleting edge(2,3) leaves a graph isomorphic to Moser.
The isomorphism is checked by all permutations. The rigidity proof is
essential and correct: each diamond has forced opposite tips at distance
sqrt(3), and the closing edge fixes the two long tip vectors' inner
product at5/2. Their two relative angular choices differ by reflection.
Permuting a diamond's two wings does not change the point set. Thus all
injective realizations of that graph are congruent to the fixed seven
points. The exact complete-pair test excludes every additional unit edge,
so the extra edge in380278 is impossible.

The same rigidity justifies applying the external-contact theorem in the
reverse direction. An old boundary point is exterior to the private
spindle because the two point sets are disjoint. It therefore has at most
two private neighbours. A mere graph isomorphism would not have supplied
this reciprocal metric bound; the proof does not omit that step.

## 3. The list arguments have no palette cutoff

The wheel proof is the standard direct even-cycle list argument after
colouring the hub. Both theta arguments were checked against their exact
edge sets. In120444, deleting0 leaves the three paths1-2,1-6-2 and
1-5-3-4-2. The two ends were not adjacent to0 and retain three colours.
After choosing1's colour outside the two-list at6, the remaining path
6-2-4-3-5 has two colours at6, at least one at5, and at least two at every
internal vertex. Greedy colouring from5 succeeds. The analogous check
for120568 uses ends0,4, the paths0-4,0-5-3-4 and0-6-2-4, and colours
6-2-4-3-5 towards5. All earlier incident colours are removed before these
list-size bounds are used.

The Moser bad-list classification is also an iff, as required. For
three-element lists, if a colour at the degree-four centre is absent
from any wing list, using it leaves a connected graph with lists at least
its degrees and a strict surplus at that wing. Rooting a spanning tree
at the surplus vertex gives a direct greedy extension. Hence every wing
list equals the centre's three-set S. Its adjacent wing pairs then use
exactly S minus the centre colour. If either tip list differs from S,
choose a centre colour leaving that tip at least two choices; the other
tip has at least one, so their closing edge can be coloured. Uniform S
at all seven vertices is precisely the remaining failure.

For lists larger than three, an alleged failure would make every choice
of three-subsets a failing assignment, and therefore uniform. Changing
just one three-subset in a larger list makes the selected assignment
nonuniform and supplies a colouring. This establishes the stated lemma
for arbitrary palettes and all list sizes at least three.

The C++ search over four million normalized five-palette assignments
corroborates this proof. Fixing vertex0's list to{0,1,2} is legitimate by
global palette permutation; all10^6 choices for the remaining lists are
exhausted for each graph. Its recursive search branches over every
available colour and verifies each returned witness. Only the Moser
graph fails, only on the uniform normalized assignment. These controls
are not used to extrapolate from five colours to arbitrary palettes.

## 4. Abstract bad lists fail the geometric capacity test

With five named boundary colours, a uniform three-set S at every private
vertex requires exactly two old contacts of distinct colours, namely the
two colours outside S. The old-to-private demand is therefore14 edges.
The old spindle has independence number two: its six noncentral vertices
form two triangles, while an independent triple involving the centre
would need both nonadjacent-to-centre tips, which are adjacent to each
other. Each excluded colour consequently occupies at most two old points.

There are at most four supplying boundary points. The reciprocal bound
proved above permits at most two contacts from each to the private rigid
spindle. Capacity is at most8, strictly below demand14. The proof includes
unused boundary colours, which only lower that capacity. This excludes
every uniform abstract bad list without asserting that arbitrary lists
can be realized by boundary contacts.

The mechanism can be stated independently. Suppose every uncolourable
private list assignment of sizes at least q-d is a uniform(q-d)-set,
where q>d>=1; private vertices have at most d boundary contacts; every
boundary vertex has at most r private contacts; and the boundary has
independence number at most a. A bad extension demands d|K| edges from
the d excluded colour classes, containing at most da vertices and supplying
at most rda edges. Necessarily |K|<=ra. The present contradiction is7>2*2.
The classification and rigidity are what certify this abstract inequality's
hypotheses in the actual geometric problem.

## 5. Exact replay and artifact scope

I reran geometry/verify.py once against its already frozen expected census
and artifacts. Its output agrees exactly with geometry/receipt.json.
The review receipt records current source hashes, unchanged inputs, and
matching normal/optimized parent core and list receipts. The independent
C++ censuses were reviewed and their receipts checked; they were not
rerun yet again by this reviewer.

The verifier uses exact arithmetic in Q(sqrt(3),sqrt(11)), recomputes
all35 old triples by both circumcentre and Gram-determinant identities,
and confirms that every unit circumcentre is already an old point.
The graph constructor checks injectivity and every pair. The translated
fourteen-point example therefore has exactly30 edges, including all eight
cross contacts; no additional contact is discarded.

The boundary enumeration uses canonical equality partitions and includes
all42 proper five-colour patterns:16 with four blocks and26 with five.
Their orbit weights are each120 at five colours. Every one of their
private extension counts is positive; all420 excluded-colour-pair capacity
rows are retained. The computed totals4,360,320 five-colourings and16,608
four-colourings are consistency controls, not the universal proof. The
example is four-chromatic and is correctly labelled a positive interacting
control rather than a six-chromatic witness.

See [geometry_review_receipt.json](geometry_review_receipt.json) for the
audited source/result hashes and replay comparison. No requested theorem
is left conditional on an unrun finite classifier or an unproved list
claim; the acknowledged dependence is the explicit exhaustive seven-graph
classification, with an independent implementation.
