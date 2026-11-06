# Deferred active edges give an exact register-state representation

Status: PROVED and independently audited in the stated representation scope.
This supplies a physical partial-graph state representation. It does not establish
an eventual odd-spectrum threshold or the staggered profile condition.
This uses standard path-forest dynamic programming; no priority claim is made.

## 1. State convention

Fix W boundary retained ports. Introduce-before-forget updates temporarily
use W+1 ports; this peak width is distinct from the boundary width. At a boundary, every edge having a bulk endpoint is
already fixed. Edges between retained ports are deliberately absent from the
counted digraph. A final tournament is obtained by choosing any tournament
orientation on the retained ports. Consequently those active edges can be
rewritten between stages: they were not used in the maintained counts.

Count unordered spanning directed path forests in the partial digraph. A path
component records only its start and end, each either a port label or the
symbol * for bulk. Keep at most p bulk-start components. Active vertices that
are internal to paths are recognized by their absence from endpoint positions;
the active label set is part of finite control. Active singleton paths have
endpoint pair (u,u). Identical (*,*) components are counted as a multiset.
Every kept forest has at most W+p components, so the state dimension is finite.

For ordinary Hamilton-path readouts p=1 suffices. The formal pattern dimensions
are 28 for W=2 and 123 for W=3. These include unreachable formal patterns;
all actually maintained coordinates count literal forests.

## 2. Three exact nonnegative transfers

**Introduce u.** Every bulk vertex beats u, and every incidence between u and
old ports remains deferred. Thus either u is a singleton path or it is
appended to exactly one component whose end is bulk. Sum these possibilities.
Deleting u uniquely recovers its antecedent forest and choice. This is an
exact nonnegative integer linear map.

**Process a deferred edge a->b.** Both endpoints are active. Preserve each
old forest (the new edge is unused). Additionally, for each distinct pair
of components whose endpoints are respectively a and b, join the component
ending at a to the component starting at b. Components must be distinct,
since otherwise a directed cycle would result. Removing that edge uniquely
recovers the two antecedent components. This is another exact nonnegative
integer linear map. Processing is used only within forgetting or terminal-completion batches,
with each physical edge processed exactly once when its batch is reached.
One could instead track processed active edges in finite control; that more
general instruction convention is not needed here.

**Forget u.** First choose and process every still-deferred incidence between
u and the surviving active ports. Then replace endpoint label u by * and
discard forests with more than p bulk starts. This discarding is safe: future
vertices lose to every bulk vertex, and future processed edges have only
active endpoints. A bulk-start component can never acquire a predecessor or
cease to be a bulk-start component. Thus an omitted forest can never contribute
to a final forest having at most p paths. Regrouping the remaining patterns
is an exact nonnegative integer linear map.

Any relabeling of retained ports is a permutation of coordinates. Terminal
completion processes all remaining deferred edges in a chosen tournament
orientation. Summing one-component forest coordinates then gives exactly H.
For p=2, keep a designated root r active through completion. Let h=H,
s count one-component paths starting at r, and L count ordered spanning
covers (P,Q), where P is nonempty and ends at r and Q may be empty.
Read L by summing one- and two-component patterns with a component ending
at r. That component is unique and determines P, so no extra ordering
factor is needed. The usual physical profile is B=h+L-s. The standalone
implementation returns (h,s,L,B); its L,B values require p>=2. The ordinary
h,s readouts need only p>=1.

## 3. Physical meaning and equivalence of final realizable classes

A boundary partial digraph is a tournament except for its absent active-active
edges. Introductions and forgetting fix every new bulk incidence exactly once.
Terminal completion therefore always produces a literal tournament. Equally,
one may complete a boundary to a tournament at every stage, then discard its
active-active orientations before the next update. The next displayed
completion can reverse these edges, so the displayed tournament sequence is
not required to be induced-monotone.

This does not reinterpret the old fully committed forest state as a state for
an edited graph. Its counts differ because active-active edges have been absent
throughout. That distinction is essential: endpoint data alone cannot erase
an edge already used internally in a counted path.

Every completed finite word also has an ordinary induced-register realization
with the same introduction and forgetting schedule. Choose each active-active
edge's eventual orientation in advance, and install it when its later endpoint
is introduced. If the earlier endpoint has already been forgotten, that edge
is forward in both constructions. These choices reproduce every final edge.
Thus deferred processing does not enlarge the final realizable tournament
class relative to the ordinary register construction of the same width.

What changes is the state representation and the compatibility required of
intermediate readouts. Two continuations of one partial-graph state may choose
different orientations of its still-deferred edges. Their completed prefixes
need not be the same induced tournament. Therefore an obstruction that assumes
one fixed completed physical prefix needs its hypotheses checked before it is
applied to such a branching readout. No claim is made to escape every
induced-birth obstruction, to strictly enlarge O/Z/J profile images, or to
establish arithmetic coverage.

## 4. Reproduction and scope

The standalone verifier regenerates formal patterns and exact transfers without
cached matrices or an inherited seed bank. It compares them against literal
subset path-cover counting and checks terminal completions. Its tests include
both p=1 and p=2, while the ordinary count readout needs only p=1.

Run `python3 verify.py` and `python3 -O verify.py`. Their receipts must agree.
Checks use explicit exceptions rather than Python assertions. The finite
controls test the implementation; Sections 1-3 establish the all-order lemma.
Exploratory digit searches are excluded from this packet and are not proof
dependencies. No all-base numerical occupancy or eventual threshold follows.
