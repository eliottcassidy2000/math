# The minimal singleton carrier and the price of an interior selector

2026-10-28. Local identifiers M1--M5 refer to this file. **PROVED** denotes
the arguments below, not an independent-agent review or a reserved canon ID.
Finite checks are labelled **FINITE-EXACT**. No priority claim is made.

## M1. The precise language and the residual blue state

Use [S1--S3 of the inherited packet](../singleton_pair_reconstruction_20261028/THEORY.md).
An element of X_n is a tournament on positions 0,...,n-1 containing the forward
Hamilton path. Its free entries are x_ij, j>=i+2. Native complement tau flips
these entries; it is not global converse. J_L prepends a universal source and
J_R appends a universal sink. The observations are ordinary isomorphism class
q and b, the converse/reversal fixed flag, suppressed at n<3. Re-marking means
a fresh uniform Hamilton path, or the entire counted multiset of such paths.
It does not mean a named legality query or a selected interior attachment.

Let E(T) be the isomorphism class with ordered first and last vertices marked.
The inherited sufficient carrier retained E(T), E(tau T), and a boundary
record Psi of O(n) bits. Here that last record is replaced by B(T):

* At n=1,2 there is just the corresponding tiny state.
* flat(c), c in {0,1}, means every free entry is c.
* center means blue and nonflat.
* shift(d,c), d nonzero, means the unique feasible reflection offset is d,
  with boundary polarity c. For d>0 the first d free rows are all the SAME
  constant c and the remaining suffix is blue. For d<0 the last -d free
  columns are all the SAME constant c and the remaining prefix is blue.
* dead means none of these conditions holds.

Empty and one-/two-vertex intervals are reflection invariant when checking
these conditions; suppression concerns the observed flag only. Nonflat live
offsets satisfy 1<=|d|<=n-3. The fact that the offset is unique is proved next.
Only flat and center have observed b=1 at n>=3.

### Why a homogeneous boundary is forced

Consider any deterministic word in tau,J_L,J_R. Suppose it has inserted l
vertices on the left and r on the right, and write d=r-l. Up to the total
tau parity, the old relation survives on a contiguous interval. Reflection
sends an old position i to old-relative position n-1+d-i.

For d>0, the first d old vertices reflect to the right exterior. Each exterior
vertex has one common free-entry value toward all old vertices, except for a
possible consecutive pair, which has no free coordinate. Let these values be
c_i. For two affected old vertices i,j<=n-3, compare the reflected cross-pairs
(i,reflection(j)) and (j,reflection(i)). These are nonconsecutive and force
c_i=c_j. Reflected old-old pairs then force every free entry of each affected
row to equal this same c: when both reflected vertices are exterior, the
later right insertion supplies the value. The unaffected suffix must itself
be blue. If d>=n-2 this forces all free entries to be the same, hence flat.
For d<0 the argument is reversed and uses terminal columns. At d=0 the
restriction to the entire old interval must be blue.

Conversely the stated row/suffix condition is sufficient: append d sinks
when c=1. For c=0, use tau J_R tau at each insertion. The column/prefix case
uses the corresponding left insertions. Thus the feasibility test is exact;
arbitrary exterior words may add further restrictions, but cannot introduce
a new feasible offset.

### Uniqueness of the nonflat offset

Fix polarity 1. Peel the maximal initial string of universal sources and
the maximal final string of universal sinks, of lengths s,t. Because the
identity is a Hamilton path, a universal source or sink must be the first
or last position. Unless the tournament is flat(1), a nonempty nontransitive
core remains. A blue suffix after removing d sources is possible exactly
when the remaining source and sink strings have equal lengths and the core
is blue. Hence its only possible offset is d=s-t. The same formula covers
negative d by removing sinks. Polarity 0 reduces to this argument by tau.

Two positive feasible offsets necessarily use the same first-row bit; two
negative ones use the same last-column bit. A positive and a negative one
use the same bit x_(0,n-1). A blue nonflat state with a constant boundary of
either polarity has equal peeled lengths and therefore cannot have a second,
nonzero offset. These facts establish uniqueness across both polarities.

### Closed transitions and the exact minimum

Tau swaps c with 1-c in flat and shift, and fixes center and dead. Let e=+1
for J_L and e=-1 for J_R. At n>=3:

| state | after the indicated singleton insertion |
|---|---|
| flat(1) | flat(1) |
| flat(0) or center | shift(e,1) |
| shift(d,1) | shift(d+e,1), or center if d+e=0 |
| shift(d,0) | dead |
| dead | dead |

At order1 insertion gives tiny(2), and at order2 it gives flat(1). These rules
follow from the homogeneous-boundary test: a new free row/column has polarity1;
the required difference r-l changes by e. A polarity0 live boundary conflicts
with that new row/column, whereas a center acquires exactly one unit of offset.
The argument includes both necessity and sufficiency, not merely a recurrence
suggested by a census.

Distinct states have distinct future blue languages. A shift(d,c) is exposed
by |d| pure insertions of polarity c on its required side; another nonflat
state with different d or c fails that test. Center is separated from dead
and shift immediately. One insertion of the chosen polarity distinguishes
center from flat, and distinguishes the two flat states. Therefore B is the
coarsest deterministic blue carrier, separately at each size.

Every claimed state is realized. Flat states are immediate. For n>=4, all
free bits1 except x_(0,n-1)=0 give a nonflat center. For n>=5, all free bits0
except x_(0,2)=x_(1,n-1)=1 give mixed boundaries and a nonblue dead state.
For each 1<=d<=n-3, prepend d universal sources to flat(0) on n-d vertices;
this realizes shift(d,1). Append sinks for negative d, and apply tau for
polarity0. At n=4 the only missing kind is dead; the eight presentations
give seven states directly. At n=3 there are only the two flats.

Consequently the exact numbers of blue residual states are

    1, 1, 2, 7, and 4n-8 for n>=5.                      (M1)

This is O(log n) bits of residual information, not O(log n) bits for the
whole tournament. Standalone B is not claimed to close fresh reset, whose
law depends on q. **FINITE-EXACT:** all presentations through7 and all three
literal one-step updates agree, giving1,1,2,7,12,16,20 states.

## M2. An endpoint decoder using four singleton insertions

Define the allowed deterministic context, composed right to left,

    D = tau J_R^2 J_L^2 tau.                           (M2)

It adds positions a,b on the left and c,d on the right. The old induced
tournament is unchanged. All newly introduced free bits are0 and consecutive
arcs are forward. Vertices a,b are exactly the vertices of outdegree1;
c,d are exactly those of indegree1. Every old vertex has both degrees>=2.
For n>=2 the old first vertex has its old forward neighbor and the new a as
out-neighbors; the old last vertex has its old predecessor and the new d as
in-neighbors. Every other old vertex gains the necessary two exterior
neighbors. For n=1 its unique old vertex has both degrees2 directly.

Thus the four exterior vertices are intrinsically recognizable, without
labels. Delete them to recover T. The unique old vertex beaten by a member
of the outdegree1 class is the old first vertex; the unique old vertex
beating a member of the indegree1 class is the old last vertex. Hence q(D(T))
determines E(T) at every size, including when its two ports coincide.

This extends the inherited two-vertex **strong-core** decoder
[U1](../unnamed_carrier_dissection_20261026/README.md): the two extra vertices
remove the strongness hypothesis and make recognition degree-only. Both
contexts can be realized by singleton/complement words when their ports
are the current endpoints. We do not claim four is the least possible
number of inserted vertices.

**FINITE-EXACT:** every T in X_n, 1<=n<=7, is decoded literally and after
a varying cyclic relabeling. The proof, not the finite test, supplies all n.

## M3. Exact minimal full carrier and a finite decision horizon

For precisely M1's language, set

    K_min(T) = (E(T), E(tau T), B(T)).                  (M3)

**Theorem.** Equality of K_min is exactly equality of all future observation
laws, for arbitrary finite adaptive use of the specified operations.

Sufficiency: the two endpoint-rooted types update by inherited S3. M1 gives
the blue update. For fresh re-marking, every Hamilton presentation of q(T)
has |Aut(T)| preimages among the H(T) paths, so the entire next-state law is
class-constant before quotienting. Its pushforward to K_min is therefore
well-defined. Counted reset has the same property, with multiplicities.

Necessity: if E(T) differs, the ordinary observation after D distinguishes.
If E(tau T) differs, use D after tau. If B differs, use a pure completion
from M1. Thus any carrier retaining all future observation laws must refine
K_min. This proves minimum as a quotient, not uniqueness of its encoding
or independence of its three coordinates.

Every unequal pair at size n has a deterministic separating word using at
most max(4,n-3) inserted vertices. Therefore refinement up to order

    max(n+4, 2n-3)                                   (M3 horizon)

already equals the all-future quotient at order n. Fresh resets are not
needed to separate unequal K_min states. This gives a stopping certificate
for finite refinement, without constructing the larger state spaces.

**FINITE-EXACT:** at n=7, all inherited endpoint-pair/blue cells have constant
B. The inherited exact endpoint bank has32272 cells among32768 presentations.
Since B also determines blue, there can be neither extra merging nor extra
splitting:32272 is now the exact all-future minimum. The inherited smaller
orders are discrete. The old minimality OPEN statement is superseded for
this language only; arbitrary joint-port assembly is not settled.

## M4. Weakest selectors for the seven-vertex pair: specify the cost

Keep the strong pair A=1072, B=4176 in the inherited free-bit encoding:

    A rows = (66,12,41,17,39,75,30)
    B rows = (66,20,9,51,37,71,30).

Their unique endpoint-preserving isomorphisms are respectively

    A -> B:          (0,1,4,2,3,5,6),
    tau A -> tau B:  (0,2,4,1,5,3,6).

Both B-states are dead. M3 therefore proves that zero additional selector
information cannot separate them at any later size.

**Port-arity minimum: one.** Mark any one current interior position i in
both rooted relation observations. For i=2,3,4 the first displayed map
already fails; for i=1,5 the complemented map fails. Endpoints do not help.
The current second vertex (position1) is the shallowest interior choice.
This is a minimum relative to the number of extra marked port roles, not
a claim that all possible one-bit observations are ordered by strength.

**Concrete selected reset.** Complement first, then re-mark uniformly among
Hamilton paths starting at the current second vertex, and complement again.
The eligible path counts are23 and21. Even after normalization, the ordinary
partner classes differ. The event "the resulting tournament has29 Hamilton
paths" has probabilities1/23 and0. This uses an ordinary isomorphism invariant
already determined by q. Before the first complement the same selection
has13 paths on each side and identical partner distributions, a positive
control for the alignment explanation.

**One Boolean legality test.** The endpoint-preserving word

    (0,1,3,4,2,5,6)

is a Hamilton path of A and not of B. It moves only three consecutive
positions. Exhausting all5040 permutation words finds178 distinguishing
guards and none with fewer than3 moved positions or smaller moved-position
span. Thus support3/span3 is sharp for this specified cost, and a single
Boolean output suffices. This is the named-guard menu, outside M1.

**FINITE-EXACT strengthening:** adding the position1 mark to both endpoint
types separates every presentation in X_7:32272 cells become32768. That
finite result does not extend to all sizes, as the following theorem shows.

## M5. No fixed boundary width reconstructs presentations

Fix a positive width w. Enhance the observer to retain the rooted types
of T and tau T with all first w and last w positions marked (overlap is
allowed), along with the original q,b. Keep the same deterministic operations
and fresh reset. This still fails to reconstruct presentations.

Proof: take the seven-vertex pair and add w-1 singleton vertices to each
side with operations tau J_L tau and tau J_R tau. Its order is2w+5. The two
explicit inherited maps extend by fixing every added vertex. They also fix
the original endpoints. Therefore both maps fix every queried boundary
position. After any further deterministic word, the old seven vertices
remain in an interval; their map is the original or complemented map,
according to total tau parity. All exterior vertices are fixed, and the
queried boundary positions never penetrate beyond an original endpoint.
The blue states remain dead by M1. At the first fresh reset, q is equal,
so the entire state law agrees and all subsequent observations agree.

The buffered presentations are distinct and strong. Each left insertion
has new -> old first -> ... -> old last -> new as a directed spanning cycle;
right insertion is dual. Thus the obstruction persists inside strong atoms.

In particular, the one-successor enhancement that reconstructs X_7 already
has a strong counterexample at order9 (w=2). **FINITE-EXACT:** w=1,...,6 and
all85 deterministic context patterns through three further insertions check
the literal maps, fixed boundary ports, strong initial states and dead blue
states. The preceding argument proves every width and every continuation.

This is a concrete alignment-memory lower bound. Increasing boundary detail
by a fixed amount can postpone an invisible discrepancy, but cannot expose
every interior alignment. It does not rule out an adaptive moving selector,
deletion operations, or named legality guards; those are different languages.
