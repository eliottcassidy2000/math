# Bounded root degree gives an exact finite weighted carrier

**PROVED:** For every fixed integer d≥1, adaptive insertion restricted to
roots of current indegree at most d has an exact finite nonnegative integer
word automaton for its Hamiltonian-path count. The alphabet has at most
2d letters after the first insertion. This is a constructive structural
representation; its numerical image is not classified.

**FINITE-EXACT:** After fixing the first triangle root to be 0, the reachable
degree-control graphs for d=1,2,3 have 1,6,24 vertices. The endpoint carrier
constructed below has respectively 15,460,9,672 reachable
(control, endpoint-signature) coordinates. These are explicit carrier
sizes, not ranks or minimal dimensions.

The mechanism is the endpoint-pairing transfer of
[Session 124, bounded-backedge tournaments](../session124_bounded_width_tournaments_20261031/README.md),
with a temporal set of live vertices in place of a fixed-width final
window, and with the forgotten-start/end convention reversed. Both
component pairing and multiplicity are retained. This specialization
closes under actual adaptive root choices and recovers a bounded sidecar
after the unrestricted four-coordinate closure failed in
[ADAPTIVE.md](ADAPTIVE.md).

## 1. A bounded temporal boundary

Start with C3 and add a new x at a chosen root v by v→x and x→all old
vertices other than v. Let A_d(T) be the set of vertices of current
indegree at most d. Under the cap-d rule, v∈A_d(T), and exactly

    A_d(T+x) = {x,v} ∪ {w≠v: indeg_T(w)≤d−1}.       (1.1)

Every old vertex except v gains one incoming arc; v retains its indegree
and x has indegree one. Once a vertex leaves A_d it is never a future
root, its indegree keeps increasing, and every future vertex dominates it.

For any strong tournament on n>2a+1 vertices, at most 2a vertices have
indegree at most a. Indeed, a set of m such vertices has at least
m(m−1)/2 incoming incidences internally, so m≤2a+1. Equality forces all
those indegrees to equal a internally and forbids incoming arcs from
outside, contradicting strongness. Apply this with a=d−1 to the old
tournament in (1.1). If n>2d−1, its right-hand side has size at most
2+(2d−2)=2d; if n≤2d−1, the whole child has size at most 2d. Therefore

    |A_d(T)|≤2d after every insertion.               (1.2)

For d=1 the initial triangle itself has three active vertices. Fixing its
first root to 0 gives the usual four-vertex child and thereafter (1.2)
applies. C3 symmetry makes this normalization harmless for counts and
ordinary generated classes. All resulting tournaments remain strong.

## 2. Finite control needs only the degrees

Order the active ports by their original birth labels. The control is the
tuple a=(a_0,...,a_{r−1}) of their current indegrees, each in {1,...,d}.
A letter j chooses the active port of rank j as root. Its deterministic
update is:

1. Keep a_j and increase every other a_i by one.
2. Delete entries exceeding d, retaining the order of the remaining ports.
3. Append a new final entry 1.

No old active arc is needed by this control: the next cut is always the
single old incoming root and all other old outgoing neighbors. Old arcs
are represented in the path-cover weights, not discarded from the count.
The finite degree space has at most Σ_{r=1}^{2d} d^r elements. A word in
active ranks reconstructs an actual named graph by applying these rules
to a list of surviving birth labels; it never asserts realization of a
formal numerical profile.

Starting from the root-0 four-vertex child, d=1 has the sole control (1,1).
For d=2 the complete set is

    (1,1), (1,2,1), (1,2,2,1),
    (2,1,1), (2,2,1), (2,2,2,1).

The 24 cap-three controls and every labelled transition for d=1,2,3 are
frozen in `active_port_transfer.json`. They are obtained by closure of the
explicit finite transition graph, not by stopping at a selected order.

## 3. The path-cover state and exact transitions

At a processed prefix, take an unordered spanning cover by nonempty
directed paths, allowing singleton paths. Record each component's paired
(start,end). An active endpoint retains its port rank; every inactive
endpoint is replaced by the symbol *. Keep the unordered multiset of
these pairs and the number of actual covers having that signature.

A forgotten end cannot acquire an outgoing arc to any future vertex.
Consequently two components with forgotten ends can never be joined into
one final path. Discard precisely these signatures. Each surviving
component has a distinct active end or the unique forgotten end, so
there are at most r+1 components when r ports are active.

When the new vertex x is added at root j, choose at most one old component
to enter x and at most one to leave x:

* An incoming component is allowed exactly when its end is port j.
* An outgoing component is allowed when its start is any other active
  port, or is forgotten; x dominates all those vertices.
* Choosing the same component for both roles creates a directed cycle
  and is forbidden.

Join the selected components through x, keeping all unselected ones.
The new component starts at the incoming component's start, or at x if
there is none; it ends at the outgoing component's end, or at x if there
is none. Apply the degree control's port deletion and renumbering,
replace newly forgotten endpoints by *, and discard two-forgotten-end
signatures. Different choices producing the same signature have their
multiplicities added. A one-component signature is accepted.

This transition is linear over the nonnegative integers and depends only
on the control, the endpoint signature, and the chosen root rank.
The whole finite (control, signature) set therefore supplies one fixed
nonnegative integer matrix M_j for each letter j∈{0,...,2d−1}; illegal
letters have no outgoing transition. With the actual four-vertex cover
vector u_d and the one-component indicator f_d,

    H(T_word) = u_d M_{j_1} ... M_{j_m} f_d.          (3.1)

Use the empty suffix for the four-vertex seed, and adjoin the triangle's
count 3 as the initial base case.

**Why the correspondence is exact.** A final Hamiltonian path selects
each used arc exactly when its later endpoint is born. At every prefix,
its selected arcs form a path cover. They cannot have two forgotten ends,
because neither end could ever obtain its required outgoing arc. Thus
the path determines a unique accepted transfer trajectory. Conversely,
each accepted trajectory builds a spanning acyclic directed graph with
indegree and outdegree at most one; the final one-component condition
makes it a Hamiltonian path. Retaining component pairing is what excludes
spurious cycle joins. Uniform future adjacency to forgotten vertices is
what makes their individual names unnecessary for subsequent counts.

The matrices retain unbounded integer weights. This is a finite weighted
carrier, not a finite set of numerical states or a finite-state decision
procedure for membership of arbitrary target integers.

## 4. Size bounds and exact small carriers

Here is an explicit signature bound independent of graph reachability.
Let g_r count abstract signatures on r available port labels with no
forgotten end, and f_r count those with at most one. Include the empty
signature in these upper bounds. Then

    g_0=1, f_0=2,
    g_r=3g_{r−1}+2(r−1)g_{r−2},
    f_r=3f_{r−1}+g_{r−1}+2(r−1)f_{r−2},             (4.1)

with the last terms zero at r=1. To derive this, classify the least port
as absent from the endpoints, a singleton component, the end of a
forgotten-start component, the start of the unique forgotten-end
component (the extra g term), or one end of a component paired with any
other active port (two orientations). A component with both endpoints
forgotten is the second f_0 possibility.

The f_r values for r=0,...,6 are

    2,7,28,123,582,2931,15576.

Thus Σ_{r=1}^{2d}d^r f_r is one explicit finite dimension bound before
discarding unrealizable controls and signatures. These are not the
Session 124 fixed-window reachable dimensions: arbitrary temporal
forgetting requires its own carrier audit.

Closing the actual positive transfer support from the four-vertex seed
gives:

| Degree cap d | Reachable degree controls | Reachable (control,signature) coordinates |
|---:|---:|---:|
| 1 | 1 | 15 |
| 2 | 6 | 460 |
| 3 | 24 | 9,672 |

Every counted coordinate is reachable under some legal word: positivity
prevents cancellation, and each transition is an actual cover extension.
The closure computation supplies an explicit exact automaton of these
dimensions. No assertion of linear independence, minimum realization
dimension, or minimum alphabet is made.

Indeed the cap-one carrier compresses to the exact positive three-coordinate
R/N representation in [ADAPTIVE §2](ADAPTIVE.md), whose rational dimension3
is sharp by an actual Hankel minor of determinant-4. This is a concrete
example of the distinction between reachable forest-coordinate count and
minimal linear dimension. No corresponding rank claim is made for caps2,3.

## 5. Positive target interface and the sharp cap-two obstruction

The parent lane's `root_degree_probe.py` exhausts actual named histories
under a count cutoff, using strict growth to make pruning complete. Its
finite conclusions are:

* Cap two reaches every allowed odd count through 99, but misses 127
  among allowed odd counts through 199. This is permanent absence from
  that language, since counts never decrease.
* Cap three reaches every allowed odd count through 199. No assertion
  beyond that cutoff follows.

Here "allowed" omits 7 and 21 and starts at 3. The inherited actual
history (0,0,3,3,0) has chosen-root indegrees (1,1,2,2,3) and final rows

    (154,4,1,102,14,23,55,126),    H=127.

Thus the minimum root-degree cap for this prime count in the adaptive
language is exactly three. The positive bounded-port conjecture must
move past cap two; the cap-two theorem remains an exact carrier theorem.

For cap three, (3.1) makes a precise constructive research problem:
find, for every target prime or 7p in the desired lane, a legal word in
six active-rank letters whose actual automaton output is that target.
Any matrix identity proposed for an induction must preserve the degree
control and act on reachable weighted states. The 24-control graph and
9,672-coordinate positive representation give a fixed platform for that
question; they do not answer it.

## 6. Temporal boundedness is different from static feedback boundedness

The inherited
[feedback-gap automaton](../feedback_word_hamilton_spectrum_20261020/proofs.md#5-exact-finite-state-grammar-at-arbitrary-fixed-feedback-vertex-number)
fixes a feedback vertex set whose deletion leaves a transitive tournament.
Its finite carrier is a different restriction.

After the first root-0 insertion, always choose the newborn vertex as the
next root. Every chosen root has indegree one. Consecutive triples of new
births (3,4,5), (6,7,8), ... induce directed triangles: the consecutive
arcs point forward, and the last vertex points back to the first. These
triangles and the initial (0,1,2) are vertex-disjoint. At order n this
gives floor(n/3) vertex-disjoint triangles, so both feedback vertex number
and feedback arc number are at least floor(n/3).

Hence even the cap-one temporal language has unbounded static feedback
cost. No fixed feedback-vertex bound can contain this language. The
positive bridge is to the endpoint-pairing transfer, not an unsupported
identification with a fixed-feedback automaton.

## 7. Reproduction and audit boundary

`active_port_transfer.py` imports no repository code. It stores the
complete degree graphs, closes reachable endpoint supports for d=1,2,3,
and compares the transfer against direct Hamiltonian subset DP on 15,61,25
actual named graphs in the respective small controls. It also compares
the entire endpoint-weight vector with an independent literal
permutation-and-cuts enumeration on 7,19,25 graphs respectively; covers
are deduplicated as unordered collections of actual directed paths.
These tests include forgotten components and multiplicity, not only H.
The 127 witness and a 30-vertex cap-one graph with ten disjoint triangles
are checked separately. No target-support census is added by this script.

```sh
python3 -I -B 05-knowledge/results/adaptive_metagraph_reformulation_20261103/active_port_transfer.py
python3 -I -B -O 05-knowledge/results/adaptive_metagraph_reformulation_20261103/active_port_transfer.py
```

Both runs agree with `active_port_transfer.json`. The all-order proof is
the transfer bijection and bounded-port argument; the finite graph closure
computes the displayed carrier sizes, while the literal checks audit its
implementation. Rank minimization and numerical-image induction remain
open work.
