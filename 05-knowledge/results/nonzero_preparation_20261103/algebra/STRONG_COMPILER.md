# Three temporary vertices suffice for universal all-strong native transport

**PROVED — [independent proof and exact replay audit](STRONG_AUDIT.md) passed.** The theorem below
upgrades the inherited stabilized single-arc compiler by preserving strongness
at every intermediate. It controls actual named adjacency, not Hamilton counts
along the route. It does not construct a target of a prescribed H.
The move language is the inherited stabilized native language; no restriction
to a prescribed move colour or to fixed order is imposed.

## 1. Statement, inheritance, and the hostile boundary

Let T and U be any strong tournaments on the same n>=3 physical vertices.
Let d be the number of unordered old pairs whose arcs differ. A native
Hamilton complement K_P reverses every arc except the consecutive arcs of
an actual directed Hamilton word P of the current graph.

**Theorem.** There is an explicit finite program from T to U
using vertex births, legally guarded K_P moves, and deletion of only its
temporary vertices, such that:

* every intermediate tournament is strong;
* every old physical name, including arbitrarily many marked ports, survives;
* at most three temporary vertices coexist;
* the program uses2d+4 native complements, d+4 temporary births and d+4
  temporary deletions, with peak order n+3.

These are bounds, not minimality assertions; when T=U the empty program
also suffices. The final old matrix is exactly U. H, F, and other responses
are not asserted constant at intermediate stages or equal at the endpoints
unless that equality is independently known for T and U.

The inherited supplier is the arbitrary-tournament one-helper compiler in
[METAGRAPH §5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md#5-one-temporary-vertex-compiles-any-old-arc-reversal),
with its direct pair-block guard in
[NATIVE_SIX §4](../../multiport_grammar_20261103/component/NATIVE_SIX.md#4-the-inherited-legal-compiler-instantiated-on-every-discrepancy).
It gives actual adjacency transport but no universal strongness guarantee.
The new ingredient is a strong scaffold protecting even a temporarily weak
old subgraph, together with guard endpoints making the middle native state
Hamiltonian. The canonical hostile is C3: a native complement at any of its
Hamilton paths is transitive. Strong input alone does not make an arbitrary
native complement strong.

The live concepts are physical names, native legality, strong endpoints,
temporary order, a protected cyclic scaffold, and the scalar target fibre.
The map forgets no old relation: it implements each requested discrepancy
literally. Its additional sidecar is a Hamilton guard with specified outer
vertices; that is what pays intermediate strongness.

## 2. A guarded pair reversal with a strong middle

Suppose P=(f,...,u,v,...,l) is a Hamilton word of a named tournament X,
with v not last and f!=u. Birth z whose only outgoing old arc is z→u.
Then Pz is legal: l!=u, so l→z. Apply K_(Pz). The word Q obtained by
inserting z between u,v in P is now legal:

    u→z and z→v

are the reversed old incidences at z, and every other consecutive old arc
of Q was protected by Pz. Apply K_Q and delete z. The old protected-edge
sets differ exactly at {u,v}, so the final old matrix is X with that pair
reversed. The helper's outgoing old set after the two moves is exactly
{v,l}: the helper incidences at u,v,l change once and all others change
twice. These vertices are distinct in this case.

The middle state has the directed Hamilton cycle

    f,...,u,v,...,l,z,f.                              (1)

Indeed f→z was an unprotected chord at the first move, because f differs
from u and l; it becomes z→f. Thus the middle state is strong, regardless
of any other chords. If both old endpoint graphs X and X^{uv} are strong,
the helper extensions at either end are also strong: the cuts {u} and
{v,l} are nonempty and proper in every application below. Deleting z leaves
the stipulated strong endpoint.

The same two involutions may be traversed backward. In that direction,
birth z with outgoing old set {v,l}, use Q and then Pz, and delete z.
It has exactly the same strong middle (1). No reversal is treated as a
primitive operation.

## 3. Build the scaffold without a weak source/sink intermediate

Choose any old vertex w. Birth a beaten only by w and beating every other
old vertex. Since old T is strong and n>=3, this is a strong extension.
Birth b beaten by every old vertex and beating a. This too is strong:
the previous strong graph has an arc into b, and b has its arc to a.
At this stage a loses to w and b and beats every old vertex except w.

We need to reverse only w→a. Choose an old Hamilton path R beginning at w.
It exists for every w of a strong tournament; the endpoint-preserving path
argument is recovered in
[construction/RESULT §4](../../multiport_defect_grammar_20261103/construction/RESULT.md#4-why-the-natural-auxiliary-choices-do-not-rearm-the-digit),
and can also be read in converse. A finite search for such a path is already
an effective construction; no unproved complexity bound is needed here.

In the desired graph with a→w, the word

    P=(b,a,R)

is Hamilton. Its selected pair is (a,w), first vertex b differs from a,
and w is not last because R contains all n>=3 old vertices. Apply the
inverse of §2's gadget. Its initial helper cut is {w,last(R)}, a nonempty
proper cut. The two old endpoint graphs are strong: the current one by the
two births just proved, and the desired one by the cyclic scaffold below.
The middle is strong by (1), so opening the scaffold never passes through
a universal source or sink.

We have now reached the exact relation

    a → every old vertex → b → a.                    (2)

For any nonempty old tournament, even a weak one, (2) is strong. Every old
vertex reaches b,a and then every other old vertex. This is the reason the
old subgraph may be edited without keeping its own strongness.

## 4. Reverse all old discrepancies inside the scaffold

For each differing old pair {u,v}, find a pair-block Hamilton word R on
the current old vertices, legal in one of the two orientations of that
pair. The inherited direct algorithm supplies it: take a Hamilton word on
the other old vertices; insert the pair before the first word vertex beaten
by either endpoint, choosing its last endpoint to beat that vertex. Every
earlier word vertex beats both endpoints, so the preceding join is legal.
If no such vertex exists, append the pair. Choose the internal orientation
to match the resulting word. The only potentially altered old arc is the
selected pair.

In that orientation extend the guard to

    P=(a,R,b).

Its selected pair is interior: the first vertex a is not u, and v is not
the last vertex b. Both old endpoint graphs, now including scaffold a,b,
are strong by (2). Section2 therefore supplies a strong forward or inverse
pair gadget, with its own temporary z. It changes no scaffold incidence.
Delete z before proceeding to the next discrepancy.

After d repetitions, the induced graph on the original physical names is
exactly U, while (2) remains intact. The number of old discrepancies is
computed from the actual matrices, not their numerical observers.

## 5. Remove the scaffold while retaining strongness

Now U is strong. Choose an old Hamilton word R_U beginning at the same w.
In the current scaffold graph, P=(b,a,R_U) is legal. Use the **forward**
gadget of §2 at (a,w), reversing a→w back to w→a. Both endpoint graphs
are strong, and (1) pays the middle. After its helper is deleted, delete b.
The remaining graph is U plus a, with w→a and a beating all other old
vertices; this is strong. Finally delete a, leaving the strong target U.

The order of these two final deletions matters. Deleting b before restoring
w→a would leave the universal source a, which is weak. That tempting shortcut
is the canonical hostile to a naive two-vertex scaffold argument.

There are d interior pair gadgets and two opening/closing gadgets, each
with two native moves and one temporary birth/deletion. The two scaffold
vertices add two more births/deletions. This gives2d+4 native moves and
d+4 temporary births/deletions. At most a,b,z coexist, proving the peak n+3.
Every deletion names a temporary vertex explicitly; no old name is renamed
away, contracted, or discarded.

## 6. Consequences and exact remaining obligation

If a strong target U has m>=|T| vertices and includes positions for every
old name/root, first make m−|T| permanent strong births, each with sole exit
at the retained old root. Apply the theorem at order m. Thus **strong actual
targets of sufficient order are always natively reachable with every
intermediate strong**, retaining all old names. The operation itself no
longer needs a separate per-target strongness conjecture.

For H-neutral preparation with strong endpoints, the remaining condition
is exactly an actual strong ready witness of the required H at order at
least the number of retained names. Necessity is immediate; sufficiency
uses these permanent births and the all-strong compiler. This refines the
weak-allowed reduction in [PREPARATION_REDUCTION.md](PREPARATION_REDUCTION.md):
transitive padding is no substitute for a sufficiently large strong witness.
The [retained-order obstruction](../surgery/ORDER_OBSTRUCTION.md) therefore
remains unchanged and decisive for the strong H-neutral+C2 (+2,+4) route.

More generally one may target an already-prepared child of H=3h+2 or3h+4
directly. Such a move is not H-neutral preparation followed by C2, so it is
outside that obstruction. The [direct finite recipes](../surgery/DIRECT.md)
use this different scheduling idea with smaller explicit budgets. Neither
those recipes nor this universal reachability theorem prove that a suitable
target exists for every subsequent integer. Global spectrum coverage and
unbounded two-branch target selection remain OPEN.
In particular this variable-order program does not settle fixed-order
black-core connectivity.

The theorem transports arbitrarily many marked ports by keeping the full
physical relation and every guard. No bounded response quotient is asserted
to determine these moves. Full named adjacency permits exact recounting of
all requested observers after each command. The inherited coloured endpoint
forests remain available for the operations in their separately proved scopes.

## 7. Audit and controls

The proof is analytic for all n and d. The standalone verifier
[`strong_compiler.py`](../../../../04-computation/nonzero_preparation_20261103/strong_compiler.py)
imports no prior producer or numerical receipt. It implements the direct
pair-block compiler and independently interprets each complement as an
adjacency-set operation. Every guard, literal middle Hamilton cycle, named
endpoint, helper-only deletion and strong state is checked.

The declared finite universes are10,648 arbitrary old-core/pair cases through
order5 (including weak old cores inside the scaffold),2,822 setup/cleanup
cases at every root of every strong labelled core through5, and2,316 complete
strong source/target/root cases through4. Both cleanup mistakes are hostile
controls: removing b before restoring w→a, or deleting a before b afterward,
creates a weak graph. A complete H5-to-H5 trace retains every guard, helper
cut, matrix and changing root response; its intermediate H is not constant.

Normal and optimized Python modes and the frozen receipt agree at
**424,665 always-active gates**, covering27,860 pair gadgets, of which13,770
use the inverse orientation. A separate coordinator replay matches the same
receipt byte for byte. These controls test the analytic proof; they are not
an extrapolation premise.

```text
python3 -B 04-computation/nonzero_preparation_20261103/strong_compiler.py --output /tmp/strong_compiler.normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/strong_compiler.py --output /tmp/strong_compiler.optimized.json
```

Source:11,052 bytes, SHA256
`71c64d0caf21f22d454a7589dc2e05e87d66ab6a52378e13700e16ed0c362326`.
[Receipt](../../../../04-computation/nonzero_preparation_20261103/strong_compiler.json):
10,168 bytes, SHA256
`b914205ecba7e9d805ece970b1bda3dc2ee5a67a059054f2add1d0ccf7803fff`.
The independent [audit](STRONG_AUDIT.md) records the proof's type, operation,
strongness, quantifier, and target-selection boundaries. No optimal helper
bound, fixed-order connectivity result, or Lean formalization is claimed.
