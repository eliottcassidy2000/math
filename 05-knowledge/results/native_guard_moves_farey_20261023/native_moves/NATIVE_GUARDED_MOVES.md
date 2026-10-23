# Actual guarded observer moves and balanced selected-port rectangles

Status: **PROVED** for the local transposition and selected-rectangle algebra;
**FINITE-EXACT** for the complete Hamiltonian fibres of the four stated
original order-eight sources. The original electrical potential and native
owner atlas are inherited and authenticated by byte hashes. The coupled
cell-phase theorem in [`../moves/COUPLED_PHASE_MOVES.md`](../moves/COUPLED_PHASE_MOVES.md)
has a different move cost and must not be substituted for this one.

## 1. The actual finite carrier

Use the original merged native source nodes `1238,2242,2246,2315`, with
ordinary source classes `1763,3481,3490,3676`. Their complete Hamiltonian
path counts are respectively `409,407,523,527`; all four automorphism
groups are trivial. The potential `x` is the original frozen integral lift
from `metagraph_order8_electrical_20261023/receipt.json`, and

\[
               \phi_x(T,P)=x_{[T]}-x_{[K_P(T)]}\pmod8.
\]

The original order-eight word atlas uses **forward** bits on the
nonconsecutive path pairs. We use that exact gauge and its original owner
map. Each candidate is an actual permutation of the eight source owners
whose seven successive arcs are checked directly.

A move replaces `P` by another actual Hamiltonian path `Q` of the same
source, with the same specified indexed pins and the same frozen `phi`.
Its cost is

\[
                 d(P,Q)=|\{j:P_j\ne Q_j\}|.
\]

This is the number of changed **actual owner clocks**. It is not the
number of cell phases, the number of free-word bits, or the number of
elementary metagraph edges. The operation is one simultaneous reordering;
there is no claim that its individual owner changes can be applied through
accepted intermediate permutations.

## 2. The strong-four core gives a precise primitive move

Let `P` be any actual Hamiltonian path in a tournament and exchange the
owners at clocks `i<j`, calling them `p,s`. If the result is Hamiltonian,
then `j-i>=3`. Adjacent owners cannot be exchanged because that reverses
their required arc. For span two, the intervening owner would have to
both beat and lose to either exchanged owner.

For span at least three, the exact legality condition is:

* each existing predecessor at clocks `i-1,j-1` beats **both** `p,s`;
* both `p,s` beat each existing successor at clocks `i+1,j+1`.

Omit nonexistent endpoint predecessors or successors. All other adjacent
path arcs are unchanged, so this is necessary and sufficient. It is a
local native-arc test, independent of an electrical claim.

At span three, write the four-clock block as `(p,q,r,s)`. The two paths
require

\[
                 \{p,s\}\to q\to r\to\{p,s\}.
\]

The remaining arc `p-s` can point either way; the induced tournament is
the strong four-vertex class. Thus a genuine small tournament mechanism
recurs inside a larger source: two outer owners can trade positions
around two fixed middle owners, provided the external boundary test also
passes. An electrical phase still requires its original global target
receipt; strong-four structure alone does not preserve it.

## 3. A phase-preserving move that changes the actual target class

In original source `1238`, take

\[
\begin{aligned}
 P&=(7,1,0,4,2,6,5,3),\\
 Q&=(7,1,0,5,2,6,4,3).
\end{aligned}
\]

Both are actual Hamiltonian paths. They retain endpoint owners `7,3`, the
actual indexed arc at clock four `2->6`, and the original electrical phase
`7`. Their target merged classes are **different**, `881` and `2762`, and
their exact forward words are `1192759` and `1257847`.

This is a cost-two swap at clocks three and six. Its strong-four core is
`(4,2,6,5) <-> (5,2,6,4)`, and the external boundary owners are `0,3`.
The frozen values are

\[
 x_{1238}=97\equiv1,\qquad x_{881}=170\equiv2,
 \qquad x_{2762}=226\equiv2\pmod8.
\]

Thus both phases equal `1-2=7 mod8`, using the actual potential rather
than a surrogate statistic. Exhaustive enumeration shows that these are
the only two paths with those endpoints, phase and indexed middle arc.
Their exact minimum connection budget is two. Since the source has no
nontrivial automorphism, this move does not merely traverse an invisible
automorphism lift, and the changed target classes independently confirm
that fact.

## 4. Exact coverage and the remaining large-block obstruction

The checker enumerates all **1,866** Hamiltonian paths in the four sources.
It partitions them by `(first owner,last owner,phi)` and then additionally
by every one of the five internal indexed arc clocks. All candidate
simultaneous reorderings inside each resulting fibre are included.

There are `1,014` endpoint/phase fibres, of which `485` are nontrivial, and
`562` nontrivial refinements with one additional indexed internal arc.
For each fibre, a minimum-bottleneck spanning tree supplies a constructive
upper bound, while the disconnected graph at one lower budget supplies
the exact matching lower bound. This is an exact finite criterion: if
edges are allowed precisely when their clock distance is at most `b`, the
least connected threshold is the largest edge cost of a minimum spanning
tree. Kruskal's cut argument proves necessity and sufficiency.

| Source | All Hamiltonian paths | Endpoint/phase fibres needing budget6 | Additional-arc fibres needing budget4 |
|---:|---:|---:|---:|
|1238|409|20|37|
|2242|407|19|27|
|2246|523|35|39|
|2315|527|31|44|

The existence of budget-six fibres is not just the trivial upper bound
from two pinned endpoints. For example, source1238 with endpoints `7,3`
and phase7 has exactly these four accepted paths:

```text
7 1 0 4 2 6 5 3
7 1 0 5 2 6 4 3
7 4 1 0 5 2 6 3
7 5 4 2 1 0 6 3
```

The budget-five graph has exactly two components: the first two paths and
the last two. Every cross-component pair changes all six unpinned owner
clocks. A budget-six simultaneous move is therefore necessary. The
receipt retains the full path list, all target owners, the lower-budget
components and the spanning tree for every tested fibre. No statement
about arbitrary order-eight sources or larger orders is inferred.

## 5. Two actual sources can balance a selected-port switch

The sources `S=1238` and `R=2246` have frozen values `97` and `65`, equal
modulo eight. Search only their actual endpoint/phase fibres. There are
14 common unordered target pairs admitting a move at both sources; two
remain when each source must also keep an indexed internal arc. The two
survivors are `(881,2762)` and `(2650,2843)`.

For the first target pair, use the source1238 paths from section3. At
source2246 use the reverse-direction change

```text
6 5 4 3 2 0 7 1   target2762
6 5 2 4 3 0 7 1   target881
```

It has cost three, retains endpoint owners `6,1`, and retains the actual
arc at clock five `0->7`. Every one of these four actual paths has original
phase7. Selecting one path at each source therefore changes the selected
incidences

\[
          (S,881),(R,2762)\quad\longrightarrow\quad
          (S,2762),(R,881).
\]

The selected source counts and target counts are both preserved. Include
each path's actual reversed converse partner to retain the native paired
black convention. All four underlying edges are pure black in the
original atlas, with black multiplicities `4,2,4,4` on
`(S,881),(S,2762),(R,881),(R,2762)`, respectively; their blue multiplicities
are zero. A selected converse pair contributes black weight two.

This changes the selected incidence pattern, not merely the order of two
identical slots. Every selected observer before and after is actually
realized in its fixed source with its stated pins. The full canonical
all-path metagraph is not thereby replaced or re-enumerated.

## 6. Exact Kirchhoff calculation, colours and the unpaid higher layer

For four distinct vertices `S,R,U,V`, replace weight-`w` edges `(S,U)` and
`(R,V)` by `(S,V)` and `(R,U)`. The Laplacian change applied to a frozen
potential is, in the order `(S,R,U,V)`,

\[
 \Delta Lx=w(x_U-x_V,\ x_V-x_U,\ x_S-x_R,\ x_R-x_S).          \tag{1}
\]

This follows by writing each edge's two opposite currents. Consequently
the unchanged frozen potential remains harmonic modulo `m` iff both
`w(x_U-x_V)` and `w(x_S-x_R)` vanish modulo `m`. Individual selected-port
phase preservation requires the target congruence itself, without
silently cancelling the factor `w`.

The two authenticated rectangles have the following exact debts:

| Targets `(U,V)` | Frozen target values | Unit switch defect | Converse-paired black defect (`w=2`) | Highest paired dyadic modulus retained |
|---|---|---|---|---:|
|881,2762|170,226|`(-56,56,32,-32)`|`(-112,112,64,-64)`|16|
|2650,2843|150,86|`(64,-64,32,-32)`|`(128,-128,64,-64)`|64|

The highest unit-weight moduli are eight and32, respectively. Both
rectangles retain the original requested mod-eight phase and aggregate
Kirchhoff condition. **Neither retains the full frozen modulus256.**
The second fails already modulo128. A paired aggregate equation can also
retain a higher modulus than an individual port phase; these tests must
remain separate.

For the second pair, all four edges are also pure black, with
multiplicities `4,2,4,6` on `(S,2650),(S,2843),(R,2650),(R,2843)`.
Both source moves have clock cost three and keep actual indexed arcs:

```text
source1238, fixed endpoints7/2 and arcclock1=(0,6):
7 0 6 4 3 1 5 2   target2650
7 0 6 4 1 5 3 2   target2843

source2246, fixed endpoints6/4 and arcclock1=(0,7):
6 0 7 3 2 1 5 4   target2843
6 0 7 3 1 5 2 4   target2650
```

All four phases are three modulo eight. Each selected-source move retains
its own actual ports. Adding the converse partners preserves the paired
black type. The required removals have available black weight at least
two, and additions retain nonnegative even black multiplicities in the
selected-shadow calculation.

**The next amplitude is blocked by an actual port, not by the aggregate
weight.** Doubling the second rectangle would give defect
`(256,-256,128,-128)`, passing modulus128. Its chosen aggregate removals
would fit the black weights four and six. But source1238 has only **one**
actual ordinary-source path to target2843, namely
`(7,0,6,4,1,5,3,2)`. The entire exhaustive 409-path list contains no second
copy. Any amplitude-two selection would duplicate that path, irrespective
of additional pins. The first rectangle likewise has only one path from
source1238 to target2762. Each rectangle therefore has exact maximum
distinct-observer amplitude one, in either orientation.

This is a concrete loss in the aggregate capacity test: nonnegative
modified edge weights do not create additional distinct native ports.
The two rectangles have disjoint target pairs. Each paired switch of
amplitude `t` contributes `64t` at its own target coordinates, so a combined
modulus256 repair would require both amplitudes divisible by four. Under
the proved literal bounds `|t|<=1`, only the zero combination works. This
closes this particular two-rectangle repair attempt without claiming that
other source/target blocks are exhausted.

**Scope boundary.** This is a balanced move of selected actual observers
and a verified electrical rectangle. It is not yet a global ownership
bijection reassigning every source port, nor a proof that the altered full
aggregate adjacency is a realizable native metagraph. Other ports,
capacity constraints, higher electrical layers and line-identity records
are not automatically transported. Those are the next obligations, not
consequences of a vanished four-entry defect.

## 7. Reproduction and loss ledger

The map retains the actual source, every full observer, its indexed pins,
the original forward word, original target class, frozen potential and
clock cost. Projection to `(endpoints,phi)` loses the internal path; the
receipt restores it by retaining all full permutations. Converse reverses
actual paths and transports their clocks, and every path in the finite
census is checked in the converse source as well.

```text
python -B 05-knowledge/results/native_guard_moves_farey_20261023/native_moves/check.py
python -O -B 05-knowledge/results/native_guard_moves_farey_20261023/native_moves/check.py
```

The four-source census also exhaustively checks the local transposition
criterion on every pair of owner clocks in every actual path, including
the genuine strong-four type at span three. Original file hashes and the
paired-rectangle colour/current records are retained in `receipt.json`.
No floating spectral test or optimization-sensitive assertion is used.
