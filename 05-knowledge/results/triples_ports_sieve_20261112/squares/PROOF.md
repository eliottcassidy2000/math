# Theta ports and cycle-opening births in square-sum graphs

**Status. PROVED:** the elementary all-size multi-arm coefficient formula,
Johnson-product exchange classification, endpoint-relocation criterion,
directed cycle-opening bijection and infinite arithmetic family below.
**FINITE-EXACT:** the complete reduction of
all 2,600 fixed-endpoint instances in the first 25 square-sum graphs, two
remaining theta kernels, and the identical-kernel/different-answer witness at
26. **INHERITED:** the first25 Hamiltonian census, degree forcing, reflection
transport, and the need for a connectivity sidecar. **OPEN:** an all-size
square-graph recurrence closing on a bounded set of these states, or an
implication for Chan–Pak. No claim of external priority is made.

## 1. Recovery and declared objects

Throughout, `Q_n` is the finite **undirected simple** graph on `[1,n]` with
`{a,b}` an edge exactly when `a!=b` and `a+b` is an integer square. Paths are
spanning, Hamiltonian, and unoriented unless an order is explicitly used. A
fixed endpoint set contains two distinct vertices. `Q_1` has its singleton
path by convention. No tournament orientation is added to this relation.

Closest inherited mechanisms:

- [Session32, port transfer](../../session32_square_monster_20261023/ports/PORT_TRANSFER.md)
  has the first25 census, the Q12 blocks, degree-feasible assemblies versus
  disconnected cycles, and square-color reflections.
- [Session40, square-sum packet](../../session40_marked_ports_square_sum_20261023/square-sum/README.md)
  has all paths through32, the endpoint-conditioned forcing argument, the
  endpoint ambiguity at26, and the exact birth polynomial involving two-path
  forests. Its `check_square_incidence.py` already implements degree
  propagation and forced-cycle rejection to exclude edge1–15 at25.
- [Session40, matrix-tree projection](../../session40_marked_ports_square_sum_20261023/audit-first/MATRIX_TREE_PATH_PROJECTION.md)
  proves the all-size exact degree sector of the spanning-tree polynomial.
  Its connectivity warning is used here, not rediscovered as a new theorem.
- [Session257, fixed and adaptive charts](../../session257_permutohedron_squares_20261111/SYNTHESIS.md)
  contracts the fourteen edges shared by all ten Q25 paths; its fixed chart
  freezes under elementary block moves, while pair-adapted cuts retain
  larger surgery. The present kernel is **endpoint-conditioned**, so is a
  different contraction.
- [Square crossroads, sections4–4.2](../../crossroads_family_20260926_squares.md)
  already proves signless-incidence transport and the whole-graph rigidity:
  Q14 has one alternating offset; Q15 contains the odd cycle
  `1–3–6–10–15–1`; for every `n>=15`, designated square-edge labels determine
  all integer vertex values uniquely. Selected paths still have an
  alternating offset. This is the inherited square/Sierpiński reflection
  bridge, not a new rigidity claim.

The canonical hostile example is a degree-two vertex used as a path endpoint,
such as22 in Q23. The corrected near miss is endpoint-blind forcing. The
least-used sidecar developed here is the **single cycle accompanying the
endpoint path**: a failed Hamiltonian candidate can become a successful path
at the next birth.

Concept board: endpoint demands; square reflections; theta parity; outside
component partitions; cycle-opening births; exchange graphs. The first25
counts themselves are inherited. The new scoped contribution is their exact
two-kernel reduction, the matching infinite family of port responses and
exchange graphs, and the explicit birth that changes the connectivity filter.

## 2. Sound degree reduction

Fix endpoints `{s,t}` and put `d_s=d_t=1`, `d_v=2` otherwise. Each graph edge
is forced present, forced absent, or undecided. At a vertex with `f` forced
incident edges and `u` undecided incident edges:

1. `f>d_v` or `f+u<d_v` is a contradiction.
2. `f=d_v` deletes every undecided incident edge.
3. `f+u=d_v` forces every undecided incident edge.

A cycle of forced edges is a contradiction for a spanning path. Iterate to
a fixed point. Every action preserves every Hamiltonian path with the named
endpoints: the three degree rules follow from the exact degree equation, and
a path contains no cycle. This is a sound reduction, not a completeness
claim for arbitrary input graphs.

The reproducible convention scans vertices in increasing order and checks
for a forced cycle after each full sweep. This convention matters to the
reported **reason** for rejecting an already impossible endpoint sector.
It does not change the final Hamiltonian answer.

If no edge remains undecided and no contradiction occurred, the forced
support is a forest with the required degrees. It has `n-1` edges, hence is
connected and is the unique path. Otherwise let `F` be the forced forest,
`K` the undecided graph, and

`r_v=d_v-deg_F(v)`.

The remaining problem is to select a spanning `r`-degree support in `K` and
test connectivity after adjoining `F`.

## 3. All-size multi-arm response and its exchange graph

Take distinct hubs `u,v` joined by `r>=1` internally vertex-disjoint arms
with positive edge lengths `l_1,...,l_r`. Arms share no edges. In a simple
graph at most one arm has length1. Give every internal vertex residual
demand1, and give the hubs demands `a,b`. No other undecided edges are
incident with these vertices. Let `o` arms have odd length and `e=r-o` have
even length.

**Multi-arm theorem.** The exact number of degree completions is

```math
C_{a,b}=[U^aV^b](1+UV)^o(U+V)^e
       =\binom{o}{(a+b-e)/2}\binom{e}{(a-b+e)/2}.
```

A binomial is zero when its lower argument is nonintegral, negative or
larger than its upper argument. Thus this also gives a necessary and
sufficient feasibility test, including demands0 and all boundary cases.

**Proof.** Successive edge indicators along an arm must alternate because
each internal vertex requires exactly one of its two edges. An odd arm has
endpoint indicators `(0,0)` or `(1,1)`, giving `1+UV`. An even arm has
indicators `(1,0)` or `(0,1)`, giving `U+V`. Independent arm choices multiply.
If `i` odd arms start with1 and `j` even arms start with1, then
`a=i+j`, `b=i+e-j`. Solve these two equations and choose the corresponding
arms. Every choice determines all its edge indicators uniquely. QED.

The edge-weighted version retains the entire selected support. For arm `k`,
let `O_k` and `E_k` be the products of edge weights in odd and even positions,
with its first edge in an odd position. Replace the factors by

```math
E_k+O_kUV\quad(l_k\text{ odd}),\qquad
O_kU+E_kV\quad(l_k\text{ even}).
```

The coefficient is a positive, multilinear generating polynomial of exactly
the residual degree supports. Specializing weights to1 is the information
loss that forgets which internal edges were selected.

**Exchange theorem.** In a nonempty degree sector, set
`i=(a+b-e)/2`, `j=(a-b+e)/2`. Join configurations when one selected and one
unselected arm of the same parity exchange their first-edge bits. The
resulting graph is exactly

```math
J(o,i)\mathbin{\square}J(e,j).
```

Here `J(N,k)` has the `k`-subsets of `[N]` as vertices, with adjacency when
one element is replaced; `J(0,0)` is one vertex. The exchange graph is
connected, vertex-transitive, has

```math
\binom{o}{i}\binom{e}{j}\text{ vertices},\qquad
i(o-i)+j(e-j)\text{ neighbors per vertex},
```

and diameter `min(i,o-i)+min(j,e-j)`.

**Proof.** A configuration is exactly the pair of chosen arm subsets of
sizes `i,j`. An allowed move replaces one element in exactly one subset.
This is the stated Cartesian product. Distance is the sum of the numbers
of replacements needed in the two coordinates; its maximum gives the
diameter. Independent permutations of odd and even arms act transitively.
On the original graph a move toggles the entire even cycle formed by the
two exchanged arms. Alternation and opposite first-edge bits ensure every
residual degree remains unchanged. QED.

This is an all-size classification of a precisely declared family of
residual states. It is not a classification of all square-sum residual
graphs. Reversing the hub gauge exchanges `a,b`, complements the even-arm
first bits, and preserves the coefficient and exchange graph. Arm labels
are retained; no multiplicities are lost by an isomorphism quotient.

## 4. The connectivity sidecar

Let a completion `S` have all required residual degrees and let `F` be the
forced forest. Contract **every** component of `F`, including isolated
vertices. The selected edges `S` give a multigraph on the resulting blocks.
Keep loops and parallel edges.

**Exact filter.** `F union S` is Hamiltonian with the fixed endpoints if and
only if the contracted multigraph is a tree.

**Proof.** If `F` has `c` components then `|F|=n-c`; the completed degree
sum requires `|F union S|=n-1`, so `|S|=c-1`. Contraction preserves
connectedness. A connected multigraph on `c` vertices with `c-1` edges is
a tree. The resulting full graph has degree1 at the endpoints and degree2
elsewhere, so connectedness makes it a path. QED.

This is the contracted form of the inherited degree-plus-connectivity
criterion. The requisite outside information is the partition of named
residual vertices by forced component, including any component meeting no
residual vertex. The total number of forced components is retained. For
path reconstruction retain the actual forced edges as well.

An arm exchange always preserves degrees. It preserves Hamiltonicity
**exactly when** the new component quotient is still a tree. Thus the
Hamiltonian state graph is an induced subgraph of the Johnson product and
need not inherit its connectedness or vertex-transitivity.

## 5. Complete first25 reduction

The independent computation checks every unordered endpoint pair for every
`2<=n<=25`, including sectors immediately rejected by isolated vertices or
excess leaves. There are `sum binom(n,2)=binom(26,3)=2600` such sectors.
The table reports the fixed deterministic reduction from section2.

| Order n | Degree contradiction | Forced cycle | Unique forced path | Unresolved |
|---|---:|---:|---:|---:|
| 2 through14, separately | binom(n,2) | 0 | 0 | 0 |
| 15 | 104 | 0 | 1 | 0 |
| 16 | 119 | 0 | 1 | 0 |
| 17 | 135 | 0 | 1 | 0 |
| 18 | 153 | 0 | 0 | 0 |
| 19 | 170 | 1 | 0 | 0 |
| 20 | 177 | 13 | 0 | 0 |
| 21 | 197 | 13 | 0 | 0 |
| 22 | 230 | 1 | 0 | 0 |
| 23 | 249 | 1 | 3 | 0 |
| 24 | 272 | 4 | 0 | 0 |
| 25 | 283 | 7 | 8 | 2 |

Only two sectors remain, both at25:

| Endpoints | Hubs | Ordered arms | Lengths | Hub demands |
|---|---|---|---|---|
| `{10,18}` | `1,3` | `1–3`; `1–15–10–6–3`; `1–24–12–13–3` | `(1,4,4)` | `(1,1)` |
| `{18,22}` | `3,12` | `3–13–12`; `3–1–24–12`; `3–22–14–11–5–4–12` | `(2,3,6)` | `(2,2)` |

All non-hub residual vertices have demand1. Both kernels have `(o,e)=(1,2)`
and exactly two completions. The odd arm is unselected at both ends in the
first kernel and selected at both ends in the second. The exchange graph
in each case is `J(2,1)=K_2`; the move toggles the length8 cycle made from
the two even arms.

Connectivity accepts exactly one completion in each case. The rejected
completion at `{10,18}` is a path on19 vertices plus the six-cycle on
`{2,3,13,14,22,23}`. The rejected completion at `{18,22}` is a path on17
vertices plus the eight-cycle on `{2,11,12,13,14,23,24,25}`. Hence the eight
fully forced Q25 paths plus these two paths recover the inherited count10.

The first25 nonzero path counts are therefore exactly the inherited values
`H_1=H_15=H_16=H_17=1`, `H_23=3`, `H_25=10`. All paths and full forcing traces
for the two residual cases are in [receipt.json](receipt.json).

## 6. The same kernel changes its answer at26

At Q26 with endpoints `{13,18}`, degree forcing leaves **literally the same
labelled graph**, arm paths, and residual demands as Q25 with endpoints
`{10,18}`. It now accepts both completions. This checks equality of every
kernel edge and demand, not just equality of arm lengths or isomorphism.

The forced-component partitions on the residual vertices explain the change:

| Forced component boundary | Q25 endpoints10,18 | Q26 endpoints13,18 |
|---|---|---|
| first block | `{1,6}` | `{1,6}` |
| second block | `{3,13}` | `{3,10}` |
| third block | `{12,15}` | `{12,15}` |
| fourth block, containing18 | `{24}` | `{24}` |
| isolated endpoint | `{10}` | `{13}` |

All forced components occur in this table. Their full vertex lists and
edges are retained in the receipt. In particular the old forced path
`13–23–2–14–22–3` becomes `10–26–23–2–14–22–3`; the endpoint13 is isolated
until residual edges are selected. The first candidate at25 is already a
Hamiltonian path. The second closes the forced path by selecting edge3–13
and creates the six-cycle. At26 that same edge instead joins the isolated
endpoint13, and the cycle disappears.

Thus the counterexample to a kernel-only Hamiltonian count is exact:

```text
identical residual theta + identical demands:  1 Hamiltonian completion -> 2.
```

No scalar parity polynomial, Cartan type, spectrum or arm-count invariant of
this unmarked kernel can distinguish the two inputs. The outside partition
does distinguish them. This does not contradict the multi-arm theorem,
which counts **degree** completions: it correctly gives2 in both cases.

## 7. A uniform cycle-opening endpoint relocation

Let `G` be any finite simple graph and let `a,b,v,t` be four distinct vertices.
Take a spanning edge support `H` with degree1 at `b,t` and degree2 elsewhere,
containing edge `{a,v}`. Add a new vertex `z` with allowed edges `{v,z}` and
`{z,b}`. Define

```math
R(H)=\bigl(H\setminus\{\{a,v\}\}\bigr)
       \cup\{\{v,z\},\{z,b\}\}.
```

The new support has degree1 at `a,t` and degree2 elsewhere.

**Endpoint-relocation theorem.** `R(H)` is Hamiltonian if and only if exactly
one of these alternatives holds:

1. `H` is a Hamiltonian path and its order from `b` to `t` meets `a` before
   `v` (the selected edge is traversed `a→v`).
2. `H` consists of the path from `b` to `t` and exactly one cycle, and
   `{a,v}` belongs to that cycle.

**Proof.** Every graph with the specified degrees is the endpoint path plus
zero or more disjoint cycles. If `{a,v}` lies on a cycle, deleting it opens
that cycle to a path from `a` to `v`. The new two-edge ear joins its end `v`
to endpoint `b`. This produces one spanning path precisely when no other
cycle remains. If `{a,v}` lies on the endpoint path, any other cycle remains
and is fatal. With no other cycle, orient the old path from `b` to `t`.
If it encounters `a` then `v`, deletion gives paths `b...a` and `v...t`,
and the new ear joins them as `a...b–z–v...t`. If it encounters `v` then
`a`, the new ear closes `b...v–z–b` into a cycle, leaving the path `a...t`.
These alternatives exhaust all possibilities. QED.

The operation is reversible on its declared supports: remove `z` and add
edge`{a,v}`. Consequently there is an **exact count identity**, not merely
an injection. Hamiltonian paths of the enlarged graph with endpoints
`{a,t}`, with the two selected neighbors of `z` exactly `{b,v}`, and with
`{a,v}` absent, are in bijection with the disjoint union of the two source
classes in the theorem. The old vertices and all other selected edges are
unchanged. If the enlarged graph has additional edges incident with `z`,
this identity is conditioned on its selected neighbors; it does not
discard those edges from the graph itself.

At the actual birth25→26 use

```text
a=13, b=10, v=23, t=18, z=26.
```

The square identities are `13+23=10+26=36`, `23+26=49`. The old positive
candidate meets13 before23 in its10→18 order and survives. The old negative
candidate's six-cycle contains13–23 and is opened into a spanning path.
The receipt verifies the resulting two supports equal the **entire**
endpoint sector `{13,18}` at26 by an independent vertex-path enumeration.

This is an explicit mechanism by which a discarded cycle-defect state
becomes productive. It changes the named endpoint and increases order;
it is not yet a target-preserving return or a Chan–Pak construction.

## 8. Infinitely many genuine square-graph hosts

For every integer `m>=4`, set

```math
z=m^2-1,\quad b=1,\quad v=2m+2,\quad a=m^2-2m-2.
```

Then `a,b,v` are distinct members of `[1,z-1]` and

```math
a+v=b+z=m^2,\qquad v+z=(m+1)^2.
```

Thus the endpoint relocation is a genuine support operation from
`Q_(m^2-2)` to `Q_(m^2-1)`, for every admissible fourth endpoint `t` and
every source support in one of the two theorem classes. The arithmetic
identities follow by substitution. Positivity and `v<z` follow from
`m>=4`; `a<z` is immediate. Equality `a=v` would require
`(m-2)^2=8`, impossible for an integer. Also `a>b` for `m>=4`.

A larger family, containing the observed25→26 birth, is

```math
v=2m+2+r,\quad b=r+1,\quad
a=m^2-2m-2-r,\quad z=m^2-1-r,
```

where `m>=4`, `r>=0`, `2r<m^2-2m-3`, and `a!=v`.
These inequalities give `a>b>=1` and `b,a,v<z`; the same square identities
hold. The example is `m=6,r=9`.

**Failure boundary:** infinitely many arithmetic hosts do not prove that
each host contains a qualifying source support. The operation only closes
after retaining at least the selected edge, endpoint orientation and
one-cycle sector. Further births can create more cycles in larger support
carriers, so this is not a bounded-state all-size enumeration claim.

## 9. An actual directed cycle-opening operator

The undirected path reversal in section7 cannot be carried over silently
to a tournament. There is a different exact directed operation.

Let `D` be a finite loopless digraph, fix four distinct vertices `a,b,v,t`,
and assume the allowed arc `a→v` is present in `D`. A source support has
one directed path `b→...→t` and zero or more disjoint directed cycles,
together spanning every vertex. Require the selected arc `a→v`. Adjoin a
new vertex `z`, retain every old arc, and require the allowed arcs `a→z`
and `z→b`. Other new incident arcs can be oriented arbitrarily. On supports
perform

```text
delete a→v; insert a→z→b.
```

**Directed cycle-opening theorem.** The resulting support has start `v`
and end `t`. It is a directed Hamiltonian path if and only if the source
has exactly one directed cycle and `a→v` lies on that cycle.

**Proof.** If `a→v` belongs to a cycle, deleting it opens that cycle as
`v→...→a`. The new arcs join it to the old endpoint path as
`v→...→a→z→b→...→t`. This is spanning and connected precisely when no
other cycle remains. If `a→v` lies on the endpoint path, the new arcs
instead close `b→...→a→z→b` into a directed cycle and leave the path
`v→...→t`. This cannot be Hamiltonian, and any old cycles also remain.
No segment has been reversed. QED.

**Exact inverse and count.** A directed Hamiltonian path in the enlarged
graph starting at `v`, ending at `t`, and with the selected predecessor and
successor of `z` equal to `a,b` necessarily has the order

```text
v→...→a→z→b→...→t.
```

Remove `z` and add the allowed arc `a→v`. This yields exactly one cycle
containing that arc and the remaining path `b→...→t`. The target path
cannot already select `a→v` because `v` is its start. These operations
are mutual inverses. Therefore the count of this precisely conditioned
target path sector equals the count of the named one-cycle source sector.

For a tournament the source directed cycle has length at least3, and the
required extension is always realizable by orienting the two named new
arcs and completing the other incident pairs arbitrarily. With
`h_{x,y}(D[S])` denoting the number of directed Hamiltonian paths from `x`
to `y` on the specified vertex subset, the common count is also exactly

```math
\sum_{\substack{S\subseteq V(D)\\a,v\in S;\ b,t\notin S}}
h_{v,a}(D[S])\,h_{b,t}(D[V(D)\setminus S]).
```

This follows by cutting the selected cycle edge; the cycle vertex set `S`
is unique. In a tournament the two-vertex term is zero since `a→v` forbids
`v→a`. This formula retains named ports and the whole vertex partition;
it is not a scalar recurrence on `H(D)`.

**Small tournament hostile.** Take `T=C3→TT2` with directed cycle
`1→2→3→1`, all its vertices dominating4 and5, and `4→5`. It has three
Hamiltonian paths. Add vertex6 with fixed arcs `3→6→4`, and choose its
remaining relations to1,2,5 arbitrarily. For all eight extensions the
conditioned sector with start1, end5 and selected triple`3→6→4` is exactly
the single path `1→2→3→6→4→5`. Nevertheless the total path counts are:

| Bits recording whether6 beats1,2,5 | Total H | Conditioned count |
|---|---:|---:|
| 000 | 9 | 1 |
| 001 | 3 | 1 |
| 010 | 17 | 1 |
| 011 | 5 | 1 |
| 100 | 17 | 1 |
| 101 | 5 | 1 |
| 110 | 23 | 1 |
| 111 | 5 | 1 |

Thus the transfer is a genuine all-order tournament count identity, but
it only controls its declared sector. The remaining positive sectors can
change the physical total. A proof of Chan–Pak would still need an
operation controlling those contributions and returning the required
source ports. The cycle count here is a support coordinate, not an
identification with the tournament strong-core defect from older work.

## 10. Connection ledger and decisive tests

| Source → target | Map and preserved predicate | Lost data / required sidecar | Decisive test |
|---|---|---|---|
| Endpoint-conditioned Q_n → theta response | Sound degree forcing; exact residual degree supports | Outside forced-component partition for connectedness | Identical Q25/Q26 kernel gives1 versus2 paths |
| Multi-arm supports → Johnson product | Selected odd/even arm subsets; same-parity exchanges preserve all degrees | Actual arm edge weights and global connectivity | Both Q25 kernels have K2 exchanges but one surviving vertex |
| Degree support + selected edge → next support | Replace a–v by v–z–b; exact new degrees and inverse | Old path orientation or cycle membership | At25→26 one path survives and one cycle opens |
| Directed one-cycle support → tournament path sector | Delete a→v and add a→z→b; exact named-port count | Every other path sector in the actual extension | Same conditioned count1, total H ranging over3,5,9,17,23 |
| Square-shell arithmetic → infinitely many hosts | The three exact square sums in section8 | Source support existence and interval-spanning property | m4 gives a host, not a new Hamiltonian path theorem |
| Reflection connection → square graph | Edge equation `x_u+x_v=q_e`; signless transport | Fixed labels, domain, injectivity, coverage | Whole Q15 is rigid while a selected path is not |

The arm exchange is a genuine local involution on an exact state space.
There is no claim that it satisfies a cluster-algebra exchange relation or
that its Hamiltonian subgraph is a cluster exchange graph. Similarly the
first25-prime exponent sieve and the first25 square graphs have different
objects and predicates. Their shared useful method is retaining a complete
finite boundary state; equality of the number25 supplies no mathematical map.

## 11. Reproduction and independent controls

Run in this directory:

```text
python -B verify.py --output receipt.json
python -O -B verify.py --output receipt.optimized.json
```

Both runs produce byte-identical JSON. The checks are explicit exceptions,
not removable assertions. The exact finite universes are:

- All2,600 endpoint sectors in Q2 throughQ25; all paths are independently
  enumerated by vertex DFS whose adjacency is generated by square shells,
  not by the edge-forcing implementation. The DFS uses only the proved
  isolated-vertex and mandatory-leaf endpoint filters.
- The entire Q26 path set, independently enumerated, and exact equality of
  the two endpoint kernels, full outside partitions and all completions.
- Every one of26,656 edge subsets of the54 simple, ordered three-arm graphs
  whose lengths lie in1..4, compared in all16 hub-demand sectors with the
  coefficient formula. Subsets are enumerated edge by edge, independently
  of the alternating-arm generator.
- All126 parity words of1..6 arms, all1,470 feasible degree sectors and all
  5,460 configurations, checking the Johnson-product degree, connectedness
  and diameter.
- Every spanning degree-(1,1,2,...) support on4..7 labelled vertices:
  respectively12,70,465,3507 supports. Every eligible oriented selected edge
  and endpoint orientation is tested. An explicit two-triangle-plus-edge
  support adds the multiple-cycle hostile at8. Together there are63,552
  relocation checks:22,584 retained paths,22,584 paths closing a cycle,
  16,680 single cycles opened, and1,704 cases leaving another cycle.
- All83,565 arithmetic gadgets in the section8 parameter range4<=m<=80,
  together with the general proof of the infinite family.
- All loopless directed endpoint-degree supports on4..7 vertices:
  respectively36,220,1590,12978. The subuniverse with no directed two-cycle
  is precisely the supports embeddable in tournaments; its corresponding
  counts are24,160,1140,8988. Every eligible selected arc is tested, together
  with an explicit two-directed-triangles-plus-path hostile at8. There are
  59,334 directed relocation checks, including40,950 tournament-realizable
  checks. The eight fully specified extensions of the five-vertex example
  are independently enumerated by all720 vertex orders.

Finite controls corroborate the all-size proofs; they are not their basis.
The complete receipt retains every first25 endpoint verdict, all successful
paths, the three full kernel certificates and the birth outputs.

Independent proof review of sections3,7,8,9 is recorded in
[the sibling audit](../ports/AUDIT_SQUARES.md). The present lane separately
reviewed the sibling short-mutation proof and regenerated every rank-three
and rank-four case table in [audit_mutation.py](audit_mutation.py), with
[its exact receipt](mutation_audit.json). Those local checks corroborate
the all-size restriction argument and do not replace it.
