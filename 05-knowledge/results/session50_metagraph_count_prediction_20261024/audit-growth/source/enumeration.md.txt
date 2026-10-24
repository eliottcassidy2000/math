# Exact class generation and the information needed for metagraph growth

**Status:** FINITE-EXACT census and hostile probes through the bounds below;
PROVED elementary obstruction from the explicit order-three and order-four
examples. No new general formula for A000568 is claimed.

## Inheritance and scope

The closest proved mechanisms are
[THM-781 — Hamiltonian-path inversion](../../../01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md)
and
[THM-1430 — the tiling/class/metagraph dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md).
The latter already implements class BFS and the Hamilton-path fibre method.
This is an independent implementation and a new compression audit, not a
novel generic enumeration algorithm.

The critical quotient loss is the Hamiltonian-path presentation: a tournament
class can have many tilings, and complementing only the non-path arcs can
send them into different classes. The sidecar to retain is the actual tiling,
equivalently a Hamiltonian path modulo automorphisms. The least-used relevant
sidecar here is the **distribution of outgoing targets**, instead of its
total degree.

We distinguish two graphs throughout:

- the blue/black **line multigraph**, obtained by pairing each tiling with
  its full tile complement;
- the **single-arc-flip class graph**, whose edges flip one tournament arc.

The first is the user's line construction. The second is used as an
enumeration engine and as one additional diagnostic coordinate.

## Exact Burnside-free census

Start with one transitive tournament. Flip each arc in each discovered
isomorphism class, canonicalize the result, and continue until the queue is
empty. This reaches every class: the labelled arc cube is connected and its
quotient by isomorphism is connected. Retaining one representative per class
suffices, since relabelling transports every arc flip to a flip of that
representative.

Canonicalization uses invariant score refinement followed by all permutations
inside the refined cells. This only restricts to an invariant ordering of the
cells; it does not drop possible isomorphisms. Independent minima over **all**
vertex permutations audit the resulting class representatives through order
six. A full fixed-path tiling census independently reaches exactly the same
class set through order seven.

| Order | Classes `a_n` | Self-converse `S_n` | Merged `M_n` | Paired merged nodes `P_n` | Tilings |
|---|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 0 | 1 |
| 2 | 1 | 1 | 1 | 0 | 1 |
| 3 | 2 | 2 | 2 | 0 | 2 |
| 4 | 4 | 2 | 3 | 1 | 8 |
| 5 | 12 | 8 | 10 | 2 | 64 |
| 6 | 56 | 12 | 34 | 22 | 1,024 |
| 7 | 456 | 88 | 272 | 184 | 32,768 |

Thus `a_n=M_n+P_n=2M_n-S_n`. Here `P_n` counts nodes representing
**two** converse classes; it does not count all original classes involved
in merging.

For an unmerged class `c`, the verified fibre is `H(c)/|Aut(c)|`.
The power `2^binom(n-1,2)` counts the **entire** fixed-path tiling cube,
not the fibre of every individual class. For example, at order four the
four class fibres have sizes `1,1,1,5`.

The script audits `H=fibre*|Aut|` for every class, the global fibre checksum,
converse involution, reflection-to-converse compatibility, and every line
incidence count. All labelled-orbit masses through order six also sum to
`2^binom(n,2)`.

Blue/black line totals at orders `3,4,5,6,7` are respectively
`(1,0),(2,2),(8,24),(32,480),(256,16128)`. The output gives every fibre
histogram, including the joint merged `(fibre,blue)` distribution.

## A minimal failure of tile complement on classes

Use vertices `0,1,2,3`, fixed path `0→1→2→3`, and free pairs
`(0,2),(0,3),(1,3)` in that order. Bit `1` directs the arc from the
smaller vertex to the larger one.

The tilings with integer masks `0` and `1` have adjacency row masks

```text
t=0: (2,4,9,3)
t=1: (6,4,8,3).
```

They are isomorphic by the permutation
`0↦1,1↦2,2↦3,3↦0`. But their full tile complements, masks `7` and `6`,
have sorted outdegree sequences

```text
(0,1,2,3),       (0,2,2,2).
```

The first is transitive. The second and its converse are nontransitive.
They therefore belong to different **merged** nodes. Tile complement does
not descend to a function on isomorphism classes, even after converse
merging. Order four is the smallest witness, as direct exhaustion through
order three confirms.

**Strongest survivor:** tile complement is an involution of the tiling cube,
and it induces a well-defined weighted correspondence between merged nodes.
Its matrix records how many tilings in one fibre have partners in another.
The lost coordinate is the location inside the source fibre.

## A minimal failure of an unmarked metagraph growth rule

At order three, the transitive and cyclic classes are both self-converse.
Each has one tiling, which is blue. The entire colored line multigraph is
one blue edge joining two indistinguishable vertices. Their local data are
identical:

```text
(fibre,blue,blue_other,black_other,blue_self,black_self,
 simple_line_degree,merged_arc_flip_degree)
 = (1,1,1,0,0,0,1,1).
```

Append a fourth vertex at the end of the declared Hamiltonian path.
There are four choices for its two remaining incident arcs. Let the three
merged child types be `(transitive T4, non-self-converse pair, other SC)`.
The child counts are

```text
transitive parent: (1,1,2),
cyclic parent:     (0,1,3).
```

If instead all eight incidence patterns for a new unmarked vertex are
allowed, the counts are `(4,2,2)` and `(0,2,6)`.

This is an automorphism obstruction: swapping the endpoints preserves the
whole unmarked colored source metagraph, while the growth profiles differ.
Consequently no rule intrinsic only to this unmarked graph can assign the
correct growth profile to each node. An additional mark is necessary.

**Boundary:** separate Hamiltonian-path and automorphism counts resolve this
particular ambiguity: they are `(H,|Aut|)=(1,1)` for the transitive triangle
and `(3,3)` for the cycle. This does **not** refute a rule retaining those
marks, and it does not rule out a global counting rule that never assigns
individual parent profiles.

## More marks do not make local line data complete

At order seven, two merged nodes have representative adjacency row masks

```text
U = (0,5,41,19,7,91,31),
V = (0,9,3,85,39,15,55).
```

An entry's binary bit `j` records whether that vertex points to `j`.
Both representatives have `H=19` and trivial automorphism group. Both
merged nodes have exactly the same ten diagnostic coordinates:

```text
(fibre,blue,blue_other,black_other,blue_self,black_self,
 simple_line_degree,merged_arc_flip_degree,H,|Aut|)
 = (38,0,0,38,0,0,16,20,19,1).
```

Nevertheless their distributions of black-line partners differ. `U` has
two incidences to a node with merged fibre `41`, whereas `V` has none.
The output records the entire partner distributions. There are nine
pairs with the same eight local coordinates and different partner
profiles at order seven; all nine remain collisions after separately
adding `H` and `|Aut|`. None occurs at smaller orders in the checked range.

This refutes reconstruction of line partner profiles from these ten
coordinates. It does not refute reconstruction from the complete marked
metagraph, since that graph explicitly retains the relevant adjacency.

The smaller tuple `(fibre,blue,H,|Aut|)` also fails to predict one-vertex
extension profiles at order six. There are six collision pairs in the
exact census. One pair has the common tuple `(2,0,3,3)` and representatives

```text
(0,1,3,23,39,15),       (0,1,11,19,7,31).
```

Their merged arc-flip degrees are respectively five and six, so the
larger local tuple separates this particular pair. With all eight local
coordinates plus `H,|Aut|`, no extension-profile collision occurs through
order six. That bounded absence is **not** a closed recurrence theorem.

## The resulting research target

The useful object is a system of compatible partitions, not a list of
degree constraints. A candidate compressed type must retain enough
information that both operators are constant on each type:

1. the weighted blue/black partner correspondence within an order;
2. the extension/deletion incidence between orders.

This is an exact closure requirement. Refining a type by its outgoing
target distribution repairs a failed test at one order, but it may
introduce more types at every larger order. The open target is a bounded
or recursively describable family of types with a proved closure law.
The witnesses above say precisely which proposed finite summaries already
lose necessary information.

## Reproduction and controls

```text
python 05-knowledge/results/metagraph_growth_20261020/enumeration.py
python -O 05-knowledge/results/metagraph_growth_20261020/enumeration.py
```

Both commands produced byte-identical stdout before newline normalization;
`enumeration.out` is that literal output. All validation raises explicit
exceptions and remains active with optimization. The universe is every
fixed-path tiling and every tournament class through order seven; extension
profiles use every incidence pattern from orders one through six.
Independent controls are full-permutation canonicalization through order
six, labelled orbit mass, Hamiltonian-path dynamic programming, known small
class counts, and the explicit order-four isomorphism above. No random
sampling or Burnside enumeration is used.
