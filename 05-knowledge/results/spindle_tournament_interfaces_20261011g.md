# From the four-arc square to spindle and multiport constraints

**Status: PROVED elementary interface statements; FINITE-EXACT small controls;
CITED published five-chromatic graphs. No new six-colour equality fragment.**

## 1. Inheritance and the precise four-vertex object

The closest mechanisms are the
[marked four-arc square](partial_four_tournament_20261011c.md), the
[unit-rhombus interfaces](unit_rhombus_interfaces_20261011c.md), and the
[complete boundary-colouring tables](boundary_palette_compression_20261011.md).
The earlier [field-escape recursion](unit_field_escape_recursion_20261011f.md)
already authenticates the spindle geometry and its equality mechanism. This
packet connects their different interfaces rather than identifying the objects.

The concept board is: named orientation ports, actual unit segments, fixed-palette
colour relations, complete multiport support, and geometric gluing. The cheapest
hostile is a proper four-colouring of a diamond whose missing-edge endpoints
differ. The least-used relevant sidecar is the **joint** boundary partition,
which cannot in general be reconstructed from all its pair projections.

The user's four prescribed arcs have the unique marked form

$$
 a\to b\to c\to d,\qquad a\to d,\qquad
 \text{missing pairs }ac,bd.
$$

Its underlying graph is $C_4$, not the diamond $K_4-e$, and certainly not the
seven-vertex spindle. With $x=1$ reversing the $ac$ orientation and $y=1$
reversing the $bd$ orientation from the path order, its four completions are

| $(x,y)$ | Marked-section tournament type |
|---|---|
| $(0,0)$ | transitive |
| $(1,0)$ | directed triangle over a sink |
| $(0,1)$ | source over a directed triangle |
| $(1,1)$ | strong |

This inherited type decoder is lossless only with its original partial
orientation and two named ports. Every completed tournament forgets to the same
undirected $K_4$. Neither an orientation bit nor a tournament type records
whether a Euclidean segment has unit length. No canonical tournament is asserted
on an arbitrary unit-distance graph.

## 2. The genuine square-to-diamond-to-spindle map

Write $\alpha=(1+i\sqrt3)/2$, and place the cyclic square labels at

$$
 (a,b,c,d)=(0,1,1+\alpha,\alpha).
$$

All four perimeter edges have length one. The diagonal $bd$ has length one;
$ac$ has length $\sqrt3$. Adding the actual unit segment $bd$ produces a
diamond whose marked tips are $a,c$. Adding an orientation on either missing
pair is a different operation: a completed tournament would also require the
nonunit pair $ac$ in its underlying graph.

For a palette of $k\ge3$ colours, and *specified* colours of the two tips, the
diamond has

$$
 e_k=(k-1)(k-2)\quad\hbox{if the tips agree},\qquad
 d_k=(k-2)(k-3)\quad\hbox{if they differ}.                 \tag{1}
$$

Indeed its other two vertices are adjacent and must avoid the colours used at
both tips. Thus the diamond is an equality gadget at $k=3$, whereas the
perimeter alone is not. At $k=4$, $(e_k,d_k)=(6,2)$; at $k=6$ it is $(20,12)$.
The implication does not persist when the palette changes.

Set $\beta=(5+i\sqrt{11})/6$. The two diamonds have common tip $0$ and respective
outer tips $1+\alpha$ and $\beta(1+\alpha)$. Since $|\beta|=1$ and
$|1-\beta|^2=1/3$, those outer tips are a unit distance apart. The complete
seven-point set is

$$
 0,\ 1,\ \alpha,\ 1+\alpha,\ \beta,\ \alpha\beta,\
 \beta(1+\alpha).
$$

It has exactly eleven unit edges: two five-edge diamonds and their outer-tip
closing edge. In every proper three-colouring, both outer tips would equal
the common tip, contradicting that edge. A proper four-colouring exists.
This is the ordinary Moser spindle. The script independently authenticates all
pair distances, a full non-three-colour proof, and a proper three-colouring
after deletion of each of its eleven edges.

The source, target, and losses are now explicit:

- source: the marked partial orientation on a unit rhombus;
- target: the diamond's palette-specific two-port colour relation;
- map: forget orientations, retain the coordinates, then add the actual unit
  diagonal $bd$;
- preserved data: named vertices and perimeter incidence;
- destroyed data: direction choices; they play no role in proper colouring;
- required sidecars: the metric, actual added edge, palette, and a complete
  colour relation.

This map requires extra data and an edge addition. It is not an equivalence of
a tournament with a spindle.

## 3. A geometric three-port hostile to pairwise compression

There is a small exact version of the higher-order issue. Take the induced
unit graph on

$$
 (z_0,z_1,z_2,p_0,p_1,p_2)
 =(0,1,\alpha,-\alpha,2,2\alpha).
$$

Its only edges are the triangle $z_0z_1z_2$ and the three spokes $z_i p_i$.
In three colours its boundary relation is exactly

$$
 \{(c_0,c_1,c_2):\text{the three colours are not all equal}\}.       \tag{2}
$$

If all three tips have one colour, every vertex of the inner triangle must
avoid it, which is impossible. Otherwise an assignment of three distinct
colours to that triangle avoiding the respective tip colours exists:
for a two-equal tip pattern or a three-distinct tip pattern, exactly two
permutations work. Consequently the graph has $24$ admissible labelled tip
assignments, each with two extensions, and $48$ proper labelled colourings.

For each of the three port pairs, both equality and inequality occur.
Therefore **all three pairwise projections are identical to those of three
independent vertices**, although the latter also permit the all-equal assignment.
An orientation imposed on these pairs supplies no missing colouring information.
This is a fully induced unit-distance hostile, rather than only an abstract
constraint counterexample. Its unequal port distances are retained; it is not
claimed to be one of the equilateral triples in the cited large constructions.

The positive repair is to store all four non-monochromatic equality partitions
of three labelled ports. More generally store the complete joint partition
relation, with lift witnesses, rather than only pairwise implications.

## 4. Complete cross-interface composition

Let $G_1,G_2$ be finite graphs, with ordered port sets $B_1,B_2$.
Require every vertex in $V(G_1)\cap V(G_2)$ to lie in both port sets.
Let $X$ be a set of additional edges whose endpoints lie in $B=B_1\cup B_2$.
For a fixed palette $k$, let $R_i$ be the complete set of equality partitions
on $B_i$ realized by proper $k$-colourings of $G_i$.

Then the complete relation of $G_1\cup G_2\cup X$ on $B$ is exactly

$$
 R=\{\pi:\ |\pi|\le k,\,
       \pi|_{B_i}\in R_i\ (i=1,2),\
       u,v\text{ lie in different blocks for every }uv\in X\}.    \tag{3}
$$

Necessity is restriction. For sufficiency choose a representative assignment of
the blocks of $\pi$ to colours. For each piece, a retained witness for
$\pi|_{B_i}$ can be permuted on the full palette to agree with that assignment.
The two witnesses then agree on their shared vertices. Their interiors are
disjoint, and the checked edges in $X$ account for every additional constraint,
so the combined witness is proper.

The exported **compose_interfaces** implements (3) on the inherited typed
**Table** objects. It audits relation completeness before composition; an
incomplete list of correct positive witnesses cannot justify absent states.
It returns an original-graph witness for each surviving partition. This extends
the earlier same-boundary intersection API to differing port sets and explicit
cross-edges; no new general graph-colouring theorem is claimed.

For geometric use the coordinate identification must separately match the
vertex labels and every retained edge must be a unit segment.
**geometric_graph** checks that condition. The default retains all unit edges;
an explicit edge subset is allowed but remains marked as such. A coincident
vertex cannot silently be treated as two independent interiors. Abstract
composition alone neither supplies a Euclidean embedding nor proves that a
declared cross-edge list is complete for an intended induced graph.

The spindle uses ports $(0,3)$ and $(0,6)$ with $X=\{36\}$. At three colours
the relation (3) is empty; at four and six it has the three partitions
$001,010,012$ in the ordered boundary $(0,3,6)$.
This is a direct, witness-preserving implementation of its constraint mechanism.

## 5. What the 509-vertex construction contributes

**CITED, not replayed here.** Jaan Parts' 2020
[graph-minimization paper](https://arxiv.org/html/2010.12665v1)
reports a five-chromatic graph with $509$ vertices and $2442$ edges.
Section 6 gives pieces of orders $374,136$ and sizes $1860,564$, sharing the
centre: the resulting incidence accounting has $18$ cross-edges. Section 4
distinguishes direct mono-pair spindle constructions from Type M constructions,
where auxiliary cross-edges matter. Section 5.3 also uses companion constraints.
Thus the reported 509 construction is not identified here with two diamonds
or a four-state tournament.

The cited Type M rotation is $\rho=(7+i\sqrt{15})/8$. It closes corresponding
radius-two tips because $4|1-\rho|^2=1$. Its quadratic generator is outside
$\mathbb Q(\sqrt{-3},\sqrt{-11})$: the latter's three quadratic subfields have
squareclasses $-3,-11,33$, not $-15$. Consequently a joined two-piece graph is
not confined to the old complex Moser field merely because each unrotated
piece was. Field escape itself still supplies no colouring obstruction.

Heule's earlier
[553-vertex paper](https://arxiv.org/html/1805.12181),
§§4.3–4.4, distinguishes vertex-critical from edge-critical graphs and analyzes
compatibility of key-vertex colour patterns between two pieces. That is the
relevant interface mechanism: compress each interior only to sufficient retained
boundary data, then impose actual cross constraints. The cited computations
are external results, not tests performed by this packet.

There is also an elementary palette hierarchy. From any finite graph of
chromatic number $k+1$, delete edges until a minimal non-$k$-colourable subgraph
remains. Removing any one of its edges gives a $k$-colourable graph. Every such
colouring must identify the removed edge's endpoints, or it would colour the
undeleted graph. If the starting graph is unit-distance, the retained endpoints
remain a unit apart and every surviving edge remains unit. This is an
edge-subgraph construction, not necessarily an induced graph.

Consequently the spindle yields three-colour equality ports; the cited
five-chromatic graph yields four-colour equality ports after certified thinning.
It does not thereby yield six-colour equality. Indeed a six-colourable
unit-distance fragment whose unit-separated ports are forced equal would become
exactly seven-chromatic when that unit edge is added: non-six-colourability
follows from the constraint, and seven colours suffice by recolouring one port.
Conversely an edge-minimal seven-chromatic unit graph gives such a fragment.
This is a genuine higher target, not a change of labels on the old construction.

## 6. Reproduction and scope

The companion uses exact rational arithmetic in the inherited marked Moser
field, complete finite partition enumeration, and source-bound colouring
refutations. Run from the repository root:

~~~powershell
python -B -X utf8 04-computation/experiments/spindle_tournament_interfaces_20261011g.py
python -B -O -X utf8 04-computation/experiments/spindle_tournament_interfaces_20261011g.py
~~~

The finite universes are four marked tournament completions; diamond counts
at palettes $3$ through $6$; the complete seven-point spindle and its eleven
edge deletions; all $27$ labelled three-tip assignments of the six-point toy;
and $384$ joins of two arbitrary three-vertex graphs sharing one retained
vertex, two cross-edge choices, and palettes $2,3,4$. Malformed source,
incomplete-table, missing-port, metric, and bool/float-label controls are included.

Normal and optimized runs agree on **1,307 always-active checks**.
No 509-vertex coordinate or SAT certificate is supplied or reverified, and no
new six-colour obstruction, equality fragment, or whole-plane colouring follows.
