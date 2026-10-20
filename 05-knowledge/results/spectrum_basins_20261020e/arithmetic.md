# Fibonacci products, parity-path matchings, and the marked torus return

**PROVED:** the golden SCC word language has a reflection-equivariant,
triangle-potential-preserving bijection with matchings on two paths. This
explains the squares/products and the sum of squares, including273.
A marked edge gives the exact split63+210; its reflected pair selects
a reflection-invariant105-element family, comprising54 merged types.
**CITED / explicitly marked construction:** consecutive Fibonacci counts
give Farey slopes and torus knots with surface framing equal to their product.
**FINITE-EXACT:** all7920 words through vertex order18, independent matching
enumeration, arithmetic through order200, normal and optimized Python.
No correspondence of individual tournaments with knots or black moves is claimed.

## 1. Inheritance and objects

The closest proved mechanism is
[the earlier golden SCC language](../scc_basins_20261020d/arithmetic.md), §2:
the alphabet consists of the three ordinary strong-component types

```text
A=one vertex, B=C3, C=R4,       weights 1,3,4.
F(z)=sum a_n z^n=1/(1-z-z^3-z^4).
```

Here `a_n` counts ordinary tournament isomorphism classes in this restricted
SCC language, including the empty word. SCC order is part of the object.
It does not count all tournaments, all black components, or just the active
ones. The previous file already gave the golden asymptotic and a binary
prefix codec. The new carrier explains the coefficients directly and retains
the converse involution and triangle potential.

The prior arithmetic bridge is **THM-3509**,
[reduced-fraction harmonic K4 face and Fibonacci unit-Cassini ray](../../../01-canon/theorems/THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray.md),
§7. It supplies the primitive Farey / Fibonacci / Euclid dictionary, with
the parity-content sidecar. The canonical hostile for losing geometric data
is **THM-2646**,
[braid-three modular central pullback and full-twist knot fibre](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md),
§2: the same exact modular braid image admits unknot and trefoil lifts, and
one modular fibre has unbounded knot genus. MISTAKE-225, in
`01-canon/MISTAKES.md`, separately forbids identifying tournament games,
torus dynamics, and antisymmetry without an actual map.

The concept board is: ordered SCC language; parity-path matching;
converse reflection; directed-triangle potential; marked Farey slope;
torus framing; actual Hamiltonian-path moves. Each map below declares
which of these it carries.

## 2. Exact matching bijection

Let `G_n` have vertices `0,...,n-1` and edges `{i,i+2}`. It is the disjoint
union of paths on the even and odd labels, with respectively
`ceil(n/2)` and `floor(n/2)` vertices. This `G_n` is an auxiliary undirected
graph, not the tournament metagraph.

Replace each word atom by the following matching on its consecutive block:

| Atom | Block length | Selected edges relative to its first label |
|---|---:|---|
| A |1|none|
| B |3|`{0,2}`|
| C |4|`{0,2},{1,3}`|

Concatenation gives a matching of `G_n`. For the inverse, inspect the least
remaining label `i`. If it is unmatched, remove the A block. Otherwise it
must be matched to `i+2`. The vertex `i+1` is either unmatched, giving B,
or matched to `i+3`, giving C. All removed vertices now have fixed matching
status, so no selected edge exits the removed block. Iterate. This proves
existence and uniqueness of the decomposition.

Each block matching is invariant under reflection of its own block. Therefore
reflection `i -> n-1-i` reverses the atom word. Since all three component
types are self-converse, this is exactly tournament converse on SCC words.
The map is thus an equivariant bijection before and after the converse
quotient. The proof was independently checked by the H-spectrum agent.

Each B contains one directed triangle and C contains two; no directed
triangle crosses SCCs. Consequently

```text
c3(word) = #B + 2#C = number of selected matching edges.       (1)
```

The whole word remains reconstructible from the full matching. If only
matching size is retained, the C inventory is lost: on seven vertices,
`BBA` and `CAAA` both have c3=2, but have respectively H=9 and H=5.
Here H counts labelled Hamiltonian paths in a fixed representative.
The greedy block decomposition restores this missing coordinate.

## 3. The products and sums of squares

Write `F_0=0,F_1=1`. A path with r vertices has `F_(r+1)` matchings:
its final vertex is either unmatched or belongs to the final edge. Thus

```text
a_(2k)   = F_(k+1)^2,
a_(2k+1) = F_(k+1) F_(k+2).                              (2)
```

The user's sequence is exactly the odd subsequence

```text
n:       1,3,5,7,9,11,13,15,17,19
a_n:     1,2,6,15,40,104,273,714,1870,4895.
```

There is also a direct combinatorial sum-of-squares decomposition. In a
pair of matchings on `P_k` and `P_(k+1)`, inspect the last vertex of the
longer path. If unmatched, delete it and obtain two matchings on `P_k`.
If matched, remove its final matching edge and exchange the two factors;
the same procedure now acts on `P_(k-1)` and `P_k`. Hence the product set
partitions into squares of matching spaces, and

```text
F_(k+1) F_(k+2) = sum_(i=1)^(k+1) F_i^2.                (3)
273 = 13*21 = 1+1+4+9+25+64+169.
```

The algebraic telescope is `F_i^2=F_i F_(i+1)-F_(i-1) F_i`.
The sum-of-squares partition has a terminal-edge stopping statistic; it
does not by itself select the proposed split63+210. At n13, neither63 nor210 is
even a union of entire fibres of any one of the following three retained
statistics (an exact subset-sum control is included):

```text
c3 counts:                  1,11,46,91,86,34,4
#R4 counts:                 88,130,51,4
#C3+#R4 counts:             1,21,113,129,9.
```

The initial coarse-inventory search therefore left the split unassigned.
The missing coordinate is a positional edge mark, supplied by the next
construction. The earlier negative inventory checks remain valid.

### A marked-edge partition, and a reflected105 family

At n13 the shorter parity path is

```text
1--3--5--7--9--11                              (zero-based G13 labels).
```

Mark its second edge `e={3,5}`. A matching containing e leaves an isolated
first vertex and a path on its last three vertices, hence
`1*3=3` possibilities on this P6. The other parity path is P7 with21
matchings, so

```text
contain e: 3*21=63,           avoid e: (13-3)*21=210.       (3a)
```

This is an actual partition of the273 SCC words under the bijection.
The edge mark is essential: different edges of P6 give different counts.
For example its endpoint edge leaves P4, giving `5*21=105` containing
matchings; its middle edge leaves two P2s, giving `4*21=84`. The chosen
second edge becomes `{4,6}` if the G13 labels are one-based.

Global reflection sends e to `e'={7,9}`, the fourth edge of P6. Containing
both leaves just the first and last vertices of P6 unmatched, so there
is one such short-path matching and21 full matchings. Inclusion-exclusion
gives the complete ledger

| Membership | Ordinary SCC words |
|---|---:|
|neither e nor e'|168|
|e only|42|
|e' only|42|
|both|21|
|at least one|105|

Consequently the reflected-pair union is

```text
105 = 21 both + 84 exactly one,
273 = 168 neither + 105 at least one.                    (3b)
```

The63-element part is exchanged with its reflected counterpart and
does not descend as a subset of converse-merged classes. The105-element
union does descend. A reflection-fixed member of that union must contain
both marked edges, leaving the unique short-path matching and one of
the three reflection-fixed P7 matchings. Thus it has three fixed words,
and Burnside gives `(105+3)/2=54` converse-merged types. The21-member
intersection gives12 merged types; the84-member exactly-one sector gives42.

Explicit preimages, writing A,B,C for weights1,3,4, are

```text
no selected edges anywhere:    A^13,
only e:                       A^3 B A^7,
only e':                      A^7 B A^3,
e and e' alone:                A^3 B A B A^3,
edges {2,4} and e:             A^2 C A^7.
```

The three fixed words in the union are

```text
A^3 B A B A^3,       B B A B B,       A^2 C A C A^2.
```

The first two are pure-blue rigid palindromes; the third contains R4 and
is active. Thus the union has103 active ordinary words and52 active merged
types. Even this marked105 family is not itself a105-vertex black core.

This construction supplies the requested structural rule using the
positional matching sidecar. Its21+84 split has the same numbers as the
earlier marked Fano/junction orbit split, but **no map intertwining those
actions or preserving their closure predicates has been established**.
The SCC carrier here is a matching family on two paths with a reflected
edge pair; it must not be silently identified with the105 junction tensors.

## 4. Weighted refinement and converse quotient

Let `P_r(y)` be the matching polynomial of a path with r vertices,
using y for each selected edge. Then

```text
P_0=P_1=1,       P_r=P_(r-1)+y P_(r-2),
P_r(y)=sum_j binom(r-j,j)y^j.
A_n(y)=sum_word y^c3 = P_floor(n/2)(y) P_ceil(n/2)(y),
sum_n A_n(y) z^n = 1/(1-z-yz^3-y^2z^4)
                  =1/[(1+yz^2)(1-z-yz^2)].            (4)
```

Thus uniform random golden words correspond to two independently uniform
path matchings, and c3 is the sum of their matching sizes. The previous
sum-of-squares partition retains its full potential grading:

```text
P_k(y)P_(k+1)(y) = sum_(i=0)^k y^(k-i) P_i(y)^2.        (5)
```

The factor `y^(k-i)` counts removed terminal edges. For n13,

```text
A_13(y)=1+11y+46y^2+91y^3+86y^4+34y^5+4y^6.
```

The matching model also computes the converse fixed points. A reflection-
invariant matching on a path of length r vertices has

```text
s_(2m)=F_(m+2),          s_(2m+1)=F_(m+1).              (6)
```

For even length, separate whether the central edge is present, obtaining
`F_(m+1)+F_m`. For odd length the central vertex must be unmatched: a
selected incident edge would force its reflected mate to share it.
Therefore, if `p_n` counts self-converse golden words,

```text
p_(2k)=F_(k+1),              p_(2k+1)=s_k s_(k+1),
merged golden count = (a_n+p_n)/2.                      (7)
```

At n13 this is15 fixed words and144 converse-merged words. Formula273
therefore belongs to the ordinary SCC-word layer, not the merged layer.
These numbers do not yet discard pure-blue words.

## 5. A precise Farey / torus / Pythagorean map

For odd n=2k+1, take the cardinalities of the two ordered matching factors:

```text
(p,q)=(F_(k+1),F_(k+2)).                               (8)
```

Consecutive Fibonacci numbers are coprime. The determinant of successive
vectors `(p,q)` and `(q,p+q)` has absolute value1 by Cassini. Thus these
vectors are primitive neighboring Farey slopes on a marked oriented torus.
The map preserves the product `pq=a_n` and the primitive / determinant-one
relations. It forgets every individual SCC word unless the actual pair of
matchings is retained as a sidecar. Advancing k changes the vertex order;
it is not an operation inside a fixed-order tournament metagraph.

There is also an exact upper-half-plane realization of this update. On
the ratio x=q/p a single Fibonacci step is `x -> 1+1/x`, with determinant
-1. Its square has matrix and action

```text
M=[[2,1],[1,1]],       det M=1,
z -> (2z+1)/(z+1),
Im(Mz)=Im(z)/|z+1|^2 > 0        when Im(z)>0.             (8a)
```

Thus the square preserves the upper half-plane and its boundary Farey
adjacency. Its real fixed points satisfy `z^2-z-1=0`, so are phi and
`-1/phi`. On the primitive column `(q,p)`, M gives `(2q+p,q+p)`, advancing
two Fibonacci steps and increasing tournament order n by4. The slope
ratios themselves lie on the real boundary; they are not interior points
of the half-plane. This provides a precise geometric model for the
recurrence, while retaining the change of order and torus marking. It
neither acts as a black move at fixed n nor invokes paradoxical
equidecomposition.

On the standard embedded torus in S3, the primitive curve of slope `(p,q)`
is a torus knot `T(p,q)`. Its torus-surface push-off has linking number pq.
This is the **surface framing**, with the torus embedding and its meridian-
longitude marking retained. The classical formula is used explicitly by
[Finashin and Kharlamov, Topology of real cubic fourfolds, §6, p25](https://arxiv.org/pdf/0906.1480).
It is not the minimal crossing number, genus, or the unrelated invariant
called a natural framing by Greene and Wiest.

This yields a lawful numerical bridge: the ordinary word count at order13
equals the torus-surface framing of `T(13,21)`. It is a map from *family
cardinality data* to marked knots, not an isomorphism between273 tournaments
and273 knots. The sharp first hostile is already `T(1,1)` versus `T(1,2)`:
both are unknots, while their torus framings are1 and2. Forgetting the
framing loses the count.

The same primitive pair gives Euclid's raw triple

```text
(q^2-p^2, 2pq, q^2+p^2),
content g=2 if p,q both odd, and g=1 otherwise.           (9)
```

For k>=1 the triple is nondegenerate. The first k=0 pair(1,1) gives the
degenerate triple and is excluded from the primitive Pythagorean lane.
Some exact controls are

| n | `(p,q)` | a_n | Primitive Euclid triple |
|---:|---|---:|---|
|3|`(1,2)`|2|`(3,4,5)`|
|5|`(2,3)`|6|`(5,12,13)`|
|7|`(3,5)`|15|`(8,15,17)`|
|9|`(5,8)`|40|`(39,80,89)`|
|11|`(8,13)`|104|`(105,208,233)`|
|13|`(13,21)`|273|`(136,273,305)`|

This explains exactly how the neighboring104 and105 can both occur without
equating them:104 is the SCC count and105 is the difference leg. At n13,
273 becomes the primitive doubled leg after division of the raw triple by2.

As q/p tends phi, the ratio `(q^2-p^2)/(2pq)` tends1/2. These triangles
approach a1:2 leg ratio; the isosceles limit would require q/p tending
`1+sqrt(2)`. No invariant yet sends complete graphs, bipartite graphs, and
snarks to these arithmetic roles.

## 6. The finite counting boundary

For a finite set X partitioned into finitely many disjoint pieces X_i,
piecewise bijections preserve `sum_i |X_i|=|X|`. They cannot give two
disjoint copies of X unless X is empty. This elementary obstruction applies
to the finite SCC and matching sets here. A quotient, duplication of labels,
or forgetting stabilizers changes the object being counted, and must be
recorded; it is not a finite Banach-Tarski mechanism. No infinite
equidecomposition theorem is used in the argument.

## 7. Reproduction and status

```text
python3 05-knowledge/results/spectrum_basins_20261020e/arithmetic_probe.py
python3 -O 05-knowledge/results/spectrum_basins_20261020e/arithmetic_probe.py
```

The checker independently enumerates all edge-subset matchings and all word
compositions through order18, checks both inverse and reflection, compares
all potential coefficients and the complete marked-edge partition and
its preimages, checks the formulas through order200, and tests
Cassini / Euclid arithmetic for100 successive pairs. It uses explicit
exceptions, not assertions disabled by `-O`. Both runs give receipt SHA256
`456135325d8b75ebbbbc4554d81f4a255cdb5980eab0921b97c6c83c18bf37ef`.
The receipt includes source hash and the full small inventory tables.
