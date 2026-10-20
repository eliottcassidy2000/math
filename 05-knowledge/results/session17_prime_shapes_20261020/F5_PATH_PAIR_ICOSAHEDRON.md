# From 625 path pairs to 576 nonzero-column matrices and an icosahedral return component

20 October 2026. **PROVED:** elementary field counts, explicit icosahedral
adjacency, and marked finite-set transport. **INHERITED PROVED:** the
balanced path decoder and the regular-five two-complement return theorem.
**FINITE-EXACT:** complete field, matrix, code, and action controls in
`verify_field.py`. No historical novelty or canonical unmarked path graph
is asserted.

## 1. The source has an actual five-by-five product

The established result in
`05-knowledge/results/zero_heavy_ports_20261020g/microcosm.md`, sections 2--3,
gives all Hamiltonian paths of the balanced two-return tournament
`T_(3,3)` as `X3 x X3`, where

\[
X_3=\{\infty,\varnothing,\{1\},\{2\},\{1,2\}\}.
\]

The involution kappa fixes infinity, interchanges empty and `{1,2}`, and
interchanges `{1}` and `{2}`. On path coordinates,

\[
J(I,K)=(\kappa I,\kappa K),\qquad R(I,K)=(K,I).
\]

Here `R` is actual positional anti-reversal and `J` is the separately
defined complement-of-choices involution. The standard increasing path
is `(infinity,infinity)`, the unique J-fixed point. T_(3,3) is rigid;
these nontrivial path permutations are not tournament automorphisms.

Choose the explicit marked bijection

\[
\ell(\infty)=0,\quad\ell(\varnothing)=1,\quad
\ell(\{1,2\})=4,\quad\ell(\{1\})=2,\quad\ell(\{2\})=3
\]

to F5. It intertwines kappa with negation. Hence the 25 paths become
vectors `v=(ell(I),ell(K))` in F5^2, with `J(v)=-v` and `R(x,y)=(y,x)`.
Every conversion retains the inherited path decoder, not merely its count.

## 2. What deleting a row and column really counts

An ordered pair of source paths becomes a 2-by-2 matrix with those two
vectors as its columns. There are 625 matrices. Removing pairs whose first
or second path is standard removes the zero-column row and column:

\[
625-25-25+1=576=24^2.
\]

This is **not** deleting the 25-element diagonal, which leaves 600.
Both are meaningful operations on the same ordered-pair carrier, but
have different predicates.

Among the 576 nonzero-column matrices, 96 have dependent columns:
choose a nonzero first vector in 24 ways, then a nonzero scalar in four
ways. The remaining 480 are invertible. For each fixed nonzero determinant
a, choose the first vector in 24 ways; the determinant functional has five
solutions for the second vector, all nonzero. Therefore

\[
576=96+4\cdot120,\quad |GL_2(5)|=480,\quad |SL_2(5)|=120. \tag{1}
\]

The user’s other sum `625+125+25+1=776` remains a distinct arithmetic
identity. No polynomial-size coincidences have been substituted for the
matrix carrier.

## 3. The twelve nonstandard J-orbits have a determinant graph

Let

\[
V=(\mathbf F_5^2\setminus\{0\})/\{\pm1\},\qquad |V|=12.
\]

These are exactly the nonstandard J-orbits of the balanced path set under
the chosen labels. Define distinct `[v],[w]` to be adjacent iff

\[
\det(v,w)^2=1\quad\text{in }\mathbf F_5. \tag{2}
\]

Changing either representative's sign does not change (2). For fixed
`[v]`, determinant zero gives exactly `[v]` and `[2v]`; the other ten
classes split into five with determinant square one and five with square
four. The map `[v] -> [2v]` is a fixed-point-free involution, since
`2^2=-1` in F5.

Quotienting the 576 ordered nonzero-vector pairs independently by their
two sign choices gives the exact decomposition

\[
144=12\text{ diagonal}+12\text{ antipodal}+60\text{ adjacent}
      +60\text{ distance-two}. \tag{3}
\]

The last description as distance two is justified by the graph certificate
below, not merely by the determinant value.

## 4. An explicit all-edge icosahedron certificate

Index i modulo five and put

\[
N=[(1,0)],\quad S=[(2,0)],\quad
U_i=[(i,1)],\quad L_i=[(2i-1,2)].
\]

These are twelve distinct classes. Direct determinant calculation gives
exactly

\[
N\sim U_i,\quad S\sim L_i,\quad
U_i\sim U_{i-1},U_{i+1},\quad
L_i\sim L_{i-1},L_{i+1},\quad U_i\sim L_i,L_{i+1}. \tag{4}
\]

For example `det(U_i,U_j)=i-j`, and
`det(U_i,L_j)=2(i-j)+1`, so cross-adjacency means `j=i` or `i+1`.
The other cases are equally direct. Formula (4) is the standard
five-cycle / five-cycle / two-pole icosahedral graph. In it `[2v]` is
the antipode at graph distance three; all remaining nonedges have distance
two. No uniqueness theorem for arbitrary five-regular twelve-vertex
graphs is needed.

## 5. Exact isomorphism to the previously proved tournament return component

The established result
`regular_return_gram_icosahedra_20261020/regular_return_theorems.md`,
section 5, has a twelve-state component encoded by tournament bits:

\[
N\mapsto947,\quad S\mapsto76,
\]
\[
(U_0,\ldots,U_4)\mapsto(613,217,710,341,668),
\]
\[
(L_0,\ldots,L_4)\mapsto(313,682,355,410,806).
\]

Bits encode `[i->j]`, for i<j in lexicographic order, first pair the low
bit. Those states are regular five-tournaments on one orientation-parity
sheet. Two distinct states are linked exactly when their codes differ in
six bits, equivalently when a legal two-step path-complement return joins
them. The prior theorem supplies this legality, not just a Hamming test.
Its full adjacency certificate is (4). Therefore the displayed table is
an explicit graph isomorphism from (2) to that **derived return graph**.
Antipodal multiplication by two corresponds exactly to tournament
converse (bitwise complement of all ten bits).

This is a true object-level connection between two previously different
carriers. It is not an isomorphism from the original tournament metagraph,
nor a statement that a source Hamiltonian path itself is a regular-five
tournament without the specified decoder and quotient choices.

### The symmetry actions can also be retained

SL2(5) preserves determinant. Its action on V has kernel exactly `{I,-I}`:
an element fixing every class has each coordinate vector as an eigenvector
with eigenvalue ±1, and fixing the class of their sum forces the two signs
to agree. Thus its induced group has order sixty.

Under the displayed code table, its action equals the action of the sixty
even relabellings of the five tournament vertices. This has an economical
explicit certificate. The standard generating matrices

\[
A=\begin{pmatrix}1&1\\0&1\end{pmatrix},\qquad
B=\begin{pmatrix}0&-1\\1&0\end{pmatrix}
\]

correspond respectively to the vertex permutations

\[
(0,1,2,3,4)\mapsto(2,3,4,0,1),\qquad
(0,1,2,3,4)\mapsto(3,4,2,0,1).
\]

The actions agree on all twelve displayed states by substitution.
These matrices generate SL2(5): powers of A give all upper elementary
unipotents, conjugating by B gives the lower ones, and elementary row
reduction generates every determinant-one matrix. Their induced action
has sixty elements and lies among the sixty even relabellings, so the
actions coincide. This explicitly realizes SL2(5)/{±I} as A5 in this
model; it does not rely on a bare equality of group orders.

### An actual source symmetry also transports exactly

The original positional anti-reversal R acts by `(x,y) -> (y,x)`, a
matrix of determinant minus one. Multiplying this matrix by the scalar
two gives determinant one. Under the code table, `2R` is the even vertex
permutation `(0 1)(3 4)`, and the scalar two is antipodal converse.
Consequently the **actual source anti-reversal** transports to

\[
T\longmapsto \big((0\ 1)(3\ 4)\cdot T\big)^{op}. \tag{5}
\]

The equality is checked on all twelve states, or by substitution in the
explicit table. Thus the connection preserves a pre-existing source
operation as well as the chosen determinant relation.

R fixes four classes in V: the two classes on `y=x` and the two on
`y=-x`. The former arise from genuine R-fixed (blue) source paths. The
latter arise from pairs of black source paths interchanged by both R and
J, so only their J-orbit is R-fixed. This recovers the prior product
formula's five R-fixed J-orbits after restoring the omitted standard
path. It also shows why fixed quotient objects need not lift to fixed
individual paths.

## 6. The marking is real mathematical data, not a cosmetic choice

There are eight sign-equivariant bijections from X3 to F5 sending infinity
to zero. The chosen one defines the field addition, determinant, and
resulting graph on the actual path-orbits. Negation alone does not select
these structures.

For a sharp hostile, the alternative relabelling `h(x)=x^3` also commutes
with negation and fixes zero. The classes represented by `(1,1)` and
`(2,1)` have determinant square one, but applying h coordinatewise gives
`(1,1)` and `(3,1)`, with determinant square four. Therefore the induced
edge relation is not invariant under every permitted five-state marking.
Every such marked model is abstractly an icosahedron by transport; it is
not a canonical path adjacency determined solely by J.

**Connection contract.** Source: actual Hamiltonian paths of T_(3,3), with
its decoded X3^2 coordinates, the standard path, and J. Target: one labelled
regular-five return component. Map: choose ell; remove the standard path;
quotient by J; use the displayed vector-to-code table. Preserved: the
specified determinant adjacency, its antipodal involution, the actual positional anti-reversal through
(5), and (after transport) the sixty-element action. Lost: the sign representative and
original seven-vertex path incidences. Needed sidecar: the explicit
five-state marking and inherited path decoder. Cheapest tests: all 144
quotient pairs, and the x^3 marking hostile.

## Reproduction

Run `python3 -B verify_field.py` and `python3 -B -O verify_field.py`.
Both produce `field_receipt.json` byte-for-byte. The finite universe is
all 625 matrices, all twelve vector classes, all 144 ordered class pairs,
all 120 determinant-one matrices, all sixty even five-vertex permutations,
and all 5040 seven-vertex orders for the balanced source. Its independent
run decoder yields exactly the same twenty-five legal source paths as the
literal census; all twenty-five also agree with the balanced-deletion
packet's separately constructed marked path table. No producer code is imported from the prior return packet;
its explicit theorem certificate is copied and checked as independent data.
