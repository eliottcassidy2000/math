# Native kernel index: a bounded positive and an all-order diagonal obstruction

**FINITE-EXACT n=3,...,7; PROVED / INDEPENDENTLY AUDITED elementary all-order
diagonal obstruction.** This continuation of the
[native parity index](overnight_hexagon_sep05_boolean_peck.md) tests the
actual Boolean adjacency, not a multiplicity-weighted surrogate.

## 1. Complete small kernel test

The closest mechanism is the bipartite index
`nullity(B)>=q_even-q_odd`. Its known equality through n=6 does not prove
equality at larger orders. The least-used sidecar is now the entire
even/odd crossblock, not its vertex count. The complete n=7 test has
32,768 labelled Eulerian graphs,54 isomorphism classes,27 classes on each
side and205 native triangle edges. Its27-by-27 Boolean crossblock has
determinant1076=4*269, so the adjacency is invertible.

The nullities for n=3,...,7 are0,1,1,4,0, exactly the parity indices in
this bounded universe. This is not an all-order saturation theorem.
The n=7 determinant also refutes a unimodularity or mod-two-invertibility
shortcut to proving saturation. Characteristic zero and characteristic
two are different kernel questions.

The checker generates every labelled cycle-space state, builds full
isomorphism orbits using adjacent transpositions (which generate S_n),
checks their counts against the independent Burnside formula, and compares
the whole orbit map with literal permutation enumeration through n=5.
It then checks **every labelled state's complete triangle-neighbor row**
against its orbit representative, not merely a sampled set of vertices.
This also checks all repeated destinations before Booleanization.

## 2. A four-entry obstruction to all diagonal repairs

Let M be the raw multiplicity-weighted triangle operator on Eulerian graph
classes at any n>=5 and B its Boolean support. Use these four classes,
padding with isolated vertices as necessary:

```text
even rows: C4, bowtie (two triangles meeting at one vertex);
odd columns: C5, H7 (K5 with one triangle removed).
```

The corresponding submatrices are exactly

```text
M_sub = [4(n-4)  2(n-4)]      B_sub = [1 1]
        [   4       4  ],             [1 1].
```

For C4 -> C5, choose any of four cycle edges and one outside vertex;
toggling that triangle replaces the edge by a two-edge path. For
C4 -> H7, choose one of the two opposite vertex pairs and one outside
vertex; the added triangle has no edge in common with C4. No other
triangle can yield these target types, as edge count and support size
already force the respective intersections.

For bowtie -> C5, choose the shared center and one noncentral vertex
from each triangle, giving four choices. The toggle removes the two
center edges and joins their other endpoints. For bowtie -> H7, choose
three of the four noncentral vertices, again four choices. Exactly one
old leaf-pair edge is removed and the two cross-pair edges are added.
A triangle meeting a new isolated vertex cannot preserve the target's
five-vertex support in either transition. These cases exhaust the choices.

The multiplicity block has determinant8(n-4), whereas the Boolean block
has determinant zero. Multiplying rows and columns by nonzero scalars
preserves the rank of every fixed submatrix. Consequently

```text
there are no invertible diagonal D_1,D_2 with B=D_1 M D_2
for any n>=5.                                       (1)
```

This rules out even independent row/column rescalings, not only a single
similarity or conductance normalization. Equivalently the alternating
product around this four-cycle has multiplicity ratio2, while the Boolean
ratio is1. It is an exact retained cycle-weight obstruction.

The threshold is minimal: at n=3 M=B. At n=4, in the order empty,C3,C4,
`M=[[0,4,0],[1,0,3],[0,4,0]]`; taking
`D_1=diag(1/4,1,1/4)`, `D_2=diag(1,1,1/3)` gives B. The n=4 failure
of the unmodified Fourier kernel therefore does not by itself rule out
diagonal repair; the four-cycle at n=5 supplies the missing obstruction.

## 3. What remains worth trying

The preserved predicate is the actual crossblock kernel. Vertex weights
discard cyclewise conductance products and fail by(1). The richer
complement-paired Fourier basis retains those products as a full matrix.
If its square free-pair block is invertible, a Schur complement would
correct each self-complementary weighted mode into a native Boolean
kernel vector. That invertibility is a new obligation, not a consequence
of either index equality or the existing weighted spectrum.

```bash
python3 -B 04-computation/eulerian_boolean_kernel_overnight_hexagon_sep05.py
python3 -B -O 04-computation/eulerian_boolean_kernel_overnight_hexagon_sep05.py
```

The [script](../../04-computation/eulerian_boolean_kernel_overnight_hexagon_sep05.py)
and frozen output retain the complete finite universe. No larger graph
census, general nullity formula or Laplacian-gap claim is made here.

The independent observer audit checked every all-order transition count,
the outside-vertex exclusion, two-sided diagonal implication and minimal
n=5 boundary. The root primary normal/optimized n=3,...,7 outputs match;
the n=7 enumeration is not described as independently reconstructed.
Raw LF source SHA256:
`658cad094d3651a9a7bb4ab88eadb1b35f50697614e8fd90801b82a9d5346484`;
[output](eulerian_boolean_kernel_overnight_hexagon_sep05.out) SHA256:
`599a24980a9d2aa5a89dbce093eb7a243e69c5eabaad2fa1f3137690fcb86b84`.
