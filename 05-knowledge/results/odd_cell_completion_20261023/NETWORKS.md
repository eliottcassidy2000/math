# Compatible odd-cell networks: marking, counting, and the first kernel

23 October 2026. **PROVED scoped**, with exact finite checks in `check_networks.py`.
The objects are tournaments with a supplied **named equitable partition**.
All internal cells and cross blocks retain their actual vertices. The claims
are not reconstruction from anonymous scalar margins, and are not an
unrestricted enumeration formula for A000568.

## 1. The new marking theorem

Suppose every cell has size 1, 3, 5, or 7. A non-singleton cell is internally
regular. Join two cells when their actual cross block is nonuniform; call
this the informative support. Equitability gives equal row degrees and equal
column degrees in every block. If the cell sizes a,b are coprime, a q=b p
forces the block to be uniform. Thus each informative component is made
entirely of cells of size 3, or entirely of size 5, or entirely of size 7.
Isolated non-singleton cells count as components.

**Theorem.** Choose one actual owner in one root cell of each size-three or
size-five component, and two distinct actual owners in one root cell of each
size-seven component. Directed equitable refinement, retaining the original
cell names and these marks, becomes discrete on the whole tournament.
If the numbers of components are c3,c5,c7, the added mark budget is at most

    c3+c5+2*c7.

Existing named singleton owners may already perform some of this work.
The statement is a sufficient budget; it is not the minimum for every input.
Two marks suffice for an arbitrarily large CONNECTED seven-cell component,
not for an arbitrary number of disconnected components.

Proof. One mark splits C3. The [five-cell proof](five/README.md) shows that
regular R5 has no proper nondiscrete equitable partition, so one mark splits
it. The [seven-cell proof](seven/PROOF.md) shows that the only possible proper
nondiscrete equitable partition of a regular seven-tournament is 1+3+3.
Two distinct marks therefore split it.

Once the root is discrete, follow an informative edge. In a five-cell,
nonuniformity supplies a proper refinement, which must become discrete.
In a seven-cell, a remaining 1+3+3 partition would admit a nonidentity C3
automorphism preserving each refined cell. Since every vertex in the already
discrete neighbour distinguishes exact incident bits, that automorphism
would preserve every column of the cross block. This contradicts the new
one-sided kernel lemma: a nonuniform biregular block into a regular seven-cell
has no such automorphism. The three-cell argument is immediate. Proceed
along a spanning tree. All later refinement only splits cells, so branches
and cycles cannot undo the result. QED.

The kernel lemma does NOT require distinct binary rows or columns. This is
the actual improvement over the incoming
[partial-bijection rule](../odd_cell_gluing_20261023/README.md): repeated
profiles can still determine the only permissible internal tournament chart.

For a source/target comparison, each informative edge consequently gives a
partial bijection on the internal chart domains. At most 3, 5, or 21 root
charts need be propagated per component. Check EVERY non-tree block on the
same propagated choices. All seven-cell groups have order 7, 3, or 21 in
the exact complete regular-seven census; the elementary kernel and marking
proofs do not depend on that census. A successful global comparison is a
cell/mark-preserving isomorphism, not a new Hamiltonian path.

## 2. Network Burnside formula: all shared choices at once

Here restrict each cell to the same p-vertex tournament R with
Aut(R)=C_p, acting by translation on F_p, for an odd prime p. This hypothesis
is stronger than 'regular of prime size'. Fix a connected simple support H
with v>=2 named cells. Give each edge an orientation and a fixed degree
1<=d_e<=p-1. Its allowed cross blocks are ALL binary p-by-p matrices with
row and column sums d_e. Blocks outside H are fixed uniform blocks. Their
orientations are part of the declared quotient, not choices to average.

Write A_e for the number of allowed matrices on edge e and B_e=binom(p,d_e).
Let N_H count configurations up to one shared cell rotation at each vertex.
Let S_H count those classes with nontrivial cell-preserving automorphism group,
and R_H those with trivial such group. Then

    N_H = (product_e A_e + (p-1)^v product_e B_e) / p^v,
    S_H = (p-1)^(v-1) product_e B_e / p^(v-1),
    R_H = (product_e A_e - (p-1)^(v-1) product_e B_e) / p^v.

These count ordinary cell-preserving classes. No converse or cell-permutation
quotient is included. One can multiply the component formulas for a
disconnected support; an isolated R-cell contributes one unmarked class and
its whole C_p automorphism group.

Proof. The group C_p^v acts independently on the named cells but with the
SAME rotation on every incident edge. For a group element (a_i), if exactly
one endpoint of an edge has nonzero a_i, invariance forces every row or every
column constant, contrary to its nonuniform regular degree. Connectivity
therefore leaves only the all-zero vector or vectors with every coordinate
nonzero as possible fixed-point contributors.

The zero vector fixes product A_e configurations. When both endpoint shifts
are nonzero, the shift on p-by-p positions has p orbits of length p. Each
orbit meets every row and every column once. An invariant degree-d_e matrix
chooses exactly d_e of these orbits, giving B_e fixed matrices, independently
on the edges. There are (p-1)^v all-nonzero group elements. Burnside averaging
gives N_H.

The one-sided kernel lemma makes restriction of a stabilizer to any root
C_p injective along the connected support. Thus a configuration stabilizer
has order 1 or p. Double-counting pairs (configuration, nonidentity fixing
group element) gives (p-1)^v product B_e on one hand, and p-1 times the
number of symmetric labelled configurations on the other. Such a symmetric
orbit has size p^(v-1). This gives S_H; subtraction gives R_H. QED.

For a matching edge, A_e=p! and B_e=p. This extends the incoming
[two-cell Burnside formula](../session40_marked_ports_square_sum_20261023/PRIME_CELL_BURNSIDE_FORMULAS.md)
to arbitrary finite named supports, retaining the compatibility that separate
edge-orbit counts lose. For instance two R5 matching blocks along a three-cell
path have 128 classes, not 8*8=64: the middle chart is shared.

| Internal cell / support / matching edges | Classes | C_p-symmetric | Rigid |
|---|---:|---:|---:|
| C3, path on three cells | 4 | 4 | 0 |
| C3, triangle | 16 | 12 | 4 |
| C3, four-cycle | 32 | 24 | 8 |
| R5, path on three cells | 128 | 16 | 112 |
| R5, triangle | 13888 | 80 | 13808 |

The checker independently enumerates the actual group orbits, including all
120^3=1,728,000 labelled R5 matching triangles. The orbit sizes are 25 and125,
as predicted, and no edgewise canonicalization is substituted for the shared
action. It separately counts all regular-five binary blocks by row words.
For arbitrary p with the stated internal group and available A_e, the displayed
formula is a theorem, not an extrapolation of these small tests.

### The missing alignment is an exact linear quotient, even on a tree

Fix the actual local rotation-orbit of the block on EACH edge, rather than
summing over all local types. Choose a reference block within each such
orbit. Its possible translates form the affine space

    W_e = F_p^2 / H_e,

where H_e is zero for a rigid block and a slope line for a symmetric block.
This is a statement about the group acting on blocks; it does NOT claim
that a nonlinear matching permutation is an affine map of its owners.

The full edge-choice space is a torsor for W=direct-sum W_e. A single shared
rotation at every cell acts by the linear map

    D: F_p^V -> W,    (a_i) -> ((a_i,a_j)+H_e)_e.

Consequently **W / image(D)** classifies the possible joint attachments
with these fixed local types. Two systems are isomorphic with these named
cells exactly when the difference of their edge coordinates belongs to
image(D). This gives a finite-field elimination test retaining the shared
middle rotations, instead of choosing a different witness on every edge.

Let beta=|E|-|V|+1 and r be the number of rigid edges. The kernel of D has
dimension a=1 precisely when r=0 and every cycle's slope product is one;
otherwise a=0. Therefore the alignment quotient has dimension

    delta = beta+r-1+a,

and there are exactly p^delta joint classes above the given local edge
types. To prove the kernel assertion, a rigid edge pins both endpoint
rotations to zero and connectivity propagates zero. Otherwise all edge
relations are a_j=s_e a_i, so one free root value survives precisely when
every cycle gain is one. Rank-nullity proves the dimension formula.

In particular a TREE with r>=1 rigid edges has p^(r-1) distinct joint
attachments even though it has no graph cycle. Every further rigid
attachment contributes a relative rotation that an independent edge
quotient erased. For the R5 two-edge path, the sixteen pairs of rigid
local types each have FIVE different joint classes. The other48 pairs
of local types each have one. This explains128 versus8*8 concretely.

The same quotient identifies the precise new sidecar and its minimum fixed
alphabet: delta field coordinates, or at least ceil(delta*log2(p)) bits
for a fixed-length lossless binary encoding once the local types are known.
It preserves the named-cell class. Restoring a literal owner gauge or an
actual Hamiltonian port requires the retained global chart/marks and path,
just as before. It is not the integer residual-flow quotient in RESIDUAL.md.

The independent orbit traversal additionally groups every global class by
its local edge-orbit tuple and checks this exact fibre dimension. The R5
matching triangle's512 local signatures split as follows:16 balanced
all-symmetric signatures with5 classes each,48 unbalanced all-symmetric
with1 each,192 one-rigid-edge signatures with5 each,192 two-rigid-edge
with25 each, and64 three-rigid-edge with125 each.

## 3. A sharp nine-cell failure: connected does not imply functional

Let R9=C3[C3,C3,C3], with vertices (a,x) in F3^2. Between different a use
the outer C3, and inside each a use the inner C3. It is regular of degree4.
Couple two copies with

    B((a,x),(b,y)) = 1 iff a=b.

The block is I3 tensor J3, regular of degree3. Fixing the whole source chart
still leaves independent rotations in each of the three target triples:
a kernel C3^3 of order27. The informative coarse support is connected, but
chart propagation is not functional.

The exact two-cell census gives 81 internal automorphisms of R9, 2187
compatible chart pairs, and 27 target charts above EACH fixed source chart.
Individualizing every source owner leaves a target partition 3+3+3.
Pointwise recovery of the 18-vertex system needs six actual owners, one in
each internal triple. Six suffice, and five cannot kill the independent
rotation of an unmarked triple. This distinguishes source-chart ownership
from target-chart ownership sharply.

The mechanism scales. Take m named R9 cells and ANY connected support, using
the same equality-of-outer-coordinate block on every informative edge and
uniform blocks elsewhere. Independent rotations of all 3m inner triples
remain. An internal R9 automorphism consists of one outer cyclic rotation
and independent inner rotations. Each cross block forces only equality of
the outer rotations. Therefore the full cell-preserving automorphism group
has order 3^(3m+1), and its minimum actual-owner base is exactly 3m.

For completeness, the three inner triples are intrinsic inside R9: two
vertices in one triple have three common out-neighbours, while two in
different triples have one. Every automorphism therefore permutes these
three parts by an outer C3 rotation and acts internally by C3 rotations.
This proves the stated internal group without an assumed modular-decomposition
algorithm. Marking one owner per inner triple destroys all these rotations.

Recursive subdivision exposes the independent freedoms; it does not erase
them. At the leaf level all attachments between the C3 modules are uniform,
so the small-cell theorem correctly counts 3m separate informative components.
Nine is the first larger odd cell size after the proved 1/3/5/7 carrier, and
already prevents a constant budget based only on coarse support connectedness.

## 4. Prime size by itself is not the repair

There is also a regular ELEVEN-vertex cell with an internal C3 module M and
eight fixed vertices. Let H,L be strong four-vertex tournaments, each with
scores 1,1,2,2. Orient H->M->L. In the H-to-L block put just a matching
from the two score-one H vertices to the two score-two L vertices; reverse
all other H/L arcs. Every vertex now has outdegree5, and rotating M fixes
the other eight vertices.

Between two such cells use J3 direct-sum a regular degree-three circulant
8-by-8 matrix. All eleven row and column sums are3. The target rotation of
M fixes this block, so it is a nonfunctional attachment despite prime cell
size. The checker retains the literal adjacency and block rows. The
prime-network formula assumes the FULL internal automorphism group C_p;
one cannot apply it merely because p is prime or because a rotation chart
was convenient in a different family.

## 5. Shared ports, converse, and the remaining target

For an accepted global map f, transport the actual ordered path P to f(P).
The literal identity f(K_P(T))=K_{f(P)}(f(T)) preserves the complement partner.
Converse reverses both tournament arcs and the actual path; retained owners
move with that path. Cross-block degree d becomes p-d, preserving all kernel
and propagation hypotheses. The electrical modulo-eight checks remain
separate integer-wiring checks; cell chirality or group order is not a
substitute for them.

The order31 source remains the required positive control: its six genuine
order21 residual groups must not be 'resolved' without two owners. The
other 168 rigid but nondiscrete cases show why automorphism uniqueness and
ordinary refinement are different claims. The new theorem supplies a uniform
refinement mechanism AFTER the stated marks; it does not declare every
unmarked rigid partition already discrete.

The precise next targets are now larger-cell one-sided kernels and a compact
way to authenticate the required block information. General nine-cells need
not split into the particular C3 modules above. Full block availability is
still a high-information input, and uniformly small observation banks for
unrestricted metagraph wiring remain OPEN.
