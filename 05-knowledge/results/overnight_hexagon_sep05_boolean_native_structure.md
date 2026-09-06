# Native cycle creation, integral elimination, and complement-localized kernels

Status: **PROVED**, with independent full written audits by root and
nc2_seed. Exact small controls are **FINITE-EXACT**, not an all-order
rank claim.

## 1. Inheritance and the two intrinsic operations

The closest proved mechanisms are the native parity index in the
[Peck obstruction](overnight_hexagon_sep05_boolean_peck.md) and the complete
[self-complementary trace dictionary](overnight_hexagon_sep05_boolean_selfdual.md).
The canonical hostile is the n=4 weighted Fourier zero mode that is not
Boolean-null. The corrected near miss is to infer full native rank from a
local pivot, a weighted spectrum, or a parity count. The least-used
sidecars are **connected-component multiplicities** and **primal graph
complementation**, not an anchored coordinate order.

The live board is: Boolean multiplicity loss; intrinsic component creation;
integer pivots; complement characters; weighted versus native zero modes;
the residual matrix. No primary literature result is needed for the proofs
below. Exact searches for isolated-cycle creation or a unitriangular graph
minor in current canon and results found no matching statement.

These operations answer different obligations. Creation preserves actual
Boolean adjacency information in a triangular submatrix but leaves a large
uneliminated matrix. Odd-order complementation preserves the actual graph
and identifies a character containing the forced kernel index; it does not
rule out extra kernels in its other character. Neither operation supplies
the Boolean Laplacian gap.

## 2. An all-order unitriangular creation minor

Let `Z_n` be the binary Eulerian graph space, with isomorphism classes as
vertices. Let `L` be any nonempty set of cycle lengths contained in
`{3,...,n}`, and put `k=max L`. Let `B_L` be the **Boolean union** of their
toggle adjacencies, with graph loops either retained or deleted. Let `M_L`
be the raw sum of the corresponding multiplicity operators. Write `q_m`
for the number of Eulerian graph classes on m vertices, including `q_0=1`.

Let `D_k` be the classes with at least k isolated vertices. For `F in D_k`,
use k isolates to form a disjoint k-cycle, obtaining

```text
h(F)=F disjoint-union C_k, with k fewer isolates.          (1)
```

This is a well-defined injection of isomorphism classes. Indeed, deleting
any isolated `C_k` component from `h(F)` recovers the same class F, with k
isolates restored. In particular, `|D_k|=q_(n-k)`.

Let `c_k(F)` denote the number of connected components exactly equal to
`C_k`. Order `D_k` by decreasing edge count, then decreasing `c_k`, with
arbitrary remaining ties. Use this same order to index rows `h(F)` and
columns F. Then

```text
U=(B_L(h(F),F'))_(F,F' in D_k)
is lower unitriangular.                                  (2)
```

### Proof, including multiplicity and quotient boundaries

If toggling a cycle T of length `ell<=k` in `h(F)` produces `F'`, then

```text
|E(F')|=|E(F)|+k+ell-2|E(h(F)) intersect E(T)|
       >=|E(F)|+k-ell >=|E(F)|.                           (3)
```

Equality with `|E(F)|` requires `ell=k` and deletion of an actual k-cycle
already contained in `h(F)`. Deleting an isolated `C_k` component gives F.
Any other deletion leaves all `c_k(F)+1` isolated `C_k` components intact,
so its result has strictly more such components than F. Thus every
off-diagonal entry of the displayed minor is in an earlier column; no
same-key nonisomorphic exception is possible. At least the added isolated
cycle gives the diagonal entry, which Booleanization turns into one.

More precisely, the matching raw minor is lower triangular with diagonal

```text
M_L(h(F),F)=c_k(F)+1.                                    (4)
```

Only isolated k-cycle deletions can return the class F; shorter cycles
cannot, by (3), and a noncomponent deletion changes `c_k`. This proves (4)
after retaining every repeated labelled destination. Flattening within a
layer or across layers therefore does not compromise (2). Deleting graph
loops cannot change its diagonal because `h(F)` and F have different edge
counts; it can only remove off-diagonal triangular entries.

Consequently, over **every field**,

```text
rank B_L >= q_(n-k),
no nonzero kernel vector is supported entirely on D_k.    (5)
```

Every nonzero native kernel therefore reaches a graph with fewer than k
isolates. For the triangle graph it must reach a graph with at most two
isolates, and cannot be supported entirely on graphs with at most `n-3`
edges: every nonisolated vertex of an Eulerian graph has degree at least
two, so such an edge bound implies at least three isolates.

This is not a full unimodularity theorem. At n=6, creation from one isolated
triangle has raw diagonal two and Boolean diagonal one, so even their local
characteristic-two behavior differs. Conversely the complete n=7 Boolean
crossblock has determinant `1076=4*269`: the local integral pivot does not
imply global mod-two invertibility or global determinant one.

## 3. Unconditional integer elimination, rather than a guessed basis

Specialize to triangles. Write `B=[[0,C^T],[C,0]]`, with C mapping even
classes to odd classes. Use the even members of `D_3` as pivot columns,
and their created odd classes as pivot rows. After these permutations,

```text
C=[U V],       S=Z-W U^(-1)V,
  [W Z]
ker C={(-U^(-1)V y,y): S y=0}.                           (6)
```

Here U is unitriangular over the integers. Thus `U^(-1)` and S are integer
matrices, and (6) is an unconditional, denominator-free initial elimination
step. It parametrizes the entire even-side kernel once the residual S is
known; it does not assert that S has full row rank. Applying (2) separately
to odd pivot columns also gives

```text
rank C >= max(q_even(n-3),q_odd(n-3))=q_even(n-3),
rank B >= 2 q_even(n-3).                                 (7)
```

The equality of the maximum with the even count uses the independently
proved nonnegative parity index, including the trivial orders below three.
The [Fourier repair](eulerian-boolean-fourier-repair-overnight-hexagon-sep05.md)
uses a different, potentially much larger pivot and retains an additional
transversality hypothesis. Formula (6) neither assumes nor proves that
hypothesis. Its advantage is intrinsic integral control, not pivot size.

## 4. Odd-order complementation localizes the forced native kernel

This section works over **Q or R**, not in characteristic two. The
all-field assertion in Section 2 concerns only the unitriangular minor;
it is not being transferred to the character decomposition.

For odd n, the complete graph `K_n` is Eulerian. Therefore ordinary graph
complementation induces an involution J on the **primal Eulerian classes**.
It commutes with every cycle-toggle Boolean adjacency: complementing
commutes with symmetric difference and with vertex relabelling. This
statement retains actual edge support, so no weighted-to-Boolean transfer
is involved. It fails as a primal map for even n.

Let `n=4k+1`. Complementation preserves edge parity. Every complement-fixed
Eulerian class has exactly

```text
binom(n,2)/2=k(4k+1)
```

edges, and hence edge parity k. The proved selfdual dictionary identifies
their number as `Delta=q_even-q_odd`. For `sigma in {+1,-1}`, let
`V_even^sigma,V_odd^sigma` be the J-character subspaces. Since the trace of
a permutation involution counts its fixed objects,

```text
dim V_even^sigma-dim V_odd^sigma
 = [Delta+sigma(-1)^k Delta]/2.                          (8)
```

Thus every Boolean union of **odd** cycle lengths has at least Delta
independent native zero modes that vanish on odd-edge classes and satisfy

```text
f(complement F)=(-1)^k f(F).                             (9)
```

The forced modes can be chosen complement-symmetric when `n=1 mod8` and
complement-antisymmetric when `n=5 mod8`. The other complement character has
a square crossblock. It may still have a kernel, and the rectangular
character may have additional nullity. Equations (8)--(9) localize the
**forced index**, not necessarily the entire native kernel.

There is an explicit native reduced matrix behind this count. For every
free complement pair `{F,JF}`, use the basis `e_F+sigma e_(JF)`; retain a
fixed class only for `sigma=+1`. At an odd representative O, the column of a
free even pair is

```text
B(O,F)+sigma B(O,JF),                                   (10)
```

and a fixed even column is simply `B(O,F)`. These integer entries are
`{-1,0,1}` for the antisymmetric character and `{0,1,2}` for free symmetric
columns. Solving this smaller rectangular matrix produces the guaranteed
modes without any Fourier transversality assumption. This is a concrete
character reduction, not a closed formula for its arbitrary-order rank.

At n=5 the unique native kernel is the difference between the C4 and bowtie
classes: they are complement twins with identical Boolean neighbors.
Its empty and K5 coordinates vanish, and it is antisymmetric as (9)
requires. Both support classes have at most two isolates, consistent with
the creation boundary.

### Relation to the Fourier compiler

For odd n, translation by `K_n` acts on a Fourier orbit sum by
`J Psi_H=(-1)^|H| Psi_H`; switching preserves this parity because every
odd-order cut has even size. A self-complementary switching class has its
unique Eulerian representative from the selfdual dictionary, with the same
parity k above. Thus the fixed Fourier modes also lie in character
`sigma=(-1)^k`. Their Boolean defect lies in that same character.

Consequently, the Fourier correction may be restricted to the free-pair
block in that character. This block is square: the complete even Fourier
basis has Delta fixed-mode columns in this character, so its free-pair
columns number `dim V_even^sigma-Delta=dim V_odd^sigma`.
Invertibility of the opposite character's square
block is unnecessary to construct the Delta repaired modes. It would still
be needed, together with the relevant character's rank condition, to rule
out additional native kernels. No such all-order invertibility is inferred.

For `n=4k+3`, primal complementation instead exchanges the two edge
parities. Pair each even class with its odd complement; then

```text
A(F,G)=B(JF,G)
```

is a symmetric square matrix, since B is symmetric and commutes with J.
This does not imply A is nonsingular. It explains an additional actual
symmetry of the balanced native problem without importing a scalar
weighted spectrum.

## 5. Exact controls and stopping boundary

The [checker](../../04-computation/overnight_hexagon_sep05_boolean_native_structure.py)
uses the complete Eulerian orbit sets for n=3,...,7 and **every nonempty
subset** of their permitted cycle lengths. It literally toggles every
cycle from each created row, checks both edge/component inequalities,
retains the exact raw diagonal, and tests the entire Boolean creation
minor. The triangle integer Schur complement is reconstructed and its full
even-side kernel lifted back to the original matrix. Odd-order complement
actions and character dimensions are separately reconstructed on actual
Eulerian graphs.

Hostile controls retain the raw diagonal-two/Boolean-one distinction at
n=6, the absence of an even-order primal complement map, the n=5
antisymmetric twin mode, and the nonunimodular full n=7 crossblock. The
inherited orbit engine already has independent literal-permutation and
Burnside controls; this checker adds a different component/creation
predicate, not another larger graph census.

```bash
python3 -B 04-computation/overnight_hexagon_sep05_boolean_native_structure.py
python3 -B -O 04-computation/overnight_hexagon_sep05_boolean_native_structure.py
```

The [output](overnight_hexagon_sep05_boolean_native_structure.out) records
the exact universe and all finite ranks. Full native nullity, universal
Fourier transversality, and the Boolean Laplacian gap remain **OPEN**.

Normal and optimized executions agree on all 21,295 optimization-live
gates and 57 complete cycle-length families. Source SHA256:
`ad4e32404d1f24fd568f6d55c0ad78f161c0b03b3b4954656e55b509c034a5fa`;
output SHA256:
`11fed612bf11c8d603b019849e15e553cf5f5fe21a5e8feb10a6c0aa0d8afd59`;
semantic trace:
`129349f1675c743935e121ec5717409afb12c12006889c5898d69befcbca2025`.
Root's complete written audit passed, including the arbitrary maximum-cycle
extension and the distinction between all-field pivots and characteristic-
zero complement subspaces. The independent full written nc2_seed audit
also passed the maximum-cycle extension, integer kernel lift, and restricted
Fourier character-block dimensions; no mathematical repair was requested.
