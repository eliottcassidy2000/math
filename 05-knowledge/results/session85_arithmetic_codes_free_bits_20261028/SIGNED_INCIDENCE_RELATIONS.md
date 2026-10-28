# Free parameters, torsion, and signed attachment geometry

Session85, 2026-10-28. **PROVED elementary statements; FINITE-EXACT controls;
INHERITED graph-labeling results identified explicitly.** No priority claim.

## Results first

* A graph's sum/difference reconstruction has two different resources: one
  **integer parameter per balanced signed component**, and one **parity
  integrality obstruction per unbalanced component**. A free integer is not a
  free bit; a torsion bit is not a choice between two integer solutions.
* Suppressing a two-edge square-sum path creates a **difference equation**.
  The three suppressed-edge marks in Petersen-minus-one -> K3,3 must carry
  signs and right-hand-side offsets to preserve the integral reconstruction.
  Forgetting them removes a genuine Z/2 obstruction and invents a parameter.
* For higher-order sums, arbitrary torsion occurs: the k-uniform hypergraph
  of all k-subsets of k+1 vertices has cokernel Z/k. Fano point-line incidence
  has Smith factors `(1,1,1,1,2,2,6)`, not the Smith factors of its 3x7 binary
  parity-check matrix. Thus there is no universal “one missing parity bit.”
* Fano is still graph-encodable, but not as ordinary cycles with merely
  existential auxiliary edge variables. Such projection always remains
  graphic. Equality/parity factor-node laws, or a non-coordinate linear
  decoder, supply a precise successful encoding and its missing sidecar.

## 1. Recovered inheritance and scope

The following repository sources are pinned to commit
`83b601158076669edc8729ce461f150da1c3877e`:

1. [`05-knowledge/results/unifying_relations_20261028/THEORY.md`](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/unifying_relations_20261028/THEORY.md) §D: signless graph incidence,
   integral Smith factors, the triangle half-integrality hostile, K3,3 and
   Petersen examples. The graph-only sum theorem is **inherited**, not new.
2. [`05-knowledge/results/incidence_attractors_20261028/THEORY.md`](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/incidence_attractors_20261028/THEORY.md): the Fano prime bank, Hamming code,
   guarded arithmetic domain, continuation quotient, and Petersen sidecars.
3. [Session70 graph-labeling `REPORT.md`](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/session70_attachment_guard_alignment_20261026/session70-graph-labeling/REPORT.md) §§4–5: graceful overwrite quotient
   `S_k x {±}^k -> S_k`, Petersen spoke contraction, marked-ear suppression.
4. [**THM-2761 — graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge**](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md):
   sum/difference gauge and the indispensable mirror collision factor.

The current mistakes entry corrects “seven Fano circuits” to **fourteen**:
seven lines and seven four-point complements. Nothing here revives the old
minimum-weight-versus-minimal-support error.

Concept board: signed reconstruction / arithmetic realization guards /
marked graph reductions / binary code coordinates / higher-arity torsion.
Anchor: distinguish parameters from bits. Niche: signed ear elimination.
Wildcard: arbitrary torsion in the simplex-facet hypergraph.

## 2. The complete signed graph law

Let G be a finite loopless graph, choose an orientation for writing each edge,
and assign signs `sigma_e in {+1,-1}`. Define the **edge-by-vertex** map

```text
(A_sigma x)_e = x_v - sigma_e x_u,   e written u -> v.
```

A positive signed edge is a difference; a negative signed edge is a sum.
A component is **balanced** when every cycle's sign product is +1; isolated
vertices are balanced. Let n,m,c,b denote vertex count, edge count, component
count and balanced-component count; put `u=c-b`.

**Theorem 1 (classical signed-incidence law, proof included).**

```text
Smith nonzero factors: 1 repeated (n-c), then 2 repeated u;
ker_Z A_sigma = Z^b;
coker_Z A_sigma = Z^(m-n+b) direct_sum (Z/2)^u.
```

The same kernel dimension b holds over Q and every field of odd
characteristic. Over F2, signs disappear and kernel dimension is c.
Consequently a nonempty solution fiber modulo p has p^b elements for odd
prime p, but has 2^c elements modulo 2.

**Proof.** In each component choose a spanning tree and root. Tree propagation
writes `x_v=epsilon_v t+c_v`, where epsilon is the product of tree-edge signs
and each c_v is integral whenever q is integral. A nontree edge either imposes
an equation independent of t, or imposes `2t=h` with integral h. The latter
occurs exactly when its fundamental cycle is unbalanced. The homogeneous
kernel therefore has one parameter precisely in balanced components.
Deleting the root column from the spanning-tree rows gives determinant ±1,
so the first n_i-1 Smith factors in a connected component are 1. In an
unbalanced component every full-rank minor is even, since modulo 2 the
all-one vector is killed. A spanning tree plus an unbalanced edge has
absolute determinant 2: eliminate leaves and switch signs around its single
cycle. Hence the final factor is exactly 2. Direct sums over components give
the stated factors, kernels and cokernels. Reduction modulo p gives the
field claims. QED.

The vocabulary and rank law are classical; see Thomas Zaslavsky,
[*Matrices in the Theory of Signed Simple Graphs*](https://people.math.binghamton.edu/zaslav/Oldcourses/510.S18/mts.pdf),
§§1 and 4. The source uses vertex-by-edge incidence, the transpose of ours.
This note gives its own integer-minor proof and makes no historical novelty
claim for signed graphs.

All-negative signs recover the inherited square-sum operator B. Balance is
then exactly bipartiteness. All-positive signs recover oriented differences,
and every component is balanced regardless of odd cycles in the underlying
graph. “Odd cycle kills the free parameter” is therefore a statement about
**all-negative signs**, not about every incidence operator on that graph.

### 2.1 What one anchor does, and what it does not do

For specified consistent edge data q on a connected graph, one root value
reconstructs every label along a spanning tree. In a balanced component the
root value is freely integral; in an unbalanced component an odd-sign cycle
already fixes it, and its integrality must be tested. Multiple cycle equations
must agree. The Z/2 cokernel does **not** mean two integer lifts: when one
exists in an unbalanced component, it is unique.

For a balanced component write `x_v=epsilon_v t+c_v`. Positive integer labels
require

```text
max_(epsilon_v=+1)(1-c_v) <= t <= min_(epsilon_v=-1)(c_v-1).
```

Empty-sided maxima/minima mean no bound on that side. Injectivity additionally
requires distinct c_v within each sign class, and excludes each opposite-sign
collision `2t=c_w-c_v` with epsilon_v=+1 and epsilon_w=-1. Thus even a free
parameter can have no realizable value. Bounds `0<=x_v<=M` are likewise
linear interval bounds followed by these collision exclusions.

For square-sum labeling q itself must satisfy `q_e=r_e^2`, with integral
roots r_e and positive vertex labels; square roots need not be distinct.
Triangle q=(1,1,1) is rationally compatible but forces all labels to 1/2.
Triangle labels `(6,19,30)` instead give the distinct positive edge sums
`25,36,49`. Rational reconstruction alone settles neither integrality nor
positivity nor injectivity nor host coverage.

On an odd sum-cycle with successive edge data q_0,...,q_(2r), the
root is `(q_0-q_1+...+q_(2r))/2`. If q_i=r_i^2, its integrality obstruction
is exactly `sum_i r_i=0 mod2`: an even number of odd square roots. Once all
rational compatibility equations hold in a connected nonbipartite graph,
one chosen odd cycle supplies the only remaining parity bit; it does not
supply the positivity or distinctness inequalities. Even cycles instead
require their alternating sum of q to vanish.

For graceful labeling of an m-edge graph, choose a vertex injection into
`{0,...,m}` and require absolute edge differences to be exactly `{1,...,m}`.
Given magnitudes, one must retain oriented signs before solving differences;
cycle sums then constrain those signs. On a tree, an anchor plus assigned
signed differences reconstructs labels, but range and collision guards remain.
Labels `(0,1,2)` on a bipartite-oriented path give differences `(-1,+1)`:
oriented separation passes while absolute separation fails. This is exactly
THM-2761's missing mirror factor. Arbitrary edge-sign bits are not guaranteed
independent graceful choices. Session70's hidden 2^k sign fiber is a distinct,
**operation-relative** overwrite result, not the kernel of A_sigma.

### 2.2 A genuine gauge-bit saving, with transported guards

For a fixed underlying graph, switching vertex signs by eta_v changes
`sigma_uv` to `eta_u sigma_uv eta_v`. Exactly `2^(n-c)` distinct switches act
on the edge-sign assignment: the kernel consists of a constant sign on each
component. Fixing a spanning forest can make every forest-edge sign positive,
leaving exactly `m-n+c` independent chord signs. They are cycle-sign products,
so there are `2^(m-n+c)` switching classes. This is a genuine removal of
`n-c` **representational** sign bits, proved by tree propagation; it neither
adds solutions nor eliminates the remaining cycle consistency equations.

Writing `x_v=eta_v y_v` transforms the edge equation to
`y_v-(eta_u sigma_uv eta_v)y_u=eta_v q_e`. Right-hand-side signs, vertex bounds,
and label reconstruction must be transported too. Positivity of x becomes
`eta_v y_v>0`; a positive square q can become its negative in the new equation.
Thus this gauge does not identify arbitrary positive square-sum realizations.
Reversing the written orientation is a separate convention: it changes q_e
to `-sigma_e q_e`, not sigma_e. Neither convention supplies a free physical
edge-orientation choice after absolute labels and their order are specified.

## 3. Petersen maps as exact elimination, not label identification

Use u_i=i, v_i=i+5 for i modulo 5; outer edges have step 1, inner step 2,
and spokes join u_i to v_i.

### 3.1 Spoke contraction with offsets

Write the spoke equations `x_(u_i)+x_(v_i)=s_i`. Eliminating v_i gives
`x_(v_i)=s_i-x_(u_i)`. Outer edges keep their sum data; an inner edge i--j
becomes

```text
x_(u_i)+x_(u_j)=s_i+s_j-q_(v_i v_j).
```

The reduced graph is K5 with all sum edges. Integral row/column elimination
therefore identifies the source matrix with `I_5 direct_sum B_K5`, up to
unimodular equivalence. The same nonunit Smith data survive. This algebraic
elimination is **not** an instruction to identify the two original labels.
Reconstruction also requires `1<=x_(u_i)<=s_i-1`, all ten labels distinct,
and any requested interval or host condition.

### 3.2 Ear suppression transports a sign

A degree-two path u--w--v with sum data q_1,q_2 satisfies

```text
x_u+x_w=q_1, x_w+x_v=q_2
    iff x_w=q_1-x_u and x_v-x_u=q_2-q_1.
```

Thus an even-length sum path suppresses to a **difference** edge, with an
offset. More generally a path of length l transports gain `(-1)^l`, by
repeated substitution; this is signed elimination, not plain graph topology.

Delete vertex 0 from Petersen and suppress paths
`2--1--6`, `3--4--9`, `7--5--8`. The graph is K3,3 with parts
`{2,8,9}` and `{3,6,7}`. Keep the three resulting edges **positive**, the
other six **negative**, and the appropriate q differences.

| Edge-by-vertex operator | Nonzero Smith factors | Integer cokernel |
|---|---|---|
| Petersen, all sums | 1^9,2 | Z^5 + Z/2 |
| K5, spoke-eliminated | 1^4,2 | Z^5 + Z/2 |
| Petersen-minus-0, all sums | 1^8,2 | Z^3 + Z/2 |
| K3,3, three positive marked ears | 1^5,2 | Z^3 + Z/2 |
| K3,3, incorrectly reset to all sums | 1^5 | Z^4 |

So bare suppression both invents one rational parameter and removes the
integrality obstruction. The same three marked edges occur in the inherited
Hamiltonian lifting problem, but there the mark means “this edge must be used
to visit the suppressed vertex”; here it means a sign, offset, and label
reconstruction. One mark set supports different **typed** predicates, not an
automatic equivalence. The inherited counts—six K3,3 Hamilton cycles, two
spanning lifts, and six disconnected Petersen 2-factors—are cited, not
re-enumerated by this packet.

## 4. Beyond pair incidence: arbitrary torsion and Fano reconstruction

For any integer matrix M, rational nullity counts free rational parameters,
while its nonunit Smith factors give integral compatibility obstructions.
The fact that pair-incidence only produces 2-torsion is special to its two
unit entries per row. Higher-order incidence need not share it.

### 4.1 A uniform hypergraph with Z/k obstruction

On k+1 vertices take all k-subsets as hyperedges, one omitting each vertex.
Its hyperedge-sum matrix is `M=J-I`. For k>=1,

```text
q_i=sum_(j!=i) x_j;
x_i=(sum_j q_j)/k-q_i;
Smith(M)=(1 repeated k times, k);
coker(M)=Z/k; ker_Q(M)=0.
```

**Proof.** Summing equations gives `sum q=k sum x`, so integral realization
requires `sum q=0 mod k`, and the displayed inverse proves sufficiency.
The homomorphism `q -> sum q mod k` is surjective, with kernel exactly the
image of M. Thus the cokernel is cyclic of order k, proving the Smith form.
Equivalently its determinant has absolute value k. QED.

All q_i=1 force all x_i=1/k. No one-anchor trick removes this obstruction.
For k=3, the distinct positive labels `(58,41,22,1)` give square triple sums
`(64,81,100,121)`, so the hypergraph has genuine positive square realizations
as well as hostile formal data. At k=5 or 15 the obstruction is exactly Z/5
or Z/15; these are instances of **all k**, not a special numerological law.

### 4.2 Two different Fano matrices

Let H be the 3x7 binary matrix whose columns are vectors 1,...,7 in F2^3.
Let N be the 7x7 **integer** matrix whose rows are the seven Fano lines
`{a,b,a xor b}` and whose columns are its seven points. They are different
operators, even though their combinatorics are linked.

```text
H: integral Smith (1,1,1); Q-kernel dimension 4; F2-kernel dimension 4.
N: integral Smith (1,1,1,1,2,2,6); Q-kernel 0; F2-kernel dimension 3.
```

**Proof for N.** Any line has three points and distinct lines meet once,
so `NN^T=2I+J`. Hence `det(N)^2=2^6*9`, giving absolute determinant 24.
Modulo 2, its rows span the Hamming code: each line is in ker H, the sum of
all seven lines is the all-one word, and adding it to a line gives every
weight-four codeword. These sixteen words exhaust ker H. Thus rank_F2 N=4.
Exactly three Smith factors are even; their product and divisibility force
`(1,1,1,1,2,2,6)`. QED.

For prescribed line-sums q, put Q=sum_L q_L and
`s_p=sum_(L contains p) q_L`. The unique rational reconstruction is

```text
x_p=(3 s_p-Q)/6.
```

Indeed `N^-1=(1/2)N^T-(1/6)J`. Integral labels exist iff

```text
Q=0 mod3, and s_p=Q/3 mod2 for every point p.
```

The cokernel is `(Z/2)^3 + Z/3`: three independent binary obstructions and
one ternary obstruction, not a single bit. Modulo 3 the kernel has dimension
one (the constant vector); modulo 2 it has dimension three; over Q it has
none. Positivity and distinctness remain additional inequalities. Applying
N's Smith form to H would therefore change the actual mathematical object.

## 5. What graph auxiliaries can and cannot preserve

**Theorem 2 (projection obstruction).** Let C(G) be a graph's binary
edge-cycle space. Split its edges into retained S and auxiliary T. Existential
projection onto S is the cycle space of the multigraph obtained by contracting
each connected component of `(V,T)`; retained internal edges become loops.

**Proof.** An assignment on S extends across T iff its induced vertex boundary
has even total on each T-component. Necessity follows by summing the degree
parities. Sufficiency follows by solving on a spanning tree of that component,
eliminating leaves. The component parity equations are precisely the incidence
equations of the contracted graph. QED.

Therefore adding and then forgetting ordinary edge-cycle variables cannot
produce the Fano/Hamming code on its original seven coordinates: it would
still be graphic. Fano is not graphic because it is simple of rank 3 with
seven elements, whereas a simple graphic rank-three matroid has at most six.
This statement permits arbitrary numbers of auxiliary vertices and edges;
it excludes neither factor graphs nor general decoded graph encodings.

**Successful factor-graph encoding.** The Tanner graph of H has seven variable
nodes, three check nodes and twelve incidence edges. Put equality of incident
bits at each variable node and even parity at each check node. These laws
produce exactly the sixteen Hamming words, with unique incidence decoration
for each word. If instead every graph node is given the ordinary even-degree
law, the graph has cycle rank `12-10+1=3` and only eight assignments. Thus
“the graph plus typed node laws” is exact; the bare graph-cycle law is not.

**Successful linear decoder.** Four independent graph cycles provide bits
(a,b,c,d). The decoder

```text
(x1,x2,x3,x4,x5,x6,x7)
 =(a+b+d, a+c+d, a, b+c+d, b, c, d), additions mod2,
```

is a bijection onto ker H. For example use four triangles sharing one vertex
(n=9,m=12,cycle rank=4). This preserves the abstract linear solution space,
but not literal edge coordinates or support weights: d=1 alone selects one
three-edge triangle but decodes to four selected arithmetic primes. Retaining
the decoder restores those predicates. It is not a matroid isomorphism or a
coordinate projection, so it does not contradict Theorem 2.

## 6. Verification and boundaries

Run, using Python's standard library only:

```sh
python3 -B check_signed_incidence.py > normal.json
python3 -B -O check_signed_incidence.py > optimized.json
cmp normal.json optimized.json
```

All checks use explicit exceptions, not removable Python assertions. The
finite universe is all 760 signed simple graphs through four vertices, with
Smith factors computed from determinantal divisors, ranks checked modulo
2,3,5, and all 75 underlying-graph switching quotients by both forest
normalization and explicit vertex-switch orbit enumeration; the five displayed Petersen operators; all simplex cases k=1,...,6
by minors (k=15 by determinant and inverse); Fano's complete minor-gcd Smith
calculation and 2,187 right-hand sides in `{0,1,2}^7`; all 4,096 Tanner edge
assignments; and all 729 absent/retained/hidden edge assignments on four
vertices for Theorem 2. These finite controls complement the general proofs.
No global graceful-tree existence, square-sum Hamiltonicity, prime infinitude,
or graph-minor preservation of Hamiltonicity is claimed.
