# Boundary support, bridges, and realizability in the tournament line metagraph

**2026-10-20. PROVED scoped statements below; FINITE-EXACT bridge replay
through order seven and degree data through order eight.**

This note concerns the **converse-merged tiling line multigraph**, whose
lines pair a fixed-path tournament presentation with the complement of all
its free arcs. It does not replace that operation by a single arc reversal.
Loops count twice in degrees; parallel lines remain separate edges.
Throughout, the tournament order is `n>=3`, so the tile complement is
fixed-point-free.

The user's clarification identifies the tournament metagraph, but the
symbols `D(z)` and `M_k` have not been identified uniquely in its sources.
The main session explicitly adopts `D_n(z)=sum_v z^(deg v)` and
`M_n=deg D_n`, with loops counted twice. Section 6 gives this declared
interpretation; it does not claim the notation was recovered. The
[companion asymptotic note](boundary_support_20261020_asymptotics.md)
handles the growth of `M_n/n!`.

## 1. Inheritance and the common support kernel

The relevant inheritance is:

- [THM-643, blue/black parity](../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md):
  ordinary fibres `h(T)=H(T)/|Aut(T)|` are odd; blue tilings occur precisely
  at self-converse classes, in odd positive counts.
- [THM-796, three-sorted incidence](../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md),
  especially sections 2, 3, 6, and 8: tilings are the actual half-edges;
  node/line incidence and its reversible kernel are already defined;
  two faces must agree on their overlap; bare node transition matrices
  lose the distribution inside a fibre.
- [THM-4731, sigma-twins](../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md):
  every black line has a distinct reflected twin with the same merged
  endpoints. The whole black multigraph is `2B`.
- [THM-840, continuation congruence](../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md):
  a quotient supports an operation precisely when equal observations have
  equal updated observations; its linear repair has an exact rank.
- The [growth packet](metagraph_growth_20261020/README.md) and
  [allocation carrier](weighted_core_transport_20261020/allocation.md)
  show concrete failures of node degrees, block words, and run lengths
  to determine complement partners. Clone allocation is essential.

Closest proved mechanism: THM-796's incidence span. Canonical hostile:
two presentations of one order-four class have different complement
classes. Corrected near miss: the older proposed pure-blue alphabet has
an odd-order proof gap, documented in the
[nonaffine-palindrome inheritance audit](weighted_core_transport_20261020/inheritance.md);
no all-order classification from that proposal is used here.
Least-used sidecar: the *full complement target row* of each
presentation, and the black connected-component quotient.

The board is: **boundary support / involution / fibre loss / black
components / parity / realization**. Anchor: exact bridges and wiring.
Niche: measure the information lost by node compression. Wildcard:
transport the earlier Pythagorean norm construction to actual line support.

For `n>=3`, let `X=F_2^m`, `m=binom(n-1,2)`, be the fixed-path tilings, `tau` their
all-bit complement, and `sigma` staircase reflection. Write `F_i` for
the merged node fibres, `f_i=|F_i|`, and `J` for the matrix with columns
their indicator vectors. If `P` is the permutation matrix of `tau`, then

    A = J^T P J,
    A_ij = |F_i intersect tau(F_j)|.                    (1)

For distinct nodes `A_ij` is exactly their line multiplicity; `A_ii`
counts a loop twice. Thus `A*1=f`. This is THM-796's existing kernel in
matrix form, not a new definition of the metagraph.

This has exactly the same contraction mechanism as the cubic boundary
formula in
[THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md):
two nonnegative state vectors have compatible support iff their inner
product is positive. Here the state is a tiling, and the operation is
`tau`; for a cubic pole it is a labelled boundary colouring and a wiring
permutation. The preserved predicate is **existence of compatible states**.
The map does not turn tournament classes into cubic graphs or transfer
the Petersen-minor theorem.

## 2. An exact measure of the forgotten fibre coordinate

Put `D=diag(f_i)` and `U=J D^(-1/2)`. Its columns are orthonormal. Let
`K=U^T P U` be the symmetric normalized incidence matrix. Then

    I-K^2 = U^T P (I-UU^T) P U >= 0.                   (2)

Equivalently, entirely over the rationals,

    D-A D^(-1) A = R^T R,
    R=PJ-J D^(-1) A.                                 (3)

**Proof.** `P^T=P`, `P^2=I`, and `U^TU=I`; expanding the right side of
(2) proves it. Expanding `R^TR`, with `J^TJ=D` and `J^TPJ=A`,
gives (3). Both right sides are Gram matrices. ∎

The defect is zero iff `P` preserves the space of functions constant
on fibres. This is equivalent to `tau` mapping every whole fibre onto
one whole fibre: every indicator `1_(tau F_i)` must be constant on
every fibre; the resulting map is an involution. It is exactly
THM-840's continuation condition, now with a positive semidefinite
certificate of its failure.

**Proposition 0 (all-order defect rank).** For any finite nonempty
fibres and involution pairing as in (1), let `r` be the number of
nodes, `c` the number of connected components of the full support
multigraph, and `b` the number of those components that are bipartite.
Loops make their component nonbipartite. Then

    rank(D-A D^(-1) A)=rank(I-K^2)=r-c-b.              (3a)

No assumption that the tournament metagraph is connected is needed.

**Proof.** By (3), the nullspace of the defect equals `ker R`.
For a vector `a` of node coefficients, a row of `Ra` at a tiling
in `F_i` is its complement target's coefficient minus the weighted
average of the coefficients at all neighbors of `i`. Every neighbor
with positive incidence occurs for at least one such tiling.
Thus `Ra=0` iff all neighbors of each node have one common coefficient.
Equivalently, coefficients are constant along every walk of even
length. In a connected bipartite graph the two shores are precisely
the two even-walk equivalence classes. In a connected nonbipartite
graph an odd closed walk changes path parity, so there is one such
class. Every node has positive degree because its fibre is nonempty.
Consequently the nullity is `c+b`, proving (3a). Equivalently, `K`
has one eigenvector of eigenvalue `1` per component, and one of
eigenvalue `-1` per bipartite component. ∎

The complete line censuses give the following exact values:

| n | nodes r | components c | bipartite components b | defect rank r-c-b |
|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 1 | 0 |
| 4 | 3 | 1 | 0 | 2 |
| 5 | 10 | 1 | 0 | 9 |
| 6 | 34 | 1 | 0 | 33 |
| 7 | 272 | 1 | 0 | 271 |

The script checks (3) entrywise and its rational matrix rank at
orders `3,4,5` using an independently authored tiling canonicalizer.
At orders six and seven the exact graph connectivity and bipartiteness
tests feed the proved formula (3a); no floating-point eigenvalue
threshold is used. An all-order theorem on the connectivity of the
specific tournament metagraph remains OPEN, but the defect rank
formula itself is unconditional.

There is also a literal Pythagorean chart. For two fibres, set
`u=1_(F_i)`, `v=P 1_(F_j)`, and `c=<u,v>=A_ij`. The tensors
`u tensor v` and `v tensor u` have equal squared norm `f_i f_j` and
inner product `c^2`. Their sum and difference are orthogonal, with
squared lengths

    2(f_i f_j+c^2),  2(f_i f_j-c^2).                   (4)

The triangle can be degenerate when the two indicator vectors coincide.
It is nondegenerate and isosceles exactly when the nodes have **no connecting
line**, because then `c=0`. This makes the earlier norm analogy exact
for tournament nonadjacency as well as cubic colouring obstruction.
It also exposes its limitation: the same right-angle construction
represents many different support zeros, not a unique Petersen
obstruction. The lengths in (4) need not be integers or primitive.

## 3. Classification of every bridge, with no leaf assumption

Let `G` be the merged line multigraph. Let `C_1,...,C_r` be the connected
components of its black support graph, including black-isolated nodes.
Contract each `C_i` to one vertex and retain every blue line separately;
call the resulting multigraph `Q`. Keep blue loops and parallel blue
edges created by contraction.

**Proposition 1 (complete bridge criterion).**

1. No individual black line is a bridge.
2. A blue line `e` is a bridge of `G` iff its image is a bridge of `Q`.
3. Equivalently, `e` is the sole line crossing a partition of the
   nodes into unions of black components.

**Proof.** A black line has a distinct parallel twin, so deleting it
does not separate its endpoints. A path in `G-e` projects to a walk
in `Q-e`. Conversely, a path in `Q-e` lifts to a walk in `G-e` by
joining successive blue endpoints inside their connected black
components. Therefore the endpoints of `e` remain connected in one
graph exactly when they do in the other. The cut characterization
of a bridge gives statement 3. ∎

This is an all-order, necessary-and-sufficient **algorithmic**
classification; it does not assert that the contracted graph always
has a particular shape.

There is a useful parity sidecar. Every self-converse merged node has
odd full degree; every paired node has even degree. If a bridge separates
one connected component into shores `S` and `S^c`, then

    #{self-converse nodes in S} = 1 mod 2,             (5)

and the same holds on the other shore. Sum the degrees on a shore:
internal lines contribute twice and the bridge contributes once.
After black contraction, a quotient vertex has blue degree congruent
to the number of self-converse nodes it contains.

**Non-leaf hostile.** Take six self-converse-labelled abstract nodes
`0,...,5`, blue lines

    0--2, 1--2, 2--3, 3--4, 3--5,

and two parallel black lines at each of `0--1` and `4--5`.
Every full degree is three. The unique bridge is `2--3`; both ends
have degree three. This wiring has an explicit commuting-involution
realization by section 4. Thus sigma-twins, odd fibre parity, and the
local colour rules alone cannot prove that all bridges meet leaves.
It is an abstract hostile, not a claimed tournament metagraph.

### Exact bridge census

The script reads THM-4731's committed full class and line TSVs and
checks bridges three ways: lowlink, deletion/component counting, and
black contraction. Independently generated tiling cubes through order
six agree in node count, line count, bridges, and degree distribution.

| n | merged nodes | line instances | black components | all bridges | non-leaf bridges |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 2 | 1 | 0 |
| 4 | 3 | 4 | 2 | 1 | 0 |
| 5 | 10 | 32 | 4 | 2 | 0 |
| 6 | 34 | 512 | 3 | 2 | 0 |
| 7 | 272 | 16384 | 5 | 3 | 0 |

In particular, a pure-blue node need not contribute a bridge: the
pure-blue fibre of size three at order five does not. Its lines are
not interchangeable with the size-one leaf fibres.

At order eight, the committed class TSV has three pure-blue classes,
all with fibre one. THM-4731's independently audited finite census
states that all remaining nodes lie in one black component. Those
two inherited finite facts imply exactly three leaf bridges and no
others. The order-eight line TSV is not committed, so this note does
not call that conclusion an independent replay of its wiring.

## 4. Sufficient criteria for abstract complete wiring data

Here “complete wiring data” means specifying every coloured edge
multiplicity, rather than demanding a complete simple graph.

Take a finite merged node set `V`, a subset `S` of self-converse nodes,
a blue multigraph `L` supported on `S`, and a black multigraph with
every multiplicity even, written `2B`. Loops are allowed. These
conditions are **necessary and sufficient** for realization by a
finite set `X` with commuting involutions `sigma,tau`, with:

- `tau` and `sigma*tau` fixed-point-free;
- blue flags precisely `Fix(sigma)`;
- `sigma` fixing the class label at `S` and exchanging two ordinary
  class sheets at each node outside `S`;
- the prescribed merged blue and black line multiplicities.

**Construction.** Give each blue line two flags; `tau` swaps them,
and `sigma` fixes both. For each edge of `B`, give four flags
`(end,sheet) in {0,1}^2`; `tau` changes `end`, while `sigma` changes
`sheet`. At an endpoint in `S`, both sheets carry the same ordinary
class label. Outside `S`, they carry the two converse class labels.
For a loop, the two ends remain distinct flags even though they
have the same merged node. Take the disjoint union of these orbits.
All involution laws and all multiplicities are immediate. Necessity
is the sigma-twin argument. ∎

If prescribed ordinary fibre sizes are `h_v`, with prescribed blue
counts `b_v` for `v in S`, the additional exact conditions are

    deg_L(v)=b_v,
    deg_B(v)=(h_v-b_v)/2        for v in S,
    deg_L(v)=0, deg_B(v)=h_v    for v not in S.         (6)

Here degrees count loops twice. To match the tournament parity
constraints, require all `h_v` odd and `b_v` positive odd at `S`.
To match a staircase's global involution orbit counts, require

    |X|=2^m,    |Fix(sigma)|=2^e,
    e=(m+floor((n-1)/2))/2,    m=binom(n-1,2).        (7)

With (7), the constructed action is isomorphic **as a finite
Klein-four set** to the staircase action: both have `2^(e-1)`
two-point orbits fixed by `sigma` and `(2^m-2^e)/4` free
four-point orbits. No other stabilizer occurs because `tau` and
`sigma*tau` are free.

This is a genuine constructive sufficient criterion for the
involution model. It does **not** imply that the class labels arise
from tournaments on that staircase.

### Actual fibres and all parity laws still do not suffice

Use the exact order-five class labels in the inherited TSV. Keep every
black line, but replace blue edges

    (0,11), (2,6)    by    (2,11), (0,6).             (8)

Every blue degree, every black multiplicity, every ordinary fibre
size, every node type, and both global counts are unchanged. The
construction above realizes the altered data.

Yet class `11` is the transitive tournament. Its unique tiling forces
its unique partner to be class `0`, the tribonacci partner with `H=9`
from THM-4731(6). The altered target class `2` has `H=13`. Therefore
(8) is not the actual wiring of these tournament nodes. This is a
stronger hostile than merely assigning nonexistent fibre sizes.

The missing coordinate is the concrete presentation and its
complement, not another degree parity.

## 5. A sufficient certificate for actual tournament wiring

There is an exact finite certificate, though no small degree-only
criterion is established.

For each proposed ordinary node supply a representative tournament
`T_i`. Supply a partition of the whole fixed-path tiling cube into
fibres `F_i`, with an explicit vertex permutation identifying every
presentation `t in F_i` with `T_i`. Supply isomorphism-separating
canonical codes for the representatives. Then compute each pair
`{t,tau(t)}` once and compare its endpoint labels and blue test with
the proposed line table.

These finite checks suffice: the permutations prove membership,
canonical codes prove that no distinct nodes were spuriously
separated or merged, the partition proves coverage, and the
involution pairs prove every line exactly once. Conversely, every
actual wiring supplies such a certificate.

THM-796's compatible-face pullback gives an inductive version:
generate only pairs of face presentations agreeing on the common
interior, add the apex bit, and retain the induced tournament
isomorphism witness. Agreement of node labels without agreement of
the presentations is insufficient. The earlier allocation carrier
supplies a compressed certificate for chain-inflated families, with
clone assignments retained.

This certificate is larger than the desired compact structural
criterion. Finding an operation-stable compressed replacement is
still OPEN; the hostile (8) states a concrete requirement on it.

If “complete metagraph” instead means complete simple underlying
graph, that possibility is excluded at every `n>=4`. The transitive
node has degree one, while there are at least three merged nodes:
transitive, its distinct self-converse tribonacci partner, and a
non-self-converse node. A non-self-converse class exists by taking
`TT_(n-3) -> C3`; its ordered strong-component list differs from
its converse. For `n=3` the line metagraph is `K2`.

## 6. The declared degree polynomial and its exact computation

Define, explicitly for this section,

    D_n(z)=sum_(merged nodes v) z^(deg v), M_n=deg D_n. (9)

For ordinary classes `T`, let `h(T)=H(T)/|Aut(T)|`.
If `a_n` is the ordinary class count and `s_n` the self-converse
class count, then

    D_n(z)=sum_(T self-converse) z^h(T)
          +(1/2)sum_(T non-self-converse) z^(2h(T)).
                                                               (10)

Thus its coefficients are exactly the two fibre-size histograms:
odd powers count self-converse fibres, powers congruent to two
modulo four count converse pairs, and positive powers divisible
by four have coefficient zero. In particular,

    D_n(1)=(a_n+s_n)/2,
    D_n'(1)=2^binom(n-1,2),
    D_n(-1)=(a_n-3s_n)/2.                             (11)

The exact degree is

    max(max_(SC T) h(T), 2 max_(non-SC T) h(T)),       (12)

with absent strata omitted. Formula (12) is a characterization in
terms of fibre extrema, not an all-order closed formula for those
extrema.

For example,

    D_3=2z,
    D_4=z+z^2+z^5,
    D_5=2z+z^2+z^3+z^5+2z^9+z^10+z^11+z^13.

The degrees for `n=3,...,8` are

    1, 5, 13, 86, 306, 1322.

The degree maximum is not the Hamiltonian-path maximum. The exact
calibration is:

| n | max H(T) | M_n | `(H,Aut,SC)` of degree maximizers | number of merged maximizers |
|---:|---:|---:|---|---:|
| 3 | 3 | 1 | `(3,3,yes)`, `(1,1,yes)` | 2 |
| 4 | 5 | 5 | `(5,1,yes)` | 1 |
| 5 | 15 | 13 | `(13,1,yes)` | 1 |
| 6 | 45 | 86 | `(43,1,no)` | 1 |
| 7 | 189 | 306 | `(153,1,no)` | 1 |
| 8 | 661 | 1322 | `(661,1,no)` | 2 |

The output includes the exact inherited class identifiers of every
maximizer. The ordinary maximum `H=15` at order five loses the degree
comparison because of its automorphisms; at orders six and seven the
factor of two for a converse pair also matters.

All coefficients through order eight are emitted in the companion
output, with (11), the modulo-four support rule, and the independent
degree-distribution controls checked exactly.

There is an exact all-order finite formula that does not presuppose a
list of isomorphism classes. Sum over all `2^binom(n,2)` labelled
tournaments on `[n]`; let `g(T)=|Aut(T)|` and `epsilon(T)` be the
self-converse indicator. Then

    D_n(z)=(1/n!) sum_(labelled T) g(T) *
      [epsilon(T) z^(H(T)/g(T))
       +(1-epsilon(T))/2 * z^(2H(T)/g(T))].            (13)

An ordinary class has exactly `n!/g(T)` labelled representatives,
which proves (13). Hence every coefficient is an explicit finite
integer-weighted count after division by `n!`; the exponents are
integers by the free action of `Aut(T)` on Hamiltonian paths.

For completeness, `H(T)` in (13) is given by the subset recurrence

    p({v},v)=1,
    p(S,v)=sum_(u in S\{v}, u->v) p(S\{v},u),
    H(T)=sum_v p([n],v).                              (14)

Enumerating all vertex permutations computes `g(T)` and determines
whether any permutation is an anti-automorphism, hence `epsilon(T)`.
Equations (13)-(14) determine `D_n` and `M_n` for every fixed input
`n`. They are expensive exact algorithms, not a closed evaluation
of the all-order extrema. Existing class BFS and canonicalization
from the growth packet are more efficient alternatives for the
finite tests.

## 7. Reproduction and remaining boundaries

Run from the repository root:

```text
python3 04-computation/boundary_support_20261020_metagraph.py
python3 -O 04-computation/boundary_support_20261020_metagraph.py
```

The [script](../../04-computation/boundary_support_20261020_metagraph.py)
and [literal output](../../04-computation/boundary_support_20261020_metagraph.out)
agree under normal and optimized Python. It uses NetworkX and SymPy.
The full universes and source tables are declared at the top; the
independent path-cube checks use the previously audited producer in
`metagraph_growth_20261020/enumeration.py`.

Proved here: the bridge quotient criterion; the exact quotient-loss
Gram identity and all-order graph-theoretic defect rank; an explicit abstract realization criterion and its
actual-tournament failure boundary; the conditional degree-polynomial
identities and all-order finite evaluation formulas. Finite here: all bridges through order seven, all degree
coefficients through order eight, and the concrete hostiles.

Still open: a structural all-order description of the blue quotient;
whether actual non-leaf bridges occur beyond these censuses; a compact
criterion forcing actual tournament wiring; a closed formula for
general degree extrema. The adopted `D_n,M_n` convention remains
explicit throughout; it should be replaced if the user intended
different invariants.
