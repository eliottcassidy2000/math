# Exact Boolean adjacency on disjoint-cycle sectors

Status: **PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED**.
The [independent referee](overnight10_20260906_boolean_independent_audit.md)
accepts the complete proof and independently rebuilds the native transitions,
Slater eigenvectors, and critical nullity. Its separate proved addendum treats
the two-state resource cap by an integer unimodular congruence.
No new theorem identifier is requested.
The all-order spectrum of the full Eulerian orbit graph remains **OPEN**.

## 1. Inheritance, precise consumer, and the recovered coordinate

The native graph is the simple isomorphism-class graph of **degree-even**, not
necessarily connected, simple graphs on exactly n vertices. An edge means that
some representative is changed into the other class by symmetric difference
with one triangle. Write its Boolean adjacency as B_n. Triangle multiplicities
and orbit volumes are not included in B_n.

The closest proved mechanism is
[THM-4073, exact cycle distance](../../01-canon/theorems/THM-4073-even-graph-diameter-layer-exact-cycle-distance.md),
especially its literal cycle-insertion path. The canonical hostile is
[THM-4078, Boolean noncommutation](../../01-canon/theorems/THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md).
The corrected near miss is MISTAKE-495's noncanonical fixed spanning-tree
image, retained in the audited
[Boolean consumer](overnight_hexagon_sep05_boolean_consumer.md).
The least-used relevant sidecar is the native parity index in the
[Peck note](overnight_hexagon_sep05_boolean_peck.md):
it forces zero modes but does not produce their eigenvectors or claim equality
with nullity. The
[flow note](overnight_hexagon_sep05_boolean_flow.md)
keeps the distinct orbit-uniform consumer and its actual edge conductance.

The concept board is: cycle-length multiset; isolated-vertex budget; labelled
move multiplicity; additive exterior adjacency; parity versus resonant kernel;
the diagonal potential and boundary edges lost by restriction. Targeted searches
for token graphs, additive compounds, exterior path spectra, Slater determinants,
and the proposed sector in current canon/results and the relevant mistakes found
no earlier repository statement. This is not a claim of external priority.

The source-target contract is exact. Restrict the native Boolean triangle graph
to c disjoint cycles of bounded length. Send the sorted lengths to distinct
positions on a path by adding their order indices. This preserves vertices,
Boolean adjacency, edge parity up to a fixed sign, and induced graph distance.
It discards labelled multiplicities and all edges leaving the class. A single
linear resource cap records the remaining isolated-vertex constraint. The cheap
decisive test is literal triangle XOR followed by component classification,
including repeated cycle lengths and a saturated resource boundary.

## 2. The exact sector and its resource cap

Fix integers c>=1, M>=3, n>=3c, and let S(n;c,M) be the subgraph induced by
classes whose nontrivial components are exactly c cycles, with lengths

```text
3 <= l_1 <= ... <= l_c <= M,        sum_i l_i <= n.
```

All other vertices are isolated. Put

```text
d=M-3,  N=c+d,  a_i=l_i-3,  b_i=a_i+i,
R=n-3c+c(c+1)/2.
```

**Theorem 1.** The map l -> {b_1,...,b_c} identifies S(n;c,M), with its
unweighted adjacency, with the induced subgraph of the c-token path F_c(P_N)
on precisely the configurations

```text
1 <= b_1 < ... < b_c <= N,        sum_i b_i <= R.       (1)
```

Here F_c(P_N) has c-element subsets of {1,...,N} as vertices; a move replaces
one occupied position by an adjacent unoccupied position. In particular, if
n>=cM the resource cap is inactive and

```text
S(n;c,M) is isomorphic to F_c(P_(c+M-3)).              (2)
```

This includes M=3: both graphs have just one vertex and zero adjacency.

**Proof.** A degree-even graph of maximum degree two is uniquely determined
up to isomorphism by its cycle-length multiset and its number of isolated
vertices. Consider a triangle T and let r be the number of its edges already
in the graph F. To keep degrees at most two after the toggle:

* If r=0, all three vertices must be isolated. A new C3 is created, changing c.
* If r=3, those edges form a component C3, which is deleted, changing c.
* If r=1, let uv be the existing edge. The third vertex w must be isolated,
  since it acquires two edges. Replacing uv by uw,vw inserts w into that cycle.
* If r=2, the existing edges form a path u-v-w in one cycle. Replacing them
  by uw isolates v and shortens that cycle by one. The cycle had length at
  least four: otherwise the third triangle edge would already be present.

These four cases are exhaustive. Thus a retained edge changes exactly one
cycle length by one. Every permitted shortening is realized by a consecutive
three-vertex path in that cycle. Every permitted growth is realized by an edge
of that cycle and an isolated vertex; if the target lengths still satisfy the
resource bound, at least one isolated vertex was available. No splitting or
merging of distinct cycles has been omitted.

After sorting, changing one length by one is exactly changing one a_i by one
while preserving weak order. At repeated values, choose the last equal entry
for an increase or the first equal entry for a decrease. Therefore adding i
to a_i turns precisely the allowed moves into adjacent unoccupied-position
moves of b_i. Different choices of identical cycles yield the same Boolean
edge, with no coefficient greater than one. The length sum becomes (1).
Its maximum over the full rectangle is cM, proving (2). QED.

The sufficient threshold n>=cM is the exact condition for this full rectangular
vertex set to be present. Below it, (1) is the strongest survivor of the same
map. For example, c=2,M=4 gives the path on the three states (3,3),(3,4),(4,4)
when n>=8, but gives only a two-state edge when n=7.

For any two retained length tuples l,l', the induced distance is also exact:

```text
dist_S(l,l') = sum_i |l_i-l'_i|.                     (3)
```

Every step gives the lower bound. Decrease the coordinates in increasing
index order to their componentwise minimum, then increase them in decreasing
index order to l'. This preserves sorted order; the first stage only removes
vertices and the second never exceeds the final length sum. It realizes (3)
even when the cap is active. In the full rectangle the induced diameter is
c(M-3). This is not a lower bound on the ambient graph's distance.

## 3. Complete adjacency spectrum on the full rectangle

Let A be the adjacency of P_N on the ordered orthonormal basis e_1,...,e_N.
On the exterior power define the **additive** operator

```text
dGamma_c(A)(v_1 wedge ... wedge v_c)
 = sum_i v_1 wedge ... wedge (A v_i) wedge ... wedge v_c.        (4)
```

This is not the multiplicative exterior matrix, whose eigenvalues are products.
In the ordered basis e_(b_1) wedge ... wedge e_(b_c), a term in (4) either
tries to occupy an already occupied position and vanishes, or moves one
position one step without crossing any other position. Every surviving sign
is therefore +1. Consequently (4) is exactly the Boolean adjacency in (2).

The ordinary path has orthonormal eigenvectors

```text
phi_j(b) = sqrt(2/(N+1)) sin(pi*j*b/(N+1)),
theta_j = 2 cos(pi*j/(N+1)),             1<=j<=N.
```

The sine addition formula checks the eigen-equation including the two endpoints;
the finite sine orthogonality identity gives completeness. Taking wedges of
distinct such eigenvectors now proves:

**Theorem 2.** For all c>=1, M>=3, n>=cM, the full adjacency spectrum of
S(n;c,M), counted with multiplicity, is

```text
{ 2 sum_(j in J) cos(pi*j/(c+M-2)) :
                    J subset {1,...,c+M-3}, |J|=c }.             (5)
```

An orthonormal eigenbasis has coordinates

```text
Psi_J(l) = det[phi_(j_s)(l_r-3+r)]_(r,s=1,...,c).                 (6)
```

There are exactly binom(c+M-3,c) mode subsets. Different subsets may have the
same eigenvalue; they still give independent orthogonal eigenvectors. The
spectral radius is the sum of the c largest theta_j. Complementing occupied
positions gives the graph symmetry F_c(P_N) ~= F_(N-c)(P_N), interpreted at
c=N as the one-state empty-configuration convention. No orbit-volume measure
or Fourier mode from the ambient Eulerian graph is used in this argument.

The framework is classical and explicitly **CITED**: Fabila-Monroy et al.,
[*Token Graphs*, Graphs and Combinatorics 28 (2012), 365-380](https://arxiv.org/abs/0910.4774),
defines this configuration graph; Mallory, Raz, Tamon and Zaslavsky,
[*Which Exterior Powers are Balanced?*, Electronic Journal of Combinatorics
20(2) (2013), P43, Section 3 and Fact 13, printed pp. 4-6 and 9](https://people.math.binghamton.edu/zaslav/Tpapers/xbal.pdf),
gives the additive exterior construction and the positive-sign path case.
The proof above is supplied in full. The repository contribution is the exact
native disjoint-cycle identification, its cap, and its explicitly typed consumers.

As a compact positive control, c=2,M=5 gives six classes and

```text
chi(lambda) = lambda^2 (lambda^2-1)(lambda^2-5).                  (7)
```

## 4. Explicit zero modes beyond the parity index

In token coordinates, native edge parity differs from sum_i b_i mod 2 by the
constant 3c-c(c+1)/2. Therefore for even N and odd c the two parity classes
are equinumerous: their signed difference, up to that constant sign, is

```text
[z^c] product_(j=1,...,N)(1+(-1)^j z)
 = [z^c](1-z^2)^(N/2) = 0.                                    (8)
```

**Theorem 3.** Let h>=3 be odd, c=3, M=3h-1, and n>=3M. Then the induced
native sector has equal parity classes and adjacency nullity at least h-1.

Here N=c+M-3=3h-1 is even. For each j=1,...,h-1 the distinct mode subset

```text
J_j={j,2h-j,2h+j}
```

has eigenvalue zero in (5), because, with x=pi*j/(3h),

```text
cos(x)+cos(2pi/3-x)+cos(2pi/3+x)=0.                            (9)
```

The subsets are different (their smallest member is j), so (6) gives h-1
orthogonal nonzero kernel vectors. This proves the all-height family; the
finite computation is only an independent control.

At the smallest member h=3, c=3,M=8, the graph has 56 vertices, split 28+28,
and its nullity is **exactly two**. Exact multiplication in
Z[zeta_18] shows that the only three-element zero-sum mode sets are
{1,5,7} and {2,4,8}; a separate integer characteristic-polynomial calculation
has zero multiplicity two. This fixed 56-state identification applies to
every n>=24 by Theorem 1. No minimality among all possible native subgraphs
is claimed, and no exact nullity h-1 is claimed for larger h.

This does not refute an equality conjecture for the full ambient B_n: the
sector's eigenvectors need not extend across its boundary. It provides a
literal native restricted class where parity alone misses a resonant kernel.

## 5. The lost data and exact repairs

The representative-level multiplicities in this sector are also explicit.
If m_s cycles have length s and I=n-sum_i l_i vertices are isolated, then a
growth s->s+1 has m_s*s*I labelled triangle realizations; a shortening s->s-1
has m_s*s realizations for s>=4. These counts apply when the endpoint remains
in the sector. The isomorphism-class volume is

```text
s(l) = n! / ( I! product_(s>=3) (2s)^m_s m_s! ).                 (10)
```

They satisfy detailed balance. At n=8,c=2,M=4 the directed retained
multiplicity matrix in order (3,3),(3,4),(4,4) is

```text
[[0,12,0], [4,0,3], [0,8,0]],
```

whereas the Boolean adjacency is the unweighted three-vertex path. The path
diagonalization above applies to the latter and does not silently retain (10).

There is also an exact **Laplacian** repair, rather than a free-sum inference.
Let L_path be the ordinary path Laplacian and let e(S) count path edges with
both endpoints in the occupied set S. Then for the full token graph,

```text
L_token = dGamma_c(L_path) - 2 diag(e(S)).                      (11)
```

The off-diagonal entries agree. The one-particle degree sum counts both ends
of every occupied-occupied edge, but neither is a legal token move; all other
incident edges give exactly one move. This proves the diagonal correction.
For F_2(P_4), the two characteristic polynomials, with coefficients listed
from highest to constant, are respectively

```text
L_token:             [1,-12,52,-100,84,-24,0],
dGamma_2(L_path):    [1,-18,128,-456,844,-744,224].
```

The second operator is not even singular. A further compression to the
resource cap in (1) requires its own boundary-degree correction.

Two more boundaries are decisive. First, at n=4,c=1,M=4 the sector is C3--C4,
but C3 has a triangle edge to the empty class. Its eigenvectors extended by
zero are not ambient eigenvectors. In general the principal compression of
the ambient Laplacian equals L_sector plus the diagonal count of edges leaving
the sector. Second, the triangle restriction is essential: at n=5 a four-cycle
toggle joins C3 directly to C5, which is a two-step separation in the triangle
path C3--C4--C5. Longer allowed generator layers change the operator.

## 6. Exact reproduction and next test

The standalone standard-library
[source](../../04-computation/overnight10_20260906_boolean_tokens.py) and
[output](overnight10_20260906_boolean_tokens.out)
contain no producer/repository imports, no floating point, no assert statements,
and no inferred labelled graph isomorphisms. They classify maximum-degree-two
graphs by literal BFS cycle components after each triangle XOR.

The complete literal universe is seven sectors/caps:
(n,c,M)=(6,1,6),(8,2,4),(10,2,5),(12,3,4),(24,3,8),(7,2,4),(9,2,5),
totalling **115,682 triangle toggles**. For every retained ordered pair the
script compares literal adjacency and multiplicity with the formulas above,
checks the token cap, and verifies orbit-volume detailed balance.
Independently, for all **36** pairs 1<=c<=N<=8 it compares integer matrix
characteristic polynomials with products of all mode sums in the exact
cyclotomic quotient ring. Resonant-family controls use h=3,5,7,9. Hostile gates
cover the resource cap, multiplicity, multiplicative versus additive exterior
operators, Laplacian diagonal, longer cycle layers, and ambient boundary.

```text
python 04-computation/overnight10_20260906_boolean_tokens.py
python -O 04-computation/overnight10_20260906_boolean_tokens.py

PASS 8338 always-active gates
normal and optimized output: byte-identical, LF
source SHA256:
985c110fe0e985d9c2c598e4aa35caf3f99756b4ca7cf00e487a41359bbb3356
output SHA256:
c35252acd8ab1cb4d98156570ad5c7281f476c21670c79e0f4e85fd9721ac4b3
```

The next bounded question is whether the exact resource-cap geometry in (1)
admits a replacement for the free-path diagonalization, or whether its
boundary destroys the useful resonances at the first possible cap. A separate
route is to retain the native length rank sum_i(l_i-3) and audit complementary
rank incidence maps on this sector, where the full labelled Peck-source
hostile no longer applies automatically. Neither statement is assumed here.
No full Boolean spectrum, gap, mixing law, or transfer to LRC/Smith/Laurent
problems has been inferred.

**Filing:** root integrated these audited artifacts in the tenth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
