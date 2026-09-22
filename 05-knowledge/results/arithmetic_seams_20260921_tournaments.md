# Pointed products, the 36 chord coordinates, and two cycle laws

**Status: PROVED elementary identities, cycle classifications, and their
stated consequences; CITED odd-cycle collection formula; FINITE-EXACT
replay.** These are scoped graph/arithmetic connections, not a Collatz
convergence theorem. No priority claim is made.

## Inheritance and working board

The closest mechanisms are [THM-781, Hamiltonian-path inversion](../../01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md),
which retains the chosen path and its chord bits, and [THM-410,
interval-reversal triangle count](../../01-canon/theorems/THM-410-interval-reversal-triangle-count.md).
The [earlier scaffolding audit, sections 2 and 5](collatz_mod6_20260917_scaffolding_audit.md)
already repairs a related triangular identity and counts path presentations.
Its formula subtracts T(AB-1); the present formula contains T(AB), so their
error terms differ. Neither historical formula should silently replace the other.

The canonical hostile is the transitive tournament: it has the full
undirected cycle space and no directed cycles. The corrected near miss is
counting free chords as directed triangles. The least-used sidecar is the
vertex-overlap pattern of odd directed cycles.

| Concept | Exact map / invariant | Information required |
|---|---|---|
| Triangular addition | unordered pairs in two blocks | block sizes |
| Pointed multiplication | root, two axes, interior | whether the root is retained |
| Fixed Hamiltonian path | chord bits and fundamental cycles | full orientation and chosen order |
| Directed triangles | quadratic polynomial in chord bits | pairwise interaction terms |
| Hamiltonian paths | weighted disjoint odd-cycle collections | all odd lengths and overlaps |
| Planarity / minors | underlying undirected graph | orientation is discarded |

Anchor: give the user's 9+36 arrangement a faithful graph interpretation.
Niche: turn its free chords into exact cycle-generating laws. Wildcard:
locate the Fermat numbers in one explicitly specified family.

## 1. Triangular addition and the product's missing root

Write T(t)=t(t+1)/2. For nonnegative integers A,B,

```text
T(A+B)=T(A)+T(B)+AB,
T(A+B+1)=T(A)+T(B)+(A+1)(B+1).                         (1)
```

For the first identity, take two vertex blocks of sizes A,B and a shared
root: the two rooted complete graphs supply T(A),T(B), and the cross
edges supply AB. For the second, use disjoint blocks of sizes A+1,B+1.
The first cross term is the additive cocycle
`AB+(A+B)C=BC+A(B+C)`, expressing that three-block pair counting does
not depend on the order of grouping.

Now let X,Y be pointed sets of sizes A+1,B+1 and N=(A+1)(B+1).
The product X times Y splits into the root, the A-axis, the B-axis,
and an interior of size AB. The complete graph on these N vertices has

```text
T(N-1)=T(A)+T(B)+T(AB)+AB(A+B+1).                      (2)
```

The three T terms count complete graphs on root-plus-one-block; they
share a vertex but no edges. The remaining cross-block pairs number
`A*B+A*(AB)+B*(AB)=AB(A+B+1)`. This is an edge-partition bijection,
not just an equality after expansion.

The user's proposed cross term is `A(B^2-1)+B(A^2-1)`.
Consequently its entire right-hand side is **exactly**

```text
T(A)+T(B)+T(AB)+A(B^2-1)+B(A^2-1)=T(N-2).             (3)
```

It counts the edges after deleting the distinguished product root.
The difference from (2) is `AB+A+B=N-1`, precisely the deleted root's
star. The terms of (3) are an algebraic rearrangement; its individual
T terms are not the three rooted subgraphs after that root is deleted.
Smallest positive witness to the unshifted claim: A=B=1 gives 6 versus 3.
Equation (3) also holds at A=0 or B=0 with T(-1)=0 at the empty boundary;
individual proposed cross terms can then be negative.

There is a second useful edge partition. If a=A+1,b=B+1,

```text
C(ab,2)=a*C(b,2)+b^2*C(a,2).                           (4)
```

An unordered pair in the product has either the same first coordinate
(a fibres) or two distinct first coordinates (b^2 choices per pair).
This is the counting substrate of the directed substitution law below.

## 2. What nine path edges and 36 chords really encode

Fix the directed Hamiltonian path `0->1->...->n-1`. Exactly

```text
C(n,2)-(n-1)=C(n-1,2)                                  (5)
```

arc orientations remain free. At n=10 this is nine prescribed arcs
and 36 independent binary choices. A repeated numerical label 3 does not
choose either direction of an arc; an orientation rule is still needed.
This is a conditional coordinate system with the path order specified.
Changing to an unmarked tournament forgets that order and its presentation
multiplicity; THM-781 describes the resulting fibres.

Each chord {i,j}, j>=i+2, together with the path segment i,...,j, is a
fundamental **undirected** cycle. These form a basis of the cycle space
over F2: each contains one different non-tree edge, and any even-degree
edge set can have all chords eliminated uniquely. Their lengths are
`j-i+1`, so length ell occurs n-ell+1 times, ell=3,...,n.
At n=10 the multiplicities are `8,7,6,5,4,3,2,1`, summing to 36.
Only eight basis cycles are triangles.

A fundamental cycle is directed precisely when its chord points backward.
This criterion does not classify the other directed cycles: sums in the
undirected cycle space need not be single cycles or retain orientation.
All-forward chords give a transitive tournament with no directed cycles.
Thus 36 is the cycle-space dimension, not a directed-cycle count.

## 3. Exact triangle interactions extend the matching-reversal theorem

For i<j let b_ij be 1 if j->i and 0 if i->j; path entries b_(i,i+1)=0.
For every i<j<k the indicator of a directed triangle is

```text
b_ik(1-b_ij)(1-b_jk)+(1-b_ik)b_ij*b_jk
 =b_ik-b_ik*b_ij-b_ik*b_jk+b_ij*b_jk.                  (6)
```

The two terms before expansion are its two possible cyclic orientations.
Summing (6) gives the full triangle count. When reversed edges form a
matching, interaction terms vanish and the answer is the interval-length
sum in THM-410. With arbitrary chord flips the interactions are essential.

An independent counting formula is

```text
c3(T)=C(n,3)-sum_v C(d_plus(v),2).                      (7)
```

Every transitive triple has a unique vertex beating its other two; cyclic
triples have none. For n=10, convexity of C(d,2) makes the maximum 40:
the minimum degree contribution occurs at five scores 4 and five scores 5.
This is achieved by deleting one vertex of the cyclic regular tournament
on 11 vertices. Hence even the number of directed triangles can exceed 36.

## 4. Two strong completions with the same scores and triangle count

For n>=3 define two tournaments in the same fixed-path coordinates.

* I_n: reverse every nonconsecutive chord; retain each i->i+1.
* E_n: retain the transitive orientation except reverse n-1->0.

Both contain the directed Hamiltonian cycle `0->1->...->n-1->0`, so
both are strongly connected. Both have score multiset
`{1,1,2,3,...,n-2,n-2}` and c3=n-2. Nevertheless,

```text
c_ell(I_n)=n-ell+1,
c_ell(E_n)=C(n-2,ell-2),               3<=ell<=n.       (8)
```

In I_n take the smallest vertex of a directed simple cycle. Its only
available outgoing edge is to its immediate successor. The cycle must
continue upwards one step at a time; a backward jump reaches a vertex
already visited and therefore must close at the start. Its cycles are
exactly the contiguous intervals. In E_n every cycle must use its sole
backward edge; the remaining ascending path chooses any nonempty subset
of the n-2 interior vertices. These descriptions prove (8) bijectively.

At n=10 their total simple directed-cycle counts are 36 and 255, despite
the same score sequence and eight triangles. Their cycle counts already
differ at n=5: six versus seven (both have three triangles).

## 5. Odd-cycle overlaps yield a tribonacci law and the Fermat family

**CITED input.** The Grinberg--Stanley odd-cycle formula, specialized to
a tournament T, is

```text
H(T)=sum_C 2^(number of cycles in C),                   (9)
```

where C ranges over collections of pairwise vertex-disjoint directed
odd cycles of length at least 3, including the empty collection. See
[Theorem 1.39 and equation (42)](https://arxiv.org/pdf/2307.05569).
Reversing every path identifies the tournament and converse counts in
that specialization. Formula (9), not a count of triangles alone, is the
external input to the interval-family formula below. The endpoint-family
formula also has an independent elementary path bijection.

For I_n the collections are disjoint odd intervals. At the first vertex
either leave a singleton or start an odd interval of length at least 3,
weighted 2. Thus, with h_0=h_1=h_2=1,

```text
h_n=h_(n-1)+2*sum_(odd ell>=3,ell<=n) h_(n-ell),
sum_(n>=0)h_n*x^n=(1-x^2)/(1-x-x^2-x^3),
h_n=h_(n-1)+h_(n-2)+h_(n-3) for n>=3.                 (10)
```

This gives `1,1,1,3,5,9,17,31,57,105,193` through n=10.
For E_n there is a direct bijection requiring no external theorem.
A Hamiltonian path avoiding the sole backward edge `n-1->0` is the
unique increasing listing. Any path using that edge has an increasing
segment before it and an increasing segment after it. Choose an arbitrary
subset S of the n-2 interior vertices: its path is

```text
(S in increasing order), n-1, 0, (the complement of S in increasing order).
```

Every such choice is valid, and the vertices before n-1 recover S uniquely.
There are therefore 2^(n-2) paths using the backward edge, in addition to
the increasing path. Thus, independently of (9),

```text
H(E_n)=1+2*2^(n-3)=1+2^(n-2),             n>=3.        (11)
```

The odd-cycle formula gives an independent derivation: all E_n cycles
share both endpoints, so a collection contains at most one cycle;
exactly half the interior subsets have odd size.

At n=10, therefore, H(I_10)=193 while H(E_10)=257. The same 45 edge
positions, the same distinguished nine-edge path, the same scores, and
the same triangle count do not specify a Hamiltonian-path count.
E_10 has just one directed fundamental cycle, of even length 10, yet it
has 128 directed odd cycles; restricting (9) to basis cycles would fail.

The Fermat numbers occur **exactly** on the subsequence

```text
n=2^r+2  ==>  H(E_n)=2^(2^r)+1=F_r.                    (12)
```

Orders `3,4,6,10,18,34` give `3,5,17,257,65537,4294967297`.
If 2^m+1 is prime, m must be a power of two: an odd divisor d>1 of m
would make x^d+1 divisible by x+1, with x=2^(m/d). The converse fails:
F_5=641*6700417. Nothing in the definition, strong connectivity, cycle
law, or path-count formula changes at order 34. Primality is a property
of the resulting integer, not a graph phase transition.

## 6. A genuine directed product law, with an explicit orientation rule

For tournaments U,V on a,b vertices, U[V] replaces each vertex of U by
a copy of V and orients every pair of distinct fibres according to U.
Then (4) is its edge count, and

```text
c3(U[V])=a*c3(V)+b^3*c3(U).                            (13)
```

Triangles lie within one fibre or in three distinct fibres. A triple
occupying exactly two fibres is transitive because its outside vertex
either beats or loses to both inside vertices. This proves (13).
The size identity becomes a cycle identity only after specifying that
orientation rule; arbitrary orientations of the same product set need
not obey it. The three-cycle substitution C3[C3] has 30 directed triangles
and, by exact replay, 3159 Hamiltonian paths; the tempting product
`H(C3)*H(C3)^3=81` is false. Triangle transport is not general path-count
multiplicativity. The strongest direct survivor is the lower bound
`H(U[V])>=H(U)*H(V)^a`: its right side counts exactly the paths that
visit each fibre in one contiguous block. Extra paths interleave fibres,
which is the coordinate forgotten by multiplying only internal counts.

## 7. The actual planarity predicate and the remaining boundary

For every tournament on n vertices the underlying graph is K_n,
regardless of chord directions. It is planar exactly for n<=4. K_4 has
a planar drawing, while K_5 has 10 edges, exceeding the planar simple
graph bound 3v-6=9 derived from Euler's formula. Thus every orientation
on ten vertices is nonplanar in this ordinary sense. Even the isolated
36-chord graph exceeds 3*10-6=24 edges. Deleting orientation cannot
distinguish the transitive example from either strong example above.

The [Robertson--Seymour graph-minor theorem](https://doi.org/10.1016/j.jctb.2004.08.001)
gives finite forbidden-minor characterizations for minor-closed classes
of finite undirected graphs. An arithmetic operation graph first needs
a defined graph and a demonstrated minor-closed predicate. The theorem
does not provide an arithmetic orientation rule or a Collatz descent
certificate. No directed-minor claim is being made here.

The constructive question is now precise: is there an independently
defined arithmetic pairwise observable that selects I_n, E_n, a
substitution, or another chord pattern and preserves a desired arithmetic
predicate? The numerical label 3 alone does not supply such a map.

## Reproduction

Run `python 04-computation/experiments/arithmetic_seams_20260921_tournaments.py`.
The companion script/JSON record exact identities, full fixed-path cubes
through n=7, independent triangle counts, simple-cycle classification and
Hamiltonian-path dynamic programming for the two families through n=10,
odd-cycle collection checks, substitution, and the Fermat subsequence.
Explicit checks remain active under Python -O. The general proofs above
do not infer infinite statements from those finite controls.
