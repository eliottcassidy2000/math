# Regular tournament moves: exact score fibres and the path-move boundary

**PROVED / classical:** fixed labelled score fibres are connected by directed
triangle reversals (Ryser). The elementary proof is included. **PROVED:**
regular vertices are independent in the actual Hamiltonian-path-complement
metagraph, and fixed-label odd-order path reversals have a parity obstruction
to literally compiling triangle reversal. **FINITE-EXACT:** every regular
labelled tournament and directed triangle move through order7; no n9 census
or general H-classification is claimed.

## 1. Inheritance and conventions

The repository's [March interchange-graph literature note](../gn-literature-search-2026-03-23.md),
§2, already distinguished score-preserving triangle moves from its then-
studied single-arc graph. This note additionally separates both from the
current Hamiltonian-path-complement graph. The classical result is stated
as Theorem1 in
[Guiduli, Gyarfas, Thomasse and Weidl, 2-Partition-Transitive Tournaments, JCTB72(1998), pp181-196](https://users.renyi.hu/~gyarfas/Cikkek/85_2partitiontranstourn.pdf),
which attributes it to Ryser. We are applying that theorem, not claiming
new regular-tournament connectivity.

Write d_T(v) for outdegree and H(T) for the number of directed Hamiltonian
paths in a fixed representative. Tournament labels are fixed until a
quotient is expressly taken. For a Hamiltonian path P, R_P reverses its
n-1 path edges; K_P complements R_P(T), reversing the nonpath edges of T.
Regularity at order n=2m+1 means every vertex has outdegree m.

The live concepts are score fibre; Eulerian difference; triangle compiler;
Hamiltonian-path count; actual R/K metagraph; label parity; isomorphism.

## 2. Exactly which arc reversals preserve regularity

Reversing one arc u->v lowers d(u) by1 and raises d(v) by1. Therefore a
single-arc flip of a nontrivial regular tournament is never regular.

More generally, let F be a set of directed arcs of T, reversed all at once.
Then

```text
d_new(v)-d_T(v)=in_F(v)-out_F(v).                        (1)
```

Thus the labelled score vector is preserved **iff** F is balanced at every
vertex. A finite balanced directed graph decomposes into edge-disjoint
directed cycles: start at an unused arc, follow until a vertex repeats,
remove that cycle, and continue. Balance remains after removal. A single
directed cycle, and in particular a directed triangle, therefore preserves
regularity. An arbitrary undirected cycle or arbitrary three-edge reversal
need not preserve it; directedness is essential.

If T and U have the same labelled score vector, orient their differing
edges as in T. At each vertex the outgoing differences and incoming
differences balance. Decompose this difference digraph into directed
cycles, and reverse the cycles in turn. Their edges are disjoint, so each
cycle is still directed when used. The final tournament is U.

## 3. A cycle-to-triangle compiler

Every directed k-cycle can be reversed by exactly k-2 directed triangle
reversals as a construction bound (not a claim of minimality). Label the
cycle `v1 v2 ... vk v1` and inspect the chord v1--v3.

* If `v3->v1`, reverse the directed triangle(v1,v2,v3), then recursively
  reverse the directed cycle(v1,v3,...,vk).
* If `v1->v3`, recursively reverse that shorter cycle first, then reverse
  the now-directed triangle(v1,v2,v3).

The chord is toggled twice; every original cycle edge is toggled once.
Induction proves both legality and the k-2 bound. Combined with §2, this
proves Ryser's fixed-score connectivity and in particular connects every
regular tournament to a labelled copy of the carousel tournament on
2m+1 vertices, where i beats the next m cyclic vertices.

Every score-preserving move also preserves the directed triangle count,
because

```text
c3(T)=binom(n,3)-sum_v binom(d_T(v),2).
```

Each transitive triple has exactly one source, which proves the formula.
For a regular tournament it reduces to `n(n^2-1)/24`. H is not preserved.

## 4. Exact small H changes and isomorphism

The independent census generates every labelled regular tournament by
completing rows subject to their exact required outdegrees. All valid
choices are retained; no random sampling or canonical-class filtering is
used before the labelled universe is complete. It then checks every
directed triangle reversal and identifies ordinary isomorphism classes.

| n | Labelled regular tournaments | Ordinary classes | H values |
|---:|---:|---:|---|
|1|1|1|1|
|3|2|1|3|
|5|24|1|15|
|7|2640|3|171,175,189|

At n7 every tournament has14 directed triangles, giving36960 directed
triangle moves. Their H differences, with directed-move multiplicities, are

```text
-18:3360, -4:5040, 0:20160, +4:5040, +18:3360.
```

The three ordinary class sizes are1680,720,240 for H171,175,189 respectively.
The quotient triangle graph has nonloop edges

```text
H175 -- H171 -- H189.
```

It has loops at175 and171, and none at189. Thus a triangle reversal can
either preserve or change isomorphism type; a change in H certifies the
latter. Through n5, uniqueness of the regular class forces every such move
to preserve isomorphism type. These are finite statements; H is not asserted
to classify regular tournaments at larger orders.

A reproducible hostile is the carousel with outgoing bitmasks

```text
(14,28,56,112,97,67,7),                 H=175.
```

Reverse its triangle(0,2,4), obtaining

```text
(26,28,41,112,100,67,7),                H=171.
```

Every vertex retains outdegree3 and c3 remains14, but the tournaments are
not isomorphic. A second witness from H171 to189 is in the receipt. Reversing
a legal move is again legal, so there is no universal sign for its H change.

## 5. Why this is not the Hamiltonian-path metagraph

For a regular T of order2m+1 and a directed Hamiltonian path from s to t,
internal path vertices lose one outgoing path edge and gain one. Its
endpoints do not balance:

```text
d_R = m*1 - e_s + e_t,
d_K = m*1 + e_s - e_t.                                  (2)
```

For n>=3, s and t are distinct, so neither R_P(T) nor K_P(T) is regular.
Hence the regular isomorphism classes form an **independent set** in the
actual path-complement metagraph, also after converse merging. This applies
to all path edges, not just black ones. Any actual metagraph route joining
two regular vertices must leave the regular family. Triangle-reversal
connectivity supplies a different normalization operation.

There is a further fixed-label obstruction. Set

```text
epsilon(T)=sum_(i<j) [i->j] mod2.
```

At odd n, an R_P move toggles n-1, hence an even number of arcs. Every
sequence of R moves preserves epsilon. A directed triangle reversal toggles
three arcs and changes epsilon. Thus **no directed triangle reversal can
be literally compiled into any number of R moves on the same fixed labels**.
This is stronger than failure of a short compiler.

The sign coordinate does not descend to ordinary tournament isomorphism:
relabeling by a permutation sigma changes epsilon by the permutation parity,
since each inversion reverses the orientation convention of one unordered
pair. In particular, a transposition changes it. Converse changes epsilon
by `binom(n,2) mod2`, which is m mod2 at n=2m+1. The two labelled sign sheets
therefore do not give two components in the converse-merged isomorphism
metagraph. For K moves the parity increment is instead m mod2. The
no-compiler statement in the preceding paragraph is specifically about R.

## 6. Reproduction and audit boundary

```text
python3 05-knowledge/results/spectrum_basins_20261020e/regular_probe.py
python3 -O 05-knowledge/results/spectrum_basins_20261020e/regular_probe.py
```

Every H value is computed by subset DP; each ordinary class representative
is independently checked by enumerating all vertex permutations. All
labelled regular universes are connected under the fully enumerated triangle
moves. Every single arc flip is checked as a hostile. Every Hamiltonian
path of every small regular class representative verifies (2) and the R
parity statement. Separately, all5338 directed-cycle presentations in all
tournaments of orders3,4,5 verify the literal cycle-to-triangle compiler.

The program uses explicit exceptions, so checks survive Python `-O`.
Normal and optimized runs produce the identical JSON SHA256
`b93d4e2ccc78997add8fa8ff118e38243a7fcf68a54b68823e977fb21718ddbd`.
The receipt records source hash, representatives, all quotient move counts,
and an explicit witness for every realized class-pair/H-change combination.
