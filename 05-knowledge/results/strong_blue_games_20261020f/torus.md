# A marked torus return: Euler length, winding, and Euclid's doubled leg

**PROVED:** on a positively directed p-by-q torus grid, every Euler circuit
has length 2pq and marked winding (q,p). Coprime p,q give a primitive simple
torus representative. The same parameters give Euclid's raw triple
`(q^2-p^2,2pq,q^2+p^2)` when q>p. This is a common parameter construction,
not an invariant assigning a torus knot to an abstract complementary game.
**CITED:** the standard torus-surface framing of T(q,p) is pq.
**FINITE-EXACT:** explicit Euler tours for 81 small ordered parameter pairs,
plus six named positive and hostile controls.

## 1. The inherited map and its required extra data

The earlier [golden arithmetic packet](../spectrum_basins_20261020e/arithmetic.md)
identifies the ordinary SCC-word count at odd order 2k+1 as
`p*q=F_(k+1)*F_(k+2)`. Its two matching factors provide a distinguished
ordered pair of coprime cardinalities; an individual word does not become
a knot. **THM-3509**,
[reduced-fraction harmonic K4 face and Fibonacci unit-Cassini ray](../../../01-canon/theorems/THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray.md),
supplies the inherited Farey/Euclid dictionary with its parity-content
coordinate.

The canonical hostile is **THM-2646**,
[braid-three modular central pullback and full-twist knot fibre](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md),
§2: the same modular braid image can lift to an unknot and a trefoil.
The missing integer central coordinate matters. MISTAKE-225 in
`01-canon/MISTAKES.md` also rejects identifying tournaments and torus dynamics
without a map. Here the map is explicit, and its losses are explicit too.

## 2. The positive torus-grid construction

For positive integers p,q, use the marked torus `R^2/Z^2`. Put vertices at

```text
(i/p,j/q),  0<=i<p, 0<=j<q,
```

with a positively oriented horizontal edge to `(i+1,j)` and a positively
oriented vertical edge to `(i,j+1)`, reducing indices modulo p,q. Retain
the E/N direction and the actual torus edge, including edges crossing the
chosen fundamental rectangle's boundary. For p,q>=3 its underlying graph
is a simple undirected torus grid. Smaller parameters require the natural
cellulation multigraph: p=1 gives horizontal loops, and p=2 has two distinct
horizontal torus edges between the two positions.

There are pq vertices and 2pq directed edges. Every vertex has indegree
and outdegree two, and the directed graph is strongly connected: repeated
positive steps can attain any pair of residues. Consequently the standard
edge-splicing algorithm constructs an Euler circuit. Every such circuit
uses pq horizontal and pq vertical arcs, so its lifted displacement is

```text
(pq/p, pq/q) = (q,p).                                     (1)
```

This proves both the length 2pq and the winding claim independently of the
order in which the circuit uses the edges.

If gcd(p,q)=1, the straight torus curve

```text
t |-> (qt,pt) mod Z^2,  t in R/Z,                         (2)
```

is embedded: equality at parameters t,t' makes both q(t-t') and p(t-t')
integers; Bezout's identity then makes t-t' an integer. On a standard
unknotted embedded torus in S3 it is the torus knot T(q,p), with this explicit
choice of coordinate convention. If the gcd is greater than one, (2) is
multiply covered and is not an embedded parametrized knot; dividing the
winding or selecting parallel components would be an additional operation.

The Euler circuit itself repeats vertices and is not the knot embedding.
The map is

```text
marked directed Euler circuit -> marked homology (q,p)
                             -> the chosen straight representative.
```

It preserves winding and, within this grid family, recovers the length
2pq from the ordered winding coordinates. It forgets traversal order,
repeated visits, local crossing choices, and all Hamiltonian-path counts.
It is a homology/representative construction, not an ambient-isotopy claim
about a self-intersecting walk.

## 3. A length / framing / Pythagorean dictionary

With positive orientation and the standard torus-surface push-off, the
framing of T(q,p) is pq; see
[Finashin--Kharlamov, Topology of real cubic fourfolds, §6, p25](https://arxiv.org/pdf/0906.1480).
The same pq is the grid's vertex count and half its Euler length. For q>p,
Euclid's identity is

```text
(q^2-p^2)^2 + (2pq)^2 = (q^2+p^2)^2.
```

Thus the Euler length is literally the raw doubled leg. If p,q are coprime
and of opposite parity, the triple is primitive. If both are odd, its
content is two and the primitive doubled leg becomes pq. The primitive
normalization must be retained as a sidecar.

| p,q | Euler length | Marked winding | Raw triple | Primitive triple |
| --- | --- | --- | --- | --- |
| 1,2 | 4 | (2,1) | (3,4,5) | (3,4,5) |
| 3,4 | 24 | (4,3) | (7,24,25) | (7,24,25) |
| 3,5 | 30 | (5,3) | (16,30,34) | (8,15,17) |
| 13,21 | 546 | (21,13) | (272,546,610) | (136,273,305) |

The root 3-4-5 example uses the loop/double-edge cellulation boundary; it is
not a simple subgraph in the K_n game. The p=13,q=21 row connects the earlier
273-count to an actual length and a framed knot using precisely the same
factor pair. It does not identify the 273 individual SCC words with grid
vertices unless their original pair of matching labels is explicitly used
to label those vertices; even that set bijection would not carry adjacency.

Length alone loses necessary information. The grids (p,q)=(5,12) and (6,10)
both have 60 vertices and Euler length 120, but only the first has primitive
winding. This is a cheap decisive hostile to an unmarked length-to-knot map.
Even among primitive slopes, forgetting the marked pair and retaining its
product discards the slope.

## 4. What this does and does not say about the game

The abstract K_n game in [game.md](game.md) has no chosen embedding in a
torus, no E/N edge directions, and no distinguished homology basis.
Consequently its number H of Hamiltonian paths does not automatically
supply a winding pair. For example K_n has n! ordered Hamiltonian paths
(n!/2 up to reversal when n>=2), whereas its Euler-trail existence is governed
by vertex parity: K_n has an Euler circuit when n is odd; for even n>=4 it
has more than two odd-degree vertices and no Euler trail. These are different
predicates and cannot be equated by using H as an arbitrary knot parameter.

The useful survivor is a concrete common carrier: **a marked embedded grid
with its directed Euler traversal**. Its winding, length, framing, and Euclid
triple are related by proved formulas. A game-to-grid construction preserving
the winning predicate is still absent. No transfer to Petersen minors,
Kuratowski obstructions, or tournament black connectivity follows from this
dictionary alone.

### An exact Hamiltonian/Eulerian count bridge from incoming work

The [incoming port-and-winding note](../kernels_games_detours_20261020f/walks_knots.md)
provides the missing incidence-preserving comparison. For a directed graph
D with labelled arcs, form its directed line graph L(D), whose vertices
are the arcs, with e->f exactly when head(e)=tail(f). Linear Eulerian arc
words in D are in bijection with Hamiltonian vertex words in L(D): both
are permutations of the same arc set, with the same consecutive-incidence
test. This also identifies their counts, with no quotient by rotation.
The closed variant adds the last-to-first incidence test.

For our torus grid D_(p,q), L(D) has 2pq vertices. Its Hamiltonian words
are the grid's linear Euler words; the underlying walks carry winding
(q,p) and length2pq. The Hamiltonian paths themselves have 2pq-1 edges.
The vertex labels retain each original directed arc and its marked lift
increment, which is how the winding survives. Dropping those labels loses
the knot parameters. L(D) is generally not a tournament, and no equality
with the HP count of T_k is asserted.

This is sensitive to direction: the ordinary undirected line graph of
K_(1,3) is K3 and has a Hamiltonian path, while the claw has no Euler trail.
A pairwise common endpoint need not be the endpoint at which the current
traversal arrives. The incoming note also retains integer carries on the
residue map i->i+q mod p; the residue graph alone identifies q only modulo p.
Thus the new exact line-graph bridge and our marked-grid bridge preserve
compatible data, while an unmarked edge count does not.

## 5. Explicit checks

```powershell
python3 05-knowledge/results/strong_blue_games_20261020f/torus.py
python3 -O 05-knowledge/results/strong_blue_games_20261020f/torus.py
```

[torus.json](torus.json) contains the script hash and literal edge-ID Euler
circuits for the six displayed/hostile examples. The checker independently
constructs circuits for all 81 ordered pairs `1<=p,q<=9`, verifies continuity,
one use of every directed edge, vertex spanning, marked winding, and the
Euclid identity and parity normalization. These are finite controls for the
general proof, not a knot-equivalence calculation.
