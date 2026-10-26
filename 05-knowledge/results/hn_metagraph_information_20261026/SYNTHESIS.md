# Keep the relation, its alignment, and the operation it must support

2026-10-26. This session moves the proposed six-colour obstruction through
explicit symmetry and metric tests, makes the inherited complete tournament
family construction executable, and recovers the older information programme
with exact reconstruction statements. The remaining Hadwiger–Nelson target is
an actual finite unit-distance graph with no proper six-colouring. None is
claimed here.

## 1. The four-state object needs more than three unprescribed ports

The previous pairing construction gives eight endpoint choices, split into
two four-state Fano completions. Its abstract transport group is `C2 x S4`.
But an ordinary unprescribed graph with just three exposed ports sees only
their colour-equality partition. All eight choices from three disjoint colour
pairs are injective triples and belong to one global `S6` orbit. Their
extension weights are therefore identical, for every finite graph.

With two corresponding paired frames, the complete restricted boundary table
has form `w(x,y)=f(x+y)`. On one four-state completion there are16 invariant
Boolean tables and only4 deterministic permutations: translations. This is a
necessary symmetry condition, not a graph-realization theorem for those16
tables. A larger jointly retained palette boundary can escape this three-port
restriction; its variables and feasible assignments must actually be retained.

There is also a sharp chromatic gate. Three-port fragments with chromatic
number at most4 extend every proper six-colour boundary assignment. With
independent ports and chromatic number at most5, only the all-distinct boundary
orbit can be absent. An exact six-colour NAE fragment requires its exterior
to be at least six-chromatic. The abstract sharp examples contain `K4`, whose
unit Gram matrix has rank3, and fail planar unit-distance realization.

The old prescribed unit gadgets survive their metric audit but change their
relations when unpinned. Their ordinary weights on the five three-port
partitions are

```
w_F = 6*5^25 (60,68,68,68,71),
w_G = 6*5^15 (1,1,1,1,1).
```

Every orbit survives. An exact seven-point graph consisting of two unit
triangles and one isolate also realizes an inconsistent parity-compatible
conditional frame loop. It nevertheless has86400 six-colourings:1728 of its
2187 ternary pair assignments lift. This locates the missing existential
quantifier in an actual induced unit graph, with every pair distance checked.
Details and full tables are in [six_color/RESULT.md](six_color/RESULT.md).

## 2. Whole planar frame transport removes the hidden obstruction itself

For any three distinct antipodal direction pairs in the plane, the three
direction vectors have a unique dependence with nonzero coefficients. Its
signs select exactly two opposite endpoint triangles containing the centre.
This balanced antipodal pair is a canonical point of the four-state space.
Every isometry of the entire six-point frame preserves that point, in the
appropriate local coordinates.

Consequently **every network of such whole-frame isometries has a full
endpoint section whenever its parity equations are consistent**. The earlier
two holonomies with disjoint fixed pairs cannot share one planar realization.
The theorem extends the flower collision-pair observation to generic planar
frames with eight distinct cube subset sums. All110592 marked flower triangle
networks pass the exact full-section/parity comparison.

This rules out a construction class. Partial, relational seams with genuine
interior vertices are outside the theorem and remain candidates. Their
complete six-colouring extension tables must establish the operation, rather
than identifying a colour permutation with a rigid isometry. See
[metric_transport/RESULT.md](metric_transport/RESULT.md).

## 3. A000568 has a complete structural generator, with an unavoidable new-core step

The recovered weighted-core classification, alternating attachments and
canonical construction path already supply the necessary mathematics. This
session combines them in one executable enumerator:

1. Contract maximal transitive clone chains to a pair-free core `Q`.
2. Retain positive chain weights modulo the actual action of `Aut(Q)`.
3. Generate new cores by one alternating attachment phase per old chain.
4. Use an invariant deletion orbit to keep one unrooted child, with duplicate
   acceptance treated as an error.

This terminates and produces every class exactly once at every finite order.
It is an all-order algorithm, not a completed infinite list or a new scalar
recurrence. The finite census now reaches all191536 order-nine classes,
including147095 new pair-free cores; the compressed bank retains every
representative and birth. An invariant individualization/refinement labeller
removed the degree-permutation bottleneck. The exact banks are specified in
[enumeration/RESULT.md](enumeration/RESULT.md).

The same full family records reconstruct the metagraph: enumerate valid
Hamiltonian vertex words, quotient their free automorphism action, and retain
every fixed-path tiling and literal complement line. The unmerged class fibre
has size `H(T)/|Aut(T)|`; hence

```
A(n) = sum_(fixed-path tilings t) |Aut(T_t)| / H(T_t).
```

Converse-merged nodes have count `(A(n)+SC(n))/2`. The historical fabricated
identity `SC(n)=A(n-1)` is not revived. The July blue/black parity convention
is recovered through literal lines, including loops and multiplicities.

For core cutoff `m`, the exact inflation image of all order-m tournaments is
the set with core order at most m. Its size satisfies

```
B(n,m) <= binom(n-1,m-1) A(m).
```

Thus every fixed finite core bank eventually covers a vanishing fraction.
Exact count arithmetic at100 and128 bounds the fraction arising from smaller
cores by `2^-85` and `2^-113`, respectively. This is an obstruction to stopping
core generation, not an obstruction to arithmetic counting. Recent integral
and reverse-valuation engines are recovered and checked through128; their
odd-partition cores are different objects from literal tournament cores.

## 4. Binary intransitivity is in the relation; information depends on the operation

The signed triangle circulation is±3 for a directed cycle and±1 for a
transitive triple. The recovered Session36 theorem says the complete labelled
curl retains the SCC interiors and circular component order, losing only the
initial cut. It preserves closed H, but can fail after reversing one named
arc, already at three vertices.

The new exact repair considers the full response `R_T(F)=H(T xor F)` for
every named arc subset F. Its adjacent-arc quadratic Fourier coefficients
recover the products of all adjacent orientation signs. The line graph of
the complete graph is connected, so **the response determines T up to global
converse**. Exactly `binom(n,2)-1` binary correlations along a spanning tree
encode this quotient; one oriented arc restores T. This is an information
bound, not a bound on how many H queries compute the correlations.

The aligned pair `(full curl, full response)` determines T exactly for n≥3.
Taking the two isomorphism quotients separately destroys this property:
`C3 -> vertex` and `vertex -> C3` are the first witness, at order4. Their
relative labelling is the missing information. For the actual fixed-path
metagraph operation, `(curl, marked valid Hamiltonian path)` is already
faithful; the valid start identifies the missing SCC cut.

The March information-geometry note also required corrections: the order-five
Fisher derivative at zero is165/32, and score-conditional variance1935/112 is
not degree-two Fourier energy135/8. The order-three H=1 sublevel is a six-cycle,
so connectivity does not imply simple connectivity. The old text is preserved
under a correction banner, with exact witnesses in the mistakes ledger.
See [information/RESULT.md](information/RESULT.md) and the
[independent audit](audits/INFORMATION_AUDIT.md).

## 5. Arithmetic/F1 work and the next boundary experiment

The separately requested session has published
[normseven_arithmetic_f1_20261026](../normseven_arithmetic_f1_20261026/README.md).
Its audited results distinguish the outer Gale minus-one twist over Q, give
the determinant criterion for marked Laurent lifts, retain the exceptional
characteristic-three inertia, and repair the central-boundary classification
to allow finitely eight rational nodes whose resolution lifts can fail.
Normalization, rational versus geometric isomorphism, and residual torsors
remain distinct. These arithmetic statements do not supply colouring gadgets.

The incoming [descent continuation](../normseven_descent_moments_20261026/RESOLVED_FIBRES.md)
was also merged and read. It separates the residue-field lift from the tangent
conic: with residue field Q(sqrt(13)), conics A²+B²=2C² and A²+B²=3C² give
minimum point-field degrees two and four although their good-prime counts
agree. Its own audit is internal proof review plus independent computation,
not a second agent proof audit. This gives an additional arithmetic example
of information lost by counts; no map to plane colouring is inferred.

The next six-colour experiment has a precise acceptance condition: choose
fragments with a retained palette boundary large enough to pass the symmetry
gate, compute their full unprescribed extension relations, realize all edges
at unit distance, include every cross-piece constraint used by the target
graph, and show the final joint composition empty across **all** feasible
palette states. With private interiors and all seam constraints accounted for,
the exact count is a sum of products of complete boundary weights. Passing
only a prescribed frame sector, its parity shadow, or separate marginals is
insufficient. The tournament maps give exact examples of why common labels
and operation data must survive, without asserting a map from tournaments to
unit-distance colouring obstructions.
