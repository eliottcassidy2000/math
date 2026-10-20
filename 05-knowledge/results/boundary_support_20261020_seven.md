# Seven, forty-two, and fifty-six as actual boundary operations

**Status: PROVED elementary matching, contraction, and involution statements;
FINITE-EXACT the specified eight-port kernels and two closure censuses.
OPEN a general obstruction-preserving minor reduction.**
No theorem ID or novelty claim. This packet extends the inherited construction;
it does not identify tournament classes with cubic graphs.

## 1. Inheritance and the intended objects

The user's numbers have a precise existing source:
[Hamming matching atlases, Proposition A2](hamming_matching_atlases_20261018e.md).
The maximum component charts in `H8 direct-sum I2` project to the 28 pairs
of the eight slopes `F2^3`. A partition into four maximum charts is a
perfect matching of those slopes. There are

```text
105 = 7 + 42 + 56
```

such partitions, using respectively one, two, or four distinct pair-XOR
directions. These are not three unrelated coincidences.

The boundary mechanism comes from
[THM-4116, ordered boundary gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md)
and the [previous graph/Pythagorean synthesis](graph_pythagorean_20261018_synthesis.md).
The tournament mechanism comes from
[THM-643, the blue/black involution frame](../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md)
and [THM-4731, self-lines, twins and depth](../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md).
Here `tau` in THM-643 is `kappa` in THM-4731.

The closest proved mechanism is contraction of *labelled* boundary states.
The hostile is a scalar count or quotient losing the attachment map.
The corrected near miss is a code-preserving action that changes its
partition or phase, recorded in the inherited atlas. The least-used
sidecar here is the actual perfect matching permutation before quotienting
to graph or tournament classes.

Board: cubic boundary support (Anchor); eight-port matching operators
(Niche); tournament fibre pushforwards and Pythagorean Gram geometry
(Wildcard). MISTAKE-507 blocks a transfer from shared numbers to a
colouring/flow theorem. No old conjecture from THM-643 is used as proved;
THM-4731's correction of the odd-order self-line guess remains in force.

## 2. Recovered proof of the direction census

Give the eight addresses the additive structure `V=F2^3`. For a perfect
matching `M`, assign direction `x+y` to each pair `{x,y}`. The XOR of all
four directions is zero because every address occurs once. Four nonzero
vectors with sum zero have multiplicity pattern `4`, `2+2`, or `1+1+1+1`:
patterns `3+1` and `2+1+1` would force two purportedly distinct directions
to agree.

For each nonzero `d`, the pairs `{x,x+d}` give one parallel matching,
so there are seven. For distinct `d,e`, their plane has two cosets.
Use direction `d` in one coset and `e` in the other, with two choices of
assignment. This gives `2 binomial(7,2)=42`. The remaining 56 use four
directions; their four-element set is an affine plane not containing zero.
There are seven such direction sets, each carrying eight matchings, by
the transitive linear action and the total count.

This classification depends on the retained affine addresses. An arbitrary
permutation of the eight ports can take any matching to any other and
therefore change its direction class.

## 3. The matching atlas acts on cubic colour boundaries

Let `P` be a cubic eight-pole: every internal vertex has total degree three
when semiedges are counted. Label its eight semiedges by `V`. Let
`f_P(s)` be its number of proper extensions of the ternary boundary word
`s:V->{0,1,2}`. For each matching define the equality kernel

```text
B_M(s) = product over {x,y} in M of [s(x)=s(y)].
```

Join the paired semiedges into edges, obtaining the closure `P/M`.
The exact gluing bijection gives

```text
C(P/M) = <f_P,B_M>.
```

Parallel edges are allowed. A pairing at the same vertex would create a
loop and contributes zero automatically, because its two incident
semiedges have distinct colours in every proper extension.

**Connection contract.** The source is a labelled cubic eight-pole and
one of the inherited matching partitions. The map closes paired ports.
It preserves exactly the counted colouring extensions through that
matching. It forgets no information if `f_P`, the port attachment map,
and `M` are retained. The direction class alone forgets almost everything.
Closing a pole is not a graph-minor operation in general; an additional
minor model is required before this becomes a Petersen-minor reduction.

The cubic parity lemma implies that all three colour multiplicities on
an eight-port boundary are even: sum the three nonzero `F2^2` colour
labels at all vertices. Thus the possible nonzero sizes are

```text
8; 6+2; 4+4; 4+2+2.
```

For any such word with multiplicities `n_0,n_1,n_2`, the number of
compatible matchings in the complete bank is

```text
product_c (n_c-1)!! = 105, 15, 9, or 3,
```

where `(-1)!!=1`. Pair within each monochromatic class independently.
Consequently every nonempty cubic eight-pole support admits a colourable
closure in the full bank. The bank cannot consist entirely of
obstruction-preserving rewirings.

The seven parallel kernels accept a word exactly when the word has a
nonzero translation period. They miss 1,176 of the 1,641 even-parity
ternary words:

| Multiplicity pattern | All compatible matchings | Parallel-compatible words | Parallel-blind words |
|---|---:|---:|---:|
| 8 | 105 | 3 | 0 |
| 6+2 | 15 | 168 | 0 |
| 4+4 | 9 | 42 | 168 |
| 4+2+2 | 3 | 252 | 1008 |

The table counts labelled ternary words, not words modulo colour
permutation. A concrete blind word is `(0,0,0,1,0,1,1,1)`; it is accepted
by the two-direction matching `01,24,35,67`. Parity alone does not prove
that this single word orbit is the complete support of a realizable pole.

## 4. The exact Gram quotient: 105 kernels, 91 dimensions

**PROVED.** With the ordinary counting inner product,

```text
||B_M||^2 = 3^4 = 81,
<B_M,B_N> = 3^c(M,N),
```

where `c(M,N)` is the number of connected components of the union of the
two matchings, counting a common pair as one component. The two equality
constraints together say precisely that each component is monochromatic.

Distinct matching kernels therefore always overlap: the inner products
are 3, 9, or 27. No pair of bare matching kernels yields an orthogonal
snark witness.

**FINITE-EXACT, independently checked by all ternary words.** The complete
105-by-105 Gram matrix has characteristic polynomial

```text
t^14 (t-945) (t-210)^20 (t-120)^14 (t-30)^56.
```

Its rank is 91. The three direction banks have ranks 7, 42, and 56
separately; combined ranks are

```text
one + two directions: 49;
one + four directions: 63;
two + four directions: 91.
```

Thus the 98 nonparallel probes already span all the probes.

The kernel has an elementary source. Divide the eight addresses into
two four-element sets, and antisymmetrize over the 24 bijections between
them. The resulting kernel sum is

```text
det([s(x_i)=s(y_j)]_(i,j=1..4)) = 0,
```

because this matrix factors through the three-dimensional colour space.
There are 35 unordered `4+4` partitions. Their exact relation matrix has
rank 14, so they give every dependence among the 105 kernels.

There is an especially small recovery formula. For a fixed parallel `M`,
write `S_(r,c)` for the sum of `B_N` over matchings `N` with `r` distinct
directions and `c(M,N)=c`. Then

```text
4 B_M = 2 S_(2,3) - 2 S_(2,2) - S_(4,2) + S_(4,1).
```

The four sums contain respectively 12, 6, 32, and 24 matchings. The
identity is certified for all seven choices of `M` in two independent
ways: Gram nullity and direct evaluation at every even-parity ternary
word. It therefore holds after contraction against any boundary vector.
Its negative coefficients matter: a zero on the left is a cancellation
identity, not a proof that all constituent closures are uncolourable.

This 91-dimensional quotient retains all *pairing-closure* counts. It
is not claimed to retain gluing against every possible cubic eight-pole.
The space of all colour-permutation-invariant even-parity boundary
functions has dimension 274, so unrestricted boundary data are much
larger. Realizability of those functions is an additional question.

## 5. Petersen and the prism: complete 105-closure experiments

Use the usual vertices `0,...,4` on the outer pentagon and `5,...,9`
on the inner pentagon/pentagram. Cut spokes `i--(i+5)` for `i=0,1,2,3`.
Ports `0,...,7` attach, in order, to vertices

```text
0,1,2,3,5,6,7,8.
```

Both resulting eight-poles have 300 boundary words, each with a unique
extension. Indeed a five-cycle pole has 30 boundary words and ten with
any specified colour at a marked port. Gluing one port of each of two
such cycles gives `3*10*10=300`; the remaining boundary determines the
glued colour by the `F2^2` parity equation, proving multiplicity one.

All 105 closures are audited by two independent exact engines: direct
boundary-colouring recursion and enumeration of ordered binary cycle
spaces whose union covers every edge, which counts nowhere-zero
`F2^2` flows. Bridges are checked in the actual multigraph.

| Eight-pole | One-direction closures with `C=0` | Two-direction closures with `C=0` | Four-direction closures with `C=0` |
|---|---:|---:|---:|
| Petersen, four spokes cut | 5 of 7 | 6 of 42 | 0 of 56 |
| Prism, four spokes cut | 3 of 7 | 8 of 42 | 0 of 56 |

In each complete bank there are exactly eleven uncolourable closures:
two simple bridgeless copies of Petersen and nine graphs with a bridge.
For the Petersen port gauge the two Petersen closures are parallel;
for the prism gauge they lie in the two-direction bank. The full labelled
colour-count histograms and simple-graph filters are in the receipt.

These are the same underlying barbell of two five-cycles joined by one
edge, with different port addresses. The change in the direction bank
of the Petersen closures is therefore a particularly direct hostile
to treating direction count as an unmarked graph invariant.
The prism control also shows that arbitrary rewiring can create a
colouring obstruction. Such rewiring cannot serve as an automatic
obstruction-preserving descent.

## 6. The same matching carrier acts on tournament fibres

For a finite labelled presentation set `X`, a fixed-point-free involution
`tau` is exactly a perfect matching. Let `pi:X->C` be the class map and
let `A_(c,x)=[pi(x)=c]`. If `P_tau` is the permutation matrix, then

```text
W = A P_tau A^T
```

records all class-to-class matching endpoints. Off-diagonal entries are
edge multiplicities; a diagonal entry is twice its loop multiplicity.
This is the same operation pattern as `C(P/M)=<f_P,B_M>`: retain labelled
states, apply an explicit matching, then contract the observation map.
Here the observation is a tournament class, not an edge colour.

At tournament order four there are exactly three tile bits:
`(0,2),(0,3),(1,3)`. Thus its whole tiling space is literally `F2^3`,
not merely a set that happens to have 105 partitions. The inherited
complement involution is translation by `111`, one of the seven parallel
matchings. Staircase reflection swaps bit zero with bit two and fixes
bit one. This is a labelled joint carrier for the two theories.

**FINITE-EXACT, all 105 pairings.** The actual four ordinary tournament
isomorphism-class fibres have sizes `1,1,1,5`. Pushing forward the 105
matchings gives only four uncoloured ordinary-class metagraph wirings:
the three leaves all joined to the
five-element centre, plus a centre loop (60 marked lifts); or one pair
of leaves joined, the third joined to the centre, and two centre loops
(15 marked lifts for each of three leaf pairs).

The user's converse-merged quotient is a further step: merge the source-
and sink-cycle classes, giving three node fibres of sizes `1,2,5`.
The two wirings that join the transitive singleton to one or the other
of those cycle classes coincide after merging. Thus the complete
105-bank pushes forward to three merged wirings with lift counts
`60,30,15`. The source preserves both class maps; it does not treat these
three merged nodes as the four ordinary classes.

The reflection imposes stronger requirements. Exactly nine matching
involutions commute with it: three parallel and six two-direction.
If additionally `sigma tau` is fixed-point-free, as it is for the actual
complement lines by the fixed apex tile, exactly six remain: two
parallel and four two-direction. All six project to the same coloured
line metagraph: the transitive leaf joins the centre by a blue line;
the centre has one blue loop; each of the two mutually converse leaves
joins the centre by one black line. None of the 56 four-direction
matchings passes this reflection test.
After converse merging the two black lines are parallel edges between
the merged cycle node and the centre. Their twin quotient is one edge,
exactly the fundamental black graph construction in THM-4731.

**PROVED explanation of nine and six.** Reflection has four fixed
states and two exchanged pairs. A commuting matching must pair the
four fixed states internally, in three ways. On the four other states,
it either pairs each reflection pair internally (one way), or pairs the
two reflection pairs across (two ways). This gives `3*(1+2)=9`.
The first case makes `sigma tau` have fixed states; removing it leaves
`3*2=6`.

Consequently a quotient metagraph wiring is insufficient to reconstruct
the complement operation, even after the blue/black split is retained:
six marked lifts already coincide at order four. At larger orders one
must retain the class map, actual Hamiltonian presentations, automorphism
orbits and the two commuting involutions. The 105-atlas does not by
itself supply a sufficient realizability theorem for arbitrary complete
tournament metagraph wiring.

## 7. Pythagorean scope and the next usable reduction test

Any two equal-norm vectors `u,v` give perpendicular `u+v,u-v`, with
squared legs `2(||u||^2 +/- <u,v>)`. Applied to the matching kernels,
the squared legs are `2(81 +/- 3^c)` and hypotenuse squared is 324.
Because `3^c>0`, these triangles are never exactly isosceles for distinct
matchings. They are real Gram triangles, not asserted integer triples.

The inherited four-chart Hamming atlases instead have four disjoint
indicator supports of size eight. Those indicators are orthogonal and
give isosceles triangles in *every* one of the 105 atlases. That geometry
detects partition disjointness; it does not detect a snark.

There is a more informative arithmetic connection in the concurrent
[flower-strip reduction](boundary_support_20261020_snarks.md). Its actual
colour-extension transfer satisfies, on each active three-state block,

```text
T^(2h) = I + alpha_h J,       alpha_h=(4^h-1)/3.
```

Thus the distinct eigenvalues are `1` and `1+3 alpha_h=4^h`, with square
root `2^h`. Euclid's map at `(m,n)=(2^h,1)` gives the primitive triple

```text
(3 alpha_h, 2^(h+1), 3 alpha_h+2)
    = (4^h-1, 2^(h+1), 4^h+1).
```

At `h=1,2,3` these are `(3,4,5)`, `(15,8,17)`, and `(63,16,65)`.
Here the arithmetic parameter is an actual extension multiplicity,
`alpha_h`, not a count chosen after the fact. This integer-triple chart
flattens as `h` grows, whereas every uncolourable odd flower has exact
isosceles geometry in the separate boundary-overlap chart. One must
retain which of the two maps is being used.

For a useful cubic reduction the extra test remains:

1. identify a smaller pole with the same labelled boundary support, or
   a proved containment in the direction needed for noncolourability;
2. retain a minor model for the replacement;
3. check bridges and the actual exterior boundary constraint;
4. use all needed matching probes, not just the seven translations.

The 91-dimensional closure quotient is a cheap exact filter for step 1
when the exterior is a pairing. A full support certificate is needed
before extending the claim to arbitrary exteriors.

## Reproduction

```text
python3 04-computation/boundary_support_20261020_seven.py
```

[Source](../../04-computation/boundary_support_20261020_seven.py) and
[exact receipt](../../04-computation/boundary_support_20261020_seven.out).
The universe is complete at eight ports: 105 matchings, 6,561 ternary
words, two named eight-poles, and every matching on the order-four
tournament tiling cube. No random sampling, floating point, or
unproved colouring classification enters the receipt.
