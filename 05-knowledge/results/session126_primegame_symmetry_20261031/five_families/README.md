# Session 126: five-member families, exact actions, and decisive hostiles

**Status:** PROVED elementary interfaces; CITED classical classification and
all-height square-sum result; FINITE-EXACT checks in the stated universes.
**Not proved:** that there are exactly five Fermat primes, a universal five-way
classification, or an identification of the user's metagraph with a Platonic graph.

Source inheritance was read at `8f76b84d2960f9e14ad62a154b8533c77593702f`.
This portable packet includes the A5/FUS addendum, the compact required census,
and frozen mathematical proof dependencies. `SOURCE_PINS.json` records the
original files and hashes. Publication copies replace Python assertions with
explicit gates; no source-worktree or remote edits were made by this lane.
Run `python3 verify_packet.py` to check normal/optimized replay and relocation;
`SHA256SUMS` authenticates the complete packet.

## 1. Three different meanings of “five”

| Family | Exact status and mechanism | Cheapest hostile to an identification |
|---|---|---|
| `{15,16,17,23,25}` | First five **nontrivial** square-sum Hamiltonian-path orders, not the full solution set. Four sporadic orders lie below25;25 begins the infinite tail. |26 also works. Our exhaustive search of2..26 yields precisely15,16,17,23,25,26. The vacuous order1 is excluded explicitly. |
| `{3,5,17,257,65537}` | Five **known** Fermat primes. Further examples, finiteness, and “exactly five” remain OPEN. |`F_5=2^32+1=641*6700417` kills “all Fermat numbers are prime,” but does not prove that no later prime exists. |
| Five convex Platonic solids | A complete classical classification, not a truncated search. Euler plus incidence forces the five Schläfli pairs below. |A regular17-gon is constructible, but17 cannot be a face-size of a convex Platonic solid: `(17-2)(q-2)>=15>4`. |

Gerbicz's own dated contribution to [OEIS A090461](https://oeis.org/A090461)
reports the all-height theorem: every order at least25 has a square-sum path,
and every order at least32 has a square-sum cycle. The linked original forum
was403-blocked in this audit, so the tail is **CITED**, not re-proved here.
The finite prefix was independently exhausted. Current five-known-prime status
is recorded by [OEIS A019434](https://oeis.org/A019434); the
[Boklan–Conway paper](https://arxiv.org/abs/1605.01371) is heuristic evidence,
not a proof of completeness.

For a convex regular polyhedron with face size `p` and vertex valence `q`,
`pF=2E=qV` and `V-E+F=2` give `1/p+1/q>1/2`, equivalently
`(p-2)(q-2)<4`. Since `p,q>=3`, the possibilities are exactly
`(3,3),(3,4),(3,5),(4,3),(5,3)`. Constructions establish existence;
[Euclid XIII.18](https://mathcs.clarku.edu/~djoyce/elements/bookXIII/propXIII18.html)
contains the classical exhaustion argument. Euler's equation alone is not an
existence proof, and dropping convexity changes the problem.

## 2. The genuinely shared Fermat mechanism is a 2-group

For an odd prime `p`, the same group `(Z/p)^*` acts on:

1. cyclotomic roots, by `zeta_p -> zeta_p^a`; and
2. circulant-tournament connection sets `S`, by `S -> aS`.

The first action identifies the cyclotomic Galois group; the second is the
actual multiplier operation, not a numerical analogy. The following criteria
are equivalent:

`p is Fermat prime <=> p-1 is a power of2 <=> the regular p-gon is constructible`
`<=> every skew connection set S has trivial multiplier stabilizer`.

The geometric statement is Gauss–Wantzel: more generally a regular `n`-gon is
constructible iff `n=2^a` times a squarefree product of Fermat primes. The
quadratic-tower mechanism is in [Wantzel's original paper](https://www.numdam.org/item/JMPA_1837_1_2__366_0/)
and the cyclotomic proof is given by [NPTEL's Galois-theory lecture](https://archive.nptel.ac.in/content/storage2/courses/111101001/modules/lec16/16.4.html).

For the tournament direction: a nontrivial subgroup of a cyclic2-group contains
`-1`, which cannot fix a skew set. Conversely an odd prime-order subgroup of
`(Z/p)^*` pairs its cosets under negation; choosing one coset from each pair
constructs a nontrivially stabilized skew set. This is
[Akin, Theorems3.13–3.14](https://arxiv.org/pdf/1806.11241), stated there even
for odd composite moduli. For odd `m`, universal multiplier freeness is
equivalent to `m` being a squarefree product of Fermat primes.

The inherited **THM-4644**
(`01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md`)
adds Burnside's prime-degree theorem to get the full automorphism statement:
for odd prime `p`, every circulant tournament has `Aut=C_p` iff `p` is Fermat.
Do not attribute this full-group conclusion to Akin alone: his cited text
explicitly leaves the17 case open. At prime degree a tournament's automorphism
group cannot be doubly transitive, hence lies in the affine group, reducing
the claim to its multiplier stabilizer.

**Bounded bridge to the user's first ten primes:** among the odd primes
`3,5,7,11,13,17,19,23,29`, the universal multiplier-rigidity property holds
exactly at `3,5,17`. Exhaustive enumeration checks every connection set for
each prime, including all16384 sets at29. The prime2 is outside the odd-order
circulant-tournament statement.

**Square-sum hostile:** among the first five square-sum orders only
`15=3*5`, `16=2^4`, and `17` are constructible polygon orders. `23` is a
non-Fermat prime; `25=5^2` violates squarefreeness. Thus square-sum Hamiltonicity
does not imply constructibility or multiplier rigidity. For composite15,
multiplier rigidity must not be inflated to full tournament-automorphism
rigidity; the prime-degree hypothesis above is essential.

## 3. A positive Platonic bridge on ten objects—and where it stops

The ten unordered pairs of five vertices are the edges of `K5`. Adjacency by
one common endpoint gives `J(5,2)=L(K5)`, a6-regular graph with30 edges;
disjointness gives its complement, the3-regular Petersen graph with15 edges.
Relabeling is the faithful `S5` action on these ten **positions**.

Icosahedral rotations form `A5`, acting on five inscribed cubes/true crosses;
see [Baez's geometric construction](https://math.ucr.edu/home/baez/mathematical/six.bac.html).
This yields an exact operation-preserving map:

`pair E of five labels -> subgroup C3 of cycles supported on E's complementary triple`.

There are exactly ten such subgroups. In the icosahedral realization they are
the ten opposite-face rotation axes, equivalently the ten opposite-vertex axes
of the dodecahedron. Conjugating the subgroup is exactly relabeling the pair.
Our code verifies all`60*10=600` equivariance cases. The geometric object is
an **unoriented axis**, not one of the20 individual faces. No physical metric
or metagraph adjacency is asserted by this map.

This is a real bridge from tournament **arc positions** to Platonic symmetry,
but it does **not** identify ten tournament **classes** with those positions:
vertex relabeling acts trivially after quotienting by isomorphism, and remains
trivial after further converse-merging. Faithful and trivial permutation
actions cannot be equivariantly isomorphic.

Additional exact five: `A5` has five conjugacy classes, with sizes
`1,12,12,15,20`, verified directly. This five counts conjugacy types of rotations,
not solids or Fermat primes. Its order is`60=2^2*3*5`; `S5` has order
`120=2^3*3*5`. The full icosahedral symmetry group is`A5 x C2`, not`S5`:
their centres have sizes2 and1. A tournament's own automorphism group has odd
order (an involution would reverse an arc), so it cannot itself be`A5`.

## 4. The actual colored n=5 metagraph is rigid on vertices

The inherited source is
`05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/census_n5.json`.
Its edges pair complementary six-bit **path-normal-form tilings**, not
single-arc flips. We independently reconstruct all64 tilings,12 classes,
8 self-converse classes,10 merged nodes, and all32 complement lines in its
explicit tile/reflection gauge. All match the inherited input exactly:
8 blue lines,24 black lines, including4 literal black loops.

The vertex signatures `(blue degree, black degree, blue loops, black loops)`
distinguish every node except pure-blue nodes0 and8. The third pure-blue node7
has blue degree3, unlike their degree1. Nodes0 and8 have blue neighbors5 and6,
whose black degrees/loops differ. Hence swapping0 and8 is impossible too.
**PROVED from the checked adjacency:** the color- and multiplicity-preserving
vertex automorphism group is trivial. This does not deny permutations of
indistinguishable parallel edges fixing every vertex.

Thus an `A5`/`S5` vertex-symmetry proposal fails decisively. Four black loops
already refute identification with either simple Johnson or Petersen graph.
“Pure blue” denotes a reflection/fiber condition; it does not make the three
vertices interchangeable graph-symmetry centres.

**Stronger mask-label hostile (suggested by the inheritance lane):** in this
gauge the pure-blue fibers are`{0}`,`{11,46,53}`,`{16}`. Choosing one mask from
each necessarily gives a pair at Hamming distance1 (`0` and`16`). The triple
`{7,17,55}` has distance multiset`{2,3,3}`. Every signed coordinate permutation
(bit reordering plus arbitrary bit complementation) preserves Hamming
distances. Therefore no such cube-coordinate gauge identifies this triple
with one representative from each pure-blue fiber. Arbitrary external labels
are not excluded; their meaning remains an input question.

## 5. Cellular automata: Fermat kernels, not Fermat primality

Recovered prior source:
`05-knowledge/results/gilbreath_fermat_platonic_20260926.md`, sections1–4.
Use the one-sided binary XOR convention
`b_(t+1)(i)=b_t(i) XOR b_t(i-1)`, elementary **Rule60**. Reflecting the spatial
axis gives Rule102. In this convention the time`t` kernel is

`K_t(x)=(1+x)^t mod2 = product_{i in bits(t)} (1+x^(2^i))`.

This is **not raw, uncompressed symmetric Rule90**. From a single seed, Rule90
has Laurent kernel `(x+x^-1)^t=x^-t(1+x^2)^t`; its occupied positions have the
same parity as`t`. Compressing position`2j-t` to`j` gives precisely`K_t(x)`.
The shift and rescaling must be retained. Already at time1, the rebased raw
Rule90 row is`101` (integer5), whereas the compressed/Rule60 row is`11`
(integer3). The checker verifies the parity-compression identity at every
time0..63 and this minimal hostile. None of these infinite-line kernels is
automatically the six-cell periodic CA used by the separate metagraph lane.

Lift its0/1 coefficients to integers and evaluate at2:

`K_t(2)=sum_{j subset t}2^j = product_{i in bits(t)} F_i`.

The disjoint binary supports prevent integer carries in this product, so the
identity is exact, not an illegitimate characteristic2 evaluation. The first32
rows are the subset products of the first five Fermat primes (including1 as
an arithmetic identity, not a one-sided polygon). At time32 the unchanged rule
produces composite`F5`. The kernel identity holds at every time, regardless of
primality. Direct binomial parity, iterative XOR, and the product formula agree
for all`0<=t<=63` in the test. The CA bridge does not prove Gilbreath's conjecture.

## 6. Prior correction lineage and reproduction

- `07-reflections/the-five-platonic-tournaments.md` has a maintained correction:
  no five-way identity; `Q6` has64 vertices and its polar12 vertices, not the
  six vertices of the3D octahedron. The legitimate octahedron is`L(K4)`.
- `05-knowledge/results/fivefold_boundary_20261018g.md` preserves marked
  character, lattice, and boundary maps, rather than raw cardinalities.
- `01-canon/theorems/THM-4847-fermat-primes-barriers-to-proving-exactly-five-the-function-field-analogue-and-why-tournament-counts-are-blind-while-three-state-counts-are-pepin.md`
  preserves the open Fermat question and separates binary-count blindness
  from the three-state Pépin test.

Run from this directory (standard-library Python only):

```sh
python3 check_five_families.py > checks.json
python3 check_inherited_metagraph.py > metagraph_checks.json
python3 check_a5_balance_meeting.py > a5_balance_normal.json
python3 verify_packet.py
```

The first script's single-arc-flip quotient is intentionally a **different**
operation carrier and is named accordingly; it has21 simple edges after
converse-merging, not the user's32 colored multigraph lines. The inherited
input SHA256 is recorded in `metagraph_checks.json`. All scope boundaries
above are retained; no global Fermat classification or arithmetic-to-pure-blue
bijection has been inferred.
