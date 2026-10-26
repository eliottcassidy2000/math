# Session 68: carousel classification, finite-field descent, and aligned CRT gluing

**Status:** elementary proofs below, with the normal-basis theorem **CITED**;
explicitly bounded checks **FINITE-EXACT**. This is an exact classification of
named composition families, not a classification of every tournament or a
claim of global novelty. Primary proofs were written in this lane; the
relative-phase construction and converse sharpening were supplied jointly
with the root. Root's independent field/CRT controls are separate.

## 1. Provenance, scope and random recovery

The read-only repository source pin is
`fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8`, verified before sampling.
`AGENTS.md`, navigation, active guardrails, the research protocol and the
core-paper guide were read. No repository script was imported or executed;
no repository files, branches or refs were changed by this lane. The startup
script was not run because this assignment was explicitly read-only.

Before inspecting sampled content, Python `random.Random(2026102668)` drew
six paths without replacement from each sorted stratum, in this order:

| Stratum | Eligibility at the pinned tree | Size |
|---|---|---:|
| Theorems | `01-canon/theorems/`, suffix `.md` | 3976 |
| Reflections | `07-reflections/`, suffix `.md` | 3337 |
| Results | `05-knowledge/results/`, suffix `.md` | 5580 |
| Source | `04-computation/`, suffix `.py` | 18224 |

This is uniform sampling within strata, not uniform sampling of all repository
files. The seed, draws, universe hashes and full sorted universe lists are in
`random_manifest.json` and `*_universe.txt`. Drawn files were inspected in
bounded prefixes (normally 75 lines), followed by targeted source/correction
reading only where useful. `draw_excerpts.txt` retains bounded extraction
provenance, not a proof authority. No draw was replaced after reading it.

The sampled Parker CRT diagram packet directly led to the aligned collision
and relative-phase bridge. The sampled digit/rotation note led to a check of
which coordinate system actually conjugates rotation; this supports the
normal-basis requirement, not a Collatz theorem. The random Paley and
Vandermonde material prompted strict matrix/quotient-loss checks, rather than
an inferred Paley classification.

**Inherited corrections recovered before reuse:**

* `04-computation/new_frontiers_s3b.py` prints the old good-cut generating
  function missing a factor `1/(1-x)`. The 2026-10-23 correction in
  `01-canon/MISTAKES.md` and
  `01-canon/theorems/THM-341-bivariate-GF-good-cuts.md` give
  `x B(xy)/((1-x)(1-x-x B(xy)))`; the old expression gives 1 instead of 2 at
  `n=4,k=2`. Its historical match claim was not imported.
* Sampled `THM-896-kendall-wei-perron-atlas.md` contains its own n=7
  cospectral failure and subsequent qualified addenda. Neither its title nor
  its finite n<=7 tower claim is an all-order classifier.
* `THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels.md`
  supplies a valid signed polynomial identity, not reconstruction from scores.
  The 2026-10-23 bosonic `THM-1810` correction preserves that identity while
  rejecting the conflation of the unsigned edge product with a Vandermonde
  permanent (at three all-one variables the counts are 8 versus 6).
* The sampled path-homology sources were not used to infer classification.
  MISTAKE-021 distinguishes the GLMY complex from the regular/interior-boundary
  complex. A common name or apparent even-degree vanishing is not a map.
* Sampled `THM-819-primitive-harmonic-law-for-interval-cores.md` already repairs
  its earlier pointwise confinement to confinement outside a finite null
  boundary; no pointwise version was inherited.
* The sampled Parker packet explicitly distinguishes extension fields from
  quotient rings and modular magic sums from ordinary integer magic sums.
  Those boundaries remain in force.

All other draws and their dispositions are recorded in `sample_dispositions.md`.
The incoming remote work reported at `d761e4cf0ff` includes
`tournament_composition_classification_20261026/` with canonical tree,
4p-threshold and Sylow-adapted/base-p results. Those are prior/concurrent work,
not new claims of this lane. A later incoming commit, `c9a8823436e`, adds
`05-knowledge/results/cyclictriple_composition_20261026/README.md`; its README
was read before final integration, without executing its code. It already
proves the carousel ordinary and full labelled fibres, Möbius necklace
automorphism strata, converse pairing, `Aut(H)=S3^m semidirect D_(2m)` and
`2^(2m+1)` labelled lifts, plus further metagraph results. Accordingly the
carousel results below are an independently rederived concurrent component,
not claimed new relative to that incoming work. The distinct increment of
this lane is the explicit normal-basis/Frobenius field-action dictionary,
its hostile boundaries, and the aligned CRT relative-phase law.

## 2. A prime quotient with exactly cyclic automorphisms

Let `m=2k+1>=3`. The carousel tournament `R_m` has vertex set `Z/mZ` and

    i -> j  iff  j-i mod m is in {1,...,k}.

This is the classical `T_(2k+1)` family. The definition and its place among
critical indecomposable tournaments are recorded in Belkhechine--Boudabbous,
[Section 2](https://arxiv.org/html/1007.3049v1). We do not need the full critical
classification: the following proofs establish the precise properties used.

**Lemma 2.1. `R_m` is prime.** A module is a vertex set viewed uniformly by
every outside vertex. If a proper module `M` contains at least two vertices,
choose an outside vertex `x`. Then `M` is contained in one of the two cyclic
intervals `N+(x),N-(x)`, each of length `k`. Rotate coordinates so this interval
is represented as consecutive integers. Choose `a<b` in `M`. The vertex
`z=a+k+1 mod m` lies outside that length-k interval, while `z->a` and `b->z`.
It distinguishes two vertices of `M`, a contradiction. Thus only singletons,
the empty set and the whole set are modules.

**Lemma 2.2. `Aut(R_m)=C_m`, not merely a subgroup containing `C_m`.** An arc
`i->i+d`, `1<=d<=k`, has exactly `k-d` common out-neighbours. The arcs maximizing
this statistic are precisely `i->i+1`. Any automorphism therefore preserves
the directed successor cycle and is a rotation. Every rotation is an
automorphism, proving equality.

**Boundary:** even `m` is excluded. A faithful even-order cyclic group cannot
be a tournament automorphism group: its involution swaps some pair and would
reverse that pair's arc. The existence of `F_(2^m)` for even `m` does not create
an even carousel counterpart.

## 3. Exact binary-decoration classification

Let `A=1 -> C3`, with the source at local label 0 and cyclic triple at labels
1,2,3. Write `A^0=A`, `A^1=A^op`. Their outdegree multisets are respectively
`{3,1,1,1}` and `{0,2,2,2}`, so they are not isomorphic; both have automorphism
group `C3`. On the fixed labelled set `Z/mZ x {0,1,2,3}` define

    T_epsilon = R_m[A^(epsilon_0),...,A^(epsilon_(m-1))],
    epsilon in {0,1}^m.

### Canonical recoverability of blocks

Every displayed block is a module. Suppose a proper module meets at least two
blocks. Its set `I` of met block indices is a module of `R_m`: any untouched
block sees the module uniformly. Primeness forces `I` to be all of `R_m`.
If the original module is proper, some block contains a vertex outside it.
That outside vertex sees all met vertices in every other block uniformly,
forcing its quotient vertex to dominate all other quotient vertices or lose
to all of them. Carousel regularity forbids both. Thus every proper module
is contained in one block. The blocks are exactly the maximal proper modules
and, in particular, the maximal proper strong modules.

Consequently an isomorphism induces a rotation of the block indices, and a
block can map only to a block with the same binary type. Conversely every
such rotation lifts. We obtain:

**Theorem 3.1.** `T_epsilon` and `T_eta` are isomorphic exactly when `epsilon`
and `eta` are rotations of one another. If `d|m` is the least period of
`epsilon`, then

    Aut(T_epsilon) = C3^m semidirect C_(m/d),
    |Aut(T_epsilon)| = 3^m m/d.

The quotient factor permutes the internal `C3` factors by rotation through
`d` block positions. This is a split extension: equal-type blocks have the
same labelled model and a quotient rotation lifts by the identity local map.

### Labelled cyclic-triple fibre

Dualizing a block preserves its cyclic-triple subsets. A triple meeting two
blocks in a 2+1 split is transitive; a triple meeting three blocks is cyclic
exactly when its quotient triple is cyclic. Thus every one of the `2^m`
labelled words has exactly the same labelled cyclic triples, not merely the
same count. Since `R_m` is regular,

    c3(R_m) = binom(m,3) - m binom(k,2) = m(m^2-1)/24,
    c3(T_epsilon) = m + 4^3 c3(R_m) = (8m^3-5m)/3.

This lane proves a complete classification **within this binary-decoration
family**. Exhaustion of all realizations of its common 3-uniform cyclic
hypergraph is supplied separately by the concurrent `c9a8823436e` packet,
using the cited realization/modular-decomposition theorem. That stronger
statement is not silently inferred from this construction or reproved here.

## 4. The divisor lattice is the symmetry classification

For a divisor `d|m`, let

    M_d(2) = (1/d) sum_(e|d) mu(e) 2^(d/e).

The number of length-d binary words with least period d is
`P(d)=sum_(e|d) mu(e)2^(d/e)`: every word has a unique least period, so
`2^d=sum_(r|d)P(r)` and ordinary Möbius inversion applies. A word repeated to
length m has a rotation orbit of size d. Therefore `M_d(2)` is exactly the
number of tournament isomorphism classes with automorphism order `3^m m/d`.
The total is the binary necklace number

    sum_(d|m) M_d(2) = (1/m) sum_(e|m) phi(e) 2^(m/e).

| m | Vertices | Common c3 | Least-period strata: (d, number of classes, Aut order) |
|---:|---:|---:|---|
| 3 | 12 | 67 | (1,2,81), (3,2,27) |
| 5 | 20 | 325 | (1,2,1215), (5,6,243) |
| 7 | 28 | 903 | (1,2,15309), (7,18,2187) |
| 9 | 36 | 1929 | (1,2,177147), (3,2,59049), (9,56,19683) |
| 11 | 44 | 3531 | (1,2,1948617), (11,186,177147) |

The n=12 anchor is `epsilon=000` versus `100`. Independent adjacency-based
backtracking gives automorphism orders 81 and 27; subset-DP gives Hamilton
path count 462915 for **all eight** words. We make no all-m Hamilton-path
claim in this lane merely from equal cyclic triples.

Prime m has only constant and aperiodic strata. Composite m genuinely has
intermediate stabilizers: at m=9, a nonconstant period-3 word retains a
quotient `C3`. Treating every nonconstant word as free would miss it.

### Converse and the merged quotient

Index reflection `i -> -i` is an anti-automorphism of `R_m`. Composing any
anti-automorphism with it gives an automorphism, so these are exactly the
reflections of the dihedral group. Global converse sends

    epsilon_i  to  1-epsilon_(-i).

A self-converse tournament would require equality with some rotation of this
word, equivalently `epsilon_i=1-epsilon_(c-i)` for every i. Since m is odd,
`i -> c-i` has a fixed index, an impossibility for a binary bit. Hence **none
of this family is self-converse**. Ordinary classes pair freely under
converse, preserving their least period and automorphism order. The merged
class counts are 2,4,10,30,94 for the displayed m. This is a statement about
these nodes and their converse pairing, not their complete metagraph edges.

## 5. A literal normal-basis/Frobenius bridge

Choose a normal basis `beta_i=alpha^(2^i)`, `0<=i<m`, of `F_(2^m)/F2`.
Existence is the classical normal-basis theorem; see Keith Conrad,
[Theorem 3.7 and its proof](https://kconrad.math.uconn.edu/blurbs/galoistheory/linearchar.pdf).
Define the **basis-dependent** bijection

    epsilon -> x_epsilon = sum_i epsilon_i beta_i.

Moving the word one position cyclically corresponds exactly to squaring the
field element. Thus the dictionary preserves the action, not just counts:

| Tournament-decoration datum | Field datum |
|---|---|
| cyclic rotation | Frobenius `x -> x^2` |
| ordinary tournament isomorphism class in this family | Frobenius orbit |
| least word period d | `[F2(x):F2]=d` |
| quotient automorphism group `C_(m/d)` | `Gal(F_(2^m)/F2(x))` |
| minimal-period necklaces of period d | monic irreducible polynomials over F2 of degree d |

The last row follows because a field element's distinct Frobenius conjugates
are exactly the roots of its minimal polynomial. Degrees occurring in
`F_(2^m)` are precisely divisors of m. This is a noncanonical bijection of
classified objects, with an exact action/stabilizer correspondence; it is not
an identification of arbitrary tournament substitution with field addition
or multiplication. The normal basis and the fixed labelled quotient are the
necessary coordinate sidecars. Passing from an element to its minimal
polynomial forgets the root/rotation phase.

### Explicit F8 witness, and hostile tests

Let `a^3+a+1=0` over F2 and take `alpha=1+a`. The normal basis is

    (alpha,alpha^2,alpha^4) = (1+a, 1+a^2, 1+a+a^2),

represented by polynomial bit masks `(3,5,7)`. Their sum is 1 and they are
linearly independent. The four necklace classes correspond to

| Word representative | Element type | Minimal polynomial | Aut order |
|---|---|---|---:|
| 000 | 0 | X | 81 |
| 111 | 1 | X+1 | 81 |
| 100 | alpha | X^3+X^2+1 | 27 |
| 110 | alpha+alpha^2 | X^3+X+1 | 27 |

These are orbit representatives, not an ordering claim about roots.

* **Bad-basis hostile:** `(a,a^2,a^4)` has bit masks `(2,4,6)` and sum 0.
  It is not a normal basis; its encoding collapses 000 and 111.
  A polynomial basis `(1,a,a^2)` is a basis but does not intertwine cyclic
  rotation with Frobenius: 100 would map to 1, despite word period 3.
* **Trace/norm hostile:** 111 and 100 both have field trace 1 and norm 1,
  yet their degrees are 1 and 3 and their tournament automorphism orders
  are 81 and 27. Norm or even the trace--norm pair cannot replace descent degree.
* **Reflection hostile:** normal-basis index reflection J is F2-linear but
  not generally a field automorphism. Here `J(a^2)=a^2+a` while
  `J(a)^2=a^2`. Converse corresponds to `x -> 1+J(x)`, not a field conjugation.
  In general J conjugates Frobenius to its inverse; an actual finite-field
  automorphism would commute with Frobenius, impossible for m odd >=3.
* **Hierarchy hostile (root's independent packet):** on nine leaves,
  `C3 wr C3` has order 81, not the cyclic order-9 Galois group. For its
  9-cycle `(0,3,6,1,4,7,2,5,8)`, the words 000000001 and 001001011 both
  have cyclic period 9 but have full hierarchical stabilizer orders 9 and 1.
  Thus equal normal-basis field degree does not classify the recursive
  ternary-tree family. The missing coordinate is its larger hierarchical
  permutation action. This control is independently checked in root's
  `check_field_bridge.py`; it is not a new check claimed by this lane.

## 6. CRT means aligned gluing, not independent local quotients

The random draw
`05-knowledge/results/session55_parker_fields_integer_lifts_20261025/diagrams/README.md`
led to its `crt-wiring/THEOREM.md` and `CLASSIFICATION.md`. Only the elementary
aligned collision identity is used, not their full Parker classification.

For coprime moduli r,s and fixed cell labels, CRT gives

    x_i=x_j mod rs iff (x_i=x_j mod r and x_i=x_j mod s).

Thus the collision relation of the joint array is the intersection of the
two labelled local relations. The explicit row-major witness modulo 35 is

    0   9   1
    16 15  14
    29 21  30.

All entries are squares modulo 35; all eight line sums are 10 modulo 35;
all nine entries are distinct. The reductions modulo 5 and 7 each have
9 unordered colliding cell pairs, with disjoint collision sets. Neither
modulus has nine distinct square residues (only 3 and 4, respectively).
Consequently the two Boolean local Parker flags cannot reconstruct the
joint existence result. Ordinary row sums are 10,45,80, so this is **not**
an integer magic square of squares.

The exact tournament version retains the same missing coordinate: relative
phase. Fix two individual cyclic word orbits, of least periods d,e dividing
m. Their product has d*e labelled pairs. Under simultaneous rotation, every
pair has least period lcm(d,e): a rotation fixes the pair exactly when it
fixes both components. Therefore

    joint stabilizer order = m/lcm(d,e),
    number of joint rotation orbits = gcd(d,e).

Equivalently, phases `(r mod d,s mod e)` are identified by simultaneous
translation; their remaining invariant is `s-r mod gcd(d,e)`. When d,e are
coprime, CRT gives one gluing orbit; when they are not coprime, independent
necklace classes forget exactly this finite phase choice. Example: d=e=3
has three possible gluings, while d=3,e=5 inside m=15 has one.

### Realization as genuine tournament substitution

The four distinct block types

    B_(epsilon,eta) = C3[A^epsilon,A^eta,1]

have 9 vertices, 18 identical labelled cyclic triples, and automorphism
group `C3 x C3`. The unique singleton maximal strong block fixes every
quotient rotation; hence the two four-vertex positions cannot be exchanged,
and the four binary pairs really are four different isomorphism types.
Composing `R_m` with these four colours turns the joint-colouring law above
into a theorem about actual tournaments, with automorphism order

    9^m m/lcm(d,e).

The same pair of individually classified binary necklaces can therefore
produce gcd(d,e) distinct joint tournament types. No claim about all possible
four-colour tournaments is implied.

Under the normal-basis dictionary the gluing also concerns Frobenius orbits
of pairs, so its preserved target is the diagonal action and its stabilizer.
The classical finite-field tensor identity
`F_(2^d) tensor_F2 F_(2^e) = product_(gcd(d,e)) F_(2^lcm(d,e))`
is a compatible interpretation: factorization over F_(2^d) partitions the
e conjugate roots into gcd(d,e) orbits of length e/gcd(d,e). It does not
remove the phase choice or make the tensor product a single field when the
gcd exceeds one. This interpretation is a proved algebraic corollary of the
same Frobenius orbit calculation, not an additional numerical experiment.

## 7. Verification boundary and reproducibility

`check.py` is dependency-free and uses explicit `require` checks, not Python
assertions. It imports no repository implementation. Run:

    python3 -B check.py > receipt.json
    python3 -O -B check.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The byte-identical receipts verify:

* exhaustive proper-module checks for R_m at odd m=3,5,7,9,11;
* independent full permutation Aut enumeration for R_m at m=3,5,7;
* every binary word and every rotation orbit at those five lengths, matching
  Möbius strata and total orbit mass;
* equality of the actual labelled cyclic-triple sets for every expanded word
  at m=3,5,7 (not every word at m=9,11);
* all eight n=12 expanded automorphism orders by adjacency backtracking,
  Hamilton counts by subset DP, and canonical strong-module blocks by full
  subset enumeration;
* converse pairing on every necklace through m=11;
* all eight explicit F8 encodings, Frobenius equivariance, degrees, traces,
  norms, bad normal basis and nonmultiplicative-reflection hostile;
* the CRT35 square/magic/collision witness and ordinary-row-sum hostile;
* all four actual B blocks' automorphism orders and canonical modules;
* 33 explicit period-pair/phase tests for m=3,5,9,15.

These checks are corroboration, not the proofs of the all-odd-m statements.
No unbounded search or external mathematical software was executed.
`random_manifest.json` and `replay_random.py` reproduce the draw at the pinned
Git tree. `SHA256SUMS` identifies the frozen authored report, scripts and receipts.

**Open frontier:** this does not classify arbitrary prime tournament quotients,
all cyclic-triple realization fibres beyond the named concurrent carousel
classification, arbitrary metagraph incidences, or all hierarchical automorphism
representations. A normal-basis field model is exact for the named cyclic
quotient action and ceases to be complete when the acting group grows.
