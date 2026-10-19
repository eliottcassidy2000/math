# Ten labels: Hamming replacements, marked Mathieu actions, and joint boundary information

**Status: PROVED elementary classification and compatibility lemmas;
FINITE-EXACT labelled incidence and group censuses; OPEN geometric realization
and universal Collatz coverage.** Proofs are self-audited, supported by independent
finite constructions. No new theorem ID or historical-priority claim is made.
Session `frontier_petersen_oct18`, continuation `20261018n`.

The main result is a complete ten-label extension of the
[previous matching–Hamming bridge](cardinality_bridges_20261018m.md):

\[
2295=945+1350=9!!+\binom{10}{2}\,30.
\]

The missing 1350 structures each have a **unique distinguished pair** and
an extended Hamming code on the remaining eight labels. A local replacement
connects each matching to ten missing structures and each missing structure
back to seven matchings, preserving a specified four-dimensional code of
observations. The cube matrix in the requested sporadic-sequence papers gives
an explicit second route to the same Hamming code; its all-order norm and
integer-clock test are in the [companion note](cube_unitary_clock_packets_20261018n.md).

## 1. Inheritance and working board

- Closest proved mechanism: the even-subset/complement symplectic model and
  the length-eight Hamming classification in the previous bridge, Sections
  3–5. Binary self-dual length-ten codes having two equivalence types is
  classical; [Harada's classification database](https://www.math.is.tohoku.ac.jp/~mharada/research/codes/sd2.htm)
  independently records this. The proof below is included for transparency.
- Hostile: all single-port patterns can be covered even when joint controller
  compatibility fails. For geometry, the
  [coupled-support census](coupled_support_branches_20261018c.md) already has
  full extensions of all 570 marked boundary patterns; failures of every
  fixed-pair chart do not contradict those witnesses.
- Corrected near miss: [recurrence reanchoring](reanchoring_kernel_20261018c.md)
  preserves the relevant sequence only with its affine carry correction.
  Likewise a new code chart must preserve the selected observations, rather
  than merely dimension or cardinality.
- Least-used sidecar: a marked pair inside the matching, plus the entire bank
  of observation rows. The earlier
  [Mathieu/Golay construction](procgen_sources_20260923_interplay.md#3-the-moonshine-chain)
  supplies explicit generators, not permission to forget this marking.

Anchor: ten-label joint coverage. Niche: boundary relations beyond fixed
binary charts. Wildcard: the signed cube and its norm-preserving operators.
Board: labelled self-dual codes; cubic intersection parity; marked pair;
Mathieu stabilizers; boundary support with fibres; integer-clock Gram forms.

The S44 theorem called Witt rigidity, [THM-4710, Virasoro singular-support
generators](../../01-canon/theorems/THM-4710-witt-rigidity-of-the-degree-n-singular-support-generators-of-virasoro-minimal-models.md),
concerns the Witt Lie action. It is not a theorem about Witt designs or Mathieu
groups. The shared name was checked and does not supply a dependency here.

## 2. Complete classification, with an elementary proof

Let E be the even-weight subspace of F2^10 and let 1 be the all-ones vector.
The quotient V=E/<1> is an eight-dimensional nondegenerate alternating space,
with pairing given by intersection parity. Its Lagrangians lift bijectively
to binary self-dual codes C of length ten and dimension five. A matching gives
the span of its five pair indicators, denoted I2^5.

**Proposition 1 (PROVED).** Every such C contains a weight-two word.

Suppose otherwise. Since 1 belongs to C, complementation excludes weight eight
as well. The remaining weights are 0,4,6,10. On the self-orthogonal code C,
q(c)=wt(c)/2 mod 2 is linear: the correction in the weight-of-sum identity is
the even intersection size. Since q(1)=1, its kernel has dimension four.
All fifteen nonzero words in that kernel have weight four, so its total weight
is 60. But each active coordinate of a four-dimensional binary code is one
in exactly eight words. Thus the total weight must be divisible by eight,
contradicting 60. QED.

Peeling a weight-two pair is exact: orthogonality forces its two coordinates
to agree in every codeword. Hence C=I2 direct-sum C8. The inherited length-eight
classification says C8 is I2^4 or H8. For completeness, in the no-weight-two
case C8 has weights 0,4,8. A basis consisting of 1 and three other rows gives
eight distinct three-bit column addresses: repeated addresses would give a
weight-two word in the dual. Thus C8 is precisely the affine-function code
on F2^3, H8=RM(1,3).

In I2+H8 the sole weight-two word identifies the distinguished pair. There
are 45 choices of pair and 30 labelled H8 codes on its complement, giving
1350. Matching codes have five weight-two words and number 945. The two
families are disjoint. Their sum equals the independent symplectic count
3*5*9*17=2295. The exhaustive span enumeration has levels
`1,255,5355,11475,2295`.

This proves that ten labels introduce no new indecomposable binary self-dual
code. It does not say the ten-label *marked* structures are equivalent.

## 3. A local replacement preserving observations

Start with a matching, retain one pair P, and write the four other pair
indicators as e0,e1,e2,e3. Preserve the code

\[
D=\langle P,e_0+e_1,e_0+e_2,e_0+e_3\rangle.
\]

**Proposition 2 (PROVED).** Exactly three self-dual codes extend D: the
original matching and two codes I2(P)+H8. Each agrees on every observation
in D, including the actual labelled pair P.

Indeed D has dimension four and contains 1, so D-perp/D is a symplectic plane
with three isotropic lines. On the remaining eight coordinates, the preserved
even-pair subcode has dimension three. One extension adjoins an individual
pair. In the other two, every new word chooses one coordinate from each of
the four pairs and has weight four. They are H8 codes. QED.

There are five choices of P and two replacements, so each matching has ten
distinct Hamming neighbours. Conversely, in the affine address model of H8,
the seven nonzero directions a pair x with x+a. The affine functions constant
on these pairs form the preserved three-dimensional subcode. These are all
the matching neighbours: hyperplanes of H8 containing 1 correspond precisely
to the seven nonzero directions. Thus

\[
945\cdot10=1350\cdot7=9450.
\]

For fixed P the 30 H8 codes form a graph: an edge joins the two H8 extensions
associated with one matching on the eight other labels. There are 105 edges.
**FINITE-EXACT:** this graph is bipartite with parts 15 and 15, degree seven,
and is the point–plane incidence graph of PG(3,2). The checker goes beyond
the parameter match: the complements of its fifteen incidence rows, with
zero adjoined, form a four-dimensional binary space. It constructs the
fifteen nonzero point and dual-plane coordinates and checks all 225 incidence
predicates as `point dot plane = 0`.

The 15+15 structure therefore has an actual incidence map. Its labels and
choice of one part as points remain part of the representation.

## 4. The first lost information is cubic

For a bank of binary observation rows r1,...,rk, label column i by its joint
signature sigma_i in F2^k, with multiplicity m_sigma.

**Proposition 3 (PROVED, all even lengths).** A single matching respects all
rows iff every m_sigma is even. Its number is

\[
\prod_\sigma(m_\sigma-1)!!,
\]

with (-1)!!=1. Pair equal signatures independently; any compatible matching
must do exactly this. This is both a test and a constructive decoder.

For two even rows the four signature counts are even iff their intersection
is even. Thus every pair of rows of a self-dual code is matching-compatible.
For three even rows with pairwise even intersections, all eight signature
counts have the same parity, namely the triple-intersection parity. This
follows by Boolean inclusion–exclusion modulo two: the moments of degrees
zero, one and two vanish, leaving the degree-three moment.

The explicit rows, as bit masks on labels 0,...,9, give:

| Successive row bank | Matching extensions | All self-dual extensions |
|---|---:|---:|
| 170 | 45 | 135 |
| 170,204 | 3 | 15 |
| 170,204,240 | 0 | 3 |
| 170,204,240,255 | 0 | 1 |

Their three weight-four rows have pairwise intersections of size two and
triple intersection {7}. The last row resolves the ambiguity between three
possible distinguished pairs. Every pair fits, but there is no common
matching for the triple. This is an exact obstruction to inferring joint
coverage from pairwise tests, not a counterexample to Collatz convergence.

## 5. Eight visible ports and two hidden coordinates

Keep ports 0,...,7 ordered and hide coordinates 8,9. Projection has 16 words
with fibre size two precisely when the hidden pair belongs to the code;
otherwise it has 32 words, each with one lift.

| Family | 16-word projections | 32-word projections |
|---|---:|---:|
| Matching | 105 | 840 |
| Hamming plus pair | 30 | 1320 |

For a specified eight-bit word of weight w, the numbers of supporting codes
are:

| w | Matching | Missing | All |
|---:|---:|---:|---:|
| 0,8 | 945 | 1350 | 2295 |
| 1,7 | 210 | 60 | 270 |
| 2,6 | 135 | 120 | 255 |
| 3,5 | 90 | 180 | 270 |
| 4 | 81 | 174 | 255 |

These counts are PROVED by extension incidence or signature pairing and
independently enumerated. For nontrivial even port weight, the two even
hidden completions are orthogonal: 135+135-15=255. For odd weight they pair
nontrivially and have disjoint extension sets: 270. Extreme words have a
zero or all-ones completion and lie in every projection.

**Every individual port word already occurs in a matching projection.** The
missing structures improve joint relations and fibre information; they add
no new word to this union. This is the decisive hostile against calling the
1350 extra codes new source coverage without specifying a joint predicate.

## 6. Recovering the Mathieu route on ten actual labels

The September construction already checked the ternary Golay code, its
S(5,6,12) hexads, M12, and the outer-twisted S6 hexad stabilizer. The recent
[factorization/vertex-algebra note](factorization_vertex_supersymmetric_bridge_20261018.md)
adds the corrected distinction between moonshine traces and physical K3
symmetries. We reuse its mathematical objects rather than equating either
notion of symmetry with a geometric port permutation.

Fix points 10 and 11 individually in the explicit twelve-point M12 action.
Its induced ten-point group M10 has order 720 and is sharply three-transitive.
This is **FINITE-EXACT** from the supplied generators, checking all ordered
triples. Keeping the two fixed points is essential.

The same action partitions the labelled codes as follows:

| Codes | M10 orbit sizes |
|---|---|
| Matching | 36,45,144,180,180,360 |
| Hamming plus pair | 90,180,180,180,720 |

Derive the Witt design twice: retain hexads containing {10,11} and remove
that pair. The resulting 30 blocks on ten points form S(3,4,10). This is
**not** one Hamming code: its binary span has dimension nine, and some block
intersections have odd size.

There is nevertheless a useful marked map. For a pair Q among the ten
points, the four design blocks containing Q partition the other eight
points into four pairs after Q is removed. Adjoin Q to get a matching.
All 45 resulting matchings are distinct. Apply Proposition 2 at the original
Q. The resulting 90 Hamming-plus-pair codes form exactly the orbit of size
90 in the table. Thus the design selects a small, invariant collection of
replacements **with the marked pair retained**.

No claim identifies M10 with the six-letter S6 action merely because both
have order 720. Nor may these eleven abstract orbits replace a search over
fixed geometric ports without an authenticated action on those ports.

## 7. Exactly four binary charts: a useful six-colour obstruction

Fix a ternary label for each vertex in the six-colour encoding c=3b+4t mod 6.
Edges within one ternary label require opposite binary shades; other edges
are automatic. After eliminating hidden vertices and imposing binary pins,
a nonempty boundary support is an affine binary relation whose independent
directions flip disjoint visible connected components. This is the exact
fixed-label criterion already inherited from the reanchoring note.

**Proposition 4 (PROVED).** If such an affine chart is contained in
C=I2+H8 on ten visible coordinates, it contains at most eight words. Exactly
four charts suffice to cover C, and four are necessary.

For the lower bound, each free visible component indicator belongs to C.
There is only one weight-two word. Every other nonzero word has weight at
least four. Disjoint supports therefore allow at most three free components
(2+4+4), hence at most eight words. C has 32 words.

For the upper bound, use H8=span(170,204,240,255) on the first eight labels,
adjoin 768, and take the four cosets

\[
r+\langle15,240,768\rangle,\qquad r\in\{0,170,204,102\}.
\]

Each coset is a signed equality relation on blocks {0,1,2,3}, {4,5,6,7},
{8,9}. They are disjoint and cover C. The checker replays the original ten
coordinates against the explicit equality/disequality predicates. QED.

This proves a precise limitation of fixed-ternary chart methods and supplies
a minimal four-chart target. It does **not** realize the selector between
the charts by unit distances. A six-colour forcing context must authenticate
that selector and its shared ports, then show incompatibility with a second
context. A nonempty Hamming boundary relation alone forces no seventh colour.

## 8. The next boundary: twelve labels require a new primitive object

The ten-label peeling rule must not be extrapolated indefinitely. There is
a concrete next construction, and the same counting argument closes twelve
labels as well.

For a matching of twelve points with pair indicators e0,...,e5, take the
five-dimensional space E of even unions of pairs. Adjoin a transversal t
choosing one coordinate from each pair. Its weight is six, and its inner
product with every even union is zero, so E+<t> is self-dual. Transversals
give exactly two codes, according to the parity of swaps within pairs.
Their weight enumerator is

\[
1+15z^4+32z^6+15z^8+z^{12}.
\]

The weight-four words are exactly unions of two pairs. Two coordinates have
the same incidence with all these words iff they are the original pair.
Thus the code recovers its matching uniquely, and there are 2*11!!=20790
distinct such codes. They have no weight-two word.

**Proposition 5 (PROVED, classical classification recovered).** The twelve-
label self-dual codes consist of

\[
75735=10395+44550+20790:
\]

matchings; H8 plus two disjoint repetition pairs; and the primitive
transversal codes just constructed. Indeed the second count is
`choose(12,4)*3*30=44550`, with the two weight-two words recovering the pairs.
The three families have respectively six, two and zero weight-two words,
so are disjoint. Their proved injections saturate the independent total
`3*5*9*17*33=75735`. No further family can remain. QED.

The same even-pair/transversal operation works for any even number k of
pairs: its enumerator is
`sum_(j even) choose(k,j) z^(2j) + 2^(k-1) z^k`.
For odd k the transversal has odd weight and fails self-orthogonality. At
k=4 the new transversal words themselves have weight four, so they mix with
the old pair-union words and the matching has seven possible decodings.
At k>=6 even, weight-four incidence recovers the matching uniquely. This
explains the exceptional Hamming multiplicity through a collision of
weight layers. It does not classify all codes at higher lengths.

The receipt checks representatives k=4,5,6,8, the twelve-label injection
counts and the unique tetrad decoder. It does not enumerate all 75735
length-twelve codes; completeness there follows from the injection proof
and the mass formula.

## 9. Research transfers and next decisive tests

| Source → target | Map and preserved predicate | Lost information / required sidecar | Next decisive test |
|---|---|---|---|
| Matching → Hamming controller | Extend fixed D; preserve all rows in D and pair P | A row outside D may change; retain full bank | Classify which actual source/port observables lie in D before switching |
| M12 → ten-label search | Pointwise two-point stabilizer; retain derived incidence | Physical coordinates are not supplied | Test a candidate's ordered port tensor, not only its abstract graph orbit |
| Fixed binary chart → Hamming relation | Four affine cosets, original ten ports | Chart selector and ternary labels | Seek one geometric selector realizing all four charts with no extra boundary words |
| Signed cube → ROOT packet | Exact clocked Gram identity in companion note | Source-specific amplitudes and erasure costs | Compute clock curls of candidate reset packets before a depth census |

The ROOT source is still the actual Eisenstein charge A_N(s), with first-arrival
seed, native depth, moving complement and integer edge clock retained, as in
[THM-4690, first-arrival ladder and clock decompositions](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md)
and [the source-energy note](root_eisenstein_first_arrival_20261018j.md).
The cube supplies an independently controlled model; the missing obligation
is a source-preserving realization of actual ROOT packets within that model.
Universal entry and the uniform source-energy target remain OPEN.

## 10. Reproduction and audit boundaries

```bash
python3 04-computation/experiments/ten_label_joint_carriers_20261018n.py --output /tmp/ten-label.json
python3 -O 04-computation/experiments/ten_label_joint_carriers_20261018n.py --output /tmp/ten-label-O.json
cmp /tmp/ten-label.json /tmp/ten-label-O.json
```

[Script](../../04-computation/experiments/ten_label_joint_carriers_20261018n.py),
[receipt](ten_label_joint_carriers_20261018n.json). All 2295 labelled codes;
262144 ordered even-row pairs; all 32768 ordered triples in a representative
H8+I2 code; all 9450 local replacements; all marked projection fibres; complete
M10 orbits. The finite universe has no hidden arithmetic or geometric filter.
The signature proof covers arbitrary row banks, while only the declared
finite banks are enumerated. Standard-library implementation; checks remain
active under optimization. Exact counts and proof statements are separate
from the OPEN transfer proposals above.
