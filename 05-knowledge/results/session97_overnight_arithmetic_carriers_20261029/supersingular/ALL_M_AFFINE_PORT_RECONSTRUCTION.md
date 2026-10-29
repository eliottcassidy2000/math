# Binary four-port reconstruction at every dimension: an exact growing-memory law

**Session97, 2026-10-29. Status:** PROVED all-m finite-geometric theorem;
FINITE-EXACT exhaustive m4/m5 controls. Classical affine-geometry and
Reed–Muller facts are credited, not claimed as new discoveries. The contribution
here is their explicit use as a sharp observer threshold and closed carrier
for the repo's marked-source/joint-port programme.

This is a new file; the earlier finite length32 and lattice checkpoints are
unchanged. It generalizes their six-code affine sector, not their exceptional
two-code sector or the rank32 lattice conclusion.

## 1. The object, observer and exact response formula

Let m≥4 and let X be a fixed set of N=2^m **named** points equipped with an
affine F2^m structure. Let B be the family of all affine subspaces of
codimension2, each of size2^(m-2). For a nonempty named subset S⊆X define

`R(S)=#{B in B : S subset B}`.

If the affine span of S has dimension d, then

`R(S)=GaussianBinomial(m-d,2;2)`.                         (1)

We set the Gaussian coefficient to zero when its upper argument is less
than2. Formula (1) includes subsets spanning a space too large to lie in a
codimension-two flat. The empty subset is different:
`R(empty)=4*GaussianBinomial(m,2;2)`.

**Proof.** Translate one point of S to zero. Any containing flat must pass
through it and is determined by an `(m-2)`-dimensional linear subspace
containing the d-dimensional direction span. Modding out by that span leaves
an `(m-2-d)`-subspace of a vector space of dimension m-d. The Gaussian count,
using duality, is `[m-d choose2]_2`. This counts labelled subspaces, not
isomorphism classes or source automorphism orbits.

## 2. A sharp fourth-order threshold

Over F2, every two distinct points span a line, and every three distinct
points span an affine plane: a binary line has only two points. Therefore
all named single-, double-, and triple-port responses are constant:

```
|S|=1: R(S)=[m choose2]_2;
|S|=2: R(S)=[m-1 choose2]_2;
|S|=3: R(S)=[m-2 choose2]_2.
```

Four distinct points span either a plane or a three-space. Their two possible
responses are respectively

`A_m=[m-2 choose2]_2`,  `B_m=[m-3 choose2]_2`,

and `A_m>B_m` for every m≥4. Thus the weighted four-port oracle recognizes
**exactly the affine planes**.

Mark one point o. For distinct a,b different from o, there is a unique fourth
point c for which `R({o,a,b,c})=A_m`. Define `a+b=c`, set `a+a=o`, and set
`a+o=o+a=a`. In any coordinate chart this is the translated XOR operation
`a+b = a XOR b XOR o`; consequently it reconstructs the entire vector-space
addition based at o. Different choices of o recover the same unbased affine
structure. No chart or group law has to be guessed from a count.

Conversely, an affine chart computes every response, of any subset size,
by taking an affine rank and applying (1). The full four-port table is
therefore a **complete carrier for all future incidence queries of this
specified kind**. Observations through order3 carry no information at all
about which labelled affine structure is present.

This is not a theorem about arbitrary future operations on tournaments,
lattices, trees or arithmetic states. Such applications require a typed map
of their operations into these incidence queries.

### Constructive query bound

A full four-port table is not needed to reconstruct a chart. Start from o,
choose a first basis point, and maintain the affine span already labelled
by binary words. Choose a new point b outside it. For each old nonzero x,
find the unique new point b+x by scanning candidate fourth points with the
four-port oracle. Label the new coset, then repeat. There are
`sum_(k=1)^(m-1)(2^k-1)=N-m-1` such searches, each with at most N candidates.
Thus `O(N²)` oracle calls suffice, and the recovered chart answers all later
incidence queries exactly. This is a simple upper bound, not a claim of
optimal query complexity.

## 3. Exact growing alignment memory

The affine automorphism group is

`AGL(m,2)=F2^m semidirect GL(m,2)`,
`|AGL(m,2)|=2^m product_(i=0)^(m-1)(2^m-2^i)`.

The number of distinct affine structures on a fixed named N-point set is

`A(m)=N! / [2^m product_(i=0)^(m-1)(2^m-2^i)]`.           (2)

**Why the stabilizer is exactly this group.** Affine maps preserve the
four-point plane relation. Conversely, any permutation preserving that
relation can be composed with a translation to fix o. The reconstruction
above then forces it to preserve addition, so it is linear. Hence no larger
permutation stabilizer is hidden by the observer bank.

If o is fixed in advance, the equivalent based-space count is
`(N-1)!/|GL(m,2)|`, the same integer. Counting choices of o *again* would
incorrectly multiply (2) by N. A coordinate basis is another redundant chart
choice, not extra geometry visible to the responses.

Every one of these A(m) structures has exactly the same responses of order
at most3, but its fourth-order response table is distinct. Therefore any
finite fixed-length carrier capable of answering all future named incidence
queries needs at least

`ceil(log2 A(m))` bits

in the worst case. A faithful finite-state encoding can achieve this bound
in principle by indexing the structures; no efficient ranking algorithm is
claimed. An explicit point-to-bit-vector chart is a more convenient
`O(m2^m)`-bit carrier, with translation/basis redundancy understood.
Stirling's estimate and the group-order formula give

`log2 A(m)=m2^m-(log2 e)2^m-m²+O(m)=Theta(m2^m)`.

This is exact growing alignment memory invisible to every three-port
summary. It does not assert that a particular yes/no decision requires that
much information: a decision certificate can still be short.

| m | Flats | Point/pair/triple responses | Plane/independent-four responses | A(m) | Fixed-length lower bound |
|---:|---:|---|---|---:|---:|
| 4 | 140 | 35,7,1 | 1,0 | 64864800 | 26 bits |
| 5 | 620 | 155,35,7 | 7,1 | 822336494953469303808000000 | 90 bits |

The 90-bit count concerns **all** labelled affine geometries on32 points.
The earlier cyclic-Frobenius constraint restricted us to six such geometries
inside an eight-kernel bank, so its much smaller source-choice count is not
in conflict with this number.

## 4. Why m≥4 is a necessary boundary

For m3, codimension-two flats are binary affine lines, hence simply **all
unordered pairs** of the eight points. There are28 flats and the response
formula is7,1,0,0 on subsets of sizes1,2,3,4. The entire incidence observer
bank—not merely its first three orders—is unchanged by every permutation.

There are nevertheless `8!/|AGL(3,2)|=30` labelled affine structures. They
are all invisible to this particular family of flats. The response stabilizer
is S8, strictly larger than AGL(3,2), so reconstruction fails. This supplies
a sharp hostile to extending the theorem to m3 by blindly extrapolating the
Gaussian formula. The formal source-geometry count of30 would need five
bits to retain the geometry, but the incidence responses alone contain zero
bits about that choice. Nothing in this statement contradicts the earlier
Fano reconstruction, which used a different observer family.

## 5. The Reed–Muller interpretation, with its limit

The binary Reed–Muller code RM(2,m) evaluates Boolean polynomials of degree
at most2 on F2^m. Its minimum nonzero weight is2^(m-2), and its minimum-word
supports are exactly the codimension-two affine flats. These are classical
facts; the standard code/affine-space connection is covered in
[Elkies's coding-theory lectures](https://people.math.harvard.edu/~elkies/M256.13/index.html), November7.
A short quadratic-form proof makes the needed scope explicit.

The polar bilinear form of a quadratic Boolean polynomial has even rank2r.
After an invertible linear coordinate change, its quadratic part is
`x1*x2 + x3*x4 + ... + x_(2r-1)*x_(2r)`.
If the remaining linear term is nonzero on the radical, the word is balanced.
Otherwise the character sum has magnitude2^(m-r), so its weight is
`2^(m-1) plus or minus 2^(m-r-1)`. A nonconstant affine polynomial is also
balanced. The smallest possible nonzero weight is attained exactly at r1
with the smaller sign. Completing the two-variable product shows that these
words are products of two independent affine linear forms. Their supports
are exactly affine codimension-two flats. This proves the claimed minimum
and its complete support classification.

Thus the all-m incidence theorem is a uniform statement about observations
of minimum words in RM(2,m). But the special Construction-A conclusion of
the previous packet does **not** extend automatically: RM(2,m) has dimension
`1+m+C(m,2)`, which equals half its length only for m5 among m≥4. Its
self-duality and the resulting rank32 even-unimodular bridge were specific
to that dimension. Neither the memory theorem nor the Gaussian recurrence
requires self-duality.

## 6. Other alphabets: triples can already see geometry

Over F_q with q>2, three distinct points may be collinear. For codimension-two
flats in F_q^m, m≥3, the same proof gives

`R(S)=[m-d choose2]_q`.

A collinear triple gives `[m-1 choose2]_q`; a noncollinear triple gives
`[m-2 choose2]_q`, and these differ. Thus the special need for a fourth port
comes from **binary lines having two points**, not a universal significance
of the number four.

This is only a collinearity detector statement. We do not assert here that
it identifies scalar labels or a preferred field chart. Field automorphisms
preserve affine incidence, so semilinear ambiguity must be retained in any
full reconstruction over a nonprime field. The binary all-m theorem avoids
that ambiguity because F2 has no nontrivial field automorphism.

## Reproduction and publication boundary

`check_all_m_affine_ports.py` exhausts all codimension-two flats for m4 and
m5; checks (1) on every named subset of size1 through4; reconstructs addition
at every possible origin; exhausts all2048/65536 RM(2,m) codewords to recover
their minimum supports; and tests a nonaffine label-swap hostile. It also
checks the sharp m3 boundary and all2925 triples in AG(3,3).
The script makes142843 explicit checks, identical under normal Python and
`-O`. This finite audit supplements, rather than replaces, the all-m proof.

No existing file in the earlier supersingular packet was changed for this
extension. The stopping boundary is this proof, code and independent audit;
no additional research branch is opened here.
