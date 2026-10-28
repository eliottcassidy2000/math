# A four-bit Fano–Paley sidecar, and an exact arithmetic transport

**Status.** PROVED elementary field/code maps and operation laws; FINITE-EXACT
complete order-seven classification and its consequences, independently
cross-checked within the producer by exact cover, permutation orbits, all line
orientations, and two Hamilton-path algorithms. No claim about unbounded
integer dynamics, idoneality, or a canonical natural-number orientation.

## Result first

Fix a labelled Fano plane F on seven points. Of the 30 labelled STS(7), exactly
eight are block-disjoint from F. Every disjoint pair F,G is the complete set
of cyclic triples of exactly two tournaments, which are opposite Paley7
orientations. Thus a fixed Fano plane admits **16 compatible Paley orientations**:
three bits choose its disjoint partner, and one bit chooses global reversal.
“Bits” here describes a finite choice space (or entropy under the uniform law),
not extra randomness manufactured by the geometry. Bit coordinates require a
chosen chart.

There is also an actual finite arithmetic map, not just matching counts. Match
the seven prime-bank columns to the **normals of the seven Fano lines**, choose
one baseline Paley tournament, and flip a line's three arcs when its matching
prime is selected. For all 128 subsets:

    sigma(product of selected primes) is square
        iff the resulting tournament is Paley7.                 (1)

The 16 positive cases have H=189; the 112 negative cases have H=171. On the
positive locus, squarefree arithmetic parity mu=lambda is exactly the chosen
Fano-preserving GL(3,2)-chirality sign. It is not the ordinary isomorphism class:
all sixteen Paley tournaments are isomorphic after unrestricted relabelling.

## Inheritance and what this note adds

Source snapshot: `83b601158076669edc8729ce461f150da1c3877e`. The pinned
public sources below supply the inherited mechanisms and conventions.

* [arithmetic_braids2_20260917_fano_code.md](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/arithmetic_braids2_20260917_fano_code.md), sections 1–2, already proves the incidence-image Hamming
  code and affine orientation locus. Its Cayley–Dickson baseline mask is 96
  in its own convention; this note instead uses a declared Paley baseline,
  so **relative** mask zero is Paley. Those absolute masks must not be mixed.
* [arithmetic_braids2_20260917_squarefree_symmetry.md](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/arithmetic_braids2_20260917_squarefree_symmetry.md), section 5, already owns the 16/112 split, H values,
  and two GL(3,2) orbits of eight exchanged by reversal.
* [collatz_mod6_20260922_paley_fano_octonion_design.md](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md), S6–S7, already proves the two cyclic Fano planes
  `dev{0,1,3}` and `dev{0,1,5}` and checks their unique decomposition. Its
  S3 and S14–S18 also own the Hamilton counts and orientation census.
* [incidence_attractors_20261028/THEORY.md](https://github.com/eliottcassidy2000/math/blob/83b601158076669edc8729ce461f150da1c3877e/05-knowledge/results/incidence_attractors_20261028/THEORY.md), A3, owns the exact prime bank, F8 chart,
  Frobenius witness 70 -> 11837, and the 168-versus-21 warning.
* Session85's [ARITHMETIC_FREE_BITS.md](ARITHMETIC_FREE_BITS.md), “A different target: one global
  bit recoverable from any Fano line,” owns the parity-access theorem. Here
  it is transported to **Paley chirality**, with the dual coordinate type explicit.

The contribution here is the assembled **typed arithmetic/line-flip map**,
its exact subset-product operation, the disjoint-partner quotient, the
four-bit normal form, and the resulting chirality decoder. The inherited
counts and Hamming locus are credited rather than claimed as discoveries.

**Concept board:** primal prime columns / dual line normals / affine orientation
masks / disjoint STS partners / arithmetic and chirality parity. The canonical
hostile is any single line flip: it preserves regularity and all seven
original Fano cycles but leaves Paley (H drops to 171). The least-used sidecar
is the explicit polarity, followed by the orientation origin.

## 1. A declared field chart, not a canonical arithmetic frame

Use V=F2^3, with nonzero vectors represented by bit integers 1,...,7. The
original Fano lines are the triples {x,y,x+y}. Identify V with
F8=F2[alpha]/(alpha^3+alpha+1), in the polynomial basis (1,alpha,alpha^2).
Then

    (alpha^i)_(i=0,...,6) = (1,2,4,3,6,7,5).

Write x=alpha^i, y=alpha^j, and orient x -> y iff
j-i mod7 belongs to {1,2,4}. Call this T0. The lines of F become
`dev{0,1,3}` in exponent coordinates. Its other seven cyclic triples are

    G0 = {{x^-1,y^-1,z^-1} : {x,y,z} in F},

which become `dev{0,1,5}`. Inversion is exponent negation, so reverses every
arrow and swaps F with G0. Frobenius is exponent doubling, hence preserves
both the XOR incidence and T0. These assertions can be verified from the
seven displayed field elements; the inherited difference-family proof
establishes the same two-plane decomposition.

The 14 cycles constitute a 2-(7,3,2) design: every unordered pair lies in
exactly two cyclic triples, one from each plane. There are 21 arcs and 21
transitive triples; neither 21 is the number of cyclic triples.

## 2. Complete partner and orientation fibres

The checker independently generates every STS(7) by exact cover of all 21
unordered pairs with triples, and as S7-images of F. The resulting sets agree
and have size 30. With F fixed, the complete intersection histogram is

| Number of shared blocks | 0 | 1 | 3 | 7 |
|---|---:|---:|---:|---:|
| Number of STS(7) | 8 | 14 | 7 | 1 |

The 5040 relabellings of T0 give 240 distinct labelled Paley tournaments.
Exactly 16 have every block of F cyclic. For each, the seven remaining cyclic
triples form one of the eight disjoint partners. Every partner occurs twice,
with the two tournaments related by reversing all arrows.

This is not just a search within Paley. Each pair of points lies in one Fano
line, so any tournament in which every F-line is cyclic is specified by
seven binary line orientations: exactly 128 possibilities. Exhausting these
shows that the disjoint-partner completions are precisely the above 16.
Each of the eight full cyclic-triple fibres has exactly two members among
**all 128 orientations**, not only among the known Paley relabellings.
The checker also tests all **120 unordered disjoint STS pairs**, exhausting
all 128 line orientations for each first plane, and finds exactly the two
opposite Paley completions. This is a complete finite proof of the stated
quantifiers, with no assumption that a sought completion is Paley.

For each Paley tournament, precisely two of the thirty STS lie in its cyclic
triple set. Consequently there are 120 distinct unions F union G, each with
two orientations. Globally, choosing a marked Fano component and a Paley
orientation gives 30*16=480 states; forgetting the mark leaves 240 tournaments.

## 3. The four-bit normal form: three partner bits plus chirality

Choose the standard dot-product polarity and label each Fano line by its
nonzero normal:

    L_a = {x != 0 : a dot x = 0},  a != 0.

Point labels x and line labels a belong to **different coordinate spaces**.
The dot product explicitly supplies their pairing. Starting at T0, let T_c
reverse all three arcs of line L_a whenever c(a)=1. Lines partition the
unordered pairs, so this defines all 128 Fano-line orientations bijectively.

The compatible Paley locus is exactly the relative Hamming code

    C = {c : XOR_(a with c(a)=1) a = 0}
      = {c_(S,v)(a)=S+a dot v : S in F2, v in F2^3}.             (2)

The first equality with the Paley locus is the inherited finite-exact
classification, replayed independently over the 128 orientations here.
The second is elementary: the displayed affine evaluations satisfy the
three parity checks; they form a four-dimensional space, equal to the
kernel since the seven nonzero check columns have rank three.

Since all seven line flips reverse all 21 arcs,

    T_(S+1,v) = opposite(T_(S,v)).                              (3)

The partner is independent of S, and the eight v give eight distinct partners
G_v=CyclicTriples(T_(0,v)) minus F. Thus the partner space is precisely
C/<all-ones>, a three-dimensional quotient **after choosing the origin and
polarity**. Equation (2) provides a literal 3+1-bit reconstruction algorithm.

The group of Fano automorphisms is GL(3,2), of order 7*6*4=168: a line
permutation preserves x+y as the third point and extends by fixing zero.
The checker finds Aut(T0) of size 21, contained in GL(3,2), and verifies

    GL(3,2).T0 = {T_(0,v) : v in F2^3},
    opposite orbit = {T_(1,v) : v in F2^3}.                    (4)

The first relative code words have weights 0 or 4; the second have weights
3 or 7. Therefore S=sum_a c(a) is a chosen chirality bit. Both group orbits
have size 168/21=8. The group acts transitively on the eight partners, also
with stabilizer 21. In particular incidence alone admits no **automorphism-
equivariant** choice of partner or compatible orientation: there is no fixed
point. An arbitrary numeric ordering can of course break symmetry; it is an
extra convention, not an intrinsic reconstruction.

## 4. Exact transport from squarefree arithmetic

Retain the inherited bank, where the parity coordinates of sigma(p)=p+1
are the valuations at (2,3,5):

| nonzero column a | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| prime p_a | 7 | 2 | 5 | 19 | 89 | 59 | 29 |
| sigma(p_a) | 8 | 3 | 6 | 20 | 90 | 60 | 30 |

For a subset mask c, let N_c=product_(c(a)=1) p_a. Distinct primes mean every
subset is legal, every product is squarefree, and all divisor sums involve
only the three indicated parity coordinates. Consequently

    sigma(N_c) square iff XOR_(c(a)=1) a=0 iff c in C.

Now use **the same a as a dual line normal**, not as a primal tournament
vertex, and send N_c -> T_c. This defines a bijection between the 128 bank
products and the 128 cyclic orientations of the fixed Fano lines. Combining
with (2) proves (1). The checker independently sums the complete divisor
list of each of the 128 products, compares the multiplicative formula,
tests squareness, and verifies the tournament classification.

The exact source operation is squarefree product

    N star M = N*M/gcd(N,M)^2.

It corresponds to symmetric difference of selected primes. On tournaments
the corresponding operation is XOR of their **relative line-flip masks**:

    T_c odot_(T0) T_d = T_(c+d),
    Phi(N star M)=Phi(N) odot_(T0) Phi(M).                     (5)

The target starts as a torsor; choosing T0 supplies its identity. This is
not a claim that tournaments have a natural multiplication, nor that ordinary
integer multiplication is preserved. In particular repeated prime factors
are removed by star. All 128^2=16384 ordered operation pairs are checked.

The seven weight-three and seven weight-four nonzero codewords are exactly
the **14 minimal nonempty square-sigma selections**, and become the minimal
sets of line-cycle reversals taking T0 to a Paley tournament. The full
seven-prime selection is the fifteenth nonzero codeword, but is not minimal;
it reverses the entire tournament. The weight-four patterns are complements
of Fano lines in the **normal-label geometry**, not single primal lines.

On the square-sigma locus, squarefreeness gives

    mu(N_c)=lambda(N_c)=(-1)^wt(c)=(-1)^S.                     (6)

By (4), this sign is exactly the chosen GL(3,2)-chirality. The matching is
explicit, origin-dependent, and loses its distinction under arbitrary S7
relabelling. It does not assert that mu or lambda is an invariant of an
unmarked tournament isomorphism class.

### Chirality can be decoded without recovering the partner

For any normal-label Fano line {a,b,a+b}, equation (2) gives

    c(a)+c(b)+c(a+b)=S.

Thus those three prime-selection bits, equivalently their matching line-flip
bits, recover chirality even though they need not recover the full four-bit
state. More generally an observed coordinate set recovers S iff it contains
a Fano line. Indeed a linear recovery combination for S is an odd subset
with XOR zero. In this length-seven code its support has size three or seven;
in either case the observed set contains a line. Conversely the displayed
triple formula suffices. If no such combination exists, linear-algebra
separation gives two codewords with identical observations and different S.

The complete recovery census, replaying Session85's parity-access result,
is 7 recoverable three-sets, 28 four-sets, all 21 five-sets, all 7 six-sets,
and the full seven-set; no set of size at most two works. Not every triple
or four-set works. This uses lines among the **coordinate normals**; it does
not claim three arbitrary tournament arcs determine chirality.

## 5. Counts and boundaries that must remain separate

* H(T0)=189 is the number of directed Hamiltonian paths, checked by subset
  dynamic programming and all 5040 orders. |Aut(T0)|=21 counts vertex
  permutations preserving its arrows. Their ratio 9 counts free Aut-orbits
  on Hamiltonian paths: fixing an ordered Hamiltonian path fixes every vertex.
  This does not identify H with 7 or 21. The inherited all-order forbidden-H
  theorem is not needed or reproved here.
* The bank's Frobenius permutation maps 70=2*5*7 to 11837=19*89*7, and
  sigma 144 to 14400. It preserves the **Boolean square-sigma predicate**,
  not the numeric divisor sum. No idoneality or aliquot-dynamics transfer
  follows. This is the inherited counterexample, not a new arithmetic claim.
  Under the chosen dot-product polarity, primal Frobenius permutes labels
  `(1,4,5,6,7,2,3)`, but its induced action on line normals is the inverse-
  transpose permutation `(1,6,7,2,3,4,5)`. Thus the arithmetic-to-line map is
  **not asserted equivariant for the same primal Frobenius permutation**;
  these are distinct predicate-preserving transports.
* The field multiplication, polarity and baseline are actual sidecars. A
  primal arithmetic Hamming code by itself contains none of these choices.
* Regularity, seven original cyclic Fano lines, or the numbers 7 and 21 do
  not characterize Paley: a single line flip is an explicit hostile.
* Scope is exactly seven points and the declared seven-prime bank. No
  general disjoint-STS completion theorem or universal arithmetic encoding
  is inferred from this census.

## Reproduction and completeness

    python3 check_paley.py > normal.json
    python3 -O check_paley.py > optimized.json
    cmp normal.json optimized.json

Standard library only; no `assert` is used, so optimized Python does not erase
the checks. Both runs pass with byte-identical deterministic JSON. The script
checks all 35 triples, all 5040 permutations, all 30 STS by two constructions,
all 128 fixed-plane orientations, all 120 disjoint-plane pairs with every
line orientation, all 128 arithmetic sources, all 16384 operation pairs,
and all 128 possible observed coordinate sets. Hamilton counts for all 128
line orientations are compared by two independent algorithms. These are
complete finite universes, not random samples. All verification checks remain
active under optimized Python.
