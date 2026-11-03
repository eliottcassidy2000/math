# Session 143: five-vertex tournament spectra and the repaired 5:4 bridge

**Status.** PROVED: Burnside counts, odd automorphism obstruction, unique regular
class, exact regular spectrum and cyclotomic algebra, four-object affine action.
FINITE-EXACT: full twelve-class table and completeness of the stated sidecars.
CITED: occurrence of `5:4 × HS:2` in the ATLAS Baby Monster subgroup list.
NOT ESTABLISHED: an explicit Baby Monster generator/intertwiner construction,
any special octonion or magic-square correspondence, or novelty.

## 1. The counting proposal survives

There are exactly **12** isomorphism classes of tournaments on five vertices;
**8** are self-converse (self-dual), leaving **two** pairs of distinct dual
classes. Hence there are **10** classes after identifying converse pairs.
The four labeled circulant orientations below must not be confused with these
12 unlabeled classes or with two distinct regular classes.

This agrees with Brendan McKay's [primary tournament catalogue](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html),
which lists 12 tournaments, 8 self-converse tournaments, and 1 regular tournament
at order five. The accompanying `catalogue_crosscheck.py` also independently
canonicalizes the complete numeric records from all three source catalogues
and matches our representatives exactly. These are established counts; no
priority claim is made.

### Independent Burnside proof

An involution cannot preserve a tournament: any transposed pair would have its
unique directed edge reversed. Every finite tournament automorphism group
therefore has odd order, by Cauchy's theorem. More directly, an automorphism
cannot have an even vertex-cycle, because half a turn on that cycle reverses
an antipodal pair. In `S5` only the following types can fix a tournament:

| permutation type | number | fixed labeled tournaments |
|---|---:|---:|
| `1^5` | 1 | 1024 |
| `3,1,1` | 20 | 16 |
| `5` | 24 | 4 |

For a 3-cycle the ten unordered edges have four orbits, permitting four
independent orientation bits; a 5-cycle has two such orbits. Thus

`(1024 + 20*16 + 24*4)/120 = 12`.

For the converse-twisted action, an anti-automorphism can only have cycle type
`2,2,1`. A permutation of odd order cannot be anti; type `2,1,1,1` fixes three
vertices; type `4,1` squares to an impossible involutory automorphism; type
`3,2` cubes to an anti-automorphism fixing three vertices. There are 15
permutations of type `2,2,1`. Each fixes 64 tournaments under relabeling followed
by converse: its two swapped-pair edges can be oriented freely and its other
eight edges form four freely orientable pairs. Consequently

`(1440 + 15*64)/240 = 10`.

If `s` of the 12 isomorphism classes are self-dual, the converse quotient has
`(12+s)/2` classes, so `s=8`.

## 2. Exact classification and the spectral information loss

The encoding fixes vertices `0,1,2,3,4`. The ten bits of a mask correspond, least
significant first, to the pairs `01,02,03,04,12,13,14,23,24,34`; bit one means
`i -> j`. Each representative is the minimum mask over all 120 relabelings.
`c_k` denotes the number of directed simple cycles of length `k`, counted once
up to cyclic rotation. SCC words are ordered from source to sink.

| mask | sorted outdegrees | converse mask | Aut order | anti count | `(c3,c4,c5)` | ordered SCC sizes |
|---:|---|---:|---:|---:|---|---|
| 0 | 0,1,2,3,4 | 0 | 1 | 1 | 0,0,0 | 1,1,1,1,1 |
| 2 | 1,1,1,3,4 | 11 | 3 | 0 | 1,0,0 | 1,1,3 |
| 4 | 1,1,2,2,4 | 9 | 1 | 0 | 2,1,0 | 1,4 |
| 5 | 0,2,2,2,4 | 5 | 3 | 3 | 1,0,0 | 1,3,1 |
| 8 | 1,1,2,3,3 | 8 | 1 | 1 | 3,3,1 | 5 |
| 9 | 0,2,2,3,3 | 4 | 1 | 0 | 2,1,0 | 4,1 |
| 10 | 1,1,2,3,3 | 10 | 1 | 1 | 3,2,1 | 5 |
| 11 | 0,1,3,3,3 | 2 | 3 | 0 | 1,0,0 | 3,1,1 |
| 12 | 1,2,2,2,3 | 12 | 1 | 1 | 4,4,1 | 5 |
| 40 | 1,2,2,2,3 | 40 | 3 | 3 | 4,3,3 | 5 |
| 41 | 1,2,2,2,3 | 41 | 1 | 1 | 4,3,2 | 5 |
| 76 | 2,2,2,2,2 | 76 | 5 | 5 | 5,5,2 | 5 |

All adjacency characteristic polynomials are given exactly by

`chi_A(x) = x^5 - c3*x^2 - c4*x - c5`.

**Proof of the formula.** In the determinant expansion, a nonidentity
contributing permutation must use a directed cycle of length at least three.
Two disjoint such cycles cannot fit on five vertices. One cycle of length
`k` contributes `-x^(5-k)`, proving the formula.

Thus there are **nine**, not two, adjacency spectra. Precisely the following
isomorphism classes collide spectrally:

- `2,5,11`: `x^5-x^2 = x^2(x-1)(x^2+x+1)`;
- `4,9`: `x^5-2x^2-x = x(x^4-2x-1)`.

The first collision still leaves two inequivalent objects after quotienting by
converse: an end-position triangle and a middle-position triangle. This is a
specific obstruction to a spectrum-only classification of the ten classes.

**Mechanism, not just a counterexample.** A tournament's SCC condensation is a
linear order, and its adjacency matrix is block triangular in that order.
Its characteristic polynomial is the product of SCC characteristic polynomials,
which forgets the order of the blocks. The triangle's position is therefore
invisible to the spectrum.

**Repaired complete invariants, order-five scope only (FINITE-EXACT).**

1. Adjacency spectrum plus sorted score sequence separates all twelve classes.
2. Adjacency spectrum plus the ordered SCC-size word separates all twelve.
3. Adjacency spectrum plus that word up to reversal separates all ten converse
   classes.

No extension of these completeness claims to larger orders is asserted.
The checker's JSON includes every adjacency matrix, automorphism and
anti-automorphism. All automorphism groups here are cyclic of order 1, 3 or 5.
For each self-dual class, the extended group generated by automorphisms and
anti-automorphisms is dihedral of order twice the automorphism-group order
(`D2=C2`, `D6=S3`, or `D10`); each anti-automorphism is an involution.

## 3. The one regular tournament and its exact spectrum

Regularity forces every outdegree to be two. Fix a vertex `0`, order its
outneighbors `a->b`, and order its inneighbors `c->d`. Vertex `b` must beat
both `c,d`; vertex `c` already beats `0,d`, so it must lose to `a,b`; vertex
`d` must then beat `a`. All arcs are now forced and give the carousel order
`0,a,b,c,d`. This proves uniqueness without exhaustive computation.

On `F5`, put `i->j` iff `j-i` belongs to `{1,2}`. Let
`P_ij=1` iff `j-i=1`, so `P^5=I` and `A=P+P^2`. If `zeta` is a primitive
fifth root of unity, Fourier eigenvectors give

`Spec(A) = {2} union {zeta^k+zeta^(2k): k=1,2,3,4}`.

Equivalently,

`Spec(A) = {2, -1/2 +/- (i/2)sqrt(5+2sqrt(5)),`
`                 -1/2 +/- (i/2)sqrt(5-2sqrt(5))}`.

Here each `+/-` supplies one complex-conjugate pair, for five values total.
The characteristic polynomial is

`(x-2)(x^4+2x^3+4x^2+3x+1) = x^5-5x^2-5x-2`.

For `q(x)=x^4+2x^3+4x^2+3x+1`, the exact identity

`16*q((y-1)/2) = y^4+10y^2+5`

gives the radical expression. The two imaginary magnitudes are within the
**same** adjacency spectrum, not the spectra of two distinct regular classes.
For skew adjacency `K=A-A^T`, the spectrum is instead
`{0, +/- i sqrt(5+2sqrt(5)), +/- i sqrt(5-2sqrt(5))}`.
The distinction between adjacency and skew adjacency is essential.

Every tournament adjacency matrix has `tr(A)=tr(A^2)=0`. In particular,
any proposed eigenvalue multiset with nonzero total sum is impossible even
before this classification is applied.

Fixing a vertex in the regular tournament fixes its two ordered in- and
out-neighborhoods pointwise; hence an automorphism is determined by one vertex
image. The five translations exist, so `Aut(T)=C5`. The reflections
`x->b-x` reverse all edges, giving five anti-automorphisms and extended group
`C5:C2=D10` of order ten. A group `5:4` has even order and cannot preserve a
tournament faithfully as its automorphism group.

## 4. Where 5:4 really acts: a family, and the cyclotomic algebra

Fix the translation group on the five labeled vertices. A translation-invariant
tournament chooses one residue from each of `{1,4}` and `{2,3}`, so there are
exactly four connection sets. Multiplication by two cycles them as follows:

`{1,2} -> {2,4} -> {3,4} -> {1,3} -> {1,2}`.

Their adjacency matrices are
`P+P^2`, `P^2+P^4`, `P^3+P^4`, `P+P^3`.
All four are isomorphic and cospectral. Two steps are negation/converse.

The affine group `x->a*x+b`, `a in F5^*, b in F5`, is `AGL(1,5)=5:4`,
order twenty. It acts transitively on this **four-object family**. Its kernel
on the family is the translation `C5`; its quotient `C4` rotates the family.
The full setwise stabilizer in `S5` of the family is exactly this affine group:
a family stabilizer conjugates the common automorphism group `C5` to itself,
so is in its normalizer, and the normalizer is affine. For example, if a
permutation normalizes translations, it satisfies `f(x+1)=f(x)+a`, hence
`f(x)=a*x+b`. The subgroup preserving one tournament up to converse has
slopes `+/-1`, and is `D10`. Quotienting the four-object family by converse
leaves two states, but still only **one** unlabeled class and **one** spectrum.

There is also an exact algebraic map, not a numerical resemblance:

`Q[A] = Q[P] ~= Q[x]/(x^5-1) ~= Q direct-sum Q(zeta5)`.

The constructive inverse is

`11P = 5A^4 - A^3 - 2A^2 - 18A - 17I`.

Also `q(A)=11J`, where `J` is all ones. Thus the rational four-dimensional
augmentation space `sum(v_i)=0` carries the cyclotomic factor. Conjugation by
the multiplier-two permutation sends `P` to `P^2`, and on the cyclotomic
factor induces the genuine Galois automorphism `zeta5 -> zeta5^2`, of order
four. Translations act trivially by conjugation on this commutative algebra.
The quartic `q` is irreducible because the displayed inverse recovers `zeta5`
from `zeta5+zeta5^2`, which therefore generates the degree-four cyclotomic field.

This action is a faithful `5:4` permutation action on the five vertices but
only a `C4` action on the four-object family or algebra by conjugation. Keeping
these two kernels separate prevents another false identification.

## 5. Baby Monster: what survives, and the next proof obligation

The [ATLAS Baby Monster entry](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/)
lists a maximal subgroup written `5:4 × HS:2`, of order 1,774,080,000.
The [current ATLAS subgroup page](https://brauer.maths.qmul.ac.uk/Atlas/v3/subgroup/BG1-max15W1)
provides a straight-line subgroup generator program. **CITED**, not an
independent construction in this session. This supplies a real subgroup-level
reason to investigate the user's `5:4` cue, but not the proposed spectrum
identification.

The precise repaired conditional interface is: given elements `c,d` in an
identified Baby Monster subgroup with

`ord(c)=5, ord(d)=4, d*c*d^-1=c^2, |<c,d>|=20`,

send `c` to translation `x->x+1` and `d` to multiplier `x->2x`. This gives an
explicit isomorphism to the affine model; its induced action on the four
orientations fixes them under `c` and rotates them under `d`. On the
cyclotomic algebra, `d` acts by `zeta5->zeta5^2`. No tournament eigenvalues are
inferred from a Baby Monster character merely because such a subgroup exists.

**Still needed before a Baby Monster representation theorem:** verified
concrete `c,d` (or explicit words and field/matrices), a specified module,
an intertwiner with the five-vertex permutation/augmentation model, and the
operator whose spectrum is being compared. None is manufactured here. The
available ATLAS program has not been executed, and its page records a file-
existence check rather than a verification by this session. The two order-five
conjugacy classes in the ATLAS table are not two spectra of regular tournaments.

## 6. Inheritance and connection contract

Pinned repository revision: `8286d710498a3b4e9b9f1ac683c6eb77716f3f4b`.

- Closest proved mechanism: `THM-3447`,
  `01-canon/theorems/THM-3447-regular-group-orbit-tournament-parity-and-involution-ties.md`:
  inverse-pair half-sets on odd regular group orbits. Our `C5` four-orientation
  family is its smallest nontrivial two-choice example.
- Hostile / corrected near miss: `THM-1858`,
  `01-canon/theorems/THM-1858-no-tournament-has-exactly-two-distinct-eigenvalues.md`:
  trace-zero and square-trace-zero exclude impossible adjacency lists. It is
  not a theorem that only two order-five tournament spectra occur.
- Sidecar mechanism: `THM-2532`,
  `01-canon/theorems/THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient.md`:
  cyclic shifts and oriented circulant algebras retain a direction that a
  converse-even quotient forgets. Here the explicit polynomial inverse
  reconstructs the chosen `P`, while an unordered spectrum does not.

**Connection contract.** Source: a chosen regular `C5` action plus an inverse
half-set. Target: a four-member family of oriented relations and its rational
circulant algebra. Map: `S -> sum_(s in S)P^s`, with affine conjugation.
Preserved: relabeling/converse action, exact spectrum, and the full chosen
operator algebra. Lost under spectrum alone: orientation, labels, and SCC
position in the general census. Sidecars: a chosen `P`/connection set or the
ordered SCC word as appropriate. Decisive tests: multiplier-two moves one
tournament instead of preserving it; middle-triangle mask5 is cospectral but
not converse-isomorphic to end-triangle mask2.

No nonassociative multiplication or squarehood predicate has been transported;
connections to octonions, magic-square-of-squares, or the named arithmetic
cycle require additional explicit maps from the other lanes.

## 7. Reproduction and hostile controls

The sole computational dependency is Python's standard library. From this
folder:

```sh
python3 check_tournaments.py > /tmp/tournaments.normal.json
python3 -O check_tournaments.py > /tmp/tournaments.optimized.json
cmp /tmp/tournaments.normal.json classification.json
cmp /tmp/tournaments.optimized.json classification.json
python3 verify_manifest.py
```

The checker performs **9,268 explicit checks**, including all 1,024 labeled
masks, all 120 vertex relabelings at each canonical representative, orbit-
stabilizer and two implementations of permutation transport, exact Newton
characteristics versus a Leibniz determinant and simple-cycle counts, both
Burnside counts by permutation type, all exact group/coset claims, the affine
family action, the rational inverse, and the radical polynomial identity.
An additional independent matrix-relabeling implementation exactly matches all
12/8/1 primary McKay catalogue records to our classes. Validation uses
exceptions, not `assert`, so `-O` does not disable checks.
The transitive nilpotent matrix, non-self-dual dual pairs, and the middle/end
triangle cospectral collision are built-in hostile controls. No floating-point
root approximation is used. Repository and external source pins are in
`source_pins.json`; computational deliverables are hashed in `manifest.json`.
