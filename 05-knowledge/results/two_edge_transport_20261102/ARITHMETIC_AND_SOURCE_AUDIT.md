# Arithmetic inspirations with their preserved data

**PROVED elementary controls / CITED literature / explicitly scoped analogies.**
The supplied text is preserved as pasted_inspiration.txt. Its suggestions
are treated as mathematical claims to test, not as instructions to assume.

## Audit of the modular and tournament claims

The five-vertex enumeration is correct:12 ordinary tournament classes,
8 self-converse classes, and2 pairs of mutually converse classes, giving
10 merged classes. This is independently reconstructed from all1024
labelled orientations. It is an application of an involution orbit count:
(number of objects + number of fixed objects)/2.

The Hecke table is not a defined calculation without a weight, character,
eigenform and normalization. Under the natural weight-two cuspidal
interpretation it is impossible: the space has dimension zero.
[Stein, Modular Forms, Proposition6.1](https://wstein.org/books/modform/stein-modform.pdf)
identifies that dimension with the genus. At p=13, the index is14,
there are two cusps and two elliptic points of each relevant type, giving
\[
 g(X_0(13))=1+14/12-2/4-2/3-2/2=0.
\]
The exact rational calculation is replayed. Atkin–Lehner W13 and a
Hecke eigenvalue a13 are also different specified operators; one cannot
substitute an unspecified W13 eigenvalue into the table. Rational Hecke
coefficients describe a coefficient field, not a collapse of all geometry
onto rational points.

Ogg's genus-zero condition concerns the quotient X_0(p)^+, not X_0(p)
itself. The [primary Duncan–Ono paper](https://uva.theopenscholar.com/files/ken-ono/files/1411.5354.pdf)
records the fifteen supersingular/Monster primes and connects appropriate
moonshine functions with supersingular j-invariants. This supplies no map
from the ten merged five-vertex tournament classes to ten “independent
paths” on X_0(13). No such path space or intertwining operator is defined
in the pasted text.

The unit-ball volume formula is
V_n=pi^(n/2)/Gamma(n/2+1), with V_(n+2)=2pi V_n/(n+2); its integer-dimensional
maximum is indeed at5. This numerical coincidence does not prove the
claimed common decay mechanism. Self-converse proportions are not
globally maximal at5: at3 the proportion is1, at4 it is1/2, and at5 it
is2/3. The claimed MRB boundary and “gravitational” prime-bias mechanism
have no specified invariant or map. They remain unsupported analogies.

### An elementary bound identifies the actual symmetry cost

**PROVED.** For n>=2, a self-converse tournament admits an anti-automorphism
that is an involution with k=floor(n/2) transpositions and at most one
fixed vertex. To see this, tournament automorphism groups have odd order:
an involution would exchange two vertices and reverse their edge.
If sigma is any anti-automorphism, its order is2m with m odd, since
sigma^2 is an automorphism. Then sigma^m is an anti-involution. Two fixed
vertices would again force a fixed edge to reverse.

There are n!/(2^k k!) such permutations. For a fixed one, the compatible
orientations have floor(n^2/4) independent pair-orbit bits, hence exactly
2^floor(n^2/4) possibilities. A union bound gives, for the labelled
self-converse fraction,
\[
 \Pr(T\text{ self-converse})\le
 \frac{n!}{2^k k!}\,
 2^{\lfloor n^2/4\rfloor-\binom n2}.
\]
If S_n,A_n are self-converse and all unlabelled class counts, then
A_n>=2^binom(n,2)/n! and S_n is at most the labelled count. Thus
\[
 \frac{S_n}{A_n}\le
 \frac{(n!)^2}{2^k k!}\,
 2^{\lfloor n^2/4\rfloor-\binom n2}
 =2^{-n^2/4+O(n\log n)}.
\]
These are upper bounds, not asymptotic equivalences. The quadratic number
of constrained edge bits drives this bound; factorial terms are lower
order in its logarithm. This differs from the unit-ball volume's
n log n scale. It also specifies what “fraction” means, which the pasted
analogy did not.

## The requested nine-cycle supplies a real symmetry test

The exact cycle
\[
 4,8,3,7,2,6,10,5,9,4
\]
is the map x -> 2+((x-2+4) mod9) on the representatives2,...,10.
It is one cycle because gcd(4,9)=1. Its return clock is9.

There is a literal tournament construction from this cycle. Order its
nine vertices by the displayed phase, and orient phase i toward phase j
when (j-i) mod9 belongs to{1,2,3,4}. Oddness removes ties. Rotation is an
automorphism, so C9 acts on the resulting tournament; its path count is
3267=27*11^2 in the exact computation. Only divisibility by9 is explained
by the rotation argument; the remaining factorization is not assigned a
new arithmetic meaning.

This transfers the selector lemma: a rotation-equivariant deterministic
map from this source to a nine-vertex H17 tournament is impossible,
because9 does not divide17. A rotation-invariant positive law supported
on actual H17 targets needs at least nine atoms. Such targets exist by
padding a six-vertex H17 example with three ordered singletons, and their
C9 orbit has exactly nine elements. A chosen phase is the needed symmetry
breaking input. The map preserves cycle rotation; it does not preserve
square-sum adjacency or magic-square equations.

For the latter distinction, start with shifted Lo Shu:
\[
 \begin{pmatrix}9&2&7\\4&6&8\\5&10&3\end{pmatrix}.
\]
Its eight magic lines sum to18. Applying the nine-cycle to its values
produces rows of sums12,21,21. Total sum and nine distinct values survive,
while the magic equations do not.

## Recovered magic-square work

The [Session55 p-adic lift proof](../session55_parker_fields_integer_lifts_20261025/padic-lifts/PROOF.md)
already gives Lo Shu transforms3L-2 and8L-7 with compatible residue-square
lifts through all sufficiently large powers of3 or2. Their fixed integer
entries need not be integer squares. Local square-root towers therefore
do not supply an integer magic square.

The [Session57 decorated biclique equivalence](../session57_square_sum_magic_paths_20261025/biclique-bridge/THEOREM.md)
retains a center C and offsets0,+/-x,+/-y,+/-(x+y),+/-(x-y). Its reflected,
parallelogram-decorated K4,4 has precisely those nine edge labels with
multiplicities. All labels being squares and all nine being distinct
are essential constraints. A bare graph, or a nine-cycle permutation
of entries, forgets them.

The transfer to this session is a concrete verification rule: preserve
shared numerical labels across paths just as the birth calculation
preserves shared orientation bits. Independent copies of a cell or an
edge can satisfy separate conditions without giving one joint object.
No integer magic square of nine distinct squares is constructed here.

## Lambda cube and octonions

I interpret “lambda cube” as Barendregt's cube of typed calculi.
The [primary paper](https://homepages.inf.ed.ac.uk/wadler/papers/barendregt/pure-type-systems.pdf)
organizes eight systems by which dependencies are permitted. This is an
inclusion poset, not an eight-dimensional multiplication algebra.
A bounded targeted search did not locate an exact prior repository
reference under lambda-cube/Barendregt; nearby branch-cube work is a
different mathematical object.

Its useful application is to specify an inverse's output type. A
tournament inverse must return n, a Boolean orientation T, and a
certificate H(T)=h. Returning a real probability vector with mean h has
a different type. Shared sources and chosen ports are dependencies that
must appear in the constructor's input and correctness statement.
This is a specification analogy, not a formal Lean construction or a
theorem identifying the lambda cube with the tournament carrier.

For octonions there is a genuine forgetful map. Basis indices lie in
F_2^3, with e_a e_b=epsilon(a,b)e_(a+b). Forgetting the sign gives XOR,
but loses the associator. The ratio of the two signed products
epsilon(a,b)epsilon(a+b,c) and
epsilon(b,c)epsilon(a,b+c) records that loss.
[Baez's primary account](https://math.ucr.edu/home/baez/octonions/node4.html)
describes this twisted group-algebra construction.

The recovered [Paley/Fano/octonion audit](../collatz_mod6_20260922_paley_fano_octonion_design.md),
S11–S18, is sharper than a seven/eight coincidence: among128 orientations
of one Fano plane,16 give the alternative composition algebra and112 do
not. It retains the actual line multiplication and exact sign convention.
The present shared lesson is to retain the operation defect, not to
identify that octonion associator with tournament re-marking. No such
identification has been proved.

All these comparisons point to a concrete rule: test a proposed quotient
against the operations it must support. They do not establish a common
Hecke spectrum, a Monster representation of the metagraph, or a proof of
the outstanding arithmetic conjectures.

