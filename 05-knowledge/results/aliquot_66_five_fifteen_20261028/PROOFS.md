# Exact bridges and their failure boundaries

Status: **PROVED elementary constructions and deductions from named classical
theorems; FINITE-EXACT complete small classifications.** This is a research
result packet, not a newly audited canon theorem. Existing Session84 results
are explicitly credited in [INHERITANCE.md](INHERITANCE.md).

## P1. Divisor sums count lattices; content separates the cyclic part

Let L(N) be the actual sublattices Λ of Z² of index N, with the ambient
basis retained. Unique column Hermite normal form writes each as

    Λ = columns of [[a,b],[0,c]],     ac=N, 0<=b<a.

There are a choices for b for each divisor a, so |L(N)|=σ(N). This counts
embedded lattices, not isomorphism types of their abstract quotient groups.
The first Smith invariant is d=gcd(a,b,c), the largest integer such that
Λ⊂dZ². Call it content. Then Λ=dΛ₀, and Λ₀ has content1 and index N/d².
Conversely each such primitive Λ₀ gives exactly one Λ. Consequently

    L(N) = disjoint union_{d²|N} d L_primitive(N/d²).       (1)

Content1 is equivalent to Z²/Λ being cyclic: its Smith invariants are
`(d,N/d)`. A cyclic quotient of order N is the kernel of a surjection
Z²→Z/N; two such surjections have the same kernel exactly when they differ
by a unit of Z/N. The primitive row vectors modulo N number
`N² product_{p|N}(1−1/p²)`, and each unit orbit has φ(N) elements. Thus

    |L_primitive(N)| = ψ(N) = N product_{p|N}(1+1/p).

For N=1 both sets have one element. Equivalently these kernels are indexed
by P¹(Z/N), with the primitive-vector/unit convention. Formula(1) proves

    σ(N) = sum_{d²|N} ψ(N/d²).                            (2)

This is classical arithmetic, proved here to make the lost coordinate
explicit. Multiplicativity also follows from primary decomposition of the
finite quotient, with actual prime labels and coprimality retained. For
squarefree N, only d=1 occurs. For three distinct primes p,q,r,

    σ(p²qr) = ψ(p²qr) + ψ(qr)
             = p(p+1)(q+1)(r+1) + (q+1)(r+1).            (3)

The last term is the complete noncyclic layer; it is not a small numerical
error. This is an exact extension of the incoming divisor-box observer.
The box counts ordinary divisors, while L(N) counts weighted choices over
those divisors. They are different carriers connected by the HNF projection
to a, whose fibre has size a.

**Where F=U+S enters.** In the repository count convention,
`F=τ(N)−2`, `U=ω(N)−[N prime]`, and
`S=2^ω(N)−1−[N squarefree]`. The inherited complete balanced shapes are
primes, prime cubes, and p²qr. Thus the mixed balanced shape has two content
layers in(3). Having two content layers is not a characterization of
balance: p², p²q, and p³ also have two layers. The atom/Boolean-subbox
counts remain necessary. No arithmetic conclusion is drawn just from the
Klein-action fold of THM-4733.

## P2. The sigma144 fibre has modular genera9 and11

The modular-curve index is μ₀(N)=ψ(N). The classical genus formula is

    g₀(N)=1+ψ(N)/12−e₂(N)/4−e₃(N)/3−c(N)/2,
    c(N)=sum_{d|N} φ(gcd(d,N/d)).                         (4)

Here e₂,e₃ count the elliptic points. These are the formulas of
[Stein, Dimension Formulas, Proposition6.1](https://wstein.org/books/modform/modform/dimension_formulas.html).
They are cited, not derived by a new modular-curve argument here.

The complete fibre, inherited from Session83 and checked by direct sums
over all1<=N<144, is

    σ^{-1}(144)={66,70,94,115,119}.

The cutoff is proved: for N>1, σ(N)>=N+1. Every source is squarefree, so
ψ=144. At66,70,94 the factor2 makes e₃=0; at115,119 the factors5,17
respectively do so. Each level also contains an odd prime3 modulo4,
making e₂=0. A squarefree level has c=2^ω(N). Hence

    g₀(66)=g₀(70)=9,
    g₀(94)=g₀(115)=g₀(119)=11.                           (5)

For N=66 or70 and any marked prime p|N, ψ(pN)=144p and
c(pN)=4(p+1). The zero elliptic-point counts persist. Therefore

    g₀(pN)=10p−1.                                       (6)

This gives genera19,29,109 at132,198,726, and19,49,69 at140,350,490.
At each of these sources the divisor count is balanced, while σ ceases to
be square. Their common modular calculation does not prove the square
obstruction; that is the independently inherited repunit-ratio theorem.

Use **X₀**, not its Fricke or full Atkin–Lehner quotient, in(4)–(6).
The supersingular/Monster lane uses specified quotients and a prime not
dividing the base level. A genus9 source curve and a genus0 quotient are
compatible statements, not conflicting calculations. In particular the
existing level lists place both66 and70 among the full Atkin–Lehner
genus-zero levels; this additional fact is inherited from that lane,
not established by(4).

## P3. Identical divisor data do not determine class inversion

Use proper classes of primitive positive binary quadratic forms of negative
discriminant D, equivalently the ideal class group of the quadratic order
of discriminant D. Classical reduction gives one representative per class
with |b|<=a<=c, taking b>=0 when |b|=a or a=c. Its inverse is represented
by (a,−b,c). Interior forms with0<|b|<a<c have distinct inverse classes;
the boundary forms b=0, |b|=a or a=c are precisely the self-inverse ones.
The complete search bound is a<=sqrt(|D|/3). These are standard reduction
facts; the [Magma class-group documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/374)
records the order/form correspondence and reduced-form enumeration route.
No Magma computation or GRH-dependent class-group algorithm is used here.

At D=−264 the complete reduced list is

    (1,0,66), (2,0,33), (3,0,22), (6,0,11),
    (5,±4,14), (7,±4,10).

There are eight classes, of which exactly four are killed by2. The
classification of finite abelian groups therefore forces C4×C2. At D=−280
the complete list is

    (1,0,70), (2,0,35), (5,0,14), (7,0,10),

so the group is C2×C2. The finite enumeration is complete by the reduction
bound, not an unexplained numerical cutoff. Both discriminants are
fundamental. In contrast, the checks at15 and1848 use orders of
discriminant−60 and−7392, of conductor2; they must not be replaced by the
maximal orders of their fields.

By classical genus theory the genus group is A/A². The incoming
[idoneal information theorem](../idoneal_information_channels_20261021/IDONEAL_INFORMATION_CHANNELS.md)
already proves that its uniform-source conditional entropy is log₂|A²|.
At66 the kernel has size2; at70 it has size1. Thus the shared observation
`(σ,τ,g₀)=(144,8,9)` does **not** determine idoneality, even in this tiny
squarefree family. This does not supply a map from divisor parity vectors
to the class group. They are different kernels over the same integer source.

## P4. An explicit non-idoneal form on the entire lifted ray

**Theorem.** Every integer N_k=72k²−6, k>=1, is non-idoneal. When
q_k=12k²−1 is prime, the source is the inherited Berggren lift N_k=6q_k,
with σ(N_k)=144k² and aliquot target72k²+6. Its quadratic-order genus
channel has a nontrivial kernel. For prime q_k that kernel has size at
least2, and the class number is at least8.

**Proof.** Set

    a=6k−1, b=4, c=12k+2.

Then

    b²−4ac = 16−4(6k−1)(12k+2)
             = 24−288k² = −4N_k.

Since a is odd, gcd(a,4,c)=1. For every k>=1,
`0<4<6k−1<12k+2`; hence both (a,4,c) and (a,−4,c) are strictly reduced
primitive forms. They are distinct canonical representatives, and thus
distinct inverse classes. The class group has an element of order greater
than2. Its square subgroup is nontrivial, which is equivalent to failure
of the idoneal condition. This argument does not need q_k prime.

If q_k is prime, N_k=2·3·q_k is squarefree and2 modulo4, so−4N_k is a
fundamental discriminant with three prime-discriminant factors. Classical
genus theory gives |A/A²|=4. As |A²|>=2, |A|>=8. Separately, elementary
multiplicativity gives σ(6q_k)=12(q_k+1)=144k². The exact scale3 Berggren
identity is the incoming one:

    (q_k−2, 2 sqrt(σ(2q_k)), q_k+4)=3 Φ(2k,1),
    Φ(m,n)=(m²−n²,2mn,m²+n²).

Its k1 root is3(3,4,5). This completes every implication claimed. QED.

**Stronger parent family and exact small boundary.** The same mechanism
excludes every `8t²−6` with integer t>=2. For t>=3 use the strictly reduced
inverse pair `(2t−1,±4,4t+2)`, of discriminant `−4(8t²−6)`.
At t2 the pair `(3,±2,9)` does the job for26. The boundary t1 gives2,
which is idoneal: discriminant−8 has only the reduced form `(1,0,2)`.
Thus within positive t the idoneal condition holds **iff t=1**.
The Berggren-compatible family is the subsequence t=3k; divisibility by3
selects the inherited square-sigma construction, not the form obstruction.

At k1 the new form is (5,4,14), one of the explicit non-self-inverse
classes at66. At k2 it is(11,4,26) at282; at k3 it is(17,4,38) at642.
The finite controls give class numbers8,8,16,32 at k=1,2,3,4. They do
not imply monotone growth: at k6 and7 the values are56 and32. The uniform
nontrivial kernel, not a monotonicity extrapolation, is the theorem.

The canonical prime-guard hostile is k5: q=299=13·23. The form theorem
still holds at N=1794, but the formal divisor-sum formula144k² no longer
computes σ(N). The two domains are stated separately. Also70 is a genuine
idoneal square-sigma source outside this ray, so the conclusion must not
be extended to every square-sigma integer.

The possible “66th idoneal number” is an ordinal classification question,
not a claim that the integer66 is idoneal. The present exclusion is
independent of the repo's separately adopted quasi-Riemann input and the
classical completeness qualifications. It does not prove RH or transfer
to exceptional Collatz trajectories.

## P5. Fifteen at two arithmetic depths: equal counts, different actions

For any prime p, P1 gives

    |L(p³)|=1+p+p²+p³=(p+1)(p²+1).

This is also the number of Lagrangian planes in a four-dimensional
symplectic space over F_p. To prove the latter count, every nonzero vector
belongs to p+1 such planes, each plane has p²−1 nonzero vectors, and the
ambient space has p⁴−1 nonzero vectors. Double counting gives the formula.
At p2 both sets have15 elements. This is an equality for all primes,
stronger than a coincidence at one number; it is still not a common action.

The lattice set uses congruence depth p³. The integral shear

    g = [[1,p],[0,1]]

is identity modulo p, yet sends the index-p³ lattice with HNF
`(a,b,c)=(p³,0,1)` to the distinct lattice `(p³,p,1)`. Any action on the
symplectic F_p carrier obtained from reduction modulo p makes g act
trivially. Consequently **there is no equivariant bijection** between
these carriers for that reduction-based action. The argument works for
every prime p. Adding a congruence-depth flag is a possible repair; an
arbitrary numbering of the two equal sets is not an operation-preserving
repair.

At p2 the lattice content split is12+3. The familiar symplectic split9+6
requires a chosen decomposition into two symplectic planes, and the
distance split1+6+8 requires a chosen Lagrangian. Equal total15 does not
identify these splittings. The older marked matching/A5/Richelot maps
remain valid with their specified actions.

There is a related weighted warning at just two steps. A chain of two
index-p inclusions of lattices has (p+1)² possibilities from a fixed
source. A cyclic index-p² endpoint has one intermediate lattice. The
single noncyclic endpoint pZ² has p+1 intermediate lattices. Thus the
endpoint projection has fibres1 and p+1, not constant fibres. At p2,
nine paths end at seven lattices with multiplicities1,1,1,1,1,1,3.
This is the elementary local mechanism behind the extra homothety term
in Hecke composition. It provides an exact example of why path-count
mass must be pushed forward with its fibre multiplicities.

## P6. The shared graph carrier and the bipartite boundary

The Petersen graph is the disjointness graph on the ten two-subsets of a
five-element set; these are also the edges of K5. It has15 edges. For a
connected graph H, a vertex-sign gauge converting every edge sum
`x_u+x_v` to a signed difference of new variables exists exactly when
H is bipartite: the signs must satisfy ε_u=−ε_v on every edge. Multiplying
these constraints around a cycle proves necessity; choosing signs on the
two parts proves sufficiency. Trees and K3,3 admit it. K5 and Petersen do
not; in Petersen an explicit odd cycle is

    {01}, {23}, {04}, {12}, {34}, {01}.

This recovers the precise operation behind the inherited graceful
sign-gauge theorem, **THM-2761 / graph-edge-sum-discriminant-codegree-
factorization-and-graceful-sign-gauge**. It changes the label domain and
does not preserve positivity, absolute differences without their mirror
collision factors, or a square-sum Hamiltonian path. It is therefore not
a proof of the graceful tree conjecture. The golden unit-distance and
finite-field realizations of Petersen in **THM-4701 / kuratowski-graphs-
the-petersen-graph-and-the-finite-unit-circle-graphs** supply further exact
maps, with metric and characteristic guards retained.

## P7. The two prime triples are different binary dependence types

For a prime p let v(p) be the vector of prime-exponent parities of p+1.
On any finite bank of distinct primes these vectors represent a binary
matroid. A squarefree source N has square sigma exactly when the columns
on its support sum to zero. This is the exponent-parity construction
already used by Session82 to build coprime square-sigma blocks; the
matroid language records which subsets carry the dependence.
The later incoming [unifying-relations packet](../unifying_relations_20261028/THEORY.md)
and [incidence packet](../incidence_attractors_20261028/THEORY.md) supply
this same interpretation and explicit graphic/Fano carriers. They are
concurrent inherited work, not a novelty claim for this section.

There is exactly one zero column: v(3)=0. Indeed p+1 square forces
`p=(h−1)(h+1)`, hence h2,p3. Therefore a squarefree square-sigma source
with exactly three prime factors has exactly one of two structures:

- It contains3: the other two columns are equal and nonzero. The support
  is a loop plus a parallel pair, and N splits as3 times a coprime
  square-sigma semiprime.
- It does not contain3: all columns are nonzero and distinct, with their
  sum zero. They form a three-element circuit. No proper nonempty
  subset has zero sum, so N cannot be split into two nontrivial coprime
  square-sigma factors.

To prove the second clause's distinctness, equal columns would force
the third to be zero. The rest is the definition of a minimal linear
dependence. This completely classifies the two dependence types at
three-prime support; it does not classify which primes realize a type.

At66 the signatures are `(3,1,3)`, giving the first type. At70 they are
`(3,6,2)`, giving the second. Thus the two triples' shared total144 is
accompanied by a **different internal factorization of the cancellation**.
The full class-group calculations go the other way:66 has a hidden genus
bit while70 has none. Neither kind of decomposability determines the other.
This is a useful hostile to identifying their two parity spaces.

**The incoming Frobenius transport loses idoneality.** In its explicitly
marked F8 model the incidence packet sends `{2,5,7}` to `{19,89,7}`.
Thus70 maps to11837, and sigma144 to14400; both are squares. But the form
`(3,2,3946)` has discriminant−47348=−4·11837, is primitive and strictly
reduced, and is distinct from its inverse. Hence11837 is non-idoneal,
whereas70 is idoneal by P3. A genuine symmetry of the arithmetic parity
incidence therefore need not preserve the class-genus kernel. The field
and basis marking selecting this Frobenius are retained; no unmarked
canonical action on natural numbers is claimed.

## Research questions arising from these proofs

1. Can the particular class `(6k−1,4,12k+2)` be described by an ideal
   composition law along the Berggren D-ray? Its non-self-inversion is
   proved, but classes live in changing discriminants. A proposed transfer
   must specify an actual correspondence before claiming a recurrence.
2. Can a congruence-depth extension of the fifteen-object symplectic carrier
   retain the integral shear in P5? The first decisive requirement is a
   nontrivial action of the first congruence kernel, not just fifteen labels.
3. Which source-preserving weighted correspondences admit a positive
   remainder after pushing through nonuniform fibres? P5 gives a fully
   solved finite local model. ROOT transport would additionally need its
   native signed weights, source-specific acceptance and integer clocks;
   the lattice count supplies none of those estimates automatically.
4. For six-colour boundary charts, which flags distinguish actions that
   agree on matching cardinalities but differ on integral or geometric
   lifts? A valid answer must preserve the complete extension relation on
   the named ports. Petersen containment or an abstract symplectic chart
   alone is insufficient.

These are **OPEN questions**, not named hypotheses promoted by finite
evidence. The proofs above resolve the stated maps and obstructions only.
