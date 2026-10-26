# A proposed unifying programme: exact transport, observable rank, and realizable alignment

Session72, 2026-10-26. The proposal is not that tournaments, arithmetic, and
geometry are secretly the same object. It is that many advances answer the
same precise questions about a specified object, operation, observation and
level. The existing `(X,M,O,ρ)` reflection and presentation-first framework
are credited in the lane reports. The abstract vocabulary is inherited;
the calculable theorems and new tests below are the substantive extension.

## 0. First separate five different kinds of equivalence

1. Equal cardinality of sets or of isomorphism types.
2. An explicit bijection of those sets.
3. A bijection compatible with chosen operations and symmetries.
4. An isomorphism of linear response spaces with transported operators.
5. An isomorphism that also preserves the required positive or geometric
   realizability constraints in a specified representation.

None of these levels can be substituted silently for the next. In particular,
the sign-even graph correspondence is initially a type count, then a signed
harmonic correspondence. The actual rank-dependent and prime-circulant codecs
are additional inherited constructions with their own preserved operations.

A basic consistency rule is useful: if O is injective on the ACTUAL state
space X, every operation on X transports through O set-theoretically by
`O F O^{-1}` on its image, including the exact domain of a partial operation.
An alleged failure after complete identification must therefore concern an
extra hidden input, a different state space, equivariance, a requested simple
composition law, or computational cost. It cannot be an information loss in
an injective map on the same X.

This resolves an apparent tension in the repo. A partner-support invariant
may identify an unmarked strong tournament in a bounded census, yet not
identify a selected Hamilton path, marked ports or a source-target alignment.
The subsequent operation may use precisely those extra choices. Likewise a
noncanonical type bijection always transports ordered join to SOME operation;
it need not turn it into ordinary graph disjoint union.

### The formal umbrella: contextual factorization

There is a precise standard theorem beneath this programme, not merely a
slogan. For finitely many declared state sorts with finite-dimensional state
spaces, multilinear operation tensors and linear output observations, let
W_s be the span of all one-hole context responses at sort s, and N_s its
annihilator. Then V_s/N_s is the smallest-dimensional linear compositional
carrier preserving those observations. An input in N_s stays unobservable
after insertion into any operation and any later context, which proves that
the operation tensors descend. Any sufficient linear carrier must kill no
more than N_s, proving minimality.

For nonnegative state-based operations, retain the exact image cones
π_s(R_+^{states_s}); a quotient basis need not make these coordinate orthants.
The declared operation tensors must include the guards and correlations one
wants to preserve. A same-sample copy, a selected port, or a postselection is
not silently supplied by a theorem about a weaker operation language.

This is standard contextual/weighted-automaton linearization, re-proved in
[the response report, section7](../session72-response/REPORT.md), and the
operation-kernel principle was already in the repo. The new theorem families
below CALCULATE these carriers and their obstructions in concrete settings;
they do not claim to invent the abstract factorization principle.

## 1. Unifying theorem I: count identities are shadows of harmonic visibility

For a finite binary vector space V, an affine V-torsor X, and a finite group G
acting affinely, choose coordinates only to write

    g·x = ρ_g x + c_g,     c_(gh)=c_g+ρ_g c_h.

Over Q, the Fourier characters form monomial G-orbits. A dual orbit contributes one
invariant line exactly when its stabilizer twist is trivial. Thus

    number of X/G states = number of admissible dual orbits.

For tournaments, V is the space of edge reversals and the affine term records
what relabelling does to the chosen reference orientation. The admissible
objects are automorphism-SIGN-even graphs. They are not degree-even graphs,
not graphs with an even number of edges, and not a vector subspace closed
under symmetric difference.

The inherited twisted Burnside/Walsh theorems already establish this count
and basis. The operation-envelope result extends the picture concretely:
source-independent, linear-G-invariant difference kernels act diagonally on these
lines; their entire observable algebra and exact invisible kernel are known.
Passing to a quotient by a G-invariant subspace removes the characters that
do not annihilate it. Kernel pushforward sums fibres, while function
conditional expectation averages them. These two operations must not be
confused.

### One operation, three exact slow modes

Retain the weighted operation “reverse one chosen arc,” including all loops
and repeated transitions. For its combinatorial Laplacian:

| State quotient | Surviving graph modes | Positive gap |
|---|---|---|
| Ordinary tournament types | sign-even graphs |4, n≥3; inherited |
| Tournament switching classes, modulo cuts | sign-even Eulerian graphs |8, n≥4 |
| Score-parity classes, modulo Eulerian edge sets | sign-even cuts |2(n−1), n≥3 |

The switching gap is the first surviving C4 mode; a concrete eigenfunction is
`tr(K^4)−n(n−1)(2n−3)` for the skew tournament sign matrix K. The score-parity
quotient is a birth–death chain on weights r≡binom(n,2) modulo2, with rates
`binom(r,2)` downward, `binom(n−r,2)` upward, and `r(n−r)` holding. Its full
spectrum is2a(n−a) for0≤a≤floor(n/2), excluding a when n=2a and a is odd.

These are consequences of an exact harmonic filter, not three unrelated
spectral coincidences. The two nontrivial quotients are not being treated as
a nested chain. These are single-arc-flip quotients, not the literal Hamilton-
path-complement metagraph, whose legal-path dependence needs separate analysis. Booleanizing edges or erasing multiplicity can even destroy
the commutation that makes the proof work.

[Full operation-envelope theorem](../session72-even/EVEN_GRAPH_OPERATION_ENVELOPE.md)
| [Score-parity specialization](../session72-even/SCORE_PARITY_QUOTIENT.md).

### Why this does not identify objects

The inherited no-go result rules out a relabelling-equivariant map on labelled
objects inducing an all-type tournament/sign-even-graph bijection for n≥5. Prime-local
symmetry equalities do not supply a simultaneous all-prime identification;
the mixed C15 symmetry control already differs32 versus31 at n=10.

There are nevertheless genuine type codecs. The SCC/Lyndon codec transports
ordered joins to conditional fusion, not disjoint union. The prime-circulant
codec is a type bijection and an AGL(1,p)-equivariant SURJECTION on
translation-invariant labelled objects, not a labelled bijection. Keeping these distinctions lets the correspondences be useful
without promising a nonexistent natural object identification.

## 2. Unifying theorem II: states, linear responses and positive responses differ

Let X be a nonempty finite set, q:X→Q onto, and M a finite transformation
monoid containing the identity. Fix one coefficient field for the linear
response ranks and reset weights. Allow resets whose ENTIRE state-transition
row depends only on q(x).
Record the deterministic signature

    e(x)=(q(m x))_(m∈M),       E=e(X).

The exact deterministic/lumpable state carrier is E. If an experiment chooses
its word in advance and reads only its terminal q, its minimum linear
response carrier for arbitrary input laws is the incidence matrix

    D[(m,a),e]=1{e_m=a}.

Its dimension is rank D, possibly strictly below |E|. Yet the pointwise
algebra generated by these indicators is all functions on E. Over the reals,
the realizable positive responses form D(R_+^E); every signature column is an
extreme ray. Thus the positive cone is a standard orthant after an invertible
coordinate change exactly when rank D=|E|. For normalized laws its polytope
has |E| vertices and dimension rank D−1.

This is a sharp separation between state minimality, linear minimality and
positive realizability. It does not apply indiscriminately to arbitrary
stochastic systems, trajectory observations, conditional tests or postselection.

### Involution/reset specialization: the missing law is exactly a cycle space

For M={id,τ}, τ²=id, the signature is(q(x),q(τx)). Regard the realized pairs
as edges E of a bipartite graph with two copies of Q. If c is its component
count, then over every field

    rank D=2|Q|−c,
    β=|E|−2|Q|+c.

The β lost directions are its alternating cycle space. A spanning forest
gives an optimal repair: retain the law mass on each nonforest edge and solve
for the remaining masses by leaf elimination. Positivity belongs to the real
cone, while modular stochastic normalization requires invertible denominators.

For ordinary-class observation of fixed-path tournament presentations:

| n | valid literal states | pair cells | terminal rank | cycle loss on pair laws |
|---:|---:|---:|---:|---:|
|4|8|7|7|0|
|5|64|50|23|27|
|6|1024|732|111|621|

At n=5, uniform laws on literal masks{15,21} and{47,5} have identical every
preselected terminal q/reset/τ response. A same-sample observation of q before
and after τ distinguishes them. The positive response polytope has50 vertices
in dimension22 and is not a simplex.

The27 is a kernel dimension on PAIR-CELL laws. The full literal-law kernel is
64−23=41: fourteen within-pair-fibre directions plus twenty-seven cycle
directions. Omitting this distinction would fabricate a new numerical law.

[Proofs, exact native witnesses and state/rank/cone boundaries](../session72-response/REPORT.md).

### A different positivity boundary: full Walsh coordinates

The full signed Walsh basis loses no type-law information, so its transported
positive cone remains simplicial. It is merely not the naive coordinate
orthant. Full tournament converse has a negative diagonal coordinate already
at order4, indexed by the three-edge sign-even star; it still preserves the
correct cone. No independent sign gauge removes that diagonal, while changing
back to point masses does.

This is NOT the same no-go as the nonsimplicial U2 cone. In the latter smaller
space no invertible basis change can turn the exact cone into an orthant.
[Explicit four-vertex matrix and inequalities](WALSH_REALIZABILITY.md).

## 3. Unifying theorem III: reversible tests turn separation into law recovery

Suppose a guarded test P_i has a lawful recovery U_i which restores the
original state on success, while failure is absorbing. The test-return word
U_i P_i is the identity on the guard's successful domain and fails otherwise.
Concatenations therefore realize products of Boolean guards as terminal-alive
responses. If the guards distinguish all valid states, these products span
all functions on them: inclusion–exclusion supplies the complementary factors
needed for point indicators.

This closes a genuine gap between “I can distinguish individual states” and
“I can reconstruct an arbitrary law.” The return mechanism is essential.
Anonymous reset can erase correlations and does not furnish such conditional
state-preserving tests merely because an involution is available.

For a named tournament re-marking, apply p and then its positional inverse.
Whenever p is legal, the inverse is legal and returns exactly the original
fixed-path presentation. Combined with the inherited all-order guard
reconstruction theorem, full named-word terminal-alive response rank is
2^M on valid initial states, M=(n−1)(n−2)/2. An explicitly totalized failure
state is a separate convention, not another valid tournament.

### A precise new open conjecture, with proved bounds

Let d_n over Q be the least r such that products of at most r whole-permutation
Hamilton guards, including the empty product, linearly span every function
on X_n. This counts test-return TESTS, not elementary relocation steps and
not tests restricted to a fixed local block.

Every guard has Boolean polynomial degree at most n−1 in the M free edge
bits. Hence

    d_n ≥ ceil((n−2)/2),        n≥3.

**Proposed sharp-depth conjecture:** equality holds for every n≥3 over Q.
This is OPEN, not an extrapolated theorem. It is verified through n=6:

- n=3,4 need one test.
- n=5 needs two. Raw guard rank is44 over Q,43 over F2;64 selected two-guard
  rows have determinant−1, so depth two works over every field.
- n=6 needs two over Q and F2, certified by1024 independent product rows
  modulo2. No every-field claim is made from an odd determinant alone.

At n=5 the raw20-dimensional rational-law kernel has a concrete Fourier
anatomy: eighteen non-linear-forest characters are absent, and the46 visible
linear-forest character columns have two further dependencies. This44 is a
rank, not the count of sign-even graph types or a proof of a bijection with
any other44-object inventory.

**First untested prediction:** n=7 should have d_7=3, with full rational rank
32768 at product order three. Orders at most two cannot suffice by the
proved degree bound. A proof or counterexample at that exact target would
advance the programme without changing the conjecture after observing data.

## 4. Unifying theorem IV: coherent alignment has a prime-layer criterion

For r≥1 cyclic factors with positive periods d_i, let L=lcm_i d_i. Their joint phase states
are(∏Z/d_i)/simultaneous translation. They are classified by pairwise residues

    r_ij=x_j−x_i mod gcd(d_i,d_j),

with antisymmetry and triangle compatibility modulo each triple gcd. There
are∏d_i/L classes; for an ambient cyclic action of order m, a positive multiple of L, their common
stabilizer has order m/L.

A sparse graph of retained pair tests identifies every ACTUAL phase state
exactly iff every NONEMPTY prime-power layer, consisting of factors whose
periods are divisible by p^a, induces a connected subgraph; singletons count
as connected. An ordinary spanning tree need not suffice:
periods(15,21,35) require all three edges, each demanded by a different prime.
Arbitrary declared sparse residues still require a lifting/consistency test;
injectivity on actual states does not declare every formal labelling realizable.

The prime layers can encode any finite set-family on a nonempty index set:
choose distinct odd primes p_j and set d_i=∏_{j:i∈S_j}p_j. The nonempty layers
are precisely the nonempty S_j. Thus this is a general simultaneous induced-
connectivity design problem, not automatically a laminar or tree problem.
Choosing an odd ambient multiple m≥3 also permits the native carousel
realization, but its size can be very large; no complexity-class reduction
or polynomial-size native encoding is claimed.

The same phase quotient appears in cyclic words, chosen-normal-basis Frobenius
coordinates, and the factors of tensor products of finite fields. It is also
realized by explicit multi-colour carousel tournaments in the odd-period
regime m≥3. This is an operation-
and-stabilizer correspondence, not an identification of tournament substitution
with field multiplication.

The cyclic hypothesis matters. In Paley7, ordered triples(0,1,3) and(0,1,5)
have identical all three named pair orbits but different triple orbits under
its order21 automorphism group. One has a common out-neighbour and the other
does not. With vertex j represented byα^j in F8, α³=α+1, their F2 ranks are2
and3. Pairwise coherence for cyclic translation coordinates cannot be promoted
to pairwise completeness for every group action.

Finally, deterministic relative phases are not probability laws. Recovering
an individual phase tuple from its pair differences is compatible with parity
mixtures whose proper marginals agree but whose global guard responses differ.

[All-arity cyclic theorem and finite/native controls](../session72-cross-domain/THEOREM-CYCLIC-ALIGNMENT.md).

## 5. What prevents a vacuous “theory of everything”

Several older results are decisive boundaries, not exceptions to be ignored:

- **Guards and carries:** affine composition is(A,b,D), with exact domain
  pullback. Equal formal centres or zero displacement do not establish legal
  partial-map commutation. Arithmetic parity/credit domains, positive length
  intervals and planar nonoverlap are different predicates.
- **Infinite rank:** THM-2810's factorial functional L(s^n)=n! has every
  finite Hankel determinant nonzero in characteristic zero. No fixed finite-
  dimensional time-homogeneous linear carrier represents its entire tower.
  Characteristic p instead has the sharp quotient K[s]/(s^p).
- **Coefficient fields:** Walsh coordinates need2 invertible; direct Boolean
  test-return certificates may remain complete over F2. These are different
  representations, not a universal “characteristic two loses everything” law.
- **Scalar norms:** Gaussian squaring preserves a precise determinant/Lorentz
  relation on a canonical raw image; a norm value alone loses placement.
- **Geometry:** flat seam transport is not nonoverlap. Penrose finite guarded
  placements are not a theorem about infinite tiling extension.

[Cross-domain proofs, recovered older mechanisms and corrected scopes](../session72-cross-domain/REPORT.md).

## 6. The research agenda this actually proposes

The unifying object is an OPERATION-RELATIVE REPRESENTATION, together with its
symmetry, response space, realization cone and admissibility domains. The
programme has concrete next obligations:

1. Resolve the sharp guard-product-depth conjecture, beginning at n=7.
2. Use the harmonic filters to design exact operators on chosen quotients;
   do not replace them by Boolean graph supports without a new proof.
3. Optimize sparse alignment tests using prime-layer connectivity, with
   realizability checked separately.
4. Compare source identification with marked-state reconstruction and efficient
   composition. The existing all-order strong partner-support classification
   question remains open beyond its finite census; it is not answered by
   counting the states or by the new terminal-rank theorem.
5. Seek representation changes that reduce computation while retaining the
   correct cone and guards, rather than simply reporting a smaller vector
   space or a matching cardinality.

These are proposed theorems and falsifiable targets with explicit inputs and
failure conditions. They do not identify Collatz, LRC, dissection optimality,
prime enumeration or other open problems with one another, and they do not
assert a shared NP-completeness mechanism without actual reductions.
