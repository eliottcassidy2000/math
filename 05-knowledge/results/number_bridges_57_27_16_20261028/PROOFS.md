# Six exact bridges, with their retained coordinates

**Status: PROVED elementary identities and constructions; self-checked with
FINITE-EXACT controls.** No independent external audit is claimed. Classical
recurrences, Schur complements, cyclic character decomposition and the golden
spectrum of a pentagon are reused, not claimed as new mathematics. The added
work is the explicit dictionary for this request, its source-preserving
elimination formulas, and its hostile controls. No major conjecture is closed.

## 1. The three square-sum components and a literal 9−8 defect

Let Q9 have vertices 1,...,9 and an edge between DISTINCT integers whose sum
is a square. Its complete edge set is

    {1,3}, {1,8}, {2,7}, {3,6}, {4,5}, {7,9}.

The components are exactly the paths

    8—1—3—6,       9—7—2,       4—5.

Write I(G;t)=sum over independent vertex sets S of t^|S|. Deleting endpoint8
of the first path gives P3; selecting8 forces deletion of its neighbour1
and leaves P2. Thus

    I(P4;t)=I(P3;t)+t I(P2;t)
            =(1+3t+t²)+t(1+2t)=1+4t+3t².

This realizes the suggested P4=P3+P2 as a **disjoint union of conditioned
configuration sets**, not as equality or disjoint union of the graphs.
The explicit relabellings are (1,3,6)→(9,7,2) in the branch excluding8,
and (3,6)→(4,5) in the branch including8. Store the branch bit to invert.
They preserve adjacency and independent sets. They do NOT preserve integer
labels or the individual square sums. Store those labels for arithmetic use.

At t=1 the counts are8=5+3. The whole graph has8·5·3=120 independent sets;
this numerical equality with5! is not asserted to identify its operations
with permutations.

Remove the central edge {1,3} of P4 instead. The result is two independent
edges, each having3 independent sets. Reinstating {1,3} excludes exactly
the one set {1,3}. Hence

    I(P4;t)=(1+2t)²−t²,       8=3²−1=2³.

The Catalan pair8,9 is therefore present in an actual operation, with the
missing configuration explicitly named. This does not use or reprove the
uniqueness theorem for consecutive perfect powers. Nor does that theorem
explain Hamiltonian extension of square-sum graphs.

There is a direct explanation of the three components. The sum9 edges form
a four-edge matching on1,...,8, leaving9 unmatched. The sum4 edge joins two
matching blocks; the sum16 edge joins a third block to9; the fourth block
stays separate. Their sizes become2+2,2+1,2. Fixed pairs2+2 and8+8 would
be loops and are excluded. The actual square targets, cutoff and distinctness
condition supply this decomposition.

### Exact generalization and its sharp boundary

For any finite simple graph G and edge uv,

    I(G−uv;t)−I(G;t)=t² I(G−(N[u]∪N[v]);t).

The lost independent sets must contain both u and v; everything else must
avoid their closed neighbourhoods. This is a bijective proof, including
weights. At t=1 the difference equals1 iff uv is a dominating edge:
N[u]∪N[v]=V(G). If any vertex remains, its singleton and the empty set
already give at least two losses.

Let D_(r,s) be two adjacent centres with r and s leaves respectively. Then

    I(D_(r,s);1)=2^(r+s)+2^r+2^s,
    I(D_(r,s)−uv;1)=(2^r+1)(2^s+1).

Their difference is always1. Among r,s≥1, the first count is a power of2
iff r=s=1. Indeed, if s>r, factor out2^r: the remaining factor
2^s+1+2^(s−r) is odd and greater than1. If s=r, factor out2^(r+1):
the remaining factor is2^(r−1)+1, which is odd greater than1 unless r=1.
Exchange r,s for the other ordering. This proves a **structural uniqueness
within double stars** for the 9−8 example, without invoking Catalan.

The eight independent sets of P4 are NOT a three-dimensional Boolean cube
under the natural operation of adding/removing one vertex. Their toggle
graph has10 edges, while the three-cube has12. Cardinality8 does not supply
three independent commuting toggles.

### The golden transfer retained at the boundary

Store two endpoint counts (A,B), for independent sets excluding/including
the endpoint. Adding one path vertex sends (A,B) to (A+B,A). This is the
Fibonacci matrix, with characteristic polynomial x²−x−1. Joining two
arbitrary rooted graphs by an edge gives AC+AD+BC for endpoint states
(A,B),(C,D). The forbidden BD term is precisely the joint occupied-port
state. Scalar totals alone lose it.

This is a model for port-preserving grammar, not a Hamiltonian-path or
six-colour extension certificate. Those problems have different boundary
states, as the inherited square-sum port work already emphasizes.

### Endpoint identification recovers the owner's 5/7/8 triple

There is also a literal abstract gluing: identify an endpoint of P3 with
an endpoint of P2, producing P4. For the original components take the
ports2 in(9,7,2) and4 in(4,5). After identification, relabel the resulting
path(9,7,[2=4],5) as(8,1,3,6). This is a graph construction; the relabelling
does not preserve the original arithmetic labels or square targets.

The independent-set boundary counts, excluding/including each port, are

    P3: (A,B)=(3,2),       P2: (C,D)=(2,1).

An independent set on the glued graph corresponds to a pair agreeing on
membership of the identified vertex. Hence among all5·3=15 pairs:

    agree = AC+BD = 6+2 = 8,
    disagree = AD+BC = 3+4 = 7.

This supplies an exact graph interpretation of(5,7,8): total states on P3,
incompatible pairs, compatible pairs. It is a proposed reconstruction of
the intended context, not evidence that this was the historical definition.
The sum–difference map of section3 then gives

    (5,7,8)→(6,10,−4)=2(A,A+B,−B).

Thus it exposes the original endpoint counts, rather than merely matching
an unrelated triple. For any boundary pair(A,B), attachment to the(2,1)
control has agreement/disagreement operator

    [[2,1],[1,2]],

whose sum and difference channels have eigenvalues3 and1. Here the total
is multiplied from5 to15 while the signed imbalance remains3−2=1.

More generally, for ANY two rooted graphs with boundary counts(A,B),(C,D),
the agreement minus disagreement is exactly(A−B)(C−D). The gluing formula
therefore retains the quantity scalar totals would erase.

For r disjoint edges, each with a marked endpoint, the parity of the number
of selected marked endpoints has counts(3^r+1)/2 and(3^r−1)/2. This is
the same two-channel operator iterated, since a single edge has counts(2,1).
At r3 there are27 configurations, divided14/13. These are precisely defined
parity classes, not an asserted interpretation of the unresolved longer lists.

The fresh Session90 tree-edge compression gives a second eight-state set
on P4: tournament completions of a prescribed orientation of its three
edges. Its three free arcs form a Boolean cube with12 single-arc-flip edges.
Our independent-set toggle graph has10. Therefore no bijection between
these eight-state sets intertwines those native single-coordinate moves.
This is a small explicit counterpart of the incoming signed flip-spectrum
obstruction: specify operations as well as cardinality.

## 2. Signed triangles, 16, and the 57/27 cubic clock

Let P_m have entry1 in position(i,i+1 mod m). For odd m,

    det(I+tP_m)=1+t^m.

Only the identity and the whole directed cycle contribute in the determinant
expansion. For m3,t=±2 this gives9 and−7. Consequently

    9−(−7)=16,       9·(−7)=−63,
    (I+2P3)(I−2P3)=I−4P3².

These are the inherited signed-triangle determinants from
[SIGNED_CORE_AND_GRAM.md](../paired_budget_selectors_matrix_20261028/SIGNED_CORE_AND_GRAM.md),
now compared across clocks. At t8=2³,

    det(I+8P3)=513=9·57=27·19,
    (2^9+1)/(2^3+1)=57=Φ18(2).

More generally D_k=1+2^(3^k) satisfies

    D_k/D_(k−1)=Φ_(2·3^k)(2),      v3(D_k)=k+1.

For the valuation formula, put x=2^(3^(k−1))≡−1 mod3 and write x=−1+3s.
Then x²−x+1=3(1−3s+3s²), exactly one factor3. Start with D0=3.
This is the already proved cyclotomic-tower mechanism, not a new claim
about primitive primes. It explains why27 and19 split the next determinant.

With ω²+ω+1=0 and N(a+bω)=a²−ab+b², the character factors are

    det(I+2P3)=3 N(1+2ω)=9,
    det(I−2P3)=−N(1−2ω)=−7,
    N(1+8ω)=N(8+ω)=N(7−ω)=57.

The useful cross-thread identity is literal:

    1+8ω=(1+2ω)(5+2ω),
    N(1+2ω)=3,       N(5+2ω)=19,
    1+8ω=ω conjugate(8+ω).

Thus the same norm57 element, up to unit and conjugation, occurs in the
cyclotomic/torus note and in the conductor40 counterexample. The latter
detects its negative norm character modulo5, since57≡2 mod5. Transporting
the element does not identify the two quotient rings or their predicates.

### Eliminating three edges really preserves ports—and changes source weights

Consider a weighted directed cycle with equations

    x_i+a_i x_(i+1)=b_i.

Keep the named positions0,3,6,... and eliminate the two intervening
coordinates in each block. The exact retained equation is

    x_i+(a_i a_(i+1) a_(i+2)) x_(i+3)
       =b_i−a_i b_(i+1)+a_i a_(i+1)b_(i+2).

Internal reconstruction is

    x_(i+2)=b_(i+2)−a_(i+2)x_(i+3),
    x_(i+1)=b_(i+1)−a_(i+1)x_(i+2).

The eliminated block is triangular with determinant1. Schur complementation
therefore preserves the determinant and the exact solution at every retained
port. This proves det(I+2P9)=det(I+8P3) with actual named ports, rather than
just coinciding numbers. Store original weights, sources and block length3
for reconstruction. Composition multiplies the integer length, not resets it.

For constant weight2, the source functional becomes (b0,b1,b2)↦b0−2b1+4b2.
Its squared Euclidean operator norm is21=1+4+16; equality occurs at
(b0,b1,b2)=(1,−2,4). After a block of length l its squared norm is
(4^l−1)/3. Determinant-preserving coarse graining is therefore not a free
energy contraction. Positive source weights d_j give the precise bound

    |sum_j (−2)^j b_j|²
       ≤ (sum_j 4^j/d_j)(sum_j d_j |b_j|²).

This Cauchy–Schwarz bound identifies the payment, not an independently
bounded Collatz weight system. A valid ROOT use still requires a native map
to this linear carrier and control of the transported source covariance.

The three-vertex matrices belong to the earlier signed-tournament core.
For m>3 a sparse cycle is not a tournament support; after elimination an
edge of weight8 is not an original ±2 edge. This construction is closed in
weighted directed networks, not proved closed under native tournament moves.

## 3. An invertible transform for the two triples

Define

    (x,y,z)=(a+c−b,b+c−a,c−a−b).

Then

    (5,7,8)↦(6,10,−4),       (0,2,10)↦(8,12,8).

Its inverse is

    a=(x−z)/2, b=(y−z)/2, c=(x+y)/2.

Its integer image consists exactly of triples with all three coordinates
of the same parity, an index-four sublattice. The fourth Hadamard output
s=a+b+c is recoverable as s=x+y−z. The complete norm identity is

    s²+x²+y²+z²=4(a²+b²+c²).

Equivalently a signed/permuted four-by-four Walsh–Hadamard matrix acts on
(a,b,c,0). The fourth output can be omitted because the input fourth
coordinate was fixed to zero. If that constraint is erased, the recovery
is no longer valid. This is an exact example of a coordinate supplied by
a known constraint, with its parity guard retained.

An exact hostile is(5,7,8,0) versus(6,8,9,1): the three retained outputs
are(6,10,−4) in both cases, while the omitted output is20 versus24.
The fixed zero coordinate, not the transform alone, supplied the recovery.

For positive triangle sides a,b,c, the strict triangle inequalities are
x>0,y>0,z<0 and Heron's formula reads

    16 area²=−(x+y−z)xyz.

Here area²=300 for(5,7,8). These are algebraic connections regardless of
whether this was the user's original intended interpretation; the original
semantic role of that tuple has not been established from the repo search.

The earlier 12/8 pair is the sum–difference pair for10 and2. Subtracting1
gives11/7; this last operation is an additional affine translation. Likewise
(9,−7) has mean1 and half-difference8. These are explicit reversible affine
maps with a retained origin, not asserted conjugacies of the dynamics.

## 4. The genuine five-cycle has a golden quotient

The corrected valuation map in Session43 uses

    N=2^5 3^7 5^3 7^4 11 13²,
    h_N(b)=max{k:b^k divides N}, b≥2.

Its actual cycle is(4,2,5,3,7);13 is a tail entering2. On the five cycle
coordinates let P be cyclic shift and set H=−(P+P^−1). Direct multiplication
gives

    H²−H−I=I+P+P²+P³+P⁴=J.

On the rational sum-zero subspace, J vanishes. Thus H obeys the Fibonacci
polynomial x²−x−1. If v=e0−e1 and E has columns(v,Hv), then

    H E=E [[0,1],[1,1]].

The columns are independent over Q since the golden polynomial has no
rational root. This is an explicit intertwiner from a Fibonacci transfer
space into the symmetrized valuation-cycle representation. It refines the
inherited golden pentagon spectrum with named cycle coordinates.

Preserved: this specified linear recurrence and the sum-zero subspace.
Destroyed: positivity of counting states and the one-way valuation action;
H uses both forward and reverse moves. The tail13 is excluded. Reversing
the full height map at2 is ambiguous between4 and13, so a tail/history
coordinate is necessary for a claim on the full arithmetic object.

## 5. The earlier {0,2,10} support belongs to the same operator language

THM-2763 records raw dual twist support{0,2,10} in C13, with actual nonzero
coefficients and an independently checked finite-field inverse transform.
An operator with that support has the form f(P13),

    f(X)=a+bX²+cX^10.

For a prime p and rational coefficients,

    det f(P_p)=f(1) Norm_(Q(ζ_p)/Q)(f(ζ_p)).

This is exactly the same cyclic-character factorization used in section2,
with a different p and a different polynomial. If f is nonzero rational,
degree≤10, then f(ζ13^j)≠0 for1≤j≤12, since the minimal polynomial has
degree12. For positive a,b,c the trivial character f(1) is nonzero too.
Thus proper rational support here has full Fourier support.

The finite-field coefficient guard cannot be dropped. For the prototype
f=1+X²+X^10, det f(P13)=159=3·53. Reducing modulo53 produces nontrivial
13th roots at which f vanishes. The checker prints them. Unit weights
are a declared prototype, not the weights of THM-2763. The actual theorem's
finite-field nonvanishing remains its own exact computation.

The common carrier is therefore a cyclic group algebra with its actual
coefficients and characteristic. No equality between a support triple,
a determinant pair, and an integer valuation cycle is being inferred.

## 6. A sixteen-element repair of the golden/Eisenstein comparison

Session76's [quantum named guards](../session76_quantum_named_guards_20261027/quantum_named_guards.md)
already records the common reduction F4 modulo2 and their different
reductions modulo3. The following sharpens the **marked-generator** issue.
Let f=X²−X−1 and g=X²+X+1. Then

    (X−2)(g−f)−2f=−2,
    Z[X]/(f,g) = F2[X]/(X²+X+1) = F4.

Indeed the displayed identity forces2=0 in the common quotient, and after
reduction modulo2 the polynomials agree and are irreducible. Thus F4 is
the universal ring in which ONE specified element satisfies both golden
and Eisenstein relations. An equality of the two distinguished generators
cannot lift to characteristic4. This does not prohibit an isomorphism of
the unmarked quotient rings with different generator images.

In fact there is an explicit repair. In R=(Z/4)[ω]/(ω²+ω+1), set

    t=2+3ω.

Then t²−t−1=0, and the basis map (1,φ)→(1,t) has determinant3, a unit
modulo4. It gives an isomorphism

    Z[φ]/4 ≅ Z[ω]/4,       φ↦2+3ω.

Both rings have16 elements. Modulo2 the map sends φ to ω, but the correction
2(1+ω) is essential modulo4. Their marked multiplicative clocks differ:

    order(ω)=3,       t³=1+2ω≠1,       order(t)=6.

This is a concrete example of a shared three-phase reduction acquiring a
second clock sheet after refinement. The extra phase is represented by an
actual ring element, not supplied by a numerical analogy. Over modulo3,
the original Eisenstein relation instead has a double root while the golden
relation is irreducible and gives F9; the old field-order8 observation has
that separate source.

The cheap hostile is to identify φ with ω unchanged modulo4: f(ω)=2+2ω≠0.
The positive repair changes the generator and retains its order6. This
does not identify global number fields, their archimedean geometry, or their
native Collatz branch rules. It suggests testing changing-generator rather
than generator-fixed transports at deeper powers of2.

### The changing-generator repair works at every binary depth

This follows by the elementary simple-root lifting argument. Fix the root
t1=ω modulo2. If t_k is a root of f modulo2^k, its four possible lifts are
t_k+2^k δ, with δ in F4. Modulo2^(k+1),

    f(t_k+2^kδ)/2^k = f(t_k)/2^k+(2t_k−1)δ mod2.

Since 2t_k−1≡1, there is exactly one correcting δ. The coefficient of ω
in every t_k is odd, so (1,φ)→(1,t_k) is an invertible basis change at
every depth. We obtain compatible isomorphisms Z[φ]/2^k≅Z[ω]/2^k, and
an isomorphism of their 2-adic completions, with the generators changed.
Each ring has4^k elements. The receipt constructs the corrections through
k16; the lifting proof covers every k.

The marked golden clock is exactly3·2^(k−1). This familiar Fibonacci
period law has a short direct proof here. Modulo2, φ has order3. Also

    φ³−1=2φ,       φ⁶−1=4φ³.

Thus v2(φ^(3·2^r)−1)=r+1 for r0,1. Once z≡1 mod4, z+1 is exactly
divisible by2, so v2(z²−1)=v2(z−1)+1. Induction proves the valuation
for all r. It supplies both the order and its sharpness modulo2^k.
By contrast ω still has order3 at every depth.

There is therefore a uniform algebraic lift but no uniform bound on the
lost clock if φ is replaced by ω. This is a precise positive construction
and a precise obstruction, rather than a declaration that the analogy
fails. One new binary clock bit is needed at each refinement; the full
ring refinement has four choices per element, so these are different
information counts.
