# Forgetting marks: arithmetic transfer, joint ports, and exact descent

Session 88, 2026-10-28. **PROVED below:** the finite weighted-fiber,
integral, modular, descent, composition, and finite-word statements.
**FINITE-EXACT:** the separately declared checker universes.
**CLASSICAL / CITED:** push-pull/groupoidification, conditional expectation,
commuting squares, two-projection theory, and finite-time gossip context.
No new general theory of these classical objects, all-order tournament
speedup, or transfer to a flagship open problem is claimed.

## 0. Inheritance, scope, and the useful extra question

Pinned source: commit `7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51`.
The exact inherited tournament theorem is
`05-knowledge/results/session79_tournament_enumeration_20261027/ROOT_SWITCH.md`,
sections 1--5: pointed tournament classes over a child U are its vertex
orbits, with block `Q_U=1 m^T`, `sum(m)=N`, and `Q^2=NQ`. Its binary rank and
even-order homology statements are **inherited**, not new in this packet.
The older root-orbit count already appears in
`05-knowledge/results/weighted_core_transport_20261020/README.md`, section 3.
The concurrent block-fill lane applies the same factorization algorithmically.

Closest mechanism: complete conditional-fiber refresh in
`05-knowledge/results/weighted_projection_20261024/PROOF.md`, sections 1--3.
Canonical hostile: root attachment Hamming weight loses which parent vertex
is beaten (the transitive/cyclic triangle collision in the inherited packet).
Corrected near miss: averaging classes uniformly instead of weighting by
actual mark multiplicity. Least-used sidecar: the *joint* table for two
forgetting operations, not either marginal list.

Concept board: (i) forgetting marks; (ii) integral averaging denominators;
(iii) joint contingency data; (iv) exact quotient descent; (v) permitted
operation words. Anchor: carry actual marked symmetry data. Niche: modular
and integral loss. Wildcard: finite exact versus asymptotic resampling.

The genuinely useful extension is an exact answer to: **when can two
individually valid forget-and-refresh operations be composed or descended
without recovering joint data?** Sections 4--6 give an iff criterion and a
finite-word obstruction. These are self-contained finite specializations of
classical theory, presented as reusable interfaces rather than novelty claims.

## 1. Finite weighted fibers and mixed marking species

Let `pi:X -> B` be a surjection of nonempty finite sets, and assign each `x in X` a
positive rational weight `w_x` (positive real weights also suffice for the
real Hilbert-space statements). Put `M_b=sum_{pi(x)=b} w_x`. Acting on **column vectors
of observables**, define

    (U h)(x)=h(pi(x)),
    (W f)(b)=sum_{pi(x)=b} w_x f(x),
    Q=UW,       WU=diag(M_b).

Thus `Q_xy=1_{pi(x)=pi(y)} w_y`. With `D_xx=M_pi(x)`,

    Q^2=DQ=QD,    Q 1=D 1,    Q^r=D^(r-1)Q  (r>=1).

No constant-fiber hypothesis is needed. Over any commutative coefficient
ring the displayed polynomial identities hold for any weights. Positivity
is used below for support recovery and Hilbert-space arguments, not here.

Over Q or R, define `P=D^(-1)Q`. Then `P^2=P`; its image is exactly the
fiber-constant functions and its kernel is the functions of weighted sum
zero on every fiber. In the inner product

    <f,g>_w=sum_x w_x f(x) conjugate(g(x)),

P is the orthogonal projection onto that image. It is conditional
expectation for the measure w. On each fiber its Markov kernel resamples
from the exact law `w_y/M_b`. It is reversible because
`w_x P_xy=w_y P_yx`. This is an operator on observables; probability mass
vectors evolve by the transpose.

Consequently

    rank_Q(Q)=|B|,
    tr(P)=sum_x w_x/M_pi(x)=|B|.

For constant mass N, `tr(Q)=N|B|`; variable masses do **not** permit division
by one universal N. Positive support identifies the actual equivalence
classes; the scalar rank only counts them.

### Literal marking construction and its groupoid weights

For each unmarked object U with finite automorphism group G_U, take a finite
G_U-set of admissible marks `A(U)`. The pointed isomorphism types over U are
the orbits `O_i`; their weights are

    m_i=|O_i|=|G_U|/|Stab_G_U(a_i)|,    M_U=|A(U)|.

Forgetting and uniformly rerandomizing the mark gives the block above.
This allows mixed object sizes, colors, ordered tuples of distinct marks,
and object-dependent admissible features. If `A(U)` is empty, U is absent
from the pointed carrier and cannot be counted by its rank. All counting
claims concern the markable objects.

The finite action-groupoid cardinalities are
`mu(U)=1/|G_U|` and `mu(U,a_i)=1/|Stab(a_i)|`, so the relative coefficient
`mu(U,a_i)/mu(U)` is exactly m_i. This is the precise link to push-pull on
finite groupoids; it is not an assertion that arbitrary arrows have
integer orbit-size weights. General finite-groupoid relative coefficients
can be rational, whereas the integral theorems below require integer m_i.
Baez--Hoffnung--Walker [R1] supplies the established groupoidification
framework; no deeper theorem from it is needed for these finite proofs.

**Not a Reynolds identification.** The Reynolds average of G_U acting on
A(U) averages *inside* each G_U-orbit. On orbit-constant observables it is
the identity. The root-switch average averages *between* those orbit
classes, with their multiplicities. The correct map is forgetting the
mark, or restricting uniform averaging on all of A(U) to orbit-constant
functions. Confusing these two averaging operations destroys the bridge.

## 2. Three distinct integral invariants of a single fiber

Suppose `m=(m_1,...,m_k)` has positive integer entries. Set
`M=sum_i m_i`, `g=gcd_i m_i`, and `Q=1_k m^T`.

**Theorem 2.1 (arithmetic separation).**

1. `SNF(Q)=diag(g,0,...,0)` and
   `coker_Z Q = Z^(k-1) direct-sum Z/gZ`.
2. For `K=ker(m^T:Z^k -> Z)` and `C=Z 1_k`,

       Z^k/(K+C) is cyclic of order M/g.

   Thus the rational splitting `Q^k=K_Q direct-sum C_Q` is an integral
   splitting iff `M/g=1` (for positive weights, iff k=1).
3. The least positive integer d such that `d(Q/M)` is integral is `M/g`.
4. For `k>=2`, the mass-defect matrix `L=M I-Q` has

       SNF(L)=diag(g, M,...,M, 0),

   where M occurs k-2 times. Thus
   `coker_Z L = Z direct-sum Z/gZ direct-sum (Z/MZ)^(k-2)`.
   At k=1, L is the 1-by-1 zero matrix instead.

*Proof.* Subtract the first row of Q from every other row, then apply the
integer Euclidean algorithm to the remaining row m; this proves (1).
The map `(m/g)^T:Z^k -> Z` is surjective with kernel K and sends `1_k` to
M/g, proving (2). Write `m_i=g a_i`, `M=g d` with gcd(a_i)=1. All entries
of Q/M are `a_i/d`, so their common denominator is exactly d: if an
integer e clears them, Bezout applied to the a_i implies d divides e.
This proves (3).

For (4), use the unimodular basis `(1_k,e_2,...,e_k)`. In this basis L has
zero first column, first row `(0,-m_2,...,-m_k)`, and remaining nonzero
block `M I_(k-1)`. Every r-by-r minor, for `1<=r<=k-1`, is zero or of
the form `+/-M^r` or `+/-m_j M^(r-1)`; all the latter possibilities needed
for their gcd occur. Their gcd is `g M^(r-1)`, since
`gcd(M,m_2,...,m_k)=g`. The determinantal divisors therefore give g followed
by k-2 copies of M and one zero. QED.

These are different losses: g governs the torsion in coker Q and the
primes where its block disappears; M/g governs exact rational averaging
and integral splitting; M governs the nonzero eigenvalue. Global cokernels
are direct sums over fibers. Listing each block's g is a direct-sum
presentation, not necessarily the globally divisibility-ordered Smith list.

### Modular reduction and the cancellation sidecar

Over F_p, a block has rank zero iff p divides **every** m_i, equivalently
p divides g; otherwise it has rank one. If p does not divide M it is a
scalar multiple of an idempotent. If p divides M then `Q^2=0`, and

    dim(ker Q / im Q)=k-2, if p does not divide g;
                       k, if p divides g.

In the first case its Jordan type is one size-two nilpotent block plus
k-2 zero blocks; k>=2 is forced. A single differential on the full X exists
only if every represented fiber has zero mass in F_p; otherwise restrict
the statement to the zero-mass sector. No unspecified topology is encoded
by this algebraic homology.

There is a useful subtlety: `P=Q/M` is p-integral exactly when
`p` does not divide `M/g`. Hence raw Q may vanish modulo p while the
*integer-lift-renormalized* P has a valid reduction. Example: m=(3),
Q mod3=0 but P=1. This cancellation cannot be reconstructed from Q mod p
alone. In the tournament binary case every m_i is odd, so g is odd:
for even N there is no such denominator rescue, and the inherited
square-zero/rank theorem applies unchanged.

Over F_p, the operator power satisfies `Q^p=M^(p-1)Q`; this comes solely
from quasi-idempotence. Entrywise Frobenius fixes the integer-reduced matrix
and is a different operation. This is not a Frobenius lift or an arithmetic
geometry theorem.

## 3. Composition: the correct measure must travel with the quotient

Let `(X,w) --f--> Y --g--> Z` be finite surjections. Put
`v_y=sum_{f(x)=y}w_x`, `u_z=sum_{g(y)=z}v_y`, and define the normalized push

    (f_bar! a)(y)=sum_{f(x)=y} w_x a(x)/v_y.

Direct expansion proves

    (g f)_bar! = g_bar! f_bar!,
    (g f)^* = f^* g^*,
    E_f E_(gf) = E_(gf) E_f = E_(gf),

where `E_f=f^* f_bar!` and g_bar! uses the actual pushed mass v.
Averaging the Y-fibers uniformly is wrong unless their v masses agree:
for X={1,2,3}, f-fibers {1,2},{3}, and observable (0,0,1), the correct
mean is 1/3, whereas the uniform mean of the two fiber means is 1/2.

Ordinary unweighted summation pushforwards compose. The unnormalized W
of section 1 uses absolute source weights; applying it again with pushed
weights v double-counts v. For a functorial weighted push use the relative
weights w_x/v_y, or keep the counting measure throughout. These are
conventions with different types, not interchangeable formulas.

### Base change with its exact hypothesis

For a Cartesian square `X'=Y' x_Y X`, let `h':X'->X` and `f':X'->Y'`
be the projections over `h:Y'->Y`. Give `(y',x)` weight w_x. Then the
fiber over y' is weight-preservingly identified with the fiber over h(y'),
so

    h^* f_bar! = f'_bar! h'^*.

Multiplying the weights in each y'-fiber by its own positive constant
also preserves this normalized identity. In a general commutative square,
the same identity holds for every observable **iff** each conditional
measure on the f'-fiber pushes under h' to the conditional measure on the
corresponding f-fiber. This follows by testing the indicators of individual
x. A Cartesian set square with arbitrary incompatible weights is not enough.
Replacing a complete reservoir by one representative usually violates this
condition, even though its support label is correct.

## 4. Joint ports: exact commutation and exact descent

Fix a finite X with positive weights w and two partitions A and B. Let C
be their **join as equivalence relations**: points are equivalent if linked
by a chain of A- or B-steps. This is the coarsest partition above both,
not the common refinement. Write `mu(S)=sum_{x in S}w_x` and E_A,E_B,E_C
for their weighted conditional expectations.

**Theorem 4.1 (joint-fiber criterion).** The following are equivalent:

1. `E_A E_B=E_B E_A`.
2. `E_A E_B=E_C` (and then also `E_B E_A=E_C`).
3. For every C-cell C_0, every A-cell A_0 inside it, and every B-cell B_0
   inside it,

       mu(A_0 intersection B_0) mu(C_0)=mu(A_0) mu(B_0).   (*)

4. E_A preserves the B-measurable functions; equivalently E_A has an
   exact induced operator on the B-quotient for **all** B-observables.

Thus every connected component of the weighted A/B contingency table
must have rank one with the specified row and column marginals. In
particular its support must be a complete bipartite graph. This is
conditional independence over the shared C-information, not unconditional
independence and not a claim that marginals determine the joint law.

*Proof.* E_A and E_B are orthogonal projections in the same weighted
inner product. If they commute, their product is the orthogonal projection
onto the intersection of their ranges. Functions constant on both A and B
are exactly C-constant, proving (1)=>(2). Conversely, if their product is
E_C then taking adjoints proves commutation. For (2)=>(3), apply the
product to the indicator of B_0: E_B fixes it, E_A gives
`mu(A_0 intersection B_0)/mu(A_0)` on A_0, while E_C gives
`mu(B_0)/mu(C_0)`. For (3)=>(2), first apply E_B and write the resulting
function as a linear combination of B-cell indicators. Condition (*)
shows that E_A sends each such indicator to its E_C value, so
`E_A E_B=E_C E_B=E_C`. Finally, invariance of ran E_B under the
self-adjoint E_A also implies invariance of its orthogonal complement;
this is exactly commutation with E_B. QED.

The same theorem is the commutative finite case of established
conditional-expectation / commuting-square theory [R2]. It is not a new
subfactor result or an identification with arbitrary Jones projections.
Here the algebra inclusions are literally
`functions on C subset functions on A,B subset functions on X`.

### Raw counting has one additional mass condition

Put `Q_A(x,y)=1_{A(x)=A(y)}w_y`. It descends to the B-quotient iff (*) holds
**and** all A-cell masses inside each C-cell are equal.

Proof: descent applied to the constant observable says the row-sum function
`mu(A(x))` is B-measurable. It is already A-measurable, so is C-measurable.
Dividing by it then reduces descent to Theorem 4.1. Conversely those two
conditions make both multiplication by the mass and E_A descend.
For A={{1,2},{3}}, B={{1,2,3}} with unit weights, E_A descends but Q_A does
not: its row sums are (2,2,1). This is why variable-mass transfer and
normalized conditional expectation must not be silently identified.

### Explicit commutator and a minimal hostile

For x,y in X,

    [E_A,E_B]_(x,y) = w_y * (
      mu(A(x) intersection B(y))/(mu(A(x))mu(B(y)))
      - mu(B(x) intersection A(y))/(mu(B(x))mu(A(y))) ).

So the joint table, with actual cell labels and weights, computes every
mixed response. Calling this order defect “holonomy” would require
additional invertible transport data; these resampling maps are not
invertible and no geometric holonomy is claimed.

The smallest noncommuting example has three unit-weight points:

    A={{1,2},{3}},     B={{1},{2,3}},     C={{1,2,3}}.

For f=1_{1},

    E_A E_B f=(1/2,1/2,0),
    E_B E_A f=(1/2,1/4,1/4).

The missing intersection {3} intersection {1} violates (*). At one or two
points all partitions are nested, so no smaller hostile exists.

Even identical marginal masses are insufficient: the positive contingency
tables `[[3,3],[3,3]]` and `[[2,4],[4,2]]` both have row and column margins
(6,6) and total 12, but only the first satisfies (*). The sidecar needed
for mixed operations is joint incidence, not just both mass lists.

## 5. Exact finite composition versus asymptotic averaging

**Theorem 5.1 (two-operation closure).** Some nonempty finite word in
E_A and E_B equals E_C as an operator on all observables iff E_A and E_B
commute. If they commute, the two-letter word already works.

*Proof.* Only necessity needs proof. Write P=E_A, R=E_B. Restrict
`S=PRP` to ran P; it is self-adjoint, positive, and bounded by identity.
Its 1-eigenspace is ran P intersection ran R. If all its other eigenvalues
were zero, then for a vector u in ran P orthogonal to the intersection,
`||Ru||^2=<u,Su>=0`; thus R preserves ran P, forcing commutation. Hence
noncommutation gives an eigenvector u with `Su=lambda u`, `0<lambda<1`,
and E_C u=0. We have `Pu=u`, `PRu=lambda u`, and Ru is nonzero. After
removing repeated adjacent letters, every nonempty word is alternating;
on u its output is a positive power of lambda times u or Ru (with exponent
possibly zero), never zero. It therefore cannot be E_C. QED.

This is a finite-dimensional consequence of classical two-projection
geometry [R3], proved here rather than imported as new research. Alternating
products still converge to E_C: the same spectral decomposition contracts
all eigenvalues strictly between zero and one. That limit is not a finite
native word. In the three-point hostile, u=(1,1,-2) has
`E_A E_B E_A u=(1/4)u`, giving an exact enduring error component.

The statement concerns fixed linear operator words and equality for every
observable. It does not rule out success for a particular initial vector,
extra operations, signed linear combinations, or different source laws.
Nor does the iff extend to three or more generators: multiple noncommuting
averages can have finite exact products.

**A sharp operation-language positive control.** In the three-point example,

    E_C = (4 E_A E_B E_A - E_A)/3.

Here E_A acts as identity on its own range, not on the whole space. Writing
S=E_A E_B E_A, an equivalent whole-space polynomial is
`E_C=(4 S^2-S)/3`, since S has eigenvalues 0,1/4,1. These formulas use
subtraction and a new denominator, not only native resampling composition.
More generally, for rational weights let S=E_A E_B E_A on the **whole**
nonempty observable space and write its minimal polynomial as `(t-1)q(t)` over Q. Self-adjointness
makes it squarefree, so q(1) is nonzero, and `q(S)/q(1)=E_C`: it is one on
the 1-eigenspace and zero on every other eigenspace. Thus exact recovery is
always possible in the rational operator algebra once those operations and
the full joint operator are allowed. This is not a claim of efficient
recovery or a short universal recurrence.

### A modular hostile even when every mass is invertible

The positive-real finite-word theorem must not be transported to F_p merely
by checking denominators. Over F3 take four points with weights (1,1,1,2)
and partitions

    A={{1,2,4},{3}},       B={{1,3},{2},{4}}.

Every point weight, A-cell mass, B-cell mass and join mass is nonzero mod3.
The reduced conditional expectations are

    P = [1 1 0 2]     R = [2 0 2 0]
        [1 1 0 2]         [0 1 0 0]
        [0 0 1 0]         [2 0 2 0]
        [1 1 0 2]         [0 0 0 1].

Direct multiplication gives `PRP=E_C`, where every E_C row is (2,2,2,1),
but `(PR)_(1,2)=1` and `(RP)_(1,2)=2`. Thus a three-letter native word can
be exactly the join projection despite noncommutation. The failed step in
the real proof is that zero squared norm forces a vector to vanish; the
nondegenerate modular bilinear form has isotropic vectors. This is a
FINITE-EXACT hostile to extending Theorem 5.1 to characteristic three,
not a counterexample to the joint-table criterion itself.

### A denominator obstruction for any finite operation alphabet

Let every allowed rational matrix have entries in `Z[1/S]`, where S is a
set of primes. Finite products and sums over that same ring remain there.
Thus an exact target matrix with an entry whose reduced denominator has a
prime outside S cannot be obtained by those operations. No spectral or
positivity hypothesis is needed for this ring obstruction.

For a weighted integer target fiber, Section 2 gives its exact common
denominator `M/g`; if any of its prime divisors is absent from the local
operation denominators, exact global refresh is impossible by a finite
word. For the three-point hostile both local kernels are dyadic, but the
target E_C has entries 1/3. This obstruction extends to an arbitrary
finite alphabet, unlike Theorem 5.1. Its converse is false: matching prime
support does not establish constructibility. Allowing an external factor
1/3 in the operation language invalidates this particular obstruction.
This is the same denominator mechanism used in finite-time symmetric
gossip; credit [R4], not a new universal consensus theorem.

### A literal arithmetic application: the missing divisor corner

Let p and q be distinct primes, N=pq, and represent its four divisors by
binary exponent pairs:

    1 <-> (0,0), p <-> (1,0), q <-> (0,1), pq <-> (1,1).

Let A forget the q-exponent (retain the p-exponent), and B forget the
p-exponent. With uniform weights their laws are independent. With the
actual divisor-value weights `w(d)=d`, the full measure is the product
`(1,p) tensor (1,q)`, so Theorem 4.1 again proves commuting expectations
and exact full-divisor averaging in two steps.

Restricting to the proper divisors {1,p,q} deletes the (1,1) corner. The
restricted partitions are A={{1,q},{p}}, B={{1,p},{q}}; their join has one
cell, but the {p}/{q} intersection is empty. Therefore (*) fails for
**every** positive retained weight law. No finite word of these two
restricted conditional expectations can give global proper-divisor
averaging. This is an explicit arithmetic-to-operator map, not a matching
integer or a claim about tournament equivalence.

**Operation scope.** Divisor-value weights define a new, declared resampling
model on proper-divisor states. This is not the aliquot iteration itself.
The sum s(pq) is directly computable by addition; the obstruction only
concerns finite words in the two specified conditional averages. Allowing
a direct weighted sum, subtraction, or the target denominator changes the
operation language and bypasses the corresponding obstruction.

With divisor-value weights the total mass is the proper-divisor sum
`s(pq)=1+p+q`; the nontrivial local masses are 1+p and 1+q. On the
A-measurable mean-zero line, S=E_A E_B E_A has eigenvalue

    lambda=pq/((1+p)(1+q)),       0<lambda<1.

Indeed write an A-observable as value a on {1,q} and b on {p}. Its S-image
has values `((1+q+pq)a+pb)/((1+p)(1+q))` and `(a+pb)/(1+p)`.
The resulting two-by-two stochastic matrix has eigenvalues 1 and lambda.
The vector taking a=1 and b=-(1+q)/p has weighted mean zero. Hence this
lambda gives an exact nonvanishing error for every finite alternating word.

| N | retained weights | total s(N) | local pair masses | new-prime obstruction |
|---:|---|---:|---|---|
| 21 | (1,3,7) | 11 | 4,8 | prime 11 absent from all local denominators |
| 35 | (1,5,7) | 13 | 6,8 | prime 13 absent |
| 14 | (1,2,7) | 10 | 3,8 | prime 5 absent |
| 6 | (1,2,3) | 6 | 3,4 | none: all target primes already occur locally |

The N=6 control still cannot finish by a finite word, by the two-projection
criterion. Thus the denominator obstruction is strictly weaker than the
positive-real two-operation obstruction. These statements make no claim
about aliquot trajectories or tournament orders 7/21. The other way an
integer can have exactly three proper divisors, N=p^3, gives the chain
{1,p,p^2}; it is not this two-prime missing-corner carrier.

### Exact boundary of the new-prime obstruction

**Proposition 5.2.** For distinct primes p<q, the proper-divisor model has
no new-prime denominator obstruction iff either `(p,q)=(2,3)` or
`1+p+q` is a power of p.

*Proof.* Let S=1+p+q and L=(1+p)(1+q). A prime dividing both S and L
must divide one of

    gcd(S,1+p)=gcd(q,1+p),
    gcd(S,1+q)=gcd(p,1+q).

Thus every such prime is p or q. Because q>p, the condition q divides p+1
requires q=p+1, and consecutive distinct primes must be 2 and 3. Outside
that special pair, every prime divisor of S that also occurs locally is
therefore p. Hence absence of a new prime forces S to be a power of p.
Conversely, if S is a power of p then p divides 1+q, so its sole prime
already occurs in a local denominator. For (2,3), S=6 and L=12 have the
same prime support. QED.

This classifies the prime obstruction, not finite-word completion: all
these missing-corner models still fail Theorem 5.1's commutation test.
The exact checker independently factors S and L for all 300 distinct
prime pairs p<q<=100. The no-new-prime cases are

    (p,q;S) = (2,3;6), (2,5;8), (2,13;16), (2,29;32),
              (2,61;64), (3,5;9), (3,23;27), (5,19;25), (7,41;49).

The proposition is proved for every prime pair; this finite list is only
the stated bounded control. No infinitude of primes `q=p^a-p-1` is claimed.

### Recovering the earlier excess-two family by an explicit radical map

The inherited Session83 theorem in
`05-knowledge/results/session83_aliquot_divisor_grammar_20261028/FIXED_EXCESS_AND_FIBERS.md`,
section 2, independently proves and credits Li--Ma--Zhang [R5],
Proposition 2.3: integers with at most two distinct prime factors and
`Delta(n)=sigma(n)-2n=2` are precisely

    n=2^b q,     b>=2,     q=2^(b+1)-3 prime.

Consequently the explicit map

    n -> rad(n)=2q

is a bijection from that classical family onto the semiprimes 2q with
q>3 whose proper-divisor resampling has **no new-prime denominator
obstruction**. Indeed their proper-divisor sum is

    s(2q)=1+2+q=2^(b+1).

Conversely Proposition 5.2 at p=2, after excluding q=3, gives
q+3=2^a with a>=3; recover b=a-1 and the unique n=2^b q.
The two-prime excess identity can also be checked directly:
`sigma(2^b q)-2^(b+1)q=2^(b+1)-1-q=2`.

The map discards the exponent b but recovers it on this restricted family
from q+3. It does **not** transport an aliquot edge or preserve excess:
`s(n)=n+2`, while `s(rad(n))=q+3` and `Delta(2q)=3-q`.
For example 20 maps to 10: the first edge is 20->22, whereas 10->8.
Nor is it an infinitude or primality theorem. All mapped resampling models
still have noncommuting expectations and cannot achieve global averaging
by a finite native two-operation word. This recovers an old arithmetic
family as the exact p=2 locus where one obstruction disappears but the
stronger operation-order obstruction survives.

### All unitary splits: abundancy as a conditional-expectation angle

**Theorem 5.3 (split-independent divisor spectrum).** Let A,B>1 be coprime
integers and N=AB. Give a divisor d its actual value as weight. Multiplication
identifies `Div(A) x Div(B)` bijectively with Div(N), preserving the product
weights `w(a,b)=ab`. Delete the unique top pair (A,B), so the retained measure
has total mass `s(N)=sigma(N)-N`.

The remaining pairs have three nonempty flags: neither coordinate top,
A-only top, B-only top. Their exact pushforward weights are

    x=s(A)s(B),       y=A s(B),       z=s(A)B,
    x+y+z=s(N).

On these three weighted states let P retain the A-top flag and let R retain
the B-top flag. Thus P has partition {{x,z},{y}} and R has partition
{{x,y},{z}}. The following hold.

1. The two mean-zero one-dimensional observable spaces have a principal
   angle theta with

       cos^2(theta)=lambda=yz/((x+y)(x+z))=N/sigma(N),
       tan^2(theta)=s(N)/N.

   In particular the angle is independent of the chosen coprime factor
   split. Perfect N gives exactly 45 degrees, abundant N gives a larger
   angle, and deficient N a smaller angle. These are statements about this
   declared weighted-observable model, not geometries of the integer itself.
2. P and R do not commute, and no finite word in them is global refresh E.
   In the enlarged rational operator algebra there is the exact formula

       E=(sigma(N) P R P - N P)/s(N).

3. On the full proper-divisor pair carrier, let P_A average at fixed a
   and P_B at fixed b, using the actual punctured product weights. The
   flag quotient intertwines these operators with P and R. On the whole
   observable space, `S=P_A P_B P_A` has precisely the nonzero eigenvalues
   1 and lambda, each with multiplicity one; its other eigenvalues are zero.
   Therefore the full global weighted expectation is

       E_full=(sigma(N) S^2-N S)/s(N).

   The three-flag formula with P_A in place of P is generally **false**
   on the full space, because P_A can have extra zero-S modes.

*Proof.* Coprimality gives unique divisor factorization, so
`sigma(N)=sigma(A)sigma(B)`. Summing the three flag regions gives the stated
x,y,z. All are positive because A,B>1. On the P-measurable plane write a
function as value u on {x,z} and v on {y}. The restriction of PRP is

    [ (x^2+z(x+y))/((x+y)(x+z))    xy/((x+y)(x+z)) ]
    [              x/(x+y)                      y/(x+y) ].

Its row sums are one, and its determinant is
`yz/((x+y)(x+z))`. Hence its eigenvalues are 1 and lambda. Equivalently,
center the indicators of the A-only and B-only flags. Their weighted
correlation is `-sqrt(yz/((x+y)(x+z)))`; the angle between their one-
dimensional spans takes its absolute value as cosine. Substituting the
factorized x,y,z cancels s(A)s(B) and gives `lambda=N/sigma(N)`.
Since `0<lambda<1`, the formula for tan squared and the abundance comparisons
follow. The missing joint top/top corner violates Theorem 4.1, and Theorem
5.1 forbids finite-word completion. On ran P the two eigenvalues are 1 and
lambda, while P and PRP both vanish on its orthogonal complement, proving
the rational three-state formula.

For the full carrier, conditioning a flag-observable at any proper a gives
the same average, independent of which proper a was chosen, and conditioning
at a=A excludes the top b=B. Thus P_A descends exactly to P, and likewise
P_B to R. This proves the intertwinings with the pullback from the flag set;
it does not say the flags reconstruct a divisor.

More strongly, start with any A-measurable function f(a). If b is proper,
P_B f is the full divisor-value average over a in Div(A); if b=B, it is
the proper-divisor average over a<A. Thus P_B f depends only on whether
b is top. Applying P_A makes the result depend only on whether a is top.
Consequently S has image in the two-dimensional A-top subspace, and its
restriction there is the same two-by-two operator just analyzed. Those two
eigenvalues are nonzero, so S has rank exactly two. Its self-adjointness
makes every remaining eigenvalue zero. The whole-space polynomial
`t(sigma(N)t-N)/s(N)` is one at t=1 and zero at t=lambda and t=0, proving
(3). QED.

**Corollary 5.4 (the precise invariant and the lost information).** For
two integers each having at least two distinct prime factors, their ordered
three-flag conditional-projection pairs are isometrically equivalent, with
the normalized constant vector distinguished, iff their abundancies
`sigma(N)/N` agree. This is independent of the chosen unitary split.

Proof: remove the common constant line. Each projection becomes a rank-one
orthogonal projection in a two-dimensional real space. Choose unit vectors
on its two ranges with positive inner product; their Gram matrix is exactly
`[[1,sqrt(lambda)],[sqrt(lambda),1]]`. Equality of lambda is necessary and
sufficient for an orthogonal map carrying the two ordered lines to the
other two; extend by the prescribed map of normalized constants. This is
the elementary rank-one case of classical two-projection classification
[R3]. For example 6 and 28 have equivalent three-flag pairs because both
are perfect. Their full divisor spaces are not identified. The isometry is
an arbitrary real orthogonal intertwiner: it need not preserve the positive
cone, individual flags, integer lattices, or the labelled divisor source.
Thus this representation remembers abundancy exactly, not the integer N.

**A finite algebra bridge, with no tower claim.** On the mean-zero flag
space set `p=P-E`, `q=R-E`. The preceding spectral formula and its symmetric
counterpart give

    p^2=p, q^2=q, pqp=lambda p, qpq=lambda q.

Set `delta=1/sqrt(lambda)=sqrt(sigma(N)/N)`, `e=delta p`, `f=delta q`.
Then `e^2=delta e`, `f^2=delta f`, `efe=e`, `fef=f`: the standard
**two-generator Temperley--Lieb relations** [R6]. For perfect N, delta is
sqrt(2). The rational unscaled relations are checked directly; the scaled
ones follow by substitution. This constructs only a finite two-projection
representation. It does not supply an infinite positive tower, a Markov
trace with the subfactor axioms, a subfactor index, or an equivalence with
any physical model. Those additional structures cannot be inferred from
these four relations.

This angle formulation uses classical projection geometry [R3]; the
explicit arithmetic-to-operator map and its split-independent spectral
value are the substantive connection. A number with at least two distinct
prime factors admits such a nontrivial unitary split; prime powers do not.
There is no claim that its divisors are reconstructed from the three flags,
or that its aliquot iteration is represented by either averaging operator.

**Full-coordinate hostile.** Take N=12 with A=4,B=3. On the full carrier,
the A-observable `f(1)=2,f(2)=-1,f(4)=0` has P_A f=f and S f=0, because
its weighted sum on proper A-divisors is `1*2+2*(-1)=0`. Its global mean
is zero, but the wrongly extended three-flag expression returns
`-N f/s(N)=-3f/4`, not zero. The S-squared formula handles this mode.

**Coprimality hostile and clone weights.** At N=36, using the noncoprime
split 6 times 6 gives full pair mass `sigma(6)^2=144`, not `sigma(36)=91`.
Deleting the top pair leaves mass 108, not s(36)=55. The flag weights
are (36,36,36), so lambda=1/4, not 36/91. Divisor multiplication now has
clone fibers: pushing the pair weights forward gives `d*c(d)` at divisor
d, where c(d) is the number of representing pairs, rather than weight d.
Merging identical products while silently resetting their mass loses
source multiplicity. Correcting each clone weight changes the product
measure and invalidates the simple product proof. Thus the unitary-split
guard is exact, not cosmetic.

**Finite controls.** The checker constructs the complete pair and divisor
carriers for N in {6,10,12,14,18,20,28,30,36,60,210,496}. It checks every
ordered nontrivial unitary split of each N, including all six splits of
30 and 60 and all fourteen of 210. It verifies weighted multiplication,
flag masses and intertwinings, both correct projector formulas, the full
rank-two trace and polynomial, the perfect/abundant/deficient angle test,
and both hostile examples. These controls do not establish perfect-number
existence beyond the examples or classify any aliquot orbit.

## 6. What transfers to the older repository, and what does not

| Source | Explicit map / preserved operation | Lost data and boundary |
|---|---|---|
| Tournament root switch | Marked type `(U,x)` to U; weight `|Aut(U)x|`; uniform resampling of the physical root | Root-specific adjacency / joint attachment data are lost under expectation. Aut(U)'s Reynolds average is a different map. |
| Complete source reservoirs | Full source configuration `(J,I)` to the actual outside set J; product weight `w(J,I)`; conditional refresh | Pushforward mass is exactly `M(J)=product_(i in J)w_i Z[k-|J|,q-ell(J)]` from `weighted_projection_20261024/PROOF.md` section 1. This is an actual map, not uniformization. Integer currents, pins, charge and legal event records remain sidecars; no all-order native supply follows. |
| Charge lattice | `Z^k -> Z`, `x -> (m/g)^T x`; kernel is the exact zero-charge lattice K | The quotient by K records one integral scalar, not the original state. The index M/g measures failure to adjoin the constant direction integrally. This alone does not transfer any GMC/LRC charge theorem. |
| Labelled prefix/LCA carrier | Finite leaves to their labelled prefix at depth d, with the original leaf weights; nested prefix partitions | Tower identity in section 3 applies literally. Longest common-prefix depth is the deepest shared prefix; bare LCA depth without prefix labels does not specify a fiber map. No Mahler/Collatz itinerary or pointwise orbit claim follows. |
| Characteristic-p arithmetic | Integral Q reduced entrywise; ordinary powers obey `Q^p=M^(p-1)Q` | This is not a geometric Frobenius action. Cancellation before reduction requires the integer lift and gcd sidecar. |

The LCA route is specifically the labelled low-bit-prefix interpretation in
`05-knowledge/results/session56_rmq_lca_connections_20261025/valuation-lca/REPORT.md`
as recovered and scope-qualified in
`05-knowledge/results/lca_ports_reverse_20261025/arithmetic/README.md`,
section 1. For a rooted weighted tree the prefix forgetting maps are
nested, hence commute; competing, nonnested ways to forget ports demand
Theorem 4.1's joint table. This connects actual operations, not a shared
word such as “tree,” “charge,” or “symmetry.”

## 7. Scope of exact verification

Run from this directory:

    python3 check_transfer.py > check-transfer.json
    python3 -O check_transfer.py > check-transfer-optimized.json
    cmp check-transfer.json check-transfer-optimized.json

The standard-library checker uses exact Fractions and explicit exceptions,
not removable asserts. It independently compares matrix commutation,
joint-table factorization, direct row-lumpability and join projections on
all pairs of set partitions through size four with every {1,2}-weight
vector, plus all size-five partition pairs for three declared positive
weight vectors. Arithmetic controls inspect **all square minors** for every
positive weight vector of lengths one through four with entries one through
four, and modular ranks at 2,3,5,7. It checks the explicit eigenvector,
finite hostile words through length 12, identical-margin controls, the F3
isotropic hostile, full/proper-divisor controls for (2,3),(2,7),(3,7),(5,7),
all 300 prime pairs through 100 for Proposition 5.2, guarded excess-two
family members with 2<=b<=10, all stated unitary-split spectra, tower
law, weighted Cartesian base change, gcd cancellation and raw-mass descent.

The finite checks support the proofs; no finite census proves the
all-parameter statements. No claim is made that the theorem supplies the
least-cost tournament enumeration, a minimal universal state encoding, or
a general lifting theorem for arbitrary native port operations.

## References and credit

Sources checked 2026-10-28. All formulas claimed as PROVED above have
self-contained proofs; these citations locate the classical machinery.

- **[R1]** John C. Baez, Alexander E. Hoffnung, Christopher D. Walker,
  *Higher-Dimensional Algebra VII: Groupoidification*, Theory and
  Applications of Categories 24 (2010), 489--553.
  [Author manuscript](https://math.ucr.edu/home/baez/hda7.pdf).
  Used for groupoids/spans and their degroupoidified operators, not an
  unverified claim that arbitrary groupoid arrows have integer weights.
- **[R2]** Remus Nicoara, *A finiteness result for commuting squares with
  large second relative commutant*, Definition 2.1 and introduction.
  [Author manuscript](https://web.math.utk.edu/~rnicoara/lattices.pdf).
  Explicitly records the conditional-expectation commutation equation and
  credits commuting squares to Popa and the basic construction to Jones.
  This packet uses only its finite commutative specialization; it does not
  assert a new Jones index, subfactor, or commuting-square classification.
- **[R3]** P. R. Halmos, *Two Subspaces*, Transactions of the American
  Mathematical Society 144 (1969), 381--389,
  [publisher record/PDF](https://www.ams.org/tran/1969-144-00/S0002-9947-1969-0251519-5/S0002-9947-1969-0251519-5.pdf).
  Classical two-projection context. The indexed publisher text was located;
  direct PDF retrieval returned 403, so no uninspected theorem statement is
  used as a proof dependency. Section 5 supplies its own elementary proof.
- **[R4]** Guodong Shi, Bo Li, Mikael Johansson, Karl Henrik Johansson,
  *Finite-Time Convergent Gossiping*, IEEE/ACM Transactions on Networking
  24(5) (2016), DOI
  [10.1109/TNET.2015.2484345](https://doi.org/10.1109/TNET.2015.2484345).
  [Author-hosted paper](https://people.kth.se/~kallej/papers/multiagent_ieeeton16finite.pdf),
  [arXiv record](https://arxiv.org/abs/1206.0992).
  Credit for finite-time gossip and the power-of-two obstruction in its
  pairwise symmetric averaging model; its model is not silently identified
  with every weighted native reservoir operation.

- **[R5]** Yan Li, Lianrong Ma, Junhua Zhang, *Odd solutions of
  sigma(n)-2n=2 have at least six distinct prime factors*, Publicationes
  Mathematicae Debrecen 82/2 (2013), 255--266, DOI
  [10.5486/PMD.2013.5125](https://doi.org/10.5486/PMD.2013.5125).
  [Publisher paper](https://publi.math.unideb.hu/paper/1759/download/),
  Proposition 2.3 checked directly. Only the two-prime excess-two
  classification is used here; no odd-number existence claim follows.

- **[R6]** H. N. V. Temperley and E. H. Lieb, *Relations between the
  percolation and colouring problem and other graph-theoretical problems
  associated with regular planar lattices: some exact results for the
  percolation problem*, Proceedings of the Royal Society A 322 (1971),
  251--280, [publisher record](https://doi.org/10.1098/rspa.1971.0067).
  Classical algebra credit; publisher metadata was located but direct
  retrieval returned 403. The displayed generator relations were checked
  in section 6.5 of the
  [Vanderbilt-hosted planar-algebra notes](https://cdn.vanderbilt.edu/vu-my/wp-content/uploads/sites/3392/2020/10/15150508/pl21.2019.pdf),
  and their validity for the present matrices is proved by direct scaling.
