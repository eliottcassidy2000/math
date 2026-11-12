# Fixed symmetry cores, finite return modes and their p-adic closures

**Status: PROVED, independently reviewed in this session.** This is a
classification of one-module transitive inflations and their p-adic closures,
not a classification of all tournament counts or a Chan–Pak proof.

## 0. Inheritance and the new question

The exact path-cover substitution kernel is inherited from
[THM-3121](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md),
and transitive substitution from
[THM-3235](../../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md).
The skeleton formula below is a specialization with a direct proof, not a new
claim to the general substitution law. The p-adic tangent idea is suggested by
[THM-4077](../../../../01-canon/theorems/THM-4077-mahler-denominator19-2adic-tangent-full-shift-isometry.md)
and the Fermat-quotient modes by
[THM-4750](../../../../01-canon/theorems/THM-4750-the-camion-function-half-speed-interpolation-and-the-camion-lift.md).
Those concern other observables. We transfer their mechanism, not their
conclusions. Session281 already proves congruence coverage with sparse physical
support; this lane gives a different exact fixed-core closure and size cost.

The new organizing separation is: a fixed modular quotient can have an
unchanging automorphism stabilizer, while the *return-mode rank* controls new
p-adic geometry. Symmetry charge, modal cancellation and prime-specific
Fermat-quotient freezing are distinct coordinates.

## 1. A finite return-mode compiler

Fix a tournament T of order n and a vertex v. Replace v by a transitive
module TT_m, with all external arcs inherited from v, and write F(m) for the
Hamiltonian path count of the resulting tournament, m>=1.

A q-skeleton is a word containing each of the n−1 other vertices once and q
copies of the symbol v, no two consecutive, with each successive arc legal in
T. Let a_q be their number. These are finite combinatorial integers. Put R=max
{q:a_q>0}. Then

    F(m) = sum_(q=1)^R a_q q! S(m,q)
         = sum_(j=1)^R b_j j^m,
    b_j = sum_(q=j)^R (-1)^(q-j) binom(q,j) a_q.          (1)

In particular b_R=a_R>0. Here S is a Stirling number of the second kind.

**Proof.** Split a Hamiltonian path at the maximal runs within TT_m. A run is
increasing, and assigning all m distinct module vertices to q nonempty ordered
runs is a surjection onto q boxes, counted by q!S(m,q). Conversely every
skeleton and assignment gives a unique path. Inclusion–exclusion for
surjections gives the second expression. The omitted j=0 term is zero for
m>=1. This is the one-variable instance of THM-3121.

Each interval between successive v symbols must contain at least two outside
vertices: one outside vertex cannot simultaneously be beaten by v and beat v.
The first and last outside vertices of these intervals must be distinct
out-neighbors and in-neighbors, respectively. Therefore

    R <= 1+min(d+(v),d−(v)) <= 1+floor((n−1)/2).          (2)

The order bound is sharp for every R. Take v and R−1 ordered pairs a_i,b_i,
with v->a_i->b_i->v and all arcs from an earlier pair to a later pair.
The word v,a_1,b_1,v,...,a_(R−1),b_(R−1),v is a legal R-skeleton.

Moreover R=1 iff v belongs to no directed cycle. In a tournament a vertex
belonging to a directed cycle belongs to a directed triangle (take a shortest
cycle through v). Two copies of v and that triangle's other vertices form a
path. Every other outside vertex can be inserted into this path without
changing its existing order: if it cannot be prepended/appended, an adjacent
change in its orientation to the path gives an insertion site. This is a
simple tournament insertion argument, and the two copies of v have identical
external orientations. Thus a_2>0. Conversely a skeleton with two v symbols
contains a directed closed walk through v, hence v belongs to a directed
cycle.

Consequences: if R=1, F is constant. If R>=2, a_2>=1, F is strictly increasing
for m>=1, and

    F(m) >= 2^m−2,    F(m) ~ a_R R^m.                  (3)

Thus every nonconstant fixed-core inflation sequence has only O(log X)
distinct values at most X. This is not merely an empirical sparsity claim.

## 2. Exact p-adic shadow, including the nonunit modes

For a prime p put (L,delta)=(p−1,1) if p is odd, and (2,3) if p=2. For
r=1,...,L define the analytic unit-mode function on Z_p

    f_r(t)=sum_(1<=j<=R, p∤j) b_j j^r (j^L)^t.        (4)

Here the p-adic power is defined by exp(t log(j^L)); j^L lies in
1+p^delta Z_p. At physical integers m=r+Lt,

    F(m)−f_r(t)=sum_(p|j) b_j j^m ∈ p^m Z_p.          (5)

The nonunit modes cannot be interpolated by simply treating j^t as a p-adic
analytic function. Equation (5), rather than a fictitious interpolation,
keeps the transient exact.

**Closure theorem.** The p-adic closure of the physical count set is exactly

    {F(m):m>=1} union union_(r=1)^L f_r(Z_p).           (6)

For an unbounded sequence of physical sizes, first take a fixed residue
r modulo L and then a p-adically convergent subsequence of t. The nonunit
terms vanish by (5), and the limit is in f_r(Z_p). A bounded size subsequence
can only give an actual F(m). Conversely every t in Z_p is approximated by
nonnegative integers t_k tending to infinity in the ordinary sense as well
as to t p-adically. Equation (5) makes F(r+Lt_k) tend to f_r(t).
This proves both inclusions and the claimed closure. No claim that a limit
itself is a finite tournament count is being made.

### 2a. Exact prime obstruction to local interior

For a nonconstant family put S={j>=2:b_j≠0} and g=gcd(S).
Then the p-adic closure in (6) has nonempty interior **iff p does not divide g**.
This statement is independent of whether the first Fermat-quotient test below
vanishes.

If p divides every active j, all unit-mode shadows are the constant b_1.
The closure is the countable set of actual F(m), together with this one limit,
and cannot contain an open p-adic ball. Conversely suppose J={j∈S:p∤j} is
nonempty, of size s. Write alpha_j=j^L and lambda_j=log(alpha_j). These are
nonzero distinct logarithms: log is injective on 1+p^delta Z_p. In

    f_r'(t)=sum_(j∈J) b_j j^r lambda_j alpha_j^t,

the coefficients are all nonzero. The s by s Vandermonde on the distinct
alpha_j shows that f_r'(t0)≠0 for some integer 0<=t0<s. Let A be its finite
p-adic valuation. All Taylor coefficients c_k of f_r(t0+z) are integral,
since v_p(lambda_j)>=delta and lambda_j^k/k! is integral. Choose N=A+1.
Then

    (f_r(t0+p^N u)−f_r(t0))/(p^N f_r'(t0))

has linear coefficient1 and every higher coefficient in pZ_p. Its difference
quotient is1 modulo p, so it is an onto isometry of Z_p by successive digit
lifting. The image contains the exact ball

    f_r(t0)+p^(A+N) Z_p.

This gives an effective test by increasing p-adic log precision until a
nonzero derivative among the finitely many samples is detected; no complexity
bound for that search is asserted. The proof is a finite-exponential inverse
function argument, not an appeal to an unproved distribution of exceptional
primes. In particular C3 inflation has interior at every odd prime even at
base2 Wieferich primes, but not at2. C5 inflation below has g=gcd(2,3)=1,
hence nonempty local interior at **every** prime despite logarithmic ordinary
integer support.

## 3. A constructive open-ball and inverse theorem

Let J={j:2<=j<=R, p∤j, b_j≠0}. If J is empty then f_r is constant. Otherwise set

    e=min_(j∈J) v_p(b_j),   h=e+delta,
    D_r=sum_(j∈J) (b_j/p^e) j^r ((j^L−1)/p^delta) mod p. (7)

If D_r≠0, then

    G_r(t)=(f_r(t)−f_r(0))/p^h

is an onto isometry Z_p -> Z_p. In particular

    f_r(Z_p)=f_r(0)+p^h Z_p.                            (8)

**Proof.** Put lambda_j=log(j^L)/p^delta. This is integral and congruent to
(j^L−1)/p^delta modulo p. For any t,u∈Z_p the quotient of their exponential
difference by p^delta(t−u) is lambda_j modulo p. The higher terms contribute
multiples of p: this follows directly from the exponential series, with
v_p(p^(delta k)/k!) and delta>=1 for odd p, delta=3 for p=2. Hence

    (G_r(t)−G_r(u))/(t−u) = D_r mod p.

Its p-adic valuation is zero, so the map is an isometry. At every finite
precision the p children of a parameter give the p different next output
digits; induction gives unique inverses. Compactness gives onto Z_p.

**Physical residue inverse.** For a>=h and an integer y congruent to f_r(0)
modulo p^h, find the unique t0 in [0,p^(a−h)) satisfying
f_r(t0)=y modulo p^a. Add multiples of p^(a−h) until
m=r+L t>=a. This does not change f_r modulo p^a; by (5), F(m)=y modulo p^a.
One may choose

    m <= max(a,r)+L p^(a−h).                            (9)

Every step is an actual transitive-module inflation. It is a residue inverse,
not an equality F(m)=y. The output tournament has n−1+m vertices. This size bound is not a claim
about total computational complexity of producing or counting the tournament.

A nonzero derivative is sufficient, not necessary, for an open p-adic image.
D_r=0 only says this first normalized tangent test does not apply.

For a single varying unit mode b_j j^m, one can bypass the first-tangent
test exactly: let kappa=v_p(j^L−1). LTE gives an onto isometry after division
by p^(v_p(b_j)+kappa), hence each phase has image

    f_r(0)+p^(v_p(b_j)+kappa) Z_p.

For triangle inflation at p=1093, the exact arithmetic
v_1093(2^1092−1)=2 gives phase balls1+2^r+1093^2 Z_1093. The first tangent
vanishes, but the deeper local chart is still onto. This is a positive control
against mistaking first-order freezing for a p-adic impossibility.

## 4. Tangent observability and symmetry charge

For p>R, use the unnormalized first tangent

    E_r=sum_(j=1)^R b_j j^r q_p(j) mod p,
    q_p(j)=(j^(p−1)−1)/p.                              (10)

The R by R matrix (j^r), r,j=1,...,R, is invertible over F_p. Therefore all
E_1,...,E_R vanish iff every b_j q_p(j) vanishes. More sharply, if s indices
have nonzero b_j q_p(j), no s consecutive phases can all have zero E_r: use
the s by s Vandermonde matrix on those indices. A single dead phase can be
cancellation; every phase being dead detects coefficient divisibility or
simultaneous Fermat-quotient freezing. This is a finite exact test, not an
assertion about the distribution of Wieferich primes.

Let G be any automorphism subgroup of T fixing v. It acts on the inflated
T_m by fixing the transitive module pointwise, and freely on Hamiltonian
paths. Hence |G| divides every F(m). If p>R, the same Vandermonde argument,
now over Z_p using m=1,...,R, yields

    min_j v_p(b_j) >= v_p(|G|).                         (11)

Thus symmetry forces modal p-adic content; content can also exceed the
symmetry bound. For p<=R the Vandermonde is not a unit matrix and (11) is NOT
asserted. These are three separate obstructions: group charge, modal content,
and the normalized tangent.

If T is prime of order>=3 and m>=2, standard modular substitution gives the
same prime quotient T, with TT_m its unique nonsingleton maximal strong
module. Consequently

    Aut(T[v→TT_m]) ≅ Aut(T)_v.                          (12)

The transitive factor is rigid. Its unique size distinguishes its quotient
vertex, and every quotient automorphism fixing v extends. Thus genuine
p-adic richness need not come from new symmetry groups at successive sizes.
Equation (12) uses the standard modular-decomposition fact, not a new claim
to the general decomposition theorem.

## 5. The sharp first dyadic birth and a completely classified example

For every core of order n<=4, (2) gives R<=2. Its only odd mode is j=1, so
its unbounded-size dyadic limit set is the singleton {b_1}. This includes a
constant family, where repeated sizes have that same limit. It refers to
limits along m→infinity, not the topological derived set of a finite set.

The prime cyclic tournament C5, arcs i->i+1,i+2 modulo5, first permits R=3.
Its root skeleton counts are (a_1,a_2,a_3)=(15,13,2), giving

    F(m)=−5+7·2^m+2·3^m.                               (13)

For p=2, the unit shadow is −5+2·3^m. Here e=1, delta=3 and both parity
charts have a unit normalized tangent. The even-size shadow is 13 modulo16,
and the odd-size shadow is 1 modulo16. For m>=4 the discarded term 7·2^m
already vanishes modulo16, so all such physical values lie in these balls.
The first three values are15,41,105, all outside the two balls. Equation (6)
therefore specializes to the EXACT closure

    closure_2 {F(m):m>=1}
      = {15,41,105} union (1+16 Z_2) union (13+16 Z_2).  (14)

Every digit cylinder within either ball has a constructive physical witness
by Section3. Nevertheless the actual count set has only O(log X) values
below X and is not all odd integers. In particular21, an odd2-adic point in
neither ball nor the finite exceptions, remains absent from this family;
7 is likewise absent. The classification does not use Chan–Pak exclusions.

The symmetry quotient in (13) is always the same prime C5; for m>=2 its
rooted stabilizer, and hence inflated automorphism group, is trivial. The
first appearance of an open dyadic part at order5 is a return-mode event,
not the appearance of a new p-group in the automorphism group.

For contrast, an ordinal copy of C3 appended to the inflated triangle has
F(m)=3(2^m+1). Its persistent C3 symmetry forces a factor3, R=2, and p=3>R.
All coefficients have3-adic content1. The two phase balls are0+9Z_3 and
6+9Z_3. This realizes the symmetry-content distinction in (11).

## 6. The unavoidable real-height bill

If R>=2, F is strictly increasing and satisfies (3). To represent all
p^(a−h) residue classes in one ball modulo p^a using this one sequence
requires at least p^(a−h) distinct physical sizes. Thus at least one necessary
representative has size m>=p^(a−h) and actual count at least

    2^(p^(a−h))−2.                                     (15)

The inverse (9) has the same exponential-in-precision scale for the size
parameter, up to a constant factor L; the count is exponential in that size.
This is a lower bound for this fixed family, not for arbitrary tournament
representations. It identifies exactly why a locally onto p-adic chart,
even with unchanged symmetry core and explicit digit inverse, is not an
absolute integer selector. A Chan–Pak bridge would need a physical mechanism
that reduces the ordinary height while retaining the selected residue and a
usable continuation interface. No such mechanism is proved here.

## Reproduction and boundaries

Run `python3 verify.py` in this directory; repeat with `python3 -O verify.py`.
The self-contained verifier compares skeleton counts with independent subset
DP, tests sharp return-rank families, finite p-adic permutations and physical
inverses (including small primes), tangent observability, and exact C5 closure
residues. Enumeration scopes are printed in RESULTS.json. It includes the
hostile dyadic triangle shadow and the base2 Wieferich prime1093. Those finite
checks support, but do not replace, the displayed all-order arguments.
