# Balanced ramified lifts: prime structure and controlled inherited spectrum

**PROVED, scoped; same-agent algebraic review.** Finite receipts are separate
from the all-depth proofs. No claim of a new unrestricted A000568 census,
general growing-orbit bound, or complement-closed carrier. The stronger
logarithm-free theorem for this particular ramified family is proved in
[LOGFREE.md](LOGFREE.md); the estimates here retain their separate
all-choices and positive-proportion conclusions.

## 1. The actual labelled carrier

Let p be an odd prime and a an integer with v_p(a²+4)=1. Then p=1 mod4,
so p>=5. Put O=Z_p[mu], mu²=a mu+1, delta=2mu-a, and

    R_k=O/p^k O,  N_k=p^(2k),  g=mu^4=(a²+1)+a(a²+2)mu.

O is a ramified quadratic discrete valuation ring, delta is a uniformizer,
and v_delta(p)=2. All symbols and reductions refer to these actual rings,
with coordinate basis (1,mu). A bank tournament has odd symbol s on the
nonzero additive group, constant on g-orbits: x points to y when s(y-x)=1.
The zero symbol is zero. An independent sign chooses one of each antipodal
orbit pair. Orientation representatives remain named ports throughout.

The inherited valuation law, also obtained directly from

    g-1=(a delta/2)[a²+2+a delta],

is v_delta(g^j-1)=1+2v_p(j). The bracket and a are units. The binomial
expansion shows that taking a pth power increases this valuation by two;
the first term has uniquely smallest valuation since p>=5. Multiplication
by an exponent prime to p does not change it. Thus |<g>|=p^k on R_k,
and an element of remaining delta-depth r has orbit length p^floor(r/2).
Summing shell sizes gives b_k=p^k-1 antipodal orbit pairs.

Under R_k -> R_(k-1), each lower pair has p child pairs; every upper orbit
has p times its lower orbit's length. The kernel contributes p-1 pairs.
This is the even-depth specialization of the incoming
[Session109 ramified carrier](../session109_pisano_symmetry_20261030/tournaments/results.md)
and [Session110 ideal-filtration correction](../session110_overnight_20261030/CHECKPOINT_01.md).

## 2. Replace coherent replication by an exact majority constraint

Fix ANY parent bank state. Orient its p children consistently with each
parent orbit O, and require

    sum_(children i of O) epsilon_i = epsilon_O.                (1)

Kernel signs are free. This gives exactly

    2^(p-1) binom(p,(p-1)/2)^(p^(k-1)-1)                        (2)

distinct labelled lifts of each fixed parent. The parent is recoverable
by the signed sum (equivalently majority). This is an observer map, not
a graph quotient: the fibres need not be modules. Converse negates all
signs and commutes with (1). No child is reused in two parent constraints.

### Theorem L1: every balanced lift is prime, for k>=2

First, any composite bank tournament has a proper nonzero ideal as a
module. Here is the argument, including the canonicality requirement.
A vertex-transitive tournament is strongly connected: its ordered strong
components would otherwise contradict transitivity. Classical modular
decomposition therefore gives a canonical partition into maximal proper
strong modules, with prime quotient. Translations preserve this partition,
so its block H through zero is an additive subgroup. The automorphism g
fixes zero and preserves H. The mu coefficient a(a²+2) in g is a unit
mod p; hence H is mu-stable, and therefore an O-ideal. Such ideals are
delta^r O/p^k, 1<=r<=2k-1. Conversely an ideal is a module exactly when
the symbol is constant on every nonzero coset of that ideal.

For r<=2k-2, all lifts of one nonzero unit parent element lie in one
coset modulo delta^r. Every one of its p child orbits meets this fibre.
Condition (1) gives both signs in that fibre, so the ideal is not a module.

For r=2k-1, use instead a parent orbit of valuation one. In the first
one-delta-step lift its remaining depth changes from 2k-3 to 2k-2:
the orbit length grows by p, so there is one intermediate orbit. In the
next step its length stays fixed and it splits into p orbits. Consequently
the p children in (1) all cover this SAME intermediate orbit. They have
both signs and meet a common nonzero coset modulo delta^(2k-1). This
excludes the final possible ideal. All proper modules are excluded. QED.

The classical input is Gallai's theorem, as stated in Boussaïri--Chergui--
Ille--Zaidi, [Theorem16](https://arxiv.org/html/1805.04607v1). The deduction
to ideals is inherited from Session110 and rechecked here. Mere failure
of coherence under reduction by p is insufficient: that session has
130560 composite counterexamples. The majority constraint defeats every
ideal level, including the last one missed by that simpler test.

## 3. Exact transport of old modes

For xi=(u,w), write chi_xi(x,y)=exp(2pi i(ux+wy)/p^k). The multiplication
matrix of g in this basis is symmetric, so dual g-orbits are the same
labelled orbit system. Let

    C_(xi,O)=2 sum_(x in O) sin(2pi <xi,x>/p^k).

The skew adjacency eigenvalue at xi is i sum_O epsilon_O C_(xi,O), up to
the harmless sign of the Fourier convention. Its norm is the maximum
absolute value of these real sums. Each nonzero g-orbit and its negative
give one absolute frequency row; the trivial character gives zero.

Every child orbit covers its parent p times. Thus (1), including the
kernel's zero contribution, proves the exact identity

    lambda_k(p eta)=p lambda_(k-1)(eta).                        (3)

This is an identity of character polynomials. Coherent replication would
give p² instead; see the previous packet's
[coherent spectral obstruction](../pinned_discrepancy_reset_tails_20261029/README.md).
Equation (3) preserves the old normalized norm lambda/sqrt(N), rather
than amplifying it. It retains the actual parent signs, not just their
spectrum or an independently resampled parent.

## 4. Exact energy of every new row

Call xi new when it is not divisible by p. Its delta-valuation v is 0
or1; in coordinates v=1 iff 2u+aw=0 mod p. Then

    sum_O C_(xi,O)^2 = p^(2v) (N_k+p)/(p+1).                   (4)

This identity concerns the fixed coefficient row, before choosing signs.

Proof. Write H=<g>, |H|=p^k, and

    S_r(t)=p^r 1_(t>=r)-p^(r-1) 1_(t>=r-1).

The coefficient-of-one pairing on O is unimodular (matrix I in the basis
(1,mu)). Therefore summing a character of valuation t over the shell
of remaining depth r gives S_r(t). Expanding each orbit's squared sine
sum and dividing by its stabilizer gives

    sum_O C_(xi,O)^2
      = sum_(r=1..2k) p^(floor(r/2)-k)
          sum_(h in H) [S_r(v+v_delta(h-1))-S_r(v)].           (5)

The subtracted term uses that h+1 is always a unit. Both antipodal
orbits are accounted for: |sum chi|²-Re((sum chi)²), summed over all
orbits, is the square sum for antipodal pairs.

For 1<=s<=2k there are p^(k-floor(s/2)) elements h with
v_delta(h-1)>=s; for s<=0 there are p^k. Insert these counts in (5).
When v=0, the even r terms vanish and the odd terms sum to

    1+(p-1) sum_(j=0..k-1) p^(2j) = (N_k+p)/(p+1).

The extra1 is the subtraction at r=1. When v=1 the odd terms vanish
after subtraction, and the even terms give p² times the same expression;
the extra p² occurs at r=2. This proves (4). QED.

In particular (4)<=p N_k, since N_k>=p². This derivation is a shell
orthogonality calculation, not an unproved independence assumption or a
generic Gaussian approximation. The exact cyclotomic checker verifies
it separately for (a,p,k)=(1,5,1),(1,5,2),(1,5,3),(3,13,1).

## 5. Two all-depth spectral conclusions

### Theorem L2: every balanced tower is quasirandom

For any choice of its signs, Cauchy--Schwarz and (4) bound each new mode by

    p sqrt((p^k-1)(N_k+p)/(p+1)) <= sqrt(p) N_k^(3/4).          (6)

Old modes obey (3). At k=1 the trivial row-sum norm bound N_1-1 is less
than sqrt(p) N_1^(3/4). Induction proves

    ||K_k|| <= sqrt(p) N_k^(3/4)                               (7)

for EVERY balanced tower from ANY base bank state. For fixed p its
normalized norm tends to zero. In particular discrepancy between two
vertex sets is bounded by this norm times the square root of their
sizes. This is the precise quasirandom property meant here.

### Theorem L3: better balanced choices always exist

Choose each group in (1) uniformly, independently between groups; choose
kernel signs independently and fairly. The expected symbol on a
nonkernel p²-fibre is its parent sign divided by p, hence constant on
that fibre. A new character has zero sum on each such fibre, and the
kernel mean is zero. Every new Fourier row therefore has mean zero.

For a group with p coefficients, its scalar range is at most
2 sum|C_i|, whose square divided by four is at most p sum C_i².
Kernel bits satisfy the same bound. Independence BETWEEN groups and
the elementary bounded-variable exponential estimate give variance
proxy at most p sum C_i² <= p² N_k. For completeness, the bounded-variable
estimate follows by differentiating log E exp(tX) twice: its tilted
variance is at most range²/4; integrate twice after subtracting the mean.
No independence within a constrained group is asserted.

There are M_k=(p-1)p^(k-1) new absolute frequency orbits. A union bound
therefore shows that at least half of all balanced lifts satisfy

    max_new |lambda| <= p sqrt(2 N_k log(4 M_k)).              (8)

The old normalized modes are fixed by (3), and the right-hand threshold
increases with k. Starting from any base, successive good choices give
a compatible prime tower with

    ||K_k|| <= p sqrt(2 N_k log(4(p-1)p^(k-1))),  k>=2.        (9)

This is existential and includes a finite rejection sampler from uniform
balanced choices with success probability at least1/2 when the bound is
tested exactly. No polynomial-bit-cost spectral algorithm, logarithm-free
bound, or number of nonisomorphic towers is claimed.

## 6. Exact finite controls and native-move boundary

For a=1,p=5 and the all-positive four-port parent at25 vertices, the
checker enumerates ALL160000 balanced625-vertex lifts. All are prime by
the ideal criterion. A best candidate is independently checked by all624
pair-module closures containing zero; translations cover all pairs.
Its new-frequency upper bound is strictly below the inherited-frequency
lower bound. Thus the optimum for this fixed-parent family is exactly

    25 sqrt(5+2sqrt(5)) = approximately76.942088429.

The parent's norm is 5 sqrt(5+2sqrt(5)): it is T5[T5], and the five-cycle
eigenvalues give the displayed radicals. This agrees with the incoming
classification of all sixteen base bank states. A15625-vertex balanced
prime completion has norm in the certified interval recorded in
`lifts.json`, approximately429.154999906. That search is not exhaustive.
Rational sine enclosures from the inherited checker certify these bounds;
floating-point eigenvalues are not used as proof.

The fixed-action bank is NOT closed under a native metagraph edge.
If T and K_P(T) share the same nonidentity labelled automorphism g, then
g preserves their agreement digraph, exactly the directed Hamilton path
P. That path has trivial automorphism group, a contradiction. This is
an application of the inherited
[pair-rigidity theorem](../session105_metagraph_recursion_20261029/metagraph/WEIGHTED_NATIVE_RECURSION.md),
not a new general rigidity theorem. Any relocation-closed extension must
retain P and the changed action embedding. Prime structure, preserved
majority ports and controlled spectrum do not themselves supply that
extension or an unrestricted A000568 recurrence.
