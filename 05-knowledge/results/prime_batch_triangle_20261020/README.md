# Smooth mass, synchronized triangles, and retained boundaries

**Status:** PROVED scoped propositions below, with the explicit same-agent
review in AUDIT.md; FINITE-EXACT controls in receipt.json. No independent-agent
audit or novelty claim. Universal Collatz entry, six-colour forcing, and a
polynomial-time A000568 algorithm remain OPEN.

## Inheritance and board

Closest mechanism: [continuation width and symmetry](../continuation_width_symmetry_20261020/README.md),
Theorem 2A (largest-prime batches), Theorems 3--6 (whole-group support,
exceptional skeleton and the slack layers below six). Its proofs, rather
than a chronological log, are the starting point. Canonical hostile: the
six-vertex pair of triangles with matching cross-orientation, whose rotations
are synchronized. Corrected near miss: bounded exceptional skeleton does
not mean bounded clone core; unlimited triangle words already have slack zero.
Least-used sidecar: the *mass* of the retained smooth multiset, alongside its
named prime powers and actual committed degree.

Anchor: sharper continuation width and the exact synchronized slack-six sector.
Niche: the complete six-colour boundary relation of a geometric triangle strip.
Wildcard: the old divisor equality and the minimum cost of odd prime symmetries.

| Live object | Retained predicate / next test |
|---|---|
| Largest-prime continuation | Future gcd response, actual degree, smooth mass |
| Slack-six words | Actual full-group orbits, not chosen generator orbits |
| Odd prime symmetry | Valuations of the actual group order |
| Divisor balance | Properness and the whole exponent box |
| Collatz block | Original source, exact exponent word and integer clock |
| Triangle strip | All extendible boundary colourings, actual unit distances |

Incoming [black descent](../black_descent_20261020c/README.md) changes cycles
while retaining actual path endpoints. Its prime-orbit block proof supplies
a useful comparison, but its triangle-count descent is not a Collatz map.
The [105/17 packet](../black_descent_20261020c/105_17.md) loses fourteen extra
coordinates after a natural incidence composition: even an exact first map
does not justify forgetting the next boundary. The proofs below test this
loss explicitly. See CONNECTIONS.md for the recovered older routes.

## 1. A stronger uniform upper bound for prime batches

Use the exact recurrence and all three series specified in the inherited
Theorem 2A: ordinary tournaments, clone cores, and fixed-point-free core mass.
Let W(N) be the peak number of occupied degree/response states with odd
lengths processed by decreasing largest prime factor. Internal batch order
is arbitrary. Logarithms in this section are natural.

**Proposition 1 (PROVED).** For every fixed c<1/2, for sufficiently large N,

    W(N) <= exp(sqrt(N) exp(-c sqrt(log N log log N))).       (1)

The same bound, after increasing the threshold N, applies to bit-operation
cost for all coefficients through N. This improves the previous upper bound
exp(O(sqrt(N)/log N)); it does not assert a sharp asymptotic.

### Restricted partitions retain more information than a coordinate count

At a stage in prime batch p, put M=floor(N/p). Old cycles a with largest
prime greater than p contribute their p-smooth part c; new cycles a in the
p-batch contribute a/p. The two coloured multisets have total mass at most
M. They determine future gcd responses; actual committed degree is separate.
Both kinds of part are p-smooth. This is the inherited kernel construction,
with its restriction on *allowed part sizes* now used rather than discarded.

Let P_p(M) count these two-coloured partitions with mass at most M and set

    D_p(s)=product_(prime q<=p) (1-q^(-s))^(-1),  1/2<=s<1.

Allowing even parts only enlarges the count. For t>0, the partition product
and e^(-x)<=x^(-s) give

    log P_p(M)
       <= t M + 2 sum_(r>=1) (1/r) sum_(m p-smooth) e^(-trm)
       <= t M + 2 zeta(1+s) D_p(s) t^(-s).

Uniformly on this s-interval, zeta(1+s)<=zeta(3/2)<3. Choosing
t=(6D_p(s)/M)^(1/(1+s)) gives

    log P_p(M) <= 12 (D_p(s) M^s)^(1/(1+s)).              (2)

The constant is deliberately loose. Also, without the smooth restriction,
the inherited two-colour partition bound gives log P_p(M)=O(sqrt(M)).

The only analytic prime input is the classical bound pi(x)=O(x/log x),
following from the [prime number theorem, DLMF 27.2.3](https://dlmf.nist.gov/27.2.E3).
Partial summation, uniformly for 0<epsilon<=1/2 and L=log p, gives

    log D_p(1-epsilon)
       <= C (1+log(1+L)+exp(epsilon L)/(epsilon L)).       (3)

For completeness, -log(1-q^(-s))<=C0 q^(-s) for s>=1/2.
Partial summation bounds the prime sum by a constant times
1+exp(epsilon L)/L+integral_(log 2)^L exp(epsilon u)/u du.
Below u=1/epsilon the integral is O(log(1+L)); above it, substitute
v=epsilon u and use integral_1^A exp(v)/v dv=O(exp(A)/A).
When A<1 the displayed exp(A)/A already absorbs the endpoint. Constants
are independent of p and epsilon. This proves (3).

### Uniform choice over every intermediate batch

Write v=log N, u=log v, B=sqrt(vu), and a=u/2-2log u>0.

* If L>=B, the unrestricted partition bound gives
  log W <= O(log N)+O(exp(v/2-B/2)).
* If L<=u, use s=1/2 in (2). Equation (3) gives
  log D_p(1/2)=O(u+sqrt(v)), so the logarithm of the right side of (2)
  is at most v/3+O(sqrt(v)+u). This is smaller than v/2-B/2.
* If u<L<B, put epsilon=a/L and s=1-epsilon. Then 0<epsilon<1/2,
  exp(epsilon L)=sqrt(v)/u^2, and (3) is o(B), uniformly in L.
  Since s/(1+s)<=1/2-epsilon/4 and log M<=v-L, (2) gives

      log(log P_p(M))
         <= v/2 - L/2 - av/(4L) + a/4 + o(B)
         <= v/2 - sqrt(av/2) + o(B)
          = v/2 - (1/2-o(1))B.                         (4)

The second inequality is AM--GM. Multiply the partition count by N+1
for actual degree. Its logarithm O(log N) is negligible compared with
the bound on log P. The final length-one stage has polynomial size.
These estimates cover stages inside batches, not only completed batches.
They imply (1) for every c<1/2. In the signed recurrences cancellation can
remove occupied states but cannot add profiles outside this set. The
inherited polynomial bit-size and transition bounds multiply cost by a
polynomial in N, absorbed by the same statement for any c<1/2. QED.

## 2. Prime batching still has a superpolynomial width obstruction

**Proposition 2 (PROVED).** For ordinary-tournament degree/response states,
the same largest-prime batching order satisfies

    W(N) >= exp(Omega(N^(1/3)/(log N)^(2/3))).             (5)

The lower bound also holds for a wordwise linear representation of the
weighted past/future operator at the specified cut. It is not a lower bound
for arbitrary algorithms, and it is not asserted after signed cancellation.

**Proof.** Choose primes P in (L,2L], and a prime q in (2L,4L]. Put
m=|P|, S=sum P and N0=(q+1)S. At the cut immediately after batch q,
the histories indexed by subsets A of P consist of lengths {qp:p in A}.
Every p in P is a future length, as is 1. They are all affordable since
past mass plus the total future bank is at most qS+S=N0. Their responses are

    R_1=|A|,    R_p=|A|+(p-1) 1_(p in A).

Thus all 2^m profiles differ. Ordinary recurrence weights are positive,
so no state is removed by cancellation. For the linear lower bound, pigeonhole
by |A| and sum A to obtain a common cell F of size at least
2^m/((m+1)(S+1)). Its actual degree is constant q sum A. Future words
{p:p in B}, B in F, have the same common degree as each other. The
cross-orbit interaction matrix, apart from positive row/column factors, is

    K(A,B)=2^(|A|^2) product_(p in A intersect B) 2^(p-1).

It is a principal submatrix of a tensor product of positive definite
matrices [[1,1],[1,2^(p-1)]], hence has full rank |F|. This proves the
wordwise dimension statement. Finally the prime number theorem gives
m=Theta(L/log L), S=Theta(L^2/log L), q=Theta(L). Choosing L a sufficiently
small constant times (N log N)^(1/3) makes N0<=N and proves (5) at every
sufficiently large horizon N. QED.

This changes the next target: a polynomial recurrence would have to change
more than this elimination order or linearly compress this exact interface.
The lower and upper exponents remain far apart.

## 3. The synchronized sector at symmetry slack six is exact

As before, f is minimum feedback, G is the full automorphism group,
s is the size of its moved support, and delta=3f-s. An atom word is an
ordered join of singleton and C3 atoms. Put A(x)=1/(1-x-x^3).

Define S_epsilon, epsilon=+1 or -1, on a_i,b_i, i in Z/3, by

    a_i -> a_(i+1),     b_i -> b_(i+epsilon),
    a_i -> b_j iff i=j.                                (6)

**Proposition 3 (PROVED).** If delta=6 and some actual three-vertex G-orbit
is not a homogeneous module, then T has the unique strong-component word

    L -> S_epsilon -> R,                               (7)

where L,R are atom words. Conversely every such word has delta=6.
The two signs give distinct isomorphism types. Consequently this sector's
ordinary class series is exactly

    2 x^6 A(x)^2.                                      (8)

**Proof.** The inherited orbit charge is
d=sum_(nontrivial orbits O) [3(|O|^2-1)/8-|O|], with d<=delta and
delta-d divisible by three. A five-orbit contributes four, a seven-orbit
eleven. At delta six a five-orbit would leave two, impossible; two five-orbits
already cost eight. Hence every nontrivial orbit has size three.
There are at most two backward arcs between distinct orbits in an optimal
order, since each triangle orbit requires one internally.

A fixed vertex sees each triangle orbit uniformly. Between two triangle
orbits A,B, invariance makes the 3-by-3 cross matrix row/column regular.
Nonuniformity means a matching or its complement. Three alternating directed
four-cycles, with no edge common to all, force at least two backward cross
arcs. Thus there is at most one nonuniform pair, and it consumes the entire
external budget. All other arcs between orbits point forward in the optimal
order. Every other triangle is a contiguous homogeneous atom, by the
inherited zero-incidence argument.

Each outside vertex sees A and B uniformly. If it saw them oppositely,
absence of backward external arcs would place all of A before it and all of
B after it, or vice versa. But either separated ordering of a nonuniform
3-by-3 pair has at least three backward cross arcs. Contradiction. Thus
A union B is homogeneous; all outside atoms lie before or after it. It is
strong because both triangles are strong and there are cross arcs both ways.

Choose A to be the side with cross outdegree one (complementing the matrix
just exchanges the two sides), and use its matching to label B. Relative
to the cyclic order on A, B has exactly the two orientations in (6).
They are nonisomorphic: inside the strong component its outdegree-two and
outdegree-three triples are intrinsic, as is their matching; whether that
matching preserves cyclic orientation is invariant. Both have f=4: internal
cost two plus the cross lower bound two, attained by the explicit orders
in receipt.json (each can also be checked directly from (6)). Their full
groups are diagonal C3, since the degree triples are invariant and the
matching forces the two rotations to agree. Hence s=6 and delta=6.
Ordered strong components add f and s and cannot permute. This proves the
converse and uniqueness of (7), and therefore (8). QED.

### Full phase and port accounting

If the outer words contain t triangles and u singletons, then

    Aut(T)=C3^(t+1),
    Z_G=((p_1^6+2p_3^2)/3) ((p_1^3+2p_3)/3)^t p_1^u.    (9)

Every individually fixed port in S_epsilon kills its *one shared* factor;
each outer triangle hit by at least one port kills its own factor. With h
outer triangles hit and b=1 if any central vertex is marked, the stabilizer
has order 3^(t+1-h-b). If u=0, the fixed-point-free fraction is (2/3)^(t+1);
otherwise it is zero. Retaining the atom address, vertex label and sign
epsilon preserves ports; retaining only six moved vertices does not.

### Incoming leaf theorem gives a separate exact boundary

The completed [deficit-two center classification](../deficit_two_centers_20261020/README.md)
and [paired-black leaf theorem](../paired_black_leaf_reduction_20261020/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md)
imply that every paired-black leaf has symmetry slack zero or three.
The non-self-converse leaves are reversal pairs of nonpalindromic atom words
(slack zero); mixed leaves are exactly W -> S4 -> W^op (slack three), since
the rigid strong four-tournament has feedback one and no moved vertices.
This is a direct intersection of the two classifications, not a numerical
analogy between their different notions of deficit.

**FINITE-EXACT central data, with all-order word consequences.** Enumerating
all 720 vertex orders of each six-vertex component gives

| Component | H | group order | path-orbit count h | blue path-orbits b | paired-black degree |
|---|---:|---:|---:|---:|---:|
| S_+ | 45 | 3 | 15 | 9 | 3 |
| S_- | 45 | 3 | 15 | 3 | 6 |

Both are self-converse: anti-isomorphisms are a_i<->b_(-i) for S_+ and
a_i<->b_i for S_-. The incoming SCC center law therefore gives the same
degrees 3 and 6 after any mirrored atom-word wrapping W -> S_epsilon -> W^op.
Nonmirrored wrappings are non-self-converse, with paired-black degree 15.
Consequently the synchronized sector, **after merging converse classes**,
has degree-resolved series

    degree 3:   x^6 A(x^2),
    degree 6:   x^6 A(x^2),
    degree 15:  x^6 (A(x)^2-A(x^2)).                    (9b)

The two orientations have identical f,s,group order and H, while their blue
counts differ. This is a concrete hostile to forgetting the matching's
orientation. The metagraph degree convention here divides literal black-line
multiplicity by two; it is not the ordinary arc-reversal graph.

## 4. The remainder of slack six has a four-endpoint grammar

**Proposition 4 (PROVED, exact recognition grammar, not a class census).**
Every remaining delta-six tournament can be generated as follows.

1. Take m=3 or 4 ordered singleton endpoints. Select two distinct pairs
   that together touch all endpoints; orient these two arcs backward.
2. Fill the m+1 gaps with arbitrary words of singleton and C3 atoms.
   All other inter-atom arcs point forward. Collapse the displayed triangle
   atoms to colour-3 quotient vertices; singletons have colour/weight 1.
3. Accept exactly when both following tests on the quotient Q pass:
   (i) its directed triangles have no common weight-one edge;
   (ii) no three colour-1 vertices induce a homogeneous directed triangle.

Every accepted expansion has delta six. Its group is exactly C3^t, acting
independently on the displayed t triangle atoms. The grammar retains actual
gap words and endpoint/port identities. Different accepted presentations
may describe the same class.

**Proof.** In the nonsynchronized sector all active orbits are homogeneous
triangles. Weighted straightening makes the quotient feedback exactly two,
with weights 3 and 1. Any backward quotient arc touching weight 3 costs at
least three. Hence an optimal quotient order has exactly two weight-one
backward arcs, involving three or four endpoints; all triangles are gap atoms.

For a tournament with at least one directed triangle, feedback at most one
is equivalent to a common edge of all directed triangles. To justify the
nontrivial direction, delete that edge: a shortest remaining directed cycle
of length at least four has a chord other than the single deleted pair and
therefore a shorter directed cycle, impossible. The remaining digraph is
acyclic. Since the original had a cycle through the deleted edge, reversing
that edge is consistent with this acyclic order. With positive vertex
weights the same assertion for cost at most one requires its endpoint
weight product to be one. In our presentation the two-backward bound is
already supplied, so test (i) is precisely weighted feedback two. The empty
triangle family fails (i): any weight-one edge is vacuously common, and the
template supplies such edges.

The coloured quotient has ordinary feedback at most two. By the inherited
whole-orbit bound, its coloured automorphism group has only singleton and
triangle orbits; a five-orbit needs feedback three. A nonuniform pair of
triangle orbits would need total feedback at least four. Thus every active
orbit is a homogeneous triangle. It must have weight one, because a triangle
of weight-three vertices alone costs nine. Test (ii) is therefore equivalent
to rigidity of the coloured quotient. Conversely any such homogeneous
weight-one triangle can rotate independently, so fails rigidity.

Distinct homogeneous directed triangles in a tournament are disjoint: if
they intersect, their two module conditions force a proper two-vertex module
in one of the triangles, impossible. Thus every automorphism of the expansion
preserves its family of displayed homogeneous triangles, unless there is an
additional homogeneous triangle of singleton vertices. Test (ii) excludes
that case. The full group is the product of the internal C3 factors with
the coloured quotient group, which is trivial. Hence s=3t and f=t+2, giving
delta=6. Conversely the actual-orbit quotient of any tournament in this
sector is coloured-rigid and passes both tests. QED.

This improves the general twelve-vertex exceptional bound at delta six to
four fixed endpoints, or one of the two six-vertex synchronized components.
A deduplicated rational series for the entire nonsynchronized sector remains
OPEN here; (8) counts only the synchronized sector.

## 5. Odd prime powers have a sharp symmetry-slack cost

For odd p>=5 define

    c_p=3(p^2-1)/8-p = (3p^2-8p-3)/8.

**Proposition 5 (PROVED).** If p^a divides |Aut(T)|, then

    delta(T) >= a c_p.                                (10)

This is sharp for every p,a. The first costs are c_5=4, c_7=11,
c_11=34, c_13=50. Prime three has zero cost: arbitrarily many quiet
triangles supply arbitrarily high powers of three at delta zero.

**Proof.** On whole-group orbits of odd sizes ell, G embeds in the product
of symmetric groups. Therefore a<=sum nu_p(ell!). For ell<p the valuation
is zero. For ell>=p, Legendre's elementary factorial count gives
nu_p(ell!)<= (ell-1)/(p-1). The ratio

    [3(ell^2-1)/8-ell]/(ell-1)=(3ell-5)/8-1/(ell-1)

is increasing for ell>=3. Thus each orbit charge is at least c_p nu_p(ell!),
and summing proves (10). For sharpness take the ordered join of a copies
of the regular circular tournament R_p, with i->j iff 1<=j-i mod p<=
(p-1)/2. Its natural order attains feedback (p^2-1)/8, matching regularity's
lower bound. Its rotation group acts transitively, so s=p. Its full group
is exactly C_p: the out-neighbourhood of each i is a transitive tournament
with intrinsic source i+1; every automorphism preserves this successor cycle.
Ordered joins multiply groups and add slack. QED.

Costs for different primes cannot simply be added using this proof: one
orbit can contribute to several factorial valuations. The next small-cost
consequence checks that issue explicitly.

**Proposition 5B (PROVED, order versus clock).** If Aut(T) contains an
element of order p^a, for any odd prime p and a>=1, then

    delta(T) >= 3(p^(2a)-1)/8-p^a.                     (10b)

This is sharp, attained by R_(p^a). A permutation of order p^a must have a
cycle of length p^a; its containing whole-group orbit has at least that
many vertices. The increasing orbit-charge function proves the lower bound.
The circular-tournament proof above applies to every odd order, including
prime powers, and its full group is the cyclic successor group. QED.

Thus group order nine can occur at slack zero as C3 x C3, whereas an actual
nine-step symmetry clock requires slack at least 21. Group order 27 can
occur at slack zero as C3^3, while a 27-step clock costs at least 246.
These objects map to the same integer group order and retain different
permutation-cycle information. At slack six all groups have exponent three;
the synchronized component reduces two independent rotations to one shared
rotation, without creating a nine-step clock.

**Corollary 5A (PROVED).** The smallest possible slack of a tournament whose
group order has the mixed balanced shape p^2 q r (three distinct primes) is
15. It is attained by the word C3,C3,R5,R7, whose group order is 315.

Below slack 15, orbits have sizes only 1,3,5,7 (size nine costs 21).
Every regular five-tournament is R5, with group C5. A regular seven-tournament
has no automorphism of order five: such an element would have one five-cycle
and two fixed points; a fixed point would then have degree 0,1,5, or 6,
contradicting regular degree three. Consequently a factor five must come from
a five-orbit and a factor seven from a separate seven-orbit, costing 4+11.
Primes at least eleven already cost 34. Three distinct odd primes must
therefore include this cost or a larger one. The displayed word has f=11,
s=18, delta=15 and |G|=3^2*5*7, proving attainment.

The connection to the inherited F=U+S is exact but scoped: apply its proper
divisor definitions to the **group order**. It classifies p,p^3,p^2qr, not
tournaments or trajectories. For example, a quiet word with t triangles has
order 3^t; for t>0 its group order is balanced precisely when t=1 or 3.
The equally lawful t=2 word fails balance. No dynamical consequence follows.

## 6. A geometric triangle word forgets six-colour restrictions in two links

Let P_ell have ell+1 ordered triangle fibres, each a K3. Between consecutive
fibres join corresponding vertices by a matching, and add no other edges.
The six boundary ports are the first and last triangles with their labels.

**Proposition 6 (PROVED).** P_ell has an induced unit-distance realization
in Q(sqrt(3))^2. For every ell>=2 and every palette size k>=4, every proper
k-colouring of its boundary triangles extends to P_ell. At six colours this
is exactly 120^2 ordered boundary assignments.
In contrast, for three colours the transfer preserves the orientation coset
S3/C3; after two links its boundary relation consists of the 18 ordered pairs
in the same coset. For ell=1 there are 12 such three-colour pairs.

**Geometry.** Take z_0=(0,0), z_1=(1,0), z_2=(1/2,sqrt(3)/2),
v=(3/5,4/5), and vertices z_i+jv. Within-fibre edges and corresponding
neighbour edges have length one. For a side vector e of the triangle,
|kv+e|=1 requires k^2+2k(v dot e)=0. For k=1 this requires v dot e=-1/2,
which fails for all six side vectors. For k=2 it requires v dot e=-1,
which also fails; for k>=3 it is impossible since |v dot e|<=1.
All vertices are distinct (the same coordinate check excludes kv=e).
Thus there are exactly the asserted unit edges, at every length.

**Colouring.** Given two ordered distinct-colour triples a,b from six colours,
an intermediate triple c must satisfy c_i not in {a_i,b_i}, with its three
entries distinct. Each allowed list has at least four colours. Choose c_0,
then c_1, then c_2 greedily; at least 4,3,2 choices remain. Five colours work
by the same argument. At four colours each list has at least two colours.
The union of all three lists has size at least three: otherwise all three
lists would be the same two-set, forcing the three distinct entries of a
into its two-element complement. Thus the three lists have distinct
representatives. An elementary justification in this three-list case:
choose different representatives for two lists; if the third is confined to
those two colours, use a third colour from one of the first two lists, and
give the freed colour to the third. Hence every pair is joined by a two-step
transfer already at four colours. Each triple also has a neighbour (cyclically
rotate its own three colours). Therefore every longer transfer is positive.
For three colours every fibre uses the same palette; a matching-compatible
permutation has no fixed point, hence is one of the two three-cycles.
Transfers stay in one orientation coset. Two steps realize all three phases,
so the relation is exactly as stated. QED.

This is an authenticated boundary no-go for this strip, not a claim about
all triangle-based gadgets. Cross-interior edges or additional palette/list
constraints are necessary to obtain a nontrivial six-colour endpoint relation
from a longer strip. The ternary orientation invariant genuinely survives at
three colours and disappears as soon as one spare colour is allowed; naming a C3 action alone
does not transfer the former obstruction to Hadwiger--Nelson.

## 7. Collatz: arbitrary block histories with exact source cylinders

The inherited affine certificate is retained: T_a(n)=(3n+1)/2^a, with a
the exact valuation at each odd step. For a nonempty word w=(a_1,...,a_m)
of positive exponents, put
K=sum a_i and

    B(w)=sum_(i=1)^m 3^(m-i) 2^(a_1+...+a_(i-1)).

The formal endpoint is (3^m n+B)/2^K. The following lemma supplies an explicit
source sidecar to the inherited sparse block grammar; the affine identity
itself is classical/inherited, not claimed as new.

**Proposition 7 (PROVED).** The positive odd sources realizing w are exactly
the positive integers in the single odd residue class

    n = (2^K-B) (3^m)^(-1) mod 2^(K+1).                 (11)

On that cylinder descent is equivalent to

    2^K>3^m  and  n > B/(2^K-3^m).                    (12)

For any fixed multiset of exponents, B is minimized by the increasing word
and maximized by the decreasing word. Consequently the decreasing-word
threshold in (12), when 2^K>3^m, certifies every larger authenticated source
for **every ordering** of that multiset. Actual words and cylinders remain
distinct; no physical trajectory is reordered.

**Proof.** The final numerator congruent to 2^K modulo 2^(K+1) says the final
quotient is odd. Reducing it successively modulo the partial powers of two
shows that each earlier quotient is odd as well: the first later carry has
exactly the current partial valuation, while all later ones have a larger
valuation and all intervening powers of three are units. This proves every
exact exponent; the reverse implication is direct composition. The residue
is odd since B is odd. Positivity of B and n gives (12).

For adjacent exponents a,b with prefix clock H and suffix length r, swapping
them changes B by

    B(...,a,b,...)-B(...,b,a,...)=2^H 3^r (2^a-2^b).

Thus bubble sorting proves both extremal statements. It retains all factors,
including the prefix integer clock. QED.

This reroute certifies a whole finite union of source cylinders from a
multiset while preserving its separate source addresses. It does not show
every integer ever enters one.

### Ordered block refinement repairs a real source-specific failure

**Proposition 8 (PROVED).** Partition an authenticated exponent word into
contiguous nonempty blocks. Replace each block only in the *upper-bound
calculation* by its decreasing rearrangement, and compose the resulting
formal affine maps on the original source n. Call the answer U_P(n).
Then

    actual endpoint <= U_refined(n) <= U_P(n)             (13)

for every splitting refinement of P. The odd-step count and integer clock
K remain exactly the same. Thus U_P(n)<n is a certificate that survives every
such refinement. At the singleton partition U_P is the actual endpoint.

**Proof.** Every within-block ordering has the same positive slope; by
Proposition 7 its intercept is bounded by the decreasing ordering. Composition
of increasing affine maps preserves the upper bound. Splitting a block
restricts its allowed permutations, so cannot increase the maximum.
This proves (13), with no replacement of the actual source or valuation word.
All permutations preserve the sums and counts defining the integer clocks.
QED.

For a fixed word of length m, the best envelope among b consecutive blocks
is computable exactly. Let D[b,j] be the minimum upper endpoint for its
first j steps, with D[0,0]=n. If U_(i,j) is the decreasing rearrangement of
the segment from step i+1 through j, then

    D[b,j]=min_(b-1<=i<j) U_(i,j)(D[b-1,i]).             (14)

Only defined predecessor entries are used. Positivity of every slope makes
the least preceding endpoint sufficient; induction proves this dynamic
program enumerates the optimum without enumerating every partition.

**FINITE-EXACT hostile and repair.** Source 27 first descends after 37 odd
steps, to 23. A single multiset envelope fails at every prefix up to its
arrival at 1. It also fails forever after that arrival: with 24 exponent-one
steps fixed, the decreasing envelope tends under appended exponent-two steps
to 2(3/2)^24-1=282421147873/8388608>27. Its initial value is also above 27;
the intervening values lie on the affine segment between the two.
Thus this is an all-future failure for that source, not an unfinished search.

Four contiguous blocks of its first 37 steps, ending at positions
6,16,32,37, give the exact upper endpoint

    7067832544468979611 / 288230376151711744 < 27.

The exact DP and exhaustive enumeration of the 630 three-block partitions
verify that fewer than four blocks cannot certify this same prefix.

A four-block *first-exit* conjecture is actually false. Source 63 first
descends to 61 after 34 odd steps, with integer clock 54. All 5,456 four-block
partitions fail, while cuts 5,14,26,30,34 certify the upper endpoint
1095707345239597/17592186044416<63. Thus five blocks are necessary and sufficient
on that prefix. This does not exclude a four-block certificate on a later
prefix of 63. The declared scan of all odd sources 3..2001 and seven additional
sources has block histogram {1:980,2:21,3:4,4:1,5:1}; it is a finite test.

Universal entry and a uniform block-budget bound remain open.
Proposition 2's rank obstruction and Proposition 6's boundary
mixing likewise caution against identifying compression with entry.

## Next productive obligations

* Narrow the gap between (1) and (5), or change the past/future interface.
* Deduplicate the four-endpoint grammar without forgetting triangle phases.
* Add geometrically authenticated cross-interior constraints to the strip;
  test their full six-colour boundary tensor before any quotient.
* Prove a source-uniform arrival rule into contracting Collatz cylinders.
  The arithmetic of individual words is exact; universal arrival is open.
