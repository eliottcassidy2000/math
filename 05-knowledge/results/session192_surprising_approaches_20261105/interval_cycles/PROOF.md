# Interval-cycle geometry, matching spectra, and directed continuation ports

Session 192, 5 November 2026. **PROVED elementary constructions and identities;
FINITE-EXACT controls. Independently reviewed.** These results address the
Hamilton-path spectrum and induced growth through connections to hard-core
interval gases, weighted path matchings, and the component-port calculus used
for square-sum graphs. No all-tournament spectrum theorem or literature-priority
claim is made.

Let H(T) count directed Hamiltonian vertex words, with H(empty)=1. Let

    Z_T(z) = sum_Gamma z^|Gamma|,

where Gamma ranges over collections of pairwise vertex-disjoint directed odd
cycles, each of length at least three and counted modulo cyclic rotation.
The inherited Grinberg--Stanley OCF is H(T)=Z_T(2). Thus Z is the independence
polynomial of the intersection graph of all directed odd cycles.

## 1. The exact interval map

For p=(p_1,...,p_k), all p_i>=2, the interval chain IC(p) has spine vertices
x_0,...,x_N, N=sum p_i. Write b_0=0, b_j=sum_(i<=j)p_i. All arcs point forward
except x_(b_j)->x_(b_(j-1)), one reversed arc per block. This is Moon's known
family; its product H=product_i(2^(p_i-1)+1) is not new here.

**Theorem 1 (cycle normal form).** Every directed cycle of IC(p) consists of:

- every block endpoint b_i,b_(i+1),...,b_j for one 0<=i<j<=k;
- an arbitrary subset S of the nonboundary interior vertices strictly between
  b_i and b_j;
- the cyclic word

      b_i, S in increasing spine order, b_j,b_(j-1),...,b_(i+1).

The choices that are directed odd cycles satisfy

    |S| congruent to j-i (mod 2).

Consequently, putting q_t=2^(p_t-1), interval [i,j] supports exactly

    w_[i,j] = (1/2) product_(i<t<=j) q_t                       (1)

odd cycles. Two such cycles meet if and only if their closed endpoint intervals
meet. The odd-cycle intersection graph is therefore an interval graph obtained
by replacing each interval [i,j] by a clique of w_[i,j] interchangeable vertices.

**Proof.** A cycle must use a reversed arc, because all other arcs increase the
spine index. Its minimum and maximum are block endpoints. Every spine cut
strictly between consecutive block endpoints has exactly one right-to-left arc,
namely that block's reversed arc. A cycle spanning that cut must use it. Thus a
cycle with extremes b_i,b_j uses every arc
b_j->b_(j-1)->...->b_i. These arcs are consecutive in the cycle: at each internal
endpoint the incoming and outgoing cycle arcs are already fixed. Its remaining
portion runs increasingly through a subset of the nonboundary interiors. The
stated word is therefore forced. Conversely every such word is a cycle, except
for the impossible two-vertex word with j=i+1 and S empty; its length is even,
so it does not enter the odd-cycle claim. The parity condition follows from its
length j-i+1+|S|. There is at least one interior vertex, so exactly half of all
subsets have the prescribed parity, proving (1). Intersecting endpoint intervals
share a mandatory boundary vertex. Disjoint endpoint intervals use disjoint
spine ranges. This proves the conflict assertion. QED.

This is an exact map of combinatorial objects and conflicts, not merely an
analogy with planar matching interfaces. It carries multiplicities of distinct
cycles on the same vertex set and does not impose noncrossing conditions on
arbitrary tournaments.

## 2. The full polynomial and a weighted matching representation

Write Z_j(z)=Z_(IC(p_1,...,p_j))(z), with Z_0=1. Weighted interval scheduling gives

    Z_j = Z_(j-1) + (z/2) sum_(i=1)^j
                        (product_(t=i)^j q_t) Z_(i-2),       (2)

where Z_(-1)=1. The last chosen interval has endpoint range [i-1,j]; previous
chosen intervals must end at most i-2, which explains the exact boundary index.
There is no orientation factor and adjacent closed intervals are incompatible.

**Theorem 2 (two-state recurrence).**

    Z_0=1,
    Z_1=1+(q_1/2)z,
    Z_j=(1+q_j)Z_(j-1)+q_j(z/2-1)Z_(j-2), j>=2.           (3)

Indeed, the interval sum in (2) equals q_j times the preceding sum plus
q_j Z_(j-2). Substituting (2) at j-1 proves (3).

At z=2 the second term disappears, recovering Moon's product. This is a precise
loss of information in the scalar observation H, not a claim that the full
polynomial determines a tournament or all of its continuation responses.

Put P_j=product_(t<=j)(1+q_t), t=z/2-1, and W_j=Z_j/P_j. Then

    W_0=1,
    W_1=1+t a_1,
    W_j=W_(j-1)+t a_j W_(j-2),

    a_1=q_1/(1+q_1),
    a_j=q_j/((1+q_(j-1))(1+q_j)), j>=2.                    (4)

Thus, for the path with vertices 0,...,k and edge weights a_1,...,a_k,

    Z_k(z)=P_k sum_(M a matching of the path)
                         (z/2-1)^|M| product_(e in M) a_e. (5)

This is an algebraic matching representation. It is not a positive, weight-
preserving bijection from odd-cycle collections to matchings at general z:
the substitution z/2-1 permits signed cancellation.

### 2.1 All-order real-rootedness, including clawed conflict graphs

**Corollary 3.** Every zero of Z_(IC(p))(z) is real and strictly negative. Its
degree is floor((k+1)/2). Its coefficients obey the usual Newton inequalities
for a real-rooted polynomial with nonnegative coefficients.

For a self-contained real-root argument, let A be the real symmetric weighted
adjacency matrix of the path, assigning edge e the entry sqrt(a_e). Expanding
its characteristic determinant by path matchings gives

    det(lambda I-A)=lambda^(k+1) M(-1/lambda^2),

where M(t)=sum_M t^|M| product a_e. All eigenvalues of A are real. Hence every
root of M is real negative: a nonreal or positive root would yield a nonreal
nonzero zero of the displayed characteristic polynomial. The affine substitution
in (5) makes every root of Z real. Since Z has nonnegative coefficients and
constant coefficient 1, none is nonnegative. Its degree is the maximum number
of disjoint closed endpoint intervals, namely floor((k+1)/2), attained by the
single-block intervals [0,1],[2,3],... . QED.

This does not follow merely from interval graphs: the independence polynomial
of the claw K_(1,3) is 1+4z+3z^2+z^3 and has nonreal roots. Nor does it follow
from claw-freeness. In IC(2,2,2,2,2), the three disjoint triangles

    (0,1,2), (4,5,6), (8,9,10)

are all met by the odd cycle (0,1,10,8,6,4,2). They form an induced claw in the
cycle-intersection graph, but

    Z=1+57z+30z^2+z^3

still has only negative real roots by Corollary 3. The repository already
records failure of real-rootedness for arbitrary tournaments at order nine
(THM-020/THM-025); our result is restricted to this explicit all-order family.

### 2.2 A scalar collision with a visible second coordinate

The strong order-nine chains IC(2,3,3) and IC(3,2,3) both have H=75, yet

    Z_(2,3,3)=1+33z+2z^2,
    Z_(3,2,3)=1+29z+4z^2.                                 (6)

The middle placement of the short block changes the number of disjoint-cycle
pairs. Under OCF weights 2^|Gamma|, the expected cycle numbers are respectively
82/75 and 90/75. More generally, (4) gives this expectation as sum_j a_j.
Therefore the product formula is not a product-independence assertion about
actual cycles, and scalar H loses arrangement-sensitive information.

### 2.3 One arc destroys the spectral conclusion

The real-rootedness theorem is not robust under arbitrary small edits. Start
with IC(2,3,3), whose polynomial is 1+33z+2z^2, and reverse the single arc
3->6. The new outgoing rows are

    (506,508,473,432,480,196,392,256,32).

Direct odd-cycle enumeration and disjoint-collection enumeration give

    Z=1+48z+11z^2+z^3, H=149.

The cubic discriminant is exactly -159431, so it has one real root and two
nonreal conjugate roots. The finite certificate retains all row data and
checks the polynomial independently of the matching recurrence. Thus the
proposed next extension by one non-interval reversed arc already crosses a
real spectral boundary at order nine. This is a new explicit control for
this architecture, not a new first counterexample to general real-rootedness:
the repository's THM-025 already disproves that broader claim.

## 3. The entire induced-deletion deck

The interval map survives arbitrary induced deletion, but the equal-parity
count (1) does not. In particular, deleting the only interior vertex of a
triangle block leaves zero interiors, rather than 'half a cycle'.

Consider any retained vertex set S. A cycle interval [i,j] remains available
only if every boundary vertex b_i,...,b_j survives. Let m_[i,j] be the number
of surviving nonboundary interior vertices of those blocks. Its exact number
of surviving odd cycles is

    w^S_[i,j] = 2^(m_[i,j]-1),                 if m_[i,j]>=1;
                 1, if m_[i,j]=0 and j-i is even;
                 0, if m_[i,j]=0 and j-i is odd.           (7)

**Theorem 4 (all induced subsets).** Z_(IC(p)[S]) is the weighted independent-
interval polynomial with weights (7), omitting every interval whose boundary
list is not entirely retained. It factors over the maximal consecutive runs
of surviving boundary indices. Surviving interiors outside those runs lie on
no cycle and do not change Z.

This is immediate from Theorem 1, whose proof works after deletion. Distinct
runs have disjoint possible cycle supports, so their collection polynomials
multiply. Empty and singleton retained sets have polynomial 1.

Each run reduces to a generalized chain with m_i>=0 surviving interiors per
block. Equivalently its part lengths are m_i+1 and may now equal one. For a run
of j blocks, write q_i=2^m_i and d_i=1 when m_i=0, otherwise d_i=0. Formula (7)
can also be written

    w_[i,j] = (product_(t=i+1)^j q_t
                + product_(t=i+1)^j (-d_t))/2.             (8)

The two terms are respectively the unrestricted subset count and the parity
correction for an entirely empty interval. Weighted interval scheduling gives
an exact algorithm using O(k^2) polynomial additions and scalar
multiplications for each boundary run. The polynomial degrees are O(k); this
is not an O(k^2) scalar-arithmetic, bit-complexity, or all-cards-at-once claim.

### 3.1 Three scalar states and the exact tribonacci crosswalk

For H=Z(2) there is a linear-time, three-integer-state recurrence. Start with
(Z,C,D)=(1,0,0). A block with m surviving interiors updates it by

    q=2^m, d=1_(m=0),
    C'=q(Z-D),
    D'=-d(Z-C),
    Z'=Z+C'+D'.                                            (9)

To derive this, let A_j and B_j be the two endpoint-interval sums from (8),
weighted by the polynomial of the earlier compatible prefix. They obey
A_j=q_j(A_(j-1)+Z_(j-2)) and
B_j=-d_j(B_(j-1)+Z_(j-2)), while Z_j-Z_(j-1)=A_j+B_j at z=2.
Substitute Z_(j-2)=Z_(j-1)-A_(j-1)-B_(j-1); set C=A and D=B.
Signed C,D are auxiliary states, not asserted counts of graph objects.

There is, however, a physical rooted interpretation of a different basis:

    h=Z, e=Z-C, B=Z-D.                                    (10)

The root is the last boundary vertex. Here e counts Hamiltonian paths ending
at it, and B counts ordered spanning two-path covers whose second path is
nonempty and starts at it; the first path may be empty. These are exactly the
inherited Session 189 sole-exit response coordinates.

For m=0, (9) becomes

    (h,e,B) -> (h-e+B, h-e, h+B),                            (11)

the inherited sole-exit tribonacci tower. The zero-interior matrix in the
(Z,C,D) basis is

    [ 1  1 -1 ]
    [ 1  0 -1 ]
    [-1  1  0 ],

and satisfies M^3=M^2+M+I. The basis change (10) is integral and invertible.
For m>=1, (9) instead becomes

    (h,e,B) -> (h+2^m B, h, h+2^m B).                      (12)

This is the inherited transitive-padding/sole-exit near-join, now derived as
the second branch of the exact induced-cycle interface. It resets D to zero.
Thus Moon's product is the repeated-reset case; deleting triangle interiors
activates the missing third coordinate and produces tribonacci propagation.

For completeness, (11) has the usual direct path proof: a new sole-exit x->r
either ends a Hamiltonian path (h-e choices) or is followed by r (B choices).
Its new end count is h-e; its new B is h+B. For (12), the m transitive interiors
may independently lie in either path of a B-cover, multiplying B by 2^m; the
new endpoint can be entered after every old Hamiltonian path. This establishes
the named combinatorial interpretation without inferring it from matrix
similarity alone.

### 3.2 Marked-deletion polynomial formulas

These formulas strengthen the scalar chain-deletion identities from Session191
by retaining every odd-cycle-collection coefficient. Number blocks 1,...,k.

If v is interior to block j with p_j>=3, deleting it replaces q_j by q_j/2,
so its polynomial loss is

    Z_T(z)-Z_(T-v)(z)=Z(q_j;z)-Z(q_j/2;z).                  (13)

If p_j=2, allow the formal recurrence (3) to take q_j=1 and call its result
Z_formal. Then

    Z_(T-v)(z)=Z_formal(z)
        -(z/2) Z_(blocks 1,...,j-2)(z)
                 Z_(blocks j+2,...,k)(z).                 (14)

An empty block list means polynomial 1. Only the one-block interval at j has
an incorrect formal weight: 1/2 instead of zero. All other intervals containing
j also contain a surviving interior in a neighboring block, so their parity
counts remain correct. Collections using the erroneous interval cannot use
either neighboring block, giving precisely the factor in (14). This also
explains why the rational intermediate polynomial ends with integral output.
At z=2, (14) recovers the Session191 triangle-deletion correction, including
both endpoint cases. It is not a new scalar minimum theorem.

If v=b_j is a boundary vertex, 0<=j<=k, all intervals containing j disappear:

    Z_(T-v)(z)=Z_(blocks 1,...,j-1)(z)
                 Z_(blocks j+2,...,k)(z).                 (15)

These formulas and Theorem 4 avoid inferring a deletion card by simply erasing
a label from an aggregated undeleted state.

### 3.3 A decisive spectrum limitation of the hereditary chain approach

Every generalized chain with at least three vertices is strong. Indeed every
vertex reaches b_0 by going to its block end and following reversed boundary
arcs; b_0 reaches b_k directly when there are at least two blocks, or through
an interior when there is one block. Descending from b_k reaches every boundary,
and forward arcs reach all interiors. The only two-vertex case is a single
reversed arc.

The inherited Busch floor is 15 for strong order six and is nondecreasing
thereafter. Therefore a generalized chain with H=11 or H=13 would have order
at most five. All compositions of total part sum at most four give exactly
H in {1,3,5,9}, as direct enumeration or (11)--(12) shows. Thus neither 11 nor
13 occurs in any generalized chain, or in any induced subtournament of an
ordinary chain: Theorem 4 factors its H into generalized-chain values, and
11 and 13 are prime. Arbitrary ordered SCC products do not repair these gaps.

This is an all-order obstruction to this proposed architecture, not to the
full tournament spectrum: strong tournaments with H=11 and H=13 already exist.
It identifies where a new construction must leave the interval-cycle geometry.

## 4. A separate bridge from square-sum component ports

The square-sum packet maintained unordered pairs of endpoints of undirected
linear forests. The correct tournament version must retain each component's
ordered pair (start,end). Let F_T(P) count spanning directed linear forests
with unordered component set P={(s_i,t_i)}. A singleton v is (v,v). Every
selected vertex has in- and out-degree at most one, and directed cycles are
forbidden. All component counts are retained.

When a new vertex z arrives, each old state has exactly these possibilities:

1. Add (z,z).
2. Choose a component (s,t) with t->z; replace it by (s,z).
3. Choose a component (s,t) with z->s; replace it by (z,t).
4. Choose distinct components (s,t),(u,v) with t->z->u; replace them by (s,v).

Each choice contributes its old coefficient once. The two components in case4
are ordered by the direction of the join. Joining a component to itself would
close a cycle and is excluded. Removing z from a resulting forest recovers
its unique previous forest and its case, proving an exact closed transfer.
The one-component entries are the full named Hamiltonian endpoint matrix.

**A physical tournament collision.** The transitive tournaments with vertex
orders (0,1,2,3) and (0,2,1,3) have outgoing rows

    A=(14,12,8,0), B=(14,8,10,0).

Their complete named Hamiltonian endpoint matrices are equal: only entry(0,3)
is one. Adjoin the same named vertex4, beating exactly {0,1}. The children are

    A+=(14,12,24,16,3), H=11;
    B+=(14,8,26,16,3), H=9.

Both children are strong, since the new vertex beats the initial spine vertex
and loses to the terminal one. Thus even the full old Hamiltonian endpoint
matrix does not determine the count after one tournament arrival. The missing
coordinate is a directed two-component joining response, not a scalar endpoint
multiplicity. This does not claim that the two old strong-component profiles
are equal or that the parents themselves are strong.

**Sharp horizon.** With r future vertices, a current forest with more than r+1
components cannot finish as one Hamiltonian path; each arrival decreases its
component number by at most one. Discarding those levels is exact. Level r+1
is necessary even inside tournaments: take the final transitive order

    old_0,new_1,old_1,...,new_r,old_r.

Its unique Hamiltonian path induces r+1 singleton components on the old
vertices. Dropping that level loses its only final path. Hence a fixed
one-/two-component profile is not closed for arbitrarily many future arrivals.
The state space is finite for a fixed graph/horizon, but no polynomial-size
bound or minimal arbitrary nonlinear encoding is asserted.

## 5. Approach ledger and next use

1. **Hard-core intervals: succeeds exactly.** The explicit cycle map identifies
   the allowable conflict graphs of chain tournaments. It supplies a full
   polynomial and an all-induced-subset evaluator, not arbitrary realization
   of interval conflict graphs by tournaments.
2. **Weighted matchings: succeeds spectrally.** A standard path determinant,
   reached by the new exact specialization, proves all-order real-rootedness
   even in a clawed family. The signed substitution blocks an unqualified
   positive combinatorial bijection.
3. **Deletion as a state transition: succeeds, with a sharp boundary.** The
   zero-interior parity mode is exactly the old tribonacci response in a
   different basis. It explains why an induced deletion can destroy the
   two-state product closure. Its hereditary architecture still excludes
   attainable 11 and13, so it cannot establish the full {7,21} conjecture.
4. **Directed continuation ports: succeeds as a general exact interface.**
   It imports the useful square-sum idea while correcting orientation and
   singleton conventions. The explicit 11/9 collision rules out the tempting
   endpoint-only compression. Its unbounded component horizon prevents a
   fixed small-state universal grammar by this direct truncation alone.

A concrete next step is to retain one non-interval reversed arc as a named
interaction with the interval states. Section2.3 proves that its full
polynomial cannot in general retain the real-rooted matching representation.
The 11/13 exclusions make those values mandatory first tests for any expanded
grammar, rather than treating more chain products as spectrum progress. No
closed transfer for this expanded architecture is claimed here.

## 6. Verification and attribution

Run in this directory:

    python3 verify.py > RESULTS.json
    python3 -O verify.py > RESULTS.optimized.json
    cmp RESULTS.json RESULTS.optimized.json

The standalone standard-library verifier uses explicit exceptions, not Python
assertions. It independently compares direct cycle DFS, literal disjoint-cycle
collections, weighted interval recurrences, weighted matching subsets, directed
Hamiltonian subset DP, and edge-subset directed forests. RESULTS.json records
the exact finite universes and examples. No numerical root-finding is used to
prove Corollary3. Auxiliary exploratory scripts are not proof dependencies.

Sources and inherited inputs:

- `01-canon/theorems/THM-002-ocf.md`: Grinberg--Stanley OCF.
- THM-4744 and THM-4745: Moon's chains, known product law, Busch minimum,
  and exact chain gluings. Moon1972 and Busch2006 are credited there; no new
  claim to their scalar theorem is made.
- [Session 191 interval-chain growth](../../session191_tournament_growth_20261104/deletion_growth/INTERVAL_CHAINS.md): marked scalar chain deletion
  formulas and minimum within the chain family.
- [Session 191 rooted growth](../../session191_tournament_growth_20261104/deletion_growth/PROOF.md): exact count-neutral deletion
  criterion and rooted-increment growth scope.
- [Session 189 guide](../../../00-navigation/S189.md): inherited sole-exit three-state law
  and near-join response. Our exact map (10) identifies, rather than renames,
  that pre-existing mechanism.
- [Session 190 square-sum ports](../../session190_constructive_links_20261104/squareports/session190-squareports/PROOF.md): undirected component-pair transfer,
  deletion-deck principle, and sharp horizon inspiration.
- [Session 190 planar interfaces](../../session190_constructive_links_20261104/planar/session190-planar/PROOF.md): full pairings versus scalar boundary counts.
- THM-020 and THM-025: real-rootedness fails in general tournaments; the
  present all-order claim is confined to ordinary interval chains.
- THM-3186: a prior, unrelated signed continuant/weighted-path expansion;
  it motivates retaining the signed specialization rather than treating
  support as a positive matching bijection.

All new all-order statements have proofs above; the finite checks do not turn
any general-tournament conjecture into a theorem. A repository search checked
the listed nearby cycle/interval/real-rootedness work, but this is not a claim
of exhaustive literature priority.
