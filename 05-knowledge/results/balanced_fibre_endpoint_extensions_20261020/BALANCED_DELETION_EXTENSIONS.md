# Boundary fibres, exact endpoint transport, and the 625/576 square

2026-10-20. Session17. Research stage only; no canonical theorem ID reserved.

**PROVED:** arbitrary interior thinning in the two-return family; general
product-fibre correction; closed endpoint-growth kernel; full endpoint
unimodularity; minimal homogeneous-twin orphan; ordered-join realization of
625 and 576; marked equivariant five-field path coordinates.
**FINITE-EXACT:** the independently implemented checks below.
**NOT CLAIMED:** arbitrary tournament deletion preserves Hamiltonian paths;
a scalar or endpoint state closes arbitrary tournament growth; determinant
adjacency or a field action is already a tournament-metagraph move.

## 0. Inheritance and changed scope

The repository was read at f3e658d7cc10, with incoming port packet first
published in ed0719d0ca3. The following are **inherited**, not new results:

- `05-knowledge/results/strong_blue_games_20261020f/double_return_counts.md`:
  two-return product decoder, square/blue-diagonal count, rigidity.
- `05-knowledge/results/zero_heavy_ports_20261020g/microcosm.md`, §§2–5:
  pointed Boolean product, commuting involutions, actual one-step deletion,
  4/2/1 fibres, scalar corrections, and subset-sum endpoint formula.
- The same packet's `profile_growth.md`, §4, and repaired
  `01-canon/theorems/THM-1900-the-insertion-response-calculus.md`:
  ordinary vertex deletion can leave an orphan, so an HP map must be proved.

The anchor is operation-preserving deletion. The niche is endpoint transport,
including its exceptional states. The wildcard is a literal square-product
carrier for 625/576 and a marked field coordinate system. The concept board
is: actual path words; pointed product; endpoint ports; blue involution;
loop-collapse of choice graphs; field marking. The method used is controlled
forgetting with a predicate-specific sidecar, not numerical pattern fitting.

A relevant older arithmetic route is
`05-knowledge/results/predecessor_prime_shapes_20261014b.md`, §1:
F counts proper nontrivial divisors, U proper prime divisors, S squarefree
proper nontrivial divisors. U's underlying set is contained in S's; F=U+S
is not itself a disjoint-set partition. None of the path-set bijections below
identifies these divisor sets without an additional explicit map.

## 1. Actual arbitrary interior deletion, not only the central pair

Let T_(a,b), a,b>=2, have vertices 0,...,a+b, all increasing arcs except
that a->0 and (a+b)->a replace the two corresponding increasing arcs.
Put A={0,...,a-1}, c=a, B={a+1,...,a+b}. Inherited unique path coordinates are

    I increasing, J increasing, c, (A minus I) increasing,
    (B minus J) increasing,

where I is any subset of A minus {0}, or I=A; and J is empty or contains
n=a+b. Represent I=A by infinity, and represent J=empty by a second
infinity. Otherwise use the finite left subset I and finite right subset
K=B minus J. Thus the path set is a product of two pointed Boolean sets.

**Theorem 1 (rectangular thinning).** Delete arbitrary r vertices from
{1,...,a-1} and s from {a+1,...,a+b-1}, with a-r,b-s>=2. Increasing
relabeling gives T_(a-r,b-s). Ordinary deletion from every source path is
an HP of the target, and the map is onto. A target state with finite/finite,
finite/infinity, infinity/finite, infinity/infinity coordinates has fibre
size respectively

    2^(r+s), 2^r, 2^s, 1.                              (1)

**Proof.** Deletion acts separately on the four increasing runs, leaves the
centre and both return endpoints in place, and restricts I and K to the
surviving interior labels. Infinity stays infinity, because its exceptional
endpoint remains. The restricted state satisfies exactly the inherited
admissibility rules in the smaller tournament. Conversely every deleted
interior bit is chosen independently when its coordinate is finite; an
infinity coordinate has exactly one possible extension. The unique path
decoder proves both surjectivity and the listed fibres. No joining arc is
assumed without checking the resulting admissible decoder. QED.

For the balanced family T_k=T_(k,k), symmetry commutes with thinning exactly
when the deleted right labels are reflections of the deleted left labels.
Arbitrary unpaired thinning still preserves the path decoder, but need not
preserve the balanced blue predicate.

A useful marking is already visible at the endpoints: left infinity iff
start=0, and right infinity iff end=a+b. Consequently (1) is an endpoint-
marked fibre law, although endpoints do not recover the complete path.

For balanced statements below, use the common-coordinate convention
`K = sigma(B minus J)` when J is nonempty, with sigma(v)=2k-v;
when J is empty, K=infinity. This reflects the raw right-shore subset
from Theorem 1, so both I and K belong to
X_k=P({1,...,k-1}) disjoint-union {infinity}.

### 1.1 All-depth correction and its exact black contribution

Write q=2^(k-1)+1 and D=2^d for d simultaneous reflected interior additions.
Each finite coordinate has D lifts; infinity has one. Therefore

    q'=D(q-1)+1,
    H'=D²H-2D(D-1)q+(D-1)²,
    Delta'=D²Delta-D(D-1)(q-1),                          (2)

where H=q² and Delta=q(q-1). These formulas remain set-theoretically valid
for any pointed q-element set under such a fibre expansion, without asserting
that arbitrary q arises from this particular T_k construction.

The proof of the correction is a disjoint classification, not interpolation:

- (q-1)(q-2) finite unequal parent pairs each have D² black lifts.
- 2(q-1) mixed parent pairs each have D black lifts, all black.
- q-1 finite diagonal parents each have D² lifts, of which D are blue and
  D(D-1) are black.
- The exceptional diagonal parent has one blue lift.

Summing gives (2). In particular a black source path can become blue after
deletion: it may have distinct new bits above equal old finite coordinates.
The valid direction is that a black *parent* has only black *lifts*.

## 2. General product fibres and the precise simple-graph boundary

For any surjection f:X'->X between finite sets with nonempty fibres of sizes
w_x, its product
f×f has w_x w_y lifts over (x,y). Coordinate swap is respected. A diagonal
parent (x,x) has w_x² lifts: w_x diagonal, w_x(w_x-1) off diagonal.
An off-diagonal parent has only off-diagonal lifts. Hence

    |X'|(|X'|-1)
      = sum_(x!=y) w_x w_y + sum_x w_x(w_x-1).            (3)

This is the abstract source of every black correction in (2). The correction
cannot be read off from the number of parent points alone when the w_x vary.

Modulo swap, off-diagonal pairs are the edges of K_X. The map to K_X is
**not** a graph homomorphism from the loopless K_(X') when some w_x>1:
an edge inside a fibre collapses to a vertex. It is a valid map only if
loops are retained, or if collapsed edges are recorded separately. The
exact decomposition is between-fibre edges plus within-fibre edges. No
orientation of these choice-state edges is intrinsic, so calling their
complete graph a tournament would add an arbitrary gauge.

## 3. The entire endpoint matrix closes under central growth

Let E_k(s,t) count HPs of T_k beginning at s and ending at t. The full path
is not determined by (s,t), and blue status is not determined by (s,t).
Nevertheless, for this **particular physical growth**, (s,t) determines
the complete multiset of lift endpoints.

Embed old labels by e(v)=v for v<k, e(k)=k+1, e(v)=v+2 for v>k.
The two added vertices are L=k and R=k+2. From an old endpoint pair read

    I=infinity iff s=0;       I=empty iff s>=k;
    K=infinity iff t=2k;      K=empty iff t<=k.            (4)

All other states are finite nonempty. Choose epsilon in {0,1} for a finite
I, or a single formal value * for infinity; choose eta similarly for K.
These are the new subset bits. Initially set s'=e(s), t'=e(t), then apply

    if I=empty:
      epsilon=1                         => s'=L;
      epsilon=0 and eta=0               => s'=R;
      otherwise                        => s'=e(s).

    if K=empty:
      eta=1                             => t'=R;
      eta=0 and epsilon=0              => t'=L;
      otherwise                        => t'=e(t).        (5)

Each allowed bit pair contributes one lift. Multiple pairs yielding the
same endpoint pair count with multiplicity. Let K_k(s,t;s',t') be this
integer kernel. Then exactly

    E_(k+1)(s',t') = sum_(s,t) E_k(s,t) K_k(s,t;s',t').   (6)

**Proof.** A nonempty I keeps its old minimum, so its start is unchanged
under e. If I was empty, adding L makes L the start; otherwise R becomes
the start precisely when it lies in the new J, equivalently eta=0. The
other cases retain the first surviving old vertex. The end argument is
symmetric with K and eta. Reading (4) from the decoder proves that no
hidden information about the parent subset is needed for (5). Every lift
is specified once by its bits. QED.

For a separately retained blue endpoint matrix, use the same kernel but
keep only epsilon=eta, including the sole (*,*) lift. This propagates blue
counts and hence black endpoint counts by subtraction. Applying that
restriction to the uncoloured E_k would be invalid: the same endpoints
can contain both blue and black paths.

### 3.1 A small state for coarse boundary strata, and formal series

The one-coordinate types (infinity, empty, nonempty) have transfer

    M = [[1,0,0], [0,1,0], [0,1,2]],

acting on a column vector of type counts. The pair-type transfer is M⊗M.
This nine-state closure records boundary strata, not labelled endpoints.
The full marked endpoint recurrence (6) has growing label sets.

With u marking the number of infinity coordinates, the stratum enumerator is

    Phi_k(u)=(2^(k-1)+u)²,
    Phi_(k+1)(u)=4 Phi_k(u/2).

The right side is a polynomial with integer coefficients: finite/finite,
exactly-one-infinity, double-infinity counts scale by 4,2,1 respectively.
The bulk, mixed boundary and exceptional series, with z marking k-2, are

    sum_(k>=2) (q_k-1)² z^(k-2) = 4/(1-4z),
    sum_(k>=2) 2(q_k-1) z^(k-2) = 4/(1-2z),
    sum_(k>=2) 1*z^(k-2)        = 1/(1-z).               (7)

These rational functions follow from the proved strata. They do not infer
a recurrence from a finite sequence. Their sum recovers the inherited
H generating function after converting from k to tournament order 2k+1.

## 4. A rank-one bulk hides a unimodular full port matrix

The incoming endpoint formula contains a rank-one k×k block plus subset
sums. Its subset sums have the following explicit coefficients. For
1<=alpha<=k-1 and 0<=beta<=k-1, put

    c(alpha,beta)=1                     if beta=alpha-1;
                  2^(beta-alpha-1)     if beta>alpha;
                  0                    otherwise.        (8)

This counts finite nonempty subsets S with min S=alpha and
max({0,...,k-1} minus S)=beta. In the first case S is the forced terminal
interval; in the second case only alpha+1,...,beta-1 are free. The case
beta=alpha is impossible, and beta<alpha-1 omits a forced complement point.

Set a_0=1 and a_i=2^(k-1-i), 1<=i<k. Then the complete endpoint matrix is
specified by disjoint supports:

    E_k(i,2k-j)=a_i a_j                    (0<=i,j<k);
    E_k(alpha,beta)=c(alpha,beta);
    E_k(2k-beta,2k-alpha)=c(alpha,beta);
    E_k(0,k)=E_k(k,2k)=E_k(k+1,k-1)=1.       (9)

Everything else is zero. In particular its support size is

    k² + 2*sum_(alpha=1)^(k-1)(k-alpha) + 3
      = 2k²-k+3.                                         (10)

**Theorem 2 (full port unimodularity).** For every k>=2,

    det E_k = permanent E_k = 1.                          (11)

Thus E_k has full rank 2k+1 over every field and is invertible over Z,
despite its large rank-one bulk.

**Proof.** Its support has exactly one perfect matching. Column 0 forces
row 1. After that row is removed, column 1 forces row 2, and successively
columns 0,...,k-2 force rows 1,...,k-1. All other entries in each such column
occur in previously removed rows, by (8). Column k-1 now forces row k+1.
Column k forces row 0. In the reflected block, column k+1 now forces row
k+2, and successively columns k+1,...,2k-1 force rows k+2,...,2k. The last
column 2k forces row k. Every forced entry equals one.

The resulting row-to-column permutation is the single cycle

    0 -> k -> 2k -> 2k-1 -> ... -> k+1
      -> k-1 -> k-2 -> ... -> 1 -> 0.

Its length 2k+1 is odd, so its sign is positive. The permanent and determinant
therefore each have a single term, equal to one. QED.

This identifies a strong failure boundary for the rank-one approximation:
the small boundary pieces retain all full-matrix rank. It does not say
that this endpoint matrix replaces the distinct-two-component tensor B_1
needed for arbitrary one-vertex insertion.

## 5. Homogeneous twins do not make arbitrary HP deletion safe

There is a minimal five-vertex counterexample. Label vertices u,v,a,b,c,
with u->v; both a,c point to both u,v; both u,v point to b; and b->a,
b->c, a->c. Thus u,v are a genuine homogeneous transitive pair.
The Hamiltonian path

    (a,v,b,c,u)

deletes to (a,b,c,u), which is not a path because b->a. In the numerical
labels (u,v,a,b,c)=(0,1,2,3,4), outgoing bit rows are

    (10,8,19,20,3).

This is minimal. An orphan after deleting an internal twin v requires
consecutive a,v,b with a->v->b->a. The other twin u satisfies a->u->b.
On at most four vertices, u would have to occur before a or after b in
that same length-four word, respectively requiring u->a or b->u, both
false. Endpoint deletions cannot create an orphan. Hence no smaller twin
example exists. The five-vertex witness shows why Theorem 1 uses its
complete product decoder, not a generic appeal to cloned vertices.

## 6. The user's 625 and 576 have a literal tournament product model

For any nonempty tournaments A,B, the ordered join A->B has all cross arcs
pointing from A to B. Every HP visits all of A before all of B, so

    HP(A->B) = HP(A) × HP(B)                              (12)

under concatenation, with a literal inverse given by restriction.

Let S=T_(3,3), with seven vertices and 25 HPs. Let P_0 be its standard
increasing HP. Then S->S has 625 HPs. Its paths whose two component paths
are both different from P_0 form exactly

    (HP(S) minus {P_0})²,

of cardinality (25-1)²=576. This is inclusion-exclusion removing a row and
a column of the 25×25 choice square:

    625 - 25 - 25 + 1 = 576.                             (13)

It is **not** deletion of the diagonal. The diagonal has 25 points, not 49.

To retain the blue predicate precisely, S is strong, rigid and selfconverse.
Its unique anti-automorphism is sigma(i)=6-i; define r(P)=sigma(reverse P).
The ordered join has exactly two strongly connected components, so any
automorphism preserves each and is trivial; any anti-automorphism swaps
them and uses sigma in each. On its two path choices the blue involution
is

    (P,Q) -> (r(Q),r(P)).

Conjugating the second coordinate by r makes this ordinary coordinate swap.
There are 25 blue paths and 600 black paths. Since r fixes P_0, the 576
bulk contains 24 blue and 552 black paths. The 49-point boundary contains
one blue and 48 black paths. This is an actual path/colour model, not a
claim that the number 576 alone determines a tournament.

Notice the parameter distinction: 25 is not 2^(k-1)+1 for an integer k.
It is H(T_3)=5². Thus 625 is not H(T_k) for any member of the original
balanced two-return family, but is H(T_3->T_3).

The arithmetic 625+125+25+1=776 is correct. No operation in this note
identifies its four summands with a naturally selected four-stratum HP set.
Such an interpretation needs its own carrier and disjoint decomposition;
the 625/576 map does not supply that missing assertion.

## 7. Marked F5 coordinates for S=T_(3,3)

The one-coordinate choice set is

    X_3={infinity, empty, {1}, {2}, {1,2}}.

Its complement involution kappa fixes infinity and swaps empty with {1,2},
and {1} with {2}. Choose a bijection ell:X_3->F5 by

    ell(infinity)=0, ell(empty)=1, ell({1,2})=4,
    ell({1})=2, ell({2})=3.                              (14)

Then ell(kappa(x))=-ell(x). Through the inherited actual path decoder,
(ell(I),ell(K)) bijects HP(S) with F5². Its exact transported symmetries are

    J=(kappa,kappa)       -> v -> -v;
    actual anti-reversal -> (x,y) -> (y,x).

The standard path maps to zero, so the 24 nonstandard paths modulo J give
the twelve classes (F5² minus {0})/{±1}. Two nonstandard component paths
in (13), independently modulo J, give 12²=144 ordered projective-sign
pairs, each with exactly four representatives in the 576 bulk.

This is a marked coordinate bridge. There are eight possible one-coordinate
bijections intertwining kappa with negation: permute the two nonzero sign
orbits and choose a sign in each. Most do not preserve the additive or
multiplicative structure of a previously chosen field labeling. Therefore
any determinant graph subsequently defined from (14) is attached to this
marking; it is not a canonical consequence of cardinality five.

Every permutation of the 25-point field set transports to a permutation
of HP(S), but S is rigid, so nontrivial field actions are not automorphisms
of S. Likewise J is a path-set permutation, not a tournament automorphism.
A determinant-defined edge relation may define a useful **new graph on
paths**; identifying it with already defined black metagraph moves requires
checking those moves and their incidence separately. The exact file
`marked_f5_paths.json` retains all 25 actual paths and coordinates.

## 8. Exact controls, scope, and provenance

Reproduce in this stage directory:

    python3 check_fibres.py > check.out
    python3 -O check_fibres.py > check.optimized.out
    cmp check.out check.optimized.out

The checker imports no repository producer code. It uses literal HP DFS,
induced-subtournament checks, exact subset decoding, direct positional-blue
checks, integer Bareiss determinants, and an independent permanent DP.
All validation failures are explicit exceptions, never removable assertions.

Frozen universe:

- All a,b=2,...,5 and every interior deletion leaving both parameters >=2:
  676 actual deletion maps, 112,896 source-path deletions, every target fibre.
- Balanced k=2,...,6: direct DFS versus decoder, every parent's endpoint
  lift multiset, aggregate full and blue endpoint matrices through k+1.
- Balanced k=2,...,8: all endpoint coefficients, support cardinality,
  determinant and permanent each equal one.
- The fourteen-vertex ordered join S->S: literal HP set equals all 625
  concatenations, direct blue predicate gives the eight counts in §6.
- Every labelled tournament on two through four vertices, every homogeneous
  pair and each twin deletion: 412 valid controls; explicit five-vertex hostile.
- 693 integer controls on (2); all 25 marked F5 paths and transported R,J.

The all-order proofs are independent of these bounded controls. The
finite-field determinant-graph and group analysis belongs to the parallel
prime-shape lane, not this packet. Literature priority and unrestricted
recurrence claims are not made. Parent publication owns repository writes.
