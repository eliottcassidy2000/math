# Parent-only unpointing projectors for tournament metagraph growth

24 October 2026. **PROVED:** the marked-extension, endpoint-fibre, unpointing,
and converse identities below, for every finite order. **FINITE-EXACT:** the
explicit universes and controls reported in `VALIDATION.md`. **OPEN:** a closed
recurrence from bare metagraph invariants or from a bounded compressed parent
type. No literature-priority or asymptotic speed claim is made.

## 1. Inheritance, target, and the missing coordinate

The inherited mechanisms are THM-781, *Hamiltonian-path inversion of merged
tournament-tiling fibres*, and THM-796, *Three-sorted recursive incidence of
tilings, complement lines, and converse-merged tournament nodes*. Their
byte-preserved texts and the October20 growth/enumeration notes are archived
as `.txt` with source paths and hashes in `source-pins.json`.

THM-796 already supplies exact word-cube growth, endpoint-face spans, their
failure of node lumpability, and the requirement to retain compatible
observer information. The October20 enumeration packet already generates
classes without Burnside and refutes specific unmarked/local summaries.
This packet does **not** claim these discoveries anew.

The next exact question is narrower: can class counts at order N=n+1 be
obtained from the complete order-n parent atlas without first classifying
order-N extensions? Yes. Retain the actual re-rooting transport, rather than
only the number of extensions. Its normalized matrix is a projector, and
its trace is the unmarked count. A refinement restricted to possible
Hamiltonian endpoints gives the literal fixed-observer cube quotient.

Live concepts:

1. Fixed-path words and their 2^(n-1) right-end lifts.
2. Incidence subsets modulo the parent automorphism action.
3. Root-change transport through order-n deletion overlaps.
4. Endpoint observer fibres and stabilizer weights.
5. Converse, which changes an endpoint into a starting point.

The anchor is exact A000568 transfer. The niche is recovery of child isotropy
from parent stabilizers. The wildcard is the idempotent/anti-idempotent
spectral description. These three use the same retained rooted incidence,
not the spectrum of the blue/black complement-line graph.

## 2. Parent data and marked extension states

Let I_n be a complete list of ordinary tournament isomorphism classes, with
one labelled representative T for each class. Write G_T=Aut(T), retaining its
permutation action, not just its order. For S subset V(T), adjoin a new vertex
star and define U(T,S) by

    v -> star exactly when v is in S.

Define the finite state set

    P_(n+1) = disjoint union over T in I_n of 2^V(T) / G_T.

A state p=(T,[S]) has stabilizer

    h_p = |(G_T)_S|.

**Lemma 1.** P_(n+1) is naturally in bijection with vertex-pointed tournament
isomorphism classes (U,v) of order n+1. Its stabilizer h_p is
|Aut(U)_v|.

**Proof.** Delete the marked vertex. An isomorphism of pointed extensions
restricts to an isomorphism of their deletions and preserves the subset of
vertices pointing into the mark. Conversely such a deletion isomorphism
extends by sending mark to mark. Once the deletion is identified with its
chosen representative T, changes of identification are exactly G_T. The
same argument for self-isomorphisms identifies the stabilizers. End proof.

The number |P_N| is a **vertex-pointed** count, not A000568(N). Dividing it
by N is also wrong when a child has a nontrivial automorphism orbit.

## 3. Constructing the transport without classifying children

For p=(T,[S]) choose U=U(T,S). For every vertex v of U:

1. Form the labelled deletion U-v, of order n.
2. Identify it with its order-n parent representative R by an isomorphism f.
3. Transport the subset {x : x -> v} by f, and take its Aut(R)-orbit.

The resulting pointed state is q=[(U,v)]. It is independent of f. Define

    D[p,q] = number of vertices v of U giving q.

This construction requires only order-n isomorphism tests and parent subset
orbits. It never asks for the isomorphism class of U as an unpointed order-N
object. The overlap isomorphisms are the exact missing cross-parent gluing
coordinate. Full child canonicalization appears only in the separate finite
audit, not in the prediction algorithm.

**Theorem 2 (unpointing projector).** With N=n+1,

    D^2 = N D,
    A000568(N) = tr(D)/N = rank(D).

The nonzero entries of every row are positive, its sum is N, and rows p and q
are identical whenever D[p,q]>0. Moreover

    |Aut(U)| = h_p D[p,p],
    D[p,q] h_q = D[q,p] h_p,
    D[p,q]/h_p = D[q,p]/h_q.

**Proof.** Fix an unpointed child U for the proof, not as algorithmic input.
Its states correspond to the vertex orbits O_1,...,O_k of Aut(U). Write
r_i=|O_i|. Re-rooting from any marked orbit to orbit O_j has exactly r_j
choices. Therefore its entire matrix block has all rows equal to
(r_1,...,r_k), and sum_j r_j=N. Distinct unpointed children have no transport
between them. Each block squares to N times itself, has rank one, and has
trace N. Summing the blocks proves the first formulas. Orbit-stabilizer gives
r_i=|Aut(U)|/h_i, proving the last formulas. End proof.

The two displayed balance identities are equivalent, not distinct claims.
The construction is a direct finite isomorphism/orbit transfer, **not** a
sum over permutations of fixed tournaments. Its proof uses elementary
orbit-stabilizer. It remains an exponential classification-strength method,
not a scalar recurrence or a newly asserted generic enumeration algorithm.

Only diagonal transport is needed for the number:

    A000568(N) = (1/N) sum_p #{v : (U_p,v) is isomorphic to (U_p,star)}.

These comparisons can all be made at order n by the deletion procedure.
The full D is useful for reproducing children and checking closure, but the
trace formula itself does not need a list of child classes or components.

## 3a. The nonnegative root-moving symmetry defect

Let p_N=|P_N| and define

    Delta_N = sum_p (D[p,p]-1).

Then

    A000568(N)=(p_N+Delta_N)/N,       Delta_N>=0.

The correction is exactly the information missed by treating every pointed
class as one of N equally distinct markings. It vanishes if and only if all
order-N children are asymmetric. Every D[p,p] is odd: a tournament
automorphism cannot be an involution, since a swapped pair would reverse its
own arc, so Cauchy's theorem makes the automorphism group odd. Its vertex
orbits are therefore odd. In particular Delta_N is even.

The parent-local part p_N is obtained by explicit subset orbits under G_T.
If desired it can instead be computed from each known parent's permutation
cycle action as sum_T (1/|G_T|) sum_(g in G_T) 2^(number of vertex cycles of g).
This optional small-group orbit formula counts only attachment subsets of an
already supplied parent; it is not a global fixed-tournament Burnside
prediction and is not used by the executable.

Only root-moving overlap tests are needed for Delta_N. They require the full
attachment alignment and parent deletion isomorphisms, not merely parent
counts, automorphism orders, or unmarked blue/black line data.

A computable partial-audit bound follows from a vertex partition of each
extension U_p whose **every cell is setwise invariant under every automorphism**.
Equivalently, it coarsens the vertex-orbit partition, as canonical invariant
color labels do. An equivariant unlabelled block partition alone is not enough:
automorphisms may permute its blocks. If b_p is the size of the cell containing
star, then 1<=D[p,p]<=b_p. Thus

    ceil(sum_p lower_p / N) <= A000568(N)
                              <= floor(sum_p upper_p / N),

where lower_p counts certified equivalent roots and upper_p counts roots not
yet excluded, initially 1 and b_p. A singleton invariant cell certifies zero
symmetry defect at that state without any re-root isomorphism test. Stable
directed degree refinement supplies the checked b_p in the finite engine.
This bound does not turn color refinement into a complete isomorphism test.

The finite Delta values at N=1,...,8 are

    0, 0, 2, 4, 12, 40, 152, 784.

For example p_8/8=6782 misses the true 6880 by Delta_8/8=98. This is a finite
exact decomposition, not a proved asymptotic estimate of the correction.

## 4. Converse and the all-root operator

Let c be the involution of P_N induced by (U,v) -> (U^op,v). In parent
coordinates it sends (T,[S]) to (T^op,[V(T)-S]), followed by the parent
identification. Let C be its permutation matrix, with (CD)[p,q]=D[c(p),q].
Converse commutes with re-rooting, so CD=DC. If SC(N) counts ordinary
self-converse child classes, then

    SC(N) = tr(CD)/N,
    merged(N) = (tr(D)+tr(CD))/(2N).

Indeed C permutes the rank-one child blocks. A self-converse block contributes
N to tr(CD); a pair of distinct converse blocks contributes zero. Counting
fixed **pointed** states of C alone would count a different object.

The all-root spectral form is especially simple. Put E=D/N and

    P_plus=(E+CE)/2,          P_minus=(E-CE)/2.

Then P_plus and P_minus are idempotent, their products in either order are
zero, and they sum to E. Their ranks are respectively

    merged(N),               (A000568(N)-SC(N))/2.

With the parent-known positive weights pi_p=1/h_p, detailed balance makes E
self-adjoint for the inner product sum_p pi_p x_p y_p. Converse preserves h_p,
so C and both P_plus/P_minus are also self-adjoint for this inner product.
Thus these are genuine mutually orthogonal weighted projectors, not merely
a formal difference of counts. They are not generally Euclidean-orthogonal
and they are not the original complement-line adjacency projectors.

## 5. The exact quotient of right-end word lifts

For a parent T define the endpoint path vector

    b_T(v) = number of directed Hamiltonian paths of T ending at v,
    L_T(S) = sum_(v in S) b_T(v).

For n=0 use the empty-parent convention L=1. An order-N Hamiltonian path
ending at star consists uniquely of a Hamiltonian path of T whose last
vertex lies in S, followed by star. Thus L_T(S) is exactly the number of
Hamiltonian paths of U(T,S) ending at star.

Call p **admissible** when L_T(S)>0, and denote this subset P_N^end. The
number of normalized fixed-path words lying over p is

    w_p = L_T(S)/h_p.                                      (1)

**Proof.** Identify the deletion of a marked endpoint with T. Normalized
words with this endpoint-pointed type correspond to Hamiltonian paths of U
ending at star modulo Aut(U)_star. The stabilizer acts freely on ordered
Hamiltonian paths, since a path names every vertex. Lemma1 identifies its
order as h_p; the path extension bijection gives L_T(S). End proof.

Consequently

    sum_(p in P_N^end) w_p = 2^binom(N-1,2),
    sum_p w_p = 2^(n-1) sum_(T in I_n) H(T)/|Aut(T)|  (n>=1).

The latter is precisely the known 2^(n-1) lift factor, now decomposed into
endpoint-pointed cells. Equation(1) preserves observer isotropy; it does not
replace it by an unweighted endpoint-marked count.

## 6. Endpoint unpointing and the correct converse operation

Let E(U) be the nonempty set of vertices that occur as the end of a directed
Hamiltonian path of U. Tournaments have Hamiltonian paths; the finite engine
also certifies positivity using the parent endpoint vectors. Write e(U)=|E(U)|.
For p,q in P_N^end set

    B[p,q] = D[p,q],
    R[p,q] = B[p,q]/e(U_p).

Equivalently, restrict D's columns to admissible states; its row sum is e(U).
No assertion that every vertex is a possible end is needed.

**Theorem 3 (native endpoint projector).**

    R^2=R,
    A000568(N)=tr(R)=rank(R).

For a fixed child U its admissible marked states are the automorphism orbits
in E(U). Every row is their common orbit-size distribution divided by e(U).
Thus each child again contributes one rank-one stochastic block and trace1.

Ordinary converse does **not** preserve P_N^end: a path end becomes a path
start. Already for a two-vertex transitive tournament, the sink is the unique
admissible end, but it becomes the source in the converse. Restricting the
pointed involution C to this state space is therefore invalid.

Instead define the anti-endpoint kernel by reversing U and then choosing a
possible end of the reversed tournament:

    K[p,q] = #{v in E(U_p^op) : [(U_p^op,v)]=q}/e(U_p^op).

It is computed as the admissible-column restriction of D[c(p),q], normalized
by that row sum; c(p) is allowed to be a nonadmissible full pointed state.
Then

    K^2=R,      RK=KR=K,
    SC(N)=tr(K),
    merged(N)=(tr(R)+tr(K))/2.

**Proof.** Both rows depend only on the underlying child class. R forgets the
mark and chooses a uniformly random Hamiltonian-eligible endpoint. K does
the same on the converse class. Applying K twice reverses twice and leaves
the same endpoint distribution as R. Composing with R on either side cannot
change the underlying class or that distribution. A self-converse child has
K's diagonal block equal to its R block, of trace1; a non-self-converse pair
has no diagonal K block. End proof.

On im(R), K is the ordinary converse permutation of unpointed children. Its
+1 eigenspace has dimension merged(N), and its -1 eigenspace has dimension
(A000568(N)-SC(N))/2. All remaining eigenvalues are zero. This is a spectrum
of the newly defined marked-incidence operator, **not** a claim that the
blue/black metagraph adjacency or Laplacian has these eigenvalues.

## 7. What is sufficient, and what remains open

Source: full parent native incidence presentations, parent automorphism
permutation actions, endpoint path vectors, and deletion-overlap transport.
Target: P_N, its admissible endpoint subset, the kernels D/R/K, and their
traces. Preserved: all child isomorphism classes, converse, pointed isotropy,
Hamiltonian observer multiplicities, and exact unpointing. Lost in scalar
traces: child identity, internal observer words, and blue/black partner
incidence. The complete source can be reconstructed from a full normalized
word presentation, but not from an unmarked colored complement-line graph.

The sufficient counting signature for an already formed endpoint extension
is its root-return multiplicity, possible-endpoint count, opposite root-return
multiplicity, and opposite endpoint count. These are computable from parent
incidence; they are **not** proved to be functions of a small parent scalar
signature. The full transporter solves gluing without pre-counting children,
but has not compressed away the hard isomorphism information.

The precise next target is a smaller parent type that predicts these diagonal
return and endpoint-support data, with a proved all-order closure law. The
October20 hostile examples remain binding: unmarked graph symmetries and
several local marks fail to determine extension profiles. This packet neither
reverses those failures nor proves that every possible compression is
impossible.

## 8. Minimal hostile controls

At N=3, the four full pointed states split into a three-state transitive
block and a one-state cyclic block. The diagonal of D is (1,1,1,3), up to
ordering. Its trace6 gives two classes; |P_3|/3=4/3 is wrong. Parent subset
stabilizers are all1, so the cyclic child's automorphism group of order3 is
created by re-root transport, not inherited as a parent automorphism.

At N=4, there are twelve vertex-pointed states, seven endpoint-pointed states,
eight fixed-path words, four ordinary classes, two self-converse classes, and
three converse-merged classes. None of these six counts may be identified
with another merely because each is a quotient of presentations.

At N=2, the pointed converse does not preserve admissible end states, as
explained above. The corrected K has trace1 and square R. The checker retains
all three controls together with complete finite class/fibre/isotropy audits.


## 9. A minimal finite failure of refinement-only root transport

Stable directed degree refinement supplies an inexpensive upper bound for a
root orbit, but it is not always exact. All pointed states through N=6 have
exact bounds. At N=7, three of the 3,040 pointed states fail, and the witness

    adjacency rows = (74,84,97,22,37,11,56)

is a regular tournament: every outdegree is3, so ordinary directed color
refinement retains one cell of size7 forever. Nevertheless its automorphism
group has order3, with vertex orbits of sizes 3,3,1; vertex6 is fixed. A
standalone all-7! checker, `refinement_hostile.py`, verifies the actual group
and orbits without importing the transfer implementation.

Using cell size in place of actual orbit size therefore contributes 21
instead of7 from this child's three pointed states. The global candidate
trace would be3206 rather than3192, wrongly giving458 instead of456 classes.
The strongest survivor is a rigorous upper bound and a way to locate the
remaining root-overlap work, not a closed exact refinement recurrence.

At N=8 the certified refinement upper bound is55068 while the exact trace
is55040; six of54,256 states have a strict bound. A singleton cell settles
53,882 states immediately. These are finite compression diagnostics, not
claims that the exceptional list or cell patterns remain bounded at larger
orders.


### Infinite family by an inherited operation

The archived October20 growth proposition proves that uniform transitive
inflation D_k(T)=T[K_k] preserves the automorphism group and multiplies the
canonical maximal transitive-module lengths by k. Here K_k is a forward
chain, not a complete undirected graph or Hadamard doubling.

Apply this inherited operation to the regular-seven witness. Label a clone
by (v,i), where 0<=i<k. Its outdegree is

    3k+(k-1-i).

Thus the first degree partition has k cells C_i={(v,i):v in V(T)}, each of
size7. From a vertex in C_i, the number of outgoing neighbors in C_j is

    3 + 1_(i<j).

These numbers depend only on i,j, not v, so no further directed color
refinement splits a cell. But Aut(D_k(T))=Aut(T)=C3, acting simultaneously on
base blocks and fixing each chain position. The k clones of base vertex6
are therefore fixed individually. Their true root orbit has size1 while
the invariant refinement cell has size7, for **every integer k>=1**.

The base tournament has no homogeneous pair: such a pair would have
outdegrees differing by1, contradicting regularity. Hence its canonical
transitive modules are singletons, and the inherited module theorem applies
without an ambiguous alternative clone decomposition. The standalone hostile
checker also verifies the clone-degree cells and homogeneous-pair components
at k=1,...,6. Those finite checks are controls; the displayed formulas and
inherited automorphism theorem prove the infinite family.

This refutes exact unpointing by stable directed degree-refinement cell sizes
at arbitrarily large orders. It does not refute stronger refinements that
retain the missing root-overlap isomorphism data.
