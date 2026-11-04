# Two-port algebra: converse separation and a strict triangle direction

**Status: PROVED statements A1--A3; FINITE-EXACT controls stated separately.**
These results do not prove universal two-port reconstruction. No historical
novelty or general convexity of the Hamiltonian polynomial is asserted.

## Inheritance and objects

[Session27, sections7--8](../../session27_repunits_lattices_20261022/session27-response/proof.md)
already proves that a fixed battery of at-most-two-port Hamilton counts
recovers H, both endpoint vectors S,E, and the ordered two-path-cover matrix
C. Specifically C_uv counts spanning ordered covers by two nonempty paths,
the first ending at u and the second starting at v. It also identifies C as
the gradient, at the tournament adjacency matrix, of

    P_n(X) = sum_(permutations pi) product_(i=1)^(n-1) X_(pi_i,pi_(i+1)).

All core vertices are named, and C_uu=0. Reversing every path and exchanging
the two ordered pieces gives `C(T^op)=C(T)^t`. The matrix gradient identity
is inherited; its global injectivity is not presumed.

The closest positive mechanism is
[C4, positive circulation of one-port collisions](../../proof_principles_20261103/continuation/PROOF.md).
The canonical hostile is the transitive three-vertex tournament: its
antisymmetric current can point against an actual arc, while some actual
arcs have zero current. The least-used sidecar below is the mod-four
divisibility of the two-colour path polynomial. The live board consists of
gradient data, path-cover mass, converse, deletion/ports, and realizable
directions in the tournament adjacency polytope.

## A1. Every odd-order tournament is separated from its labelled converse

**Theorem.** For every tournament T of odd order n>=3, C(T) is not symmetric.
Consequently the at-most-two-port response of T differs from that of T^op.
This holds for all tournaments, without regularity or transitivity assumptions.

Let q_2 be the number of ordered spanning two-path covers. Then

    q_2 = sum_(u!=v) C_uv.                                      (1)

Every cover has one first endpoint and one second start, so contributes once.
Every Hamiltonian path has exactly n-1 arcs, so marking one of its arcs gives

    sum_(u->v in T) C_uv = (n-1)H(T).                           (2)

This is also Euler's homogeneous-polynomial identity at the adjacency point.
If C were symmetric, choosing one orientation from each unordered pair in
(2) would give

    q_2 = 2(n-1)H(T).                                          (3)

Now let Q_T(2) count vertex colourings by two labelled colours with one
Hamiltonian path chosen in each nonempty colour class. Empty classes have
count one. Directly partitioning by the first colour class gives

    Q_T(2)=sum_(S subset V) H(T[S]) H(T[V\S]) = 2H(T)+q_2.      (4)

Every induced tournament has an odd Hamiltonian count by classical Redei
parity, and H(empty)=1. Pair each subset with its complement. Each pair
contributes twice an odd integer, and there are 2^(n-1) pairs. Thus

    Q_T(2) = 0 mod4,                 for every n>=2.             (5)

Equations(3)--(4) would give `Q_T(2)=2nH(T)`, which is2 modulo4 for odd n.
This contradiction proves C is asymmetric. Since converse transposes C and
the two-port battery recovers C, the response functions differ. QED.

The n>=3 odd restriction is genuine to this proof. At n=1 the converse is
identical. At n=2 the cover matrix is symmetric, although the endpoints and
one-port response already distinguish the two labelled orientations. No
all-even-order C-asymmetry theorem is claimed here.

## A2. Reversing a directed triangle strictly changes the gradient

**Theorem.** Let T contain a directed triangle on labels R={a,b,c}, and let
T' reverse precisely its three internal arcs. Put

    A(t)=(1-t)A_T+tA_(T'),     f(t)=P_n(A(t)).

Then f is a quadratic polynomial and

    f''(t) >= 12 H(T[V\R]) > 0,                              (6)

where H(empty)=1. In particular C(T') differs from C(T), so the quadratic
size Session27 two-port battery distinguishes them even in an arbitrary
ambient tournament. No module hypothesis is needed.

**Proof.** A Hamiltonian word uses at most two edges of the three-vertex
triangle. If it uses two, they form a consecutive path through all three
vertices. Relative to the original directed cycle, those two edges are
either both forward or both reversed; mixed choices would enter or leave
the middle vertex twice. Discard words having an illegal edge outside R.
The remaining terms using two triangle edges contribute

    M_+(1-t)^2 + M_- t^2,

where M_+ and M_- are their nonnegative integer counts. All other surviving
terms have degree at most one. Hence

    f''=2(M_++M_-).                                            (7)

For completeness, each directed Hamiltonian word on the complement admits
at least three insertions of a consecutive forward triangle path. The
same holds for the reversed triangle. Here is the endpoint calculation,
which avoids assuming a chosen rotation always inserts.

Write the complement word as v_1,...,v_r. Let d_j(i) be1 when v_j points to
triangle label i and0 otherwise, and let z_j=sum_i d_j(i). A path obtained
by cutting the directed triangle cycle starts at i and ends at its cyclic
predecessor rho(i). The three possible rotations give:

- 3-z_1 legal initial insertions;
- z_r legal final insertions;
- at internal gap j, exactly
  `sum_i d_j(i)(1-d_(j+1)(rho(i))) >= z_j-z_(j+1)` legal insertions.

Summing telescopes to at least3. If the complement is empty there are
exactly three rotations directly. Distinct complement words and distinct
insertions give distinct Hamiltonian words, since deleting the triangle
block recovers its complement word and location. Applying the same count
to the reversed cycle gives

    M_+ >= 3H(T[V\R]),     M_- >= 3H(T[V\R]).                  (8)

Combining(7)--(8) proves(6). The path-existence input for the complement is
the usual tournament Hamiltonian-path theorem. The constant is attained
when the triangle is an ordered-sum factor, for example T=C3 plus a
transitive singleton suffix: both counts equal three times the complement's
Hamiltonian count.

Finally, the chain rule gives

    f'(0)=<C(T), A_(T')-A_T>,
    f'(1)=<C(T'),A_(T')-A_T>.

The strict increase f'(1)-f'(0)=f'' excludes equal gradients. QED.

This is a quantitative separation theorem for this specific cyclic
three-edge reversal. It does not say the sign of a single C_uv-C_vu reveals
the arc, nor that the full gradient map is injective on all tournaments.

## A3. The parity theorem sharpens deletion-based port recovery

The independent deletion argument is as follows. Split an existing external
port into two consecutive ports x,y with x->y, giving them identical
relations to every remaining core vertex and to the other ports. For a
chosen deleted core vertex v, vary just the two edges along x->v->y and
take the signed sum over their four orientations. Every Hamiltonian word
omitting either edge cancels. A surviving word contains the consecutive
block x-v-y or y-v-x; both have character sign +1. Collapsing either block
gives the original one-port response on the core with v removed. Hence
the signed sum is twice that response. The actual x->y edge is unused in
those blocks, so its fixed orientation imposes no additional constraint.

The same argument applies with other ports present because x,y inherit
identical relations from the split port. Iterating s times uses s extra
ports and4^s terms. Therefore equality of response functions through m
ports implies equality of one-port responses on every identically labelled
deletion T-S with |S|<=m-1.

Let D be the oriented disagreement graph of two such cores, oriented as
in the first core. C4's positive weighted-circulation theorem applies to
every D-S. Every active vertex of D must have at least m outneighbors and
m inneighbors in D. To see this, first use balance on D to ensure both
kinds exist. If it has fewer than m outneighbors, delete all of them. No
inneighbor is deleted, because D is an oriented graph without two-cycles.
The remaining positive incoming edge cannot be balanced by zero outgoing
edges. The reverse argument gives the other degree bound. In particular
any nonempty disagreement has at least2m+1 active vertices.

For even n=2m this already proves m-port injectivity. For odd n=2m+1 with
m>=2, any residual disagreement must use every vertex, have every possible
edge, and have m incoming and m outgoing arcs at every vertex. Thus the
second core must be the labelled converse of the first. But m>=2 gives
access to C, and A1 rules out that remaining collision. Consequently:

    floor(n/2) ports suffice for labelled reconstruction, n>=4. (9)

For n=2 one port suffices. For n=3 two ports suffice, whereas one port
cannot distinguish opposite directed triangles. At n=1 no observation is
needed. These are sufficient port counts; no all-order optimality claim
is made. In particular the argument proves universal two-port faithfulness
through order5, but leaves orders>=6 unresolved by this argument.
The inherited exhaustive order-six computation is separate evidence.

A finite query upper bound for obtaining the deletion tables is

    sum_(s=0)^(m-1) binom(n,s) 4^s [1+(n-s)+binom(n-s,2)].      (10)

For each deletion set the bracket is the constant/singleton/pair battery
recovering the quadratic one-port response, and each entry expands into
at most4^s queries using at most m ports. This alone suffices when n<=2m.
For the improvement at odd n=2m+1>=5, append Session27's C-tomography
battery of `1+2n+3binom(n,2)` entries, each using at most two ports. Thus
an explicit sufficient reconstruction battery is bounded by(10) plus
that additional quantity in the odd boundary case. No recovery of C from
the deletion tables alone is asserted. Repeated entries are permitted;
no deduplication or query-optimality assertion is made. Exact injection
allows a finite enumeration decoder over all labelled cores of the known
order; this is not a polynomial-time reconstruction algorithm.

## Why a general variational shortcut still fails

The Hamiltonian polynomial is neither convex nor concave on the affine
tournament adjacency polytope. Along global reversal of the transitive
three-vertex tournament,

    P((1-t)A+tA^t)=1+2t-2t^2.

Along global reversal of the directed triangle the value is

    3-6t+6t^2.

The former is strictly concave and the latter strictly convex. Moreover
for C3 followed by a singleton sink, the entire converse segment has the
constant value3, despite its endpoint gradients C being different. A
scalar directional derivative or strict convexity hypothesis cannot
replace the full gradient argument. A2 succeeds because every quadratic
term along its particular triangle direction has the same curvature sign;
that mechanism has not been proved for arbitrary disagreement patterns.

## Verification and scope

`check.py` enumerates all1,099 labelled tournaments through order5. It
computes H,S,E,C by literal vertex permutations and two-path cuts, and
independently computes induced Hamiltonian counts by endpoint subset DP.
It verifies Euler's identity, the two-colour complement formula, mod-four
divisibility, converse transposition, and every odd-order separation in
that universe. For every cyclic triple it directly counts M_+,M_-, checks
the bound(8), computes the directional derivative change from the two
independent endpoint C matrices, and checks equation(7). The three
variational hostiles are computed as exact polynomials. All checks remain
active under `-O`; normal and optimized receipts agree. Universal claims
rest on the preceding proofs, not finite extrapolation.

From this directory, run `python3 check.py` and `python3 -O check.py` and
compare stdout with `receipt.json` and `receipt.optimized.json`. The optional
`--output PATH` writes UTF-8 JSON to a separate replay file. The checker
hashes its source at startup and refuses output if that source changes.

Map/loss ledger: two-port responses map exactly to the endpoint-cover
gradient; A1 uses its symmetry predicate, and A2 a particular directional
contrast. These operations preserve concrete distinctions but do not
reconstruct every named arc by sign. The deletion mechanism restores a
family of positive-cut constraints; its remaining degree residual is
eliminated only at the stated port/order thresholds. Universal two-port
faithfulness beyond the audited finite orders remains open.
