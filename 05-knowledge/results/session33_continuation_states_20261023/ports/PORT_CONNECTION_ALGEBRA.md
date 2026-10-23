# Port continuation algebra: what may be compressed, and for which observer?

Session 33, 2026-10-23.

**PROVED:** the finite-context minimal-state/linear-rank theorem; a direct
separator for every two distinct perfect port pairings; the exact criterion
for safe symmetry aggregation; and a same-response/same-determinant electrical
hostile. These are elementary specializations of established continuation and
connection-matrix ideas, not priority claims.

**FINITE-EXACT:** complete 2,4,6,8-port matching matrices, exact rational and
binary ranks, all pairwise separator certificates, and a context-relative
minimal compression of the three specified G25 blocks. The G25 quotient retains
the count 10 and, with endpoint-polynomial weights, its complete endpoint spectrum.

**CITED:** the matching-connectivity literature, with scope given below. No
asymptotic Hamiltonicity, new complexity lower bound, or directed tournament
claim is made. No Tutte determinant or electrical response is asserted to retain
Hamiltonian paths.

## 1. Inheritance and concept board

Closest proved mechanism: THM-2176,
`01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md`,
§2, in its corrected **commutative-monoid** scope: equality against every
continuation is the largest congruence inside equality of the target observable.
MISTAKE-249 warns that right contexts give only a right congruence in a general
noncommutative monoid; two-sided contexts are needed for a two-sided congruence.
The present finite response rows are
that same mechanism, with a graph-gluing operation rather than knot connected
sum. We do not transfer the metric, cocycle, or knot conclusions.

Closest exact graph result: Session32 `session32-ports/PORT_TRANSFER.md`, §§3–5:
named ports plus component pairings, prospective endpoints, and realization
multiplicities give a correct undirected Hamiltonian-path gluing law. Its
four-port hostile already proves that degrees and component counts alone fail.
Its endpoint-blind degree-two forcing near miss remains retracted.

The underused sidecar is the *family of allowed continuations*. It determines
whether a coarse state is harmless or fatal. The concept board is:

1. continuation rows (target-preserving deterministic equivalence);
2. linear connection rank (count-preserving linear compression);
3. symmetry (simultaneous action versus fixed external wiring);
4. endpoint polynomials (identity as a coefficient rather than a state);
5. Schur response (an exact compression for a different observer).

Anchor: lossless HP transfer. Niche: actual minimal G25 response rows.
Wildcard: Gaussian elimination/electrical equivalence versus Hamiltonicity.

Related THM-1975,
`01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md`,
has a VERIFIED finite directed substitution census and explicitly flags its
universal kernel derivation. It is not used as an undirected theorem here.
Uniform complete interblock adjacency can make unnamed path-cover counts
sufficient; sparse named wiring generally cannot. This is an operation
restriction, not a contradiction between the two objects.

## 2. Minimal observer theorem for a finite continuation family

Let S be a finite collection of local states, C a finite collection of permitted
external continuations, and M(s,c) the number of accepted completions of s with c.
A binary acceptance matrix is the important special case. State multiplicities
x_s produce continuation counts

    h_c = sum_s x_s M(s,c).

Define s ~ t iff M(s,c)=M(t,c) for every c in C.

**Theorem 1 (deterministic, linear, and tensor compression).**

(a) A quotient pi:S -> Q preserves every h_c for every multiplicity vector by
summing weights within pi-fibers and using one quotient response per fiber iff
M is constant on each fiber, equivalently every fiber lies in a ~-class. Thus
its unique coarsest deterministic quotient consists of distinct rows of M.

(b) Over a field F, the smallest dimension d of a linear encoding E:F^S -> F^d
from which every h_c can be recovered linearly is rank_F M.

(c) Let T(s1,...,sr,j) be any finite acceptance tensor. In coordinate i, identify
two states only if their slices agree against *all* choices of the remaining
original coordinates and j. All such coordinate quotients can be performed
simultaneously, and contracting T against multiplicity vectors equals contracting
the quotient tensor against the sums of multiplicities in its fibers.

**Proof.** For (a), sufficiency is distributing a sum by fibers. Necessity follows
by taking x to be a unit mass at either of two identified states. For (b), any
encoding/recovery gives the matrix factorization M^T=RE, so rank M<=d. A rank
factorization realizes d=rank M. For (c), replace coordinates one at a time:
the entire slice is unchanged at each replacement, even after earlier coordinates
have been replaced. The quotient tensor is therefore well defined. Distributing
its finite sum gives the claimed contraction identity. QED.

This theorem distinguishes two meanings of "smallest state": the number of
*distinct rows* and the dimension of their *linear span*. The latter permits
signed or field-valued coordinates and is not in general a positive-state
quotient. Over F_p it preserves counts modulo p, not integer counts or existence.

This is a static theorem for the specified continuation set. If states are to
be reused under future operations, that set must contain the new continuations
(or be proved closed under their composition). A quotient verified for one fixed
G25 environment is not automatically a universal graph substitution invariant.

## 3. Every distinct perfect port pairing has a hostile continuation

Let S_k be all perfect matchings on 2k named terminals. Interpret each matching
edge as one internally disjoint path fragment. A continuation R is a matching
on all but two terminals a,b. It is accepted precisely when P union R is a
single spanning path with endpoints a,b, with edge multiplicities retained.

**Theorem 2 (pairings are universally distinguishable).** For every k>=1, all
(2k-1)!! states in S_k have different continuation rows. For any P != Q and an
edge ab of P absent from Q, there is a near-perfect R missing a,b such that
Q union R is a spanning path and P union R is disconnected.

**Proof.** The Q-edges containing a and b are distinct. Order the k Q-edges so
that these are the first and last. Orient the first from a and the last toward
b; orient intervening Q-edges arbitrarily. Join the end of each ordered Q-edge
to the start of the next with an R-edge. Every terminal except a,b occurs in
exactly one R-edge, and Q union R is the displayed alternating path. In P union R,
a and b form the isolated P-edge ab, while 2k-2>0 terminals remain. Thus the
union is disconnected. The only k=1 state needs no separation. QED.

This lower bound applies to deterministic encodings of *individual states*
that must answer all these continuations. All local states here have the same
terminal degree vector and component count. It does not prohibit algebraic
compression of weighted sums, nor context-relative compression.

The abstract multigraph convention has an ordinary simple-graph realization:
subdivide each local and each continuation edge with its own private vertex.
Then duplicate terminal pairs become internally disjoint length-two paths,
not duplicated simple edges; connectivity and the path-vs-cycle distinction
remain unchanged. Endpoint and coverage requirements are preserved.

### The matching-connectivity matrix appears exactly, not analogically

Let H_k(P,Q)=1 when the edge-multiset P union Q is connected (a single alternating
cycle, with the k=1 double-edge convention). Removing any one edge ab of Q gives
a path context R=Q minus ab. Then

    [P union R is a single spanning path] = H_k(P,Q).

Indeed P union Q is 2-regular, hence a union of alternating cycles. Removing one
Q-edge breaks exactly one cycle into a path and cannot destroy any other cycle.
Consequently the full near-perfect-path connection matrix is H_k with each
column repeated k times. In particular their ranks over every field agree.

The explicit exact computations are:

| named ports 2k | matching states | distinct response rows | rank over Q | rank over F2 |
|---|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 |
| 4 | 3 | 3 | 3 | 2 |
| 6 | 15 | 15 | 15 | 4 |
| 8 | 105 | 105 | 105 | 8 |

For each matrix the certificate includes all entries, exact rational pivot rows
and columns, a nonzero full-rank minor determinant, a rational kernel basis,
and a separator for each pair of distinct rows. The verifier checks the
path/cycle closure identity for every P,Q and every removed Q-edge. There are
5,460 row-pair separators at eight ports.

A useful hostile to confusing parity with existence already occurs at four
ports. Ordering the three pairings gives H=J-I. The sum of its three rows is
(2,2,2) over the integers but (0,0,0) over F2. The collection of all three states
has two accepted completions under each context, yet its parity representation
is zero. A rank-two parity encoding does not supply exact integer counts or a
sound nonexistence test.

### Literature boundary

Cygan, Kratsch and Nederlof, *Fast Hamiltonicity checking via bases of perfect
matchings*, arXiv:1211.1506, study precisely this matching-connectivity matrix and
its F2 factorization/bases. Curticapean, Lindzey and Nederlof, *A Tight Lower Bound
for Counting Hamiltonian Cycles via Matrix Rank*, arXiv:1709.02311v2, prove its
rational rank has size 4^t/poly(t), where t is their number of terminals.
These are cited established results, not new deductions from the small table.
Our notation has t=2k. Their field-dependent asymptotic results explain why
linear rank and deterministic state count diverge, but we do not reproduce
those representation-theoretic proofs or claim a new algorithmic bound.

Primary sources checked 2026-10-23:
https://arxiv.org/abs/1211.1506
https://arxiv.org/abs/1709.02311

## 4. Symmetry: exactly when orbit merging is safe

Suppose a group G acts on S and C and M(gs,gc)=M(s,c). This says that M is
*equivariant*, not that states in an orbit have equal rows.

For fixed-context queries, orbit merging is valid iff

    M(gs,c)=M(s,c) for all g,s,c.

This is immediate from Theorem 1, and stronger than simultaneous equivariance.
It holds, for example, if G fixes every allowed context. More generally, if
only a G-invariant weighted mixture y of contexts is queried, then the scalar
response f_y(s)=sum_c M(s,c)y_c is constant on state orbits. Change variables
c=g c' and use y_(g c')=y_c'. Aggregating state weights by orbit is safe for
that restricted averaged observer. Orbit sums, not unweighted representatives,
retain counting multiplicities.

At four ports S4 acts transitively on the three pairings and simultaneously on
contexts, but H=J-I has three different rows. With a single context, its matching
state rejects and the other two accept. Merging all three fails. With uniform
context weights, every row sums to 2, so the averaged query really does admit
one state. This is the smallest nontrivial exact symmetry hostile.

## 5. G25: minimal response classes are much smaller than universal signatures

Use the exact Session32 partition:

    A={2,7,9,16,18}, B={4,5,11,12}, C={1,3,6,8,10,17},
    K={13,14,15,19,20,21,22,23,24,25},
    E(K)={13-23,14-22,15-21,24-25}.

The full named-port signatures have 5,101,149 states. For coordinate A, a
continuation is every pair of B,C states and every subset of E(K); similarly
for B,C. Thus the continuation universes have sizes 240,784; 11,920; 8,080.
There is no pruning by observed solutions in the definition.

An independent exhaustive vertex DFS from the unique graph leaf 18 enumerates
all ten full Hamiltonian paths. Each path has a unique local restriction in
the freshly reenumerated signature tables. The Session32 gluing theorem implies
that these ten entries are exactly the nonzero entries of the complete binary
acceptance tensor; every other entry is zero. This gives its sparse exact
representation without testing all 1,203,920 entries separately.

The following table uses zero-based state IDs in the attached certificate.

| block | original states | dead states | live response classes | minimum deterministic classes including dead | rational rank |
|---|---:|---:|---|---:|---:|
| A | 5 | 2 | {1,3}, {4} | 3 | 2 |
| B | 101 | 96 | {16,52}, {19}, {82}, {92} | 5 | 4 |
| C | 149 | 143 | {32,58,88}, {91}, {105}, {106} | 5 | 4 |

All live response rows within each coordinate have disjoint nonempty supports,
so their linear independence is also visible without elimination. Different
live classes cannot be merged; all dead states can. The dead class may be
omitted entirely in an implementation. With all original realization weights
one, the nontrivial live class weights are 2 in A, 2 in B, and 3 in C; the
remaining live weights are 1. Simultaneously contracting the 2-by-4-by-4 live
quotient with the sixteen core-edge choices still gives exactly ten.

These are context-relative minima, not estimates or universal minima. The
coincidence 3,5,5 including dead states does not reinstate the original three
supervertex contraction: these are distinguishable response classes within
three separate block tables, and the named external continuation remains.

### Endpoint identities can move into coefficients

Every live response class here has the same selected port set across all its
members. Therefore each fixed quotient tuple/core choice determines the core
endpoint set independently of representatives. Give each original state the
monomial product of z_v over its local prospective endpoints; sum these
monomials within a response class. The nontrivial weights become

    A {1,3}: z18 (z9+z2),       A {4}: z18,
    B {16,52}: z11+z4,
    C {32,58,88}: z10+z8+z3.

Other live B,C states have weight 1. For each accepted quotient tuple multiply
by z_v for each core endpoint v and sum. Exact contraction returns

    z18 (z2+z3+z4+z8+z9+z10+z11+z13+z22+z23).

Setting every z=1 gives 10. Taking coefficients recovers the endpoint spectrum.
The verifier checks this polynomial against the independent DFS endpoints.

This is not a claim that endpoint labels are universally dispensable: the
observer has explicitly retained them in polynomial coefficients, and the
same-selected-ports property was checked for these live classes. For arbitrary
quotients whose core endpoint output differs, the scalar response equivalence
must instead be refined to equality of polynomial-valued response rows.

The merge of A states 1 and 3 trades endpoint 9 for endpoint 2. It need not arise
from a block automorphism: vertex 9 has internal degree two while vertex 2 has
internal degree one. Thus safe compression is controlled by external behavior,
not merely by visible geometric symmetry.

## 6. An electrical response can be exact and still be the wrong observer

Consider two connected simple unit-conductance networks on vertices a,x,b,c,
with the same boundary B={a,b} and interior I={x,c}:

    P: a-x-b-c,           T: edges a-x, x-b, x-c.

P is P4 and has one unoriented Hamiltonian path. T is K1,3 and has none because
three leaves would all have to be endpoints. Both have four vertices, three
edges, one spanning tree, and cycle rank zero.

Let L be the graph Laplacian. Their interior matrices are

    L_II(P) = [[2,0],[0,1]],
    L_II(T) = [[3,-1],[-1,1]],

both with determinant 2. Direct exact elimination gives in both cases

    Lambda = L_BB - L_BI L_II^(-1) L_IB
           = (1/2) [[1,-1],[-1,1]].

For P, c is a dangling leaf attached at boundary b; for T it is a dangling leaf
attached at interior x. A dangling unit edge carries zero current at equilibrium
when no current is injected into its leaf. Hence moving that leaf changes the
Hamiltonian obligation while preserving the two-terminal response. The explicit
Schur calculation also verifies the determinant factor, not only resistance.

This gives a stronger exact failure than matching just an effective resistance:
for any same boundary quadratic modification K, the eliminated determinant
identity is

    det(L + extension(K on B)) = det(L_II) det(Lambda+K).

Thus these networks have the same such determinant observables and the same
boundary quadratic energy, including when grounding makes the determinants
nonzero, but different HP counts. No general decoder from (Lambda,det L_II)
to Hamiltonicity exists, even at four vertices. Keeping the original interior
factorization or extra degree/connectivity data is a different, richer carrier.

This does not refute matrix-tree theorems or a determinant algorithm that keeps
additional variables/state. It refutes only the proposed *unlabelled eliminated
response plus determinant* as a lossless Hamiltonian summary. The loss ledger is:

- source: a graph with interior vertex coverage obligations;
- target: its quadratic boundary response and Gaussian determinant factor;
- map: Schur elimination;
- preserved: all stated boundary quadratic/determinant responses;
- destroyed: degree obligations and path-component connectivity;
- repair: endpoint-aware path-forest transfer, not a scalar resistance correction.

## 7. Reproduction and stopping boundary

From this directory:

    python3 verify_connection.py
    python3 -O verify_connection.py

Dependencies: Python 3 and SymPy (tested with 1.14.0). No randomized arithmetic.
All rational ranks use exact DomainMatrix elimination; parity ranks use bitwise
Gaussian elimination. All checks raise exceptions explicitly, so -O removes
none. `certificate.json` records both source hashes. Normal and -O certificates
were byte-identical. `inherited_ports.py` is an unchanged copy of the Session32
verifier used only for its elementary graph/profile/DFS helpers; its main is not
executed by the new verifier. The new G25 certificate is regenerated from graph
edges, not read from the old result file.

Finite universes and controls: every pairing at 2,4,6,8 ports; every pair of such
pairings for separators and every deleted closure edge; every G25 HP by complete
DFS from its forced leaf; every local signature by exhaustive edge/port subsets;
all fixed G25 continuations via the proved sparse-tensor correspondence; and
exact P4/star Laplacian elimination. Positive controls include each explicit HP;
hostiles include unequal port pairings, wrong fixed-context orbit merging,
parity-zero positive counts, and equal electrical response with unequal HP count.

No repository or Git mutations were made. Startup/protocol/guardrail/current
routes were consulted. The coordination helper was blocked by absent
`.machine-id`; no identity was created. No theorem IDs were reserved. The scope
is intentionally bounded: larger square-sum graphs require fresh port families;
restricted-context compression does not prove a reusable all-n gluing scheme.
