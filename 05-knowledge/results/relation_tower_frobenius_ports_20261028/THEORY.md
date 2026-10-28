# Relation towers, normal ports, and attachment transport

2026-10-28. **PROVED** elementary deductions below; inherited results are
identified explicitly. Computations are **FINITE-EXACT** in their stated
universes. Same-agent proof review and independent computational paths;
no independent-agent audit or literature-priority claim.

## T1. The relation-object tower continues through 2047

Let H_r have, once each, the nonzero columns of F2^r, and let C_r=ker H_r.
Its dimension is d=2^r-1-r. Regard each nonzero member of C_r as a new
point, retaining addition of codewords. This new configuration is exactly
the nonzero-vector configuration of a d-dimensional binary space. In
particular its three-point relations are {a,b,a+b}, with a,b distinct
nonzero words. No coordinate choice is required for this assertion.
A chosen kernel basis identifies it with the usual check matrix H_d.

Consequently the rank update is r -> 2^r-1-r, and the point counts starting
at the Fano configuration are

    ranks:   3 -> 4 -> 11 -> 2036
    points:  7 -> 15 -> 2047 -> 2^2036-1.

The 2047-point stage has 2047*2046/6=698027 three-point relations and
relation-space dimension2036. This follows directly from rank-nullity and
the fact that each triple is counted six times by ordered distinct pairs.
It iterates at every rank; no prime-existence theorem is used.

The earlier [relation-gluing T3](../relation_gluing_7_15_20261028/THEOREMS.md)
is the 7->15 instance. Its 15-point five-element circuits, five-line
spreads, and six pentads are different structures; their equal small
numbers are not used to derive this recurrence. F64 and F256 below are
degrees6 and8, not the next fields of this tower (degree11 comes next).

### A finite arithmetic realization of the next stage

Take coordinate primes l=(2,3,5,7,11,13,17,19,23,29,31). For each nonzero
11-bit mask b, let d_b be the squarefree product indicated by b. The saved
bank gives a positive a_b such that p_b=d_b*a_b^2-1 is prime, for every
one of the2047 masks. The primes are distinct because p+1 has a unique
squarefree part. Their maximum is10028024506499. Every prime is backed by
a recursively checked Lucas primality certificate, not merely a probable
prime test. Search used SymPy; the certificate verification uses integers.

For any subset J of this **distinct-prime** bank, multiplicativity gives

    sigma(product_(b in J) p_b) = product_(b in J)(p_b+1),
    this is a square iff XOR_(b in J) b = 0.

Thus this finite bank realizes Hamming [2047,2036,3]. Its square-total
subset count is 2^2036. Only698027 triple relations were explicitly
enumerated; the much larger subset count follows from the proved rank.
The first15 masks reproduce the inherited fifteen-prime bank exactly.
No assertion that the same prime-search recipe works at every rank is made.

The old Hamming15 codewords are attached to new primes via their11
coefficients in the saved kernel basis. Addition of old words corresponds
to symmetric difference of old prime supports, equivalently
x star y=xy/gcd(x,y)^2. New triple square-total relations preserve this
binary dependence, not old integer magnitude, ordinary multiplication,
or aliquot iteration. The shared-prime sigma correction in the inherited
packet remains necessary for unrestricted multiplication.

## T2. Normal correlation profiles recover a relative phase

Let r=2m, E=F_(2^r), q=2^m, F(x)=x^2, and

    B(x,y)=Tr_(E/F2)(x*y^q).

The inherited [P3](../free_bits_signed_incidence_20261028/THEORY.md)
proves that B is the unique nondegenerate alternating binary form that
is F-invariant and vanishes on every F_q-line. Retaining the actual
field multiplication and subfield spread is part of the hypothesis.
Frobenius invariance alone is not enough.

Define P(x,y)_j=B(F^j x,y), 0<=j<r. If x is normal, meaning its r
conjugates are a binary basis, then y -> P(x,y) is a linear isomorphism
E -> F2^r: a zero profile means y is orthogonal to a basis and hence zero.
Also P(Fx,Fy)=P(x,y), by F-invariance. For two fixed normal Frobenius
orbits O_x,O_y, choose representatives x,y. Then

    (F^a x,F^b y) / simultaneous Frobenius
                    <-> b-a mod r
                    <-> P(x,F^(b-a)y).

The last map is injective by the preceding isomorphism. Normality is a
sufficient uniform condition on the anchor, stronger than having degree r.
For a particular pair of orbits the exact weaker condition is simply that
the r profiles are distinct. The theorem does not assert normality is
necessary for every such pair.

### Sharp failure of replacing normality by full period

In the polynomial field F64=F2[t]/(t^6+t+1), integer masks mean polynomial
coefficient bits. The triples (x,y,y')=(2,34,35) have x and y of period6,
y'=F^3(y), and P(x,y)=P(x,y')=36. The span of x's conjugates has rank5.
Because x has trivial Frobenius stabilizer, (x,y) and (x,y') are different
diagonal orbits despite equal marginal orbits and complete profiles.

In F256=F2[t]/(t^8+t^4+t^3+t+1), (2,18,19) is the corresponding
period8 hostile: y'=F^4(y), span rank7, common profile221. The missing
direction is1: all conjugates of x have trace zero, so B(F^j x,1)=0.
This explains the collision, rather than just reporting it.

### Storage is not the same as legal probe cost

For normal x,y, the XOR of the profile bits is1. Indeed the trace of a
normal element is1 (otherwise its basis conjugates have a nontrivial
sum relation), so sum_j F^j x=1 and B(1,y)=Tr(y)=1. Hence any r-1
profile bits determine the last. In the worst case r-1 individual bit
probes are necessary, even adaptively: choose the B-dual normal element
y to x. The profile and its rotations are the r unit vectors. After
r-2 negative probes two possibilities remain. The dual conjugates form
a basis because B is invariant and nondegenerate.

Thus a relative phase has only r states, but ceil(log2 r) arbitrary-bit
storage is not always attainable by that many allowed correlation-bit
queries. At r=6 the16 ordered pairs of normal orbits split8 with minimum
nonadaptive cost3 and8 with cost5. At r=8 the256 pairs split192 with
cost4,32 with cost6,32 with cost7. These distributions are finite exact;
the worst-case r-1 theorem holds at every positive even r.

The normal-basis/trace setting is classical; compare
[Zhang, Feng, Liao and Gao (2013)](https://arxiv.org/abs/1303.2283), which
studies prescribed trace self-correlations. The present pair-profile
decoder does not assume a self-dual normal basis: those need not exist
in degrees divisible by4. The proof above supplies the needed dual pair.

## T3. A smallest closed carrier for normal relative ports

Fix k>=1 named ports, each with a specified normal Frobenius orbit in E.
The metadata (field, orbit identities, names) is held fixed. The object
is the joint tuple modulo simultaneous Frobenius. Observations are all
pair profiles P(x_u,x_v). Operations are known individual Frobenius
updates, restriction, re-rooting, and gluing two tuples along named
common ports when their restrictions agree.

Choose one representative per orbit and a rooted spanning tree on the
ports. Store on each directed edge u->v the increment

    delta_uv = a_v-a_u in Z/r.

**Claim.** These k-1 increments form a closed, lossless carrier, with
exactly r^(k-1) states. No smaller state space supports all the specified
observations. Neither an absolute root phase nor all pair profiles need
be stored.

**Proof.** Set the root potential a_root=0 and sum increments to every
vertex. This reconstructs the tuple up to the omitted common phase.
Conversely every edge-increment assignment occurs and is recovered
uniquely from its joint orbit by T2. Pair profiles are computed from
a_v-a_u. If w=LCA(u,v), this difference is the sum from w to v minus
the sum from w to u. Root re-selection changes all potentials by a common
constant; a local F-update adds1 at that port and changes only incident
increments. Restriction can rebuild any tree on the surviving ports.

For gluing along a nonempty overlap, align one shared potential and check
all remaining shared potentials agree. These checks are derived from
the held increments. On agreement the two potential lists join, and a
new spanning tree can be chosen. On disagreement no joint object exists.
Gluing along an empty overlap requires a new relative phase: it cannot
be recovered for free. The same is true for adding a new port before
any relative observation has been supplied. If a cycle of measurements
is added, consistency is exactly zero sum of increments around it.

Finally any two distinct increment assignments differ in some pair
phase, so T2 supplies a distinguishing allowed observation. Every exact
carrier therefore has at least r^(k-1) states. QED.

This is minimality for the specified language and fixed metadata, not
a universal minimum for tournament attachments. Pair data is obtained
"for free" only after the connecting tree data and the injectivity
hypotheses are present. LCA is an evaluator, not an information source.
Adding absolute named-value queries requires restoring the root phase.
For nonnormal ports retain the unresolved profile fibres or additional
measurements; the two explicit hostiles forbid dropping them.

## T4. The same path-composition mechanism handles signed labels

For a tree with integer vertex labels and fixed edge equations

    x_v = s_e*x_u + w_e,     s_e in {+1,-1},

there are unique coefficients epsilon_v in {+1,-1}, t_v in Z with
x_v=epsilon_v*z+t_v, where z is the root label. This follows by composition
of the invertible affine maps z -> +/-z+w. LCA paths evaluate all other
pair transports; adding a chord requires

    (epsilon_v-s_e*epsilon_u)*z = s_e*t_u+w_e-t_v.

The coefficient is0 or+/-2. A zero coefficient requires zero right side;
a nonzero coefficient fixes an integer root iff the indicated division
is integral. This recovers the inherited signed-incidence holonomy rule.

There is an exact further feasibility reduction for bounded, distinct
integer labels. Given intervals L_v<=x_v<=U_v, intersect their pullbacks
to obtain one integer root interval I. Intersect with any chord-imposed
singleton. If two vertices have equal epsilon and t, feasibility fails.
Pairs of opposite epsilon exclude at most one integer root,
(t_w-t_v)/(epsilon_v-epsilon_w). Let D be the set of all such integral
roots. The exact feasible root set is I minus D. The solution count for
a finite I is |I|-|I intersect D|, unless a prior inconsistency occurred.
This is an iff: the interval checks cover all bounds, chord checks cover
all equations, and D covers every remaining collision.

Applications with precise limits:

- Square-sum edges x_u+x_v=k_e^2 use s=-1,w=k_e^2. On an even cycle
  there is an alternating-sum obstruction; an odd cycle fixes twice the
  root. Prescribed labels in1..N and N vertices are a permutation exactly
  when the bounds and distinctness checks pass.
- Graceful edges, after choosing a sign and a permutation of lengths
  1..|E|, use s=+1,w=+/-length. Restricting the |E|+1 distinct vertices
  to0..|E| gives exactly the graceful condition. Choosing the edge signs
  and permutation is still an existence problem; this does not settle
  the graceful tree conjecture.
- Binary arithmetic relations use additive boundary syndromes and
  intersection spaces at LCAs, as in inherited relation-gluing T1. Their
  group is linear; the signed-label group is the integer affine group;
  T3 uses a cyclic group. Equal path-composition syntax is not an
  identification of these different carriers.

Retain the named coefficient list for future gluing. Merely retaining
the feasible root count or endpoint coefficients loses interior
collisions. All-different constraints connect nonadjacent vertices and
remain necessary even though the original equation graph is a tree.

## T5. Tournament transport: all-order generation, finite recreation

The inherited [Session79 root-switch theorem](../session79_tournament_enumeration_20261027/ROOT_SWITCH.md)
already gives an all-order constructive algorithm. Given every parent
class T of order n and its actual action Aut(T) on named vertices, form
the states (T,S)/Aut(T), S a subset specifying the new vertex's out-neighbours.
For each resulting child, reroot at every vertex. Classify the deletion
parent and carry S through the **same** isomorphism, then quotient its
mask by the target parent's automorphisms. Let Q_ij count these rerootings.

Each child class contributes one block with identical rows, whose
entries are its vertex-orbit sizes. Thus

    Q^2=(n+1)Q,     A000568(n+1)=tr(Q)/(n+1)
                              = number of distinct rows.

Choosing one raw extension per distinct row produces the next complete
bank, without canonicalizing the full child. Starting with the empty
tournament recursively recreates the sequence at every finite order.
This statement is inherited, not a new formula inferred from small n.
The new implementation extends the earlier parent-only finite check
through order9 and exports every representative at that order.

There is also a useful transport identity checked independently through8:
if i is a rooted state represented by (T,S), then

    |Aut(child)| = Q_ii * |Stab_Aut(T)(S)|.

The second factor is exactly the stabilizer of the new root; Q_ii is its
vertex-orbit size. This is orbit-stabilizer, so no new child isomorphism
search is needed for the automorphism order. For future attachments the
**action**, not just this order, remains necessary. The next stage computes
that action from the retained representative.

Every tournament automorphism has odd order (an involution would exchange
an oriented pair), so every vertex-orbit size is odd. Each block therefore
has rank1 over F2, even when n+1 is even and Q^2=0 mod2. This is not a
Laplacian: rank and rank-of-square differ at even child orders.

The independent [Session67 native metagraph construction](../session67_metagraph_enumeration_information_20261026/session67-root/COUNT_FREE_ENUMERATION.md)
uses all fixed-Hamilton-path presentations and legal relocations.
Its connected components are ordinary tournament classes; quotienting
the vector space on presentations by edge differences has dimension
A000568(n) over every field. That construction had already enumerated
through8. The present root-switch census is checked against its full
Hamilton-path-orbit-size histograms through8, not merely the total counts.
These are different auxiliary graphs; literal complement lines of the
metagraph must not be mistaken for rerooting or relocation edges.

This links the two meanings of "relations become objects": a kernel
turns dependencies into a new configuration (T1); a quotient turns
equivalent presentations into classes (T5). A kernel and a quotient are
different constructions. Neither alone retains the full group action
needed for later attachments.

### The 81/27 control remains decisive

Let A=C3->1, D=1->C3. The two named12-vertex tournaments
C3[A,A,A] and C3[D,A,A] have the same67 labelled cyclic triples and
462915 Hamiltonian paths. Their automorphism orders are81 and27.
In the first, all three child types agree and the root C3 rotation
survives; in the second the distinguished child destroys that rotation.
The internal group contributes3^3 in both. The identity vertex order
is Hamiltonian in the first and fails at the changed block in the second.

Frobenius/trace recovery has no supplied intrinsic map that restores
this missing child-type bit. A minimal closed relative-field carrier
therefore does not compress the full tournament carrier. All named
legal re-marking queries already reconstruct the tournament at every
size by the inherited guard-reconstruction theorem. Unnamed restricted
attachment families remain the appropriate place to seek a smaller
carrier; dropping actual child types fails this hostile immediately.

## T6. Recovered five/seven chart transport, and its failure at fifteen

The random inheritance pull recovered a particularly relevant established
result: [regular-seven attachment transport](../odd_cell_completion_20261023/seven/PROOF.md)
and [regular-five transport](../odd_cell_completion_20261023/five/README.md).
Their partial-bijection gluing mechanism predates this session. This is
the direct tournament application of the same transport principle, with
actual chart compatibility in place of an invented field encoding.

For completeness, the elementary general statement is this. On a
connected graph of finite state domains X_v, suppose each edge relation
is a partial bijection. On a spanning tree, a root state determines at
most one state at every other vertex. Discard root states where a path
map is undefined or an off-tree edge fails. The remaining root set is
in bijection with all global compatible tuples. Hence there are at most
min_v |X_v| tuples, not a product of independent domain sizes. Pair
relations follow by path composition and inverses through the LCA,
restricted to the surviving root states. This proof is induction along
the tree followed by exact cycle checks. An empty surviving set is a
rejected comparison, not a theorem that either input object is invalid.

The inherited five-cell lemma proves the partial-bijection condition for
nonuniform equitable blocks between cells whose internal automorphism
groups are regular C5. A surviving connected symmetric support needs
at most one external actual owner; some supports already name one
intrinsically. The seven-cell lemma proves it for **every** nonuniform
equitable block between regular-seven tournaments, even with repeated
binary rows or columns. Two actual owners in one root cell suffice for
literal directed refinement throughout that connected support.

The seven proof uses its special nonidentity automorphism shapes7 and
1+3+3. In a row invariant under1+3+3, the fixed-point entry is determined
by row sum mod3; equitability on the other axis makes a nonuniform
one-sided kernel impossible. This is stronger than a raw matrix-rank
condition, and more specific than general rigidity. The F64/F256
normality test is another sufficient condition for a functional edge
decoder; neither hypothesis implies the other across these domains.

**PROVED hostile at15.** Let R=C3[T5,T5,T5], where T5 is the cyclic regular
five-tournament with out-steps1,2. R is regular of degree7. Take the15-by15
attachment block M=diag(J5,J5,J5). It is nonuniform and equitable, with
every row and column sum5. Yet its one-sided chart kernel in Aut(R) is
exactly C5^3, of order125. Each internal five-rotation fixes M, and fixing
M forces each of its three column classes to be preserved, leaving
exactly these independent rotations. Thus the seven-cell theorem does
not extend to arbitrary regular-fifteen cells.

For identity comparison across this attachment, one endpoint chart has
125 compatible charts at the other endpoint, before more constraints.
The destroyed data is the three internal cyclic phases. Retaining
only the nonuniform block's existence, marginal degrees, or a presumed
tree gain cannot restore them. One can instead retain the fibre itself
or quotient by these kernels **only if** later observations are also
kernel-invariant. Literal owner queries force them back.

This15 is the order of a composite tournament; the15 in Hamming [15,11,3]
is a set of nonzero vectors in a4-dimensional space. Here they meet
through a precise contrast in stabilizer kernels, not through an
asserted equivalence. The 81/27 pair is the same warning at a different
composition boundary: child-type and action data survive after cycle
supports have forgotten them.

## Incoming closure for aggregate ordinary composition

At close-out, commit08840182798686f53fda188add62a82fbb60a8ad supplied
[Session86's universal rigid-palette theorem](../session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md).
For prime inventories of bounded arity3..K, equality of aggregate
vertex-action cycle indices is equivalent to equality of ordinary
inflation counts for every finite nonempty palette of distinct rigid
transitive children. The single palette with child sizes
1,(K+1),...,(K+1)^(K-1) already separates different indices.

Its elementary mechanism is injectivity of power-sum substitution into
K color variables (triangular by merging partition parts), followed by
collision-free base-(K+1) encoding of their exponent vectors as sizes.
Every digit is at most K. For tagged root multisets, tags must remain
counted; deleting them changes the observer. We reviewed this proof and
replayed the supplied exact computation, agreeing with all our generated
ordinary counts through9 and prime counts through8.

This answers a different minimality question from T3: universal aggregate
ordinary substitution counts, with no request for individual class or
stabilizer distributions. It also differs from T6's exact chart comparison
and the named legal-move language that reconstructs a full tournament.
Thus there are three precise levels: aggregate cycle index; relative
port/gain data where injective; and actual joint attachment fibres and
actions where kernels survive. The81/27 pair does not contradict the
aggregate theorem because a single decoration's stabilizer is outside
its promised observations.
