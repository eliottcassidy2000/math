# Two-sided seams: exact factor selection and a renewable strong +2 bridge

**Status: PROVED IN NOTE for identities and construction; finite controls are
reported separately. Chan–Pak remains OPEN.** The raw source/sink laws are
inherited from `session194/two-sided-padding/PROOF.md`, equations (2)–(3).
Session274 supplies the framed one-sided seam and its three-mode profile.
The contribution here is an exact two-sided scalar selection criterion,
a continuation invariant visible in four physical experiments, and a
root-transport pair whose difference is exactly two at every fan length.
No novelty is asserted for general path-cover substitution or two-sided
padding itself.

## 1. Physical inputs and conventions

Fix a nonempty tournament T and root a. Write U=T−a and A=N_T^+(a).
Let H(empty)=1, and let s_i(V), e_i(V) count Hamilton paths starting/ending at i.
Define the following integer packet, with empty sums zero:

    h = H(U),   F = H(T),
    q = sum_{S subset U} H(S) H(U\S),
    R = sum_{S subset U} sum_{i in S cap A} s_i(S) H(U\S),
    L = sum_{S subset U} sum_{i in S cap A} e_i(S) H(U\S),
    s = sum_{i in A} s_i(U),   e = sum_{i in A} e_i(U).

Here q=Q_U(2) includes empty slots. In particular q=1 when U is empty.
These are ordinary positive combinatorial counts, not signed formal states.

A source seam adds x with a→x and x→u for every other old u. A sink seam
adds y with y→a and u→y for every other old u. Source seams insert after a
in a chosen Hamilton path; sink seams insert before a. Thus each operation
has a physical framed implementation. Apply m source seams and n sink seams
at the SAME named root, obtaining B_{m,n}(T,a). The source block X and sink
block Y are transitive, X→U→Y, X→Y, a→X, Y→a. The two operation types
commute up to the canonical relabelling preserving each block's order.

For m,n≥1 the result is strong, irrespective of T: each old u≠a lies on
a→X→u→Y→a; the same routes include every vertex of X,Y. For T={a},
use a→X→Y→a. This is not a claim about zero-fan boundaries.

## 2. Exact bifan response, including initial cases

Put

    beta  = −L−h,
    gamma = R−q,
    delta = F+L+2h−s−R+e.

Then, for all **m,n≥1**,

    H(B_{m,n}) = q 2^{m+n} + beta 2^m + gamma 2^n + delta.       (1)

**Proof.** The inherited source law with apex digit one is

    q'=2q, R'=R+q, L'=2L+h, s'=h, e'=e,
    F'=F+q−L−s.

Iterating m≥1 times gives

    R_m=R+(2^m−1)q,
    F_m=F+(2^m−1)(q−L−h)+h−s.

A sink seam has apex digit zero. Its law is R'=2R, e'=0,
F'=F+R+e. Hence n≥1 sink seams give
F_{m,n}=F_m+(2^n−1)R_m+e. Expansion proves (1). These equations also
hold when U is empty: the first source turns q=1 into q'=2, and s'=h=1.
They alternatively follow by assigning the new vertex to one of two
ordered path slots. No division or graph subtraction enters the construction. ∎

The omitted boundaries must use

    H(B_{0,0})=F,
    H(B_{m,0})=F+(2^m−1)(q−L−h)+h−s      (m≥1),
    H(B_{0,n})=F+(2^n−1)R+e               (n≥1).

For a singleton, (q,beta,gamma,delta)=(1,−1,−1,3), giving
(2^m−1)(2^n−1)+2. For a rooted directed triangle the packet is
(4,−2,−2,3), giving (2^{m+1}−1)(2^{n+1}−1)+2. The latter happens to
extend to m=0 or n=0; that special convenience is not universal.

## 3. Seam curvature: a measured invariant of the whole controller

For the scalar packet p=(q,beta,gamma,delta), define

    Delta = q delta−beta gamma,    kappa=Delta/q.                 (2)

Call kappa **seam curvature** only as a name for this rational invariant,
not a geometric curvature theorem. Four physical measurements at
(m,n)=(1,1),(1,2),(2,1),(2,2) recover p exactly:

    q=(H22−H21−H12+H11)/4,
    beta=(H21−H11)/2−2q,   gamma=(H12−H11)/2−2q,
    delta=H11−4q−2beta−2gamma.

In particular

    kappa=(H11 H22−H12 H21)/(H22−H21−H12+H11).                  (3)

Continuing from a bifan with m0,n0≥1 scales the response to

    (q 2^{m0+n0}, beta 2^{m0}, gamma 2^{n0}, delta).

Therefore **kappa and delta remain unchanged by all further same-root
source/sink seams**. The proof is substitution in (1), not a fitted
recurrence. The four-experiment theorem concerns scalar bifan continuations;
it does not classify rooted tournaments, whole path-cover profiles, gauges,
or contexts with a changed root.

## 4. Exact target selection is a shifted divisor problem

For a prescribed integer N, (1) is equivalent to

    (q 2^m+gamma)(q 2^n+beta)=qN−Delta.                         (4)

Moreover −q≤beta,gamma≤−1. For beta this follows from L≤q−h:
summing L over ALL possible selected vertices counts two-slot covers with
a nonempty specified slot, hence q−h; gamma follows identically from R.
The empty-U case has beta=gamma=−1. Thus both factors in (4) are positive
for m,n≥1.

Consequently a target belongs to this fixed bifan family **iff** a positive
factor d of M=qN−Delta has both

    (d−gamma)/q = a power of two at least 2,
    (M/d−beta)/q = a power of two at least 2.                    (5)

This is an exact finite inverse algorithm, including a physical witness
(m,n), not merely a congruence test. If M≤0 the target is absent. Factoring
M is not asserted fast. Trial division is sufficient for termination.

The family is sparse: for N≤X, (4) and the coefficient bounds give
q²(2^m−1)(2^n−1)≤qN−Delta, hence
2^{m+n}≤4(X−kappa)/q. A fixed seed contributes O((log X)^2)
values, as does any fixed finite seed bank. This is a scoped explanation
of why this renewable controller alone cannot cover all odd targets.
It is consistent with Session194's broader lacunarity results.

The shifted product+2 triangle example is therefore a mechanism, not a
universal factorization of N−2. In particular many N−2 are prime.

## 5. A root change supplies a physical all-length difference of two

Take the labelled seven-vertex tournament with outgoing bit rows

    T = [96,97,59,115,35,0,52].

It has H(T)=19. Its two root packets are

    root 0: (168,−104,−137,109),
    root 4: (168,−134,−107,107).

Each packet is computed directly from the explicit subset sums in §1;
`check.py` gives independent permutation controls for these small inputs.
The two packets have identical q and beta+gamma=−241, but delta differs
by exactly two. Therefore, for EVERY m≥1,

    H(B_{m,m}(T,0))=168·4^m−241·2^m+109,
    H(B_{m,m}(T,4))=168·4^m−241·2^m+107.                       (6)

These are strong tournaments on 7+2m vertices. First pairs are
(297,299), (1831,1833), (8931,8933). They are independently constructed
from the same H=19 source, not formal differences of graphs.

The root-transport transformation is bounded-width in a precise sense.
Changing root r to root t flips exactly the edges between {r,t} and the
new fan vertices. Perform the full vertex-cut switch at {r,t}, then repair
its 2(|T|−2) flipped OLD boundary edges. This realizes the transport with
one cut and ten old-edge repairs, independently of m. The old tournament
itself and all fan-internal edges are restored. This statement concerns
ordinary full cut switching; it does NOT assert that every intermediate
repair preserves a fixed Hamilton path or stays within a |Delta H|≤2
corridor. Each final tournament has its own canonical framed construction.

More generally, two packets give a constant gap along balanced fans iff
they have the same q and beta+gamma. Necessity follows because their
difference is a polynomial of degree at most two in 2^m vanishing after
subtracting the constant at infinitely many m. The constant gap is exactly
the difference of delta. This is a cheap exact criterion for finding
renewable sibling bridges rather than guessing from a finite run of counts.

## 6. Hidden continuation information and a hostile near-theorem

Two six-vertex rooted inputs have identical ENTIRE one-sided scalar fan
sequences but differ once the opposite seam is available:

    ([58,36,17,6,10,28], root 1): (84,−31,−50,41),
    ([38,16,34,7,45,10], root 5): (84,−34,−50,43).

Both have H=37 and sink-fan counts 34·2^n+5 for n≥1 (73,141,...).
Their (1,1) bifans have counts 215 and 211. Thus keeping every future
one-sided H-value does not keep the two-sided target menu. This is a
concrete continuation-loss example, not a claim that their full coloured
profiles agree.

A tempting quantitative extension, kappa≥2H(T−a), survives every rooted
labelled input through order five but is FALSE. At

    T=[44,25,34,36,45,2], root 2,
    p=(84,−33,−60,41), H(T−2)=9,

kappa=122/7=18−4/7. The failure is in promoting a small-order curvature
floor to a universal inequality. No universal sign/floor theorem for
kappa is asserted here. The exact identities do not need one.

## 7. Connection and next constructive obligation

The less obvious historical analogy is the Rule30 moving-observer theorem
`THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary.md`:
a fixed observer can have a completely controlled infinite continuation,
while moving the observer changes which information matters. Here the map
is 'fix root and read H along seam words'; the preserved predicate is the
exact bifan response, lost data are other roots' packets, and the required
sidecar is the root menu. No Rule30 theorem is used to infer a tournament
statement. Root transport is the explicit replacement for the analogy.

Equation (6) is a reusable +2 **sibling bridge**, not a path visiting every
odd count. Increasing m jumps exponentially. The next useful target is a
physical transition between different root packets which matches a high
endpoint of one bridge to a low endpoint of another while retaining the
next admissible roots. A positive solution could connect these bridges to
the low-jump metagraph corridor. Neither unbounded sibling heights,
curvature change, nor exact local gap two proves such a connection.

## 8. Balanced-seam translational transport of a whole cut family

This extends the isolated sibling bridge to actual move geometry.
Fix the deletion core U, and vary only the outgoing apex set A. In the
notation of §1 let

    c(A)=R(A)−L(A),   z(A)=e(A)−s(A).

Both charges are additive set functions of A. For m≥1, equation (1) gives

    H(B_{m,m}(U+a_A)) = F(A) + C_m(U,c(A),z(A)),
    C_m(U,c,z)=q4^m−(q+h)2^m+2h+c(2^m−1)+z.                (7)

**Translation theorem.** On every family of apex cuts with the same c,z,
balanced fans translate ALL Hamilton counts by the same integer C_m.
Every old edge change, and every shared Hamilton frame, is retained by the
extension. In particular a single-free-tile move with count difference two
remains that same single-free-tile move with difference two at every scale.

**Proof.** Substitution into (1) proves (7). Fan incidences depend only on
the named apex, not on A, so comparing two children leaves exactly the old
changed edges. Source-seam insertion immediately after the apex and
sink-seam insertion immediately before it apply equally to a shared old
frame; an old edge absent from that frame remains absent after insertion. ∎

For each vertex v of U define R_v,L_v by the §1 formulas with A={v}.
Call v a **balanced seam port** if R_v=L_v and e_v(U)=s_v(U).
Arbitrary changes of apex incidences at ANY set of balanced ports preserve
both charges, so the theorem applies to the complete cut subcube on those
ports. This is a physically realizable interface condition, not a claim
that arbitrary cut subcubes have the property.

A converse symmetry of U fixing v certifies a balanced port by path
reversal. More generally, unions of orbits of a specified converse
isomorphism have c=z=0. A direct elementary source exists at every odd
transitive order: in U=TT_{2r+1}, r≥1, its middle vertex has

    R_v=L_v=2^r,   s_v=e_v=0.

Indeed in a transitive subtournament the only starting vertex is its
smallest and the only ending vertex its largest. A subset whose smallest
vertex is v chooses any of the r later vertices; each contributes one.
The ending computation is dual. This also explains why the endpoints and
the singleton require separate treatment.

### Vertex-transitive cores supply an entire balanced cut cube

**Corollary.** Every vertex of a nonempty vertex-transitive tournament U
is a balanced seam port. Therefore the ENTIRE 2^{|U|}-vertex apex-cut cube,
with all its old edge-flip geometry and all count differences, translates
by the same C_m(U,0,0) under balanced seams. This is the labelled apex-cut
cube; it is not asserted to possess one Hamilton frame common to all cuts.
Frame preservation remains conditional on an already shared frame.

**Proof.** Automorphisms act transitively on vertices and preserve each
count. Summing single-vertex endpoint counts gives sum_v s_v=sum_v e_v=h,
so s_v=e_v=h/|U|. Summing R_v or L_v over v counts every ordered two-path
cover whose specified slot is nonempty once, hence q−h. Thus
R_v=L_v=(q−h)/|U|. This is a direct counting argument; vertex transitivity
need not supply a converse automorphism fixing each vertex. ∎

For U=C3, h=3,q=12. The eight apex cuts have counts 3 or5. The common
translation is C_m=12·4^m−15·2^m+6, giving at m=1 the two values27 and29
on the entire same cut cube. This supplies a reusable whole-family
interface, not merely a pair of prescribed inputs. No assertion is made
that larger vertex-transitive cut spectra are intervals.

### Exact lift of the 61→63 metagraph portal

Use the source from the landscape lane:

    T63=[16,5,17,7,74,31,175,63],
    T61=reverse the edge {3,4} in T63.

Retain root a=4. Deleting it gives the SAME transitive seven-vertex core,
and original vertex3 is the middle balanced port. The exact packets are

    T61 at4: (128,−18,−60,11),
    T63 at4: (128,−26,−52,13).

Thus, for every m≥1, the strong children are related by the SAME single
edge reversal {3,4} and have respective counts

    128·4^m−78·2^m+11,
    128·4^m−78·2^m+13.                                    (8)

At m=1 this is 367→369. Both old inputs admit the common frame
recorded in `controls.json`; applying the frame insertions proves that
{3,4} remains a free tile, so this is an honest low-jump framed edge at
EVERY scale. Unlike the bounded-width root-transport observation in §5,
no intermediate repair path is needed here.

**Remaining limitation.** Equation (7) preserves the width of a cut-family
interval, not necessarily its baseline or overlap with its next translated
copy. The lifted portals alone do not connect all scales. The constructive
next obligation is to connect these translated families through shared
endpoint states or to supply sufficiently wide balanced-port cut intervals
whose successive translates overlap. No such all-order overlap is claimed.
