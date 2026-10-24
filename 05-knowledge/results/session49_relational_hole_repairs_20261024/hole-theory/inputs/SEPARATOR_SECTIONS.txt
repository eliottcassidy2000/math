# Charge-complete sections across actual shared-record separators

Status: **PROVED conditional all-size compiler.** Actual native applications
require the declared full-record domains and factors to be authenticated.
Finite controls and source-specific reset witnesses are separate certificates;
no all-order native supply is inferred.

## 1. Inheritance and the missing construction

The closest proved mechanism is
`compatible_charged_helpers_20261023/compatibility/NEUTRAL_FIBRE_LIFTING.md`:
binary endpoint records give neutral upset intervals, with exact minimum
boundary cost. Its useful hostile is an electrically neutral shared-owner
component whose entire collection of literal records must change together.
The incoming
`session47_resilient_native_wiring_20261024/audit-helper-resilience-v2/AUDIT.md`
certifies twelve disjoint quadruples per charge, but only in the stated
compatible native bank. The corrected near miss is composing charge-onto
tables while forgetting their shared boundary values. The least-used sidecar
here is the exact source record, rather than only the charge of its fibre.

Concept board: complete local records; separators; physical row charge;
conditional cosets; source-specific reset work; shared-owner neutrality.
Anchor: construct compatible sections. Niche: split charge directions among
separate components. Wildcard: low separator width with arbitrarily large
mandatory native boundary cost.

## 2. Full-record factor carrier

There are finitely many labelled slots i with finite authenticated domains
D_i of FULL records. The record retains the actual path, source and target,
kernel label, normalized resource, inverse/converse partners, original local
guards, shared owner values and integer clocks. A resource fibre is not extra
capacity. Every factor is an actual allowed relation on a specified set of
slots, possibly of arity greater than two. These include all dynamic owner
equalities, capacity conflicts, background constraints and requested pins.

Fix an external core arrival a. Its incidences may change the local factor
tables, but not the meaning of a record or a resource. All factors must be
listed. In particular, a cross-component resource collision cannot be omitted
just because the individual components are internally injective. Disjoint
complete resource reservoirs are a sufficient reason that such a factor is
absent. Declared independently varying coordinates are not a substitute.

Let G be an abelian group. Slot i contributes gamma_i(r) in G. Exact integer
margins may be retained as Z coordinates. Electrical coordinates are the full
required physical current rows, including target diagonal contributions; a
smaller finite code may be used only through an authenticated equivalence.
Only finitely many charges occur, even when G is infinite.

The cost of replacing record r by t in slot i is kappa_i(r,t)>=0. It is the
actual changed-slot cost or owner-clock OCCURRENCE cost. It need not be a
function of charge or normalized resource. Additivity over labelled slots is
assumed. Counting distinct global changed coordinates is a different cost
model, requiring those coordinates to be assigned and counted once instead.

Supply a rooted tree of bags B_v with the running-intersection property:
bags containing any given slot form a connected subtree. Every factor has
its entire scope in an assigned bag. A slot domain is retained wherever the
slot appears. Set S_v=B_v intersect B_parent for a nonroot bag and S_root=empty.
Let O_v=B_v minus S_v. Each slot belongs to exactly one O_v. Thus charges and
costs are never duplicated at a separator, even when its record appears in
many bags. All external core work is retained separately unless its slots
are also explicitly part of this carrier.

## 3. Exact sections and counts

For a separator record s on S_v, let F_v^a(s,q) be the set of compatible full
assignments in the subtree at v whose separator is s and whose total charge
on slots OWNED by that subtree is q. Separator charges belong to an ancestor
and are not included again. Let N_v^a(s,q) be its cardinality.

For every local bag record t extending s and satisfying factors assigned to
v, write c_v(t)=sum_(i in O_v) gamma_i(t_i). Its child separator records are
t|S_w. The exact recurrence is

    N_v^a(s,q)
      = sum_t sum_(c_v(t)+sum_w q_w=q)
                    product_w N_w^a(t|S_w,q_w).                 (1)

Every nonzero summand retains t and its child witnesses. This reconstructs
an actual full section, not just a positive integer. At the root, N(q)>0 is
exactly feasibility at charge q. Prescribed full endpoints can be inserted
by restricting domains to their supplied records.

**Proof.** A compatible assignment determines its bag tuple and one child
assignment uniquely. Running intersection says that different child
interiors are disjoint and have no unrecorded shared slot. Every factor is
checked in its assigned bag. Conversely compatible child witnesses agreeing
with t glue to one full assignment. The owned sets partition the slots, so
their charges add exactly. These two constructions are inverse, proving the
count as well as existence. Induction proves (1) at every bag. QED.

This algorithm is exponential in bag width, not necessarily in the total
number of slots. For finite G, incremental child convolution uses at most
|G|^2 operations for each child and allowed bag tuple. Sparse finite charge
supports apply to Z margins. Reconstruction stores actual record choices;
it must not replace a marked lift with a convenient orbit representative.

## 4. A compositional sufficient condition, before global enumeration

A conditional coset certificate at (v,s) consists of a finite subgroup
H_v(s)<=G and an offset alpha_v(s), with an explicit section for every charge
in alpha_v(s)+H_v(s). In particular every Z-margin coordinate is fixed by
the offset; a finite subgroup contributes zero in those coordinates.

Choose a locally valid bag tuple t extending s. Suppose each child has such
a certificate at t|S_w. Then the parent contains the entire coset

    c_v(t)+sum_w alpha_w(t|S_w) + sum_w H_w(t|S_w).             (2)

Given a requested element of this coset, split its difference from the
offset into one element of each child subgroup, and invoke the supplied
child sections. This splitting can be tabulated by finite group addition;
it does not presume that actual slot choices are independent across an
unretained separator.

**Charge-direction composition theorem.** If, for every required boundary
record s, one supplied valid t makes the child subgroup sum contain H and
its offset equal alpha(s) modulo that sum, then it constructs a section of
alpha(s)+H. In particular, at the root a subgroup sum equal to the complete
required finite guard group Q, with the requested exact margins in its
offset, proves full Q coverage without enumerating all full assignments.

Several local tuples may instead provide several cosets (2). Their union
can certify coverage when no one tuple suffices: check the desired charges
against that finite coset union and retain a witnessing tuple for each.
This is a sufficient construction; it is not a claim that every feasible
charge set is a coset. Exact recurrence (1) handles arbitrary sets.

For example, in G=C2 x C3, one child can supply the C2 direction and another
the C3 direction, with offsets depending on the SAME retained separator.
Their sum is all G at each permitted separator value, even though neither
child is individually charge-onto. This is stronger than requiring a single
isolated full-charge reservoir.

The boundary qualification is essential. Two binary slots with charges x,y
in C2 are each individually onto. Imposing the shared equation x=y forces
total charge x+y=0. Forgetting the common separator would falsely report both
charges. Conditional on the separator, each child's charge set is a
singleton, so the theorem correctly supplies only the zero subgroup.

## 5. Exact two-arrival boundary cost

Take two external core arrivals a,b. They may impose different factor tables.
Define D_v(s,t,q,r) as the minimum total owned-slot cost between an assignment
in F_v^a(s,q) and one in F_v^b(t,r), with +infinity for an empty pair of fibres.
For locally valid bag tuples u extending s and z extending t, set

    k_v(u,z)=sum_(i in O_v) kappa_i(u_i,z_i).

Then D is computed EXACTLY by

    min_[u,z, q=c_v(u)+sum q_w, r=c_v(z)+sum r_w]
       { k_v(u,z)+sum_w D_w(u|S_w,z|S_w,q_w,r_w) }.            (3)

Retain minimizers and child backpointers. This produces two actual boundary
records attaining the minimum; no infimum or discarded marking is involved.

**Proof.** Apply the same unique decomposition as in (1) to BOTH assignments.
The owned sets partition their cost, giving the lower bound by child minima.
Conversely minimizers in the child fibres glue with the two bag tuples and
attain the displayed value. Finite domains guarantee attainment. QED.

The two-arrival table has charge pairs, so a naive finite-group convolution
may cost |G|^4 per child/bag-tuple pair. This is a finite-width constructive
compiler, not a polynomial-time claim in arbitrary width or group size.
Fixing the source full record collapses its charge/record choices and gives
the smaller source-specific min-cost section problem. A given current record
can therefore choose its CHEAPEST compatible next section, rather than pay
the diameter of an entire helper reservoir.

Minimum cost between two fibres alone does not prove bounded reachability
from every point in the first fibre. Either fix the actual current source in
(3), or separately prove within-fibre connectivity with its own budget. The
neutral-upset theorem cited above supplies one such mechanism. This caveat
is why the present result retains both the exact source and the separator.

A sharp hostile has two interior records0,1 at distance2 and three boundary
fibres A={0}, B={0,1}, C={1}. Their minimum distances are D(A,B)=D(B,C)=0
but D(A,C)=2. At interior budget1, B cannot change its record. Consecutive
zero-cost minima choose DIFFERENT records of B and cannot be concatenated.
Thus minimum fibre distance need not satisfy the triangle inequality. The
source-fixed message for current record0 correctly charges2 to reach C.

## 6. Native boundary cost need not be small at width one

Use the actual S_r=T_r -> C3 family and distinct prefix resources from
`native_packet_decomposition_20261023/workspace/ONE_BIT_WORKSPACE.md` §3 and
`compatible_charged_helpers_20261023/compatibility/NEUTRAL_FIBRE_LIFTING.md` §7.
For m slots choose m distinct prefix addresses, and pin every prefix owner.
Each slot has three actual full kernel records, obtained by rotating the
final C3. These are three lifts of ONE resource per slot, not three units.

Equate the actual first-final-C3 owner along a path of m slots. The whole
assembly has exactly three constant phases. Its factor graph is a tree;
bags contain only two slots. Nevertheless changing the shared phase changes
all m slots and exactly3m owner-clock occurrences. Current charge is constant
because every slot's target owner class is unchanged under the literal source
automorphism. This does not transport an order-eight electrical potential
to the different source order.

With phase0 pinned in the source arrival and phase1 in the target arrival,
formula (3) returns EXACTLY3m. Thus bounded separator width controls compiler
size, not event size. This is the native analogue of the lost-cost hostile,
now visible in a general multivalued separator compiler rather than only a
binary implication quotient.

## 7. Causal and physical scope

At a new core arrival, reveal its factor tables before selecting the new
section. Equations (1)--(3) retain all resulting joint constraints. Instant
physical destruction of a currently occupied resource is a different model.
Later return to a released reserved option is allowed; simultaneous reuse is
not. A requested final helper return requires that exact full record to pass
the final factors. Charges alone cannot repair an incompatible final pin.

No theorem here supplies missing native record domains, proves all local
ports extend, or turns atomic replacements into native metagraph microedges.
The gain is a compositional section certificate and a constructive exact
boundary-cost compiler for arbitrary finite full-record factors of bounded
separator width, including higher-arity, multivalued and resource predicates.
