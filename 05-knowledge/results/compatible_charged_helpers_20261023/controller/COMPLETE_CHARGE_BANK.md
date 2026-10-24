# A complete charge bank lifts every constant-margin core route

Status: **PROVED**, conditional all-size full-record lifting theorem.
The order-eight specialization is **FINITE-EXACT** plus this theorem.
Unrestricted native wiring at every tournament order remains **OPEN**.

## 1. Interface: charge coverage, compatibility, and actual capacity

Let G be a finite abelian group. Supply n core slots and h additional helper
slots, each with two authenticated full records, called U and V. These names
refer to record choices, not an identification of their owner coordinates.
Retain source frames, literal paths, all fixed shared pins, kernel lifts,
normalized reversal/inverse resources, chosen local phases, and integer clocks.
Every entire choice under consideration must be compatible with the fixed
background and other selected records, with no simultaneous resource reuse.

For this theorem assume **all** constant-margin core choices and the helper
choices used below are mutually compatible. This is a substantive hypothesis;
Section 5 explains how to keep rather than discard its failures. Suppose the
remaining retained aggregate guards, after fixing the target margin, are
exactly equality of the additive charge

    sum_(i selected at V) lambda_i in G.                    (1)

Establish this equivalence using complete physical-row currents before using
the theorem. It is not inferred from a convenient checksum. Exact integer
guards not implied by (1) must remain additional coordinates or restrictions.

Fix a helper occupancy r, 0<=r<=h. A **complete charge bank** is a supplied
section

    sigma : G -> {J subset H : |J|=r},
    sum_(j in sigma(g)) lambda_j = g.                      (2)

The h slots ALL exist at every arrival: r use V and h-r use U. The unused
alternative of a slot is not spare capacity, and a kernel lift is not another
slot. Completeness of (2) is charge coverage; compatibility and resource
authentication are separate hypotheses. Put t=min(r,h-r).

## 2. Uniform lifting theorem

Fix a core margin k and a total charge q in G. Let Omega(k,q) contain every
accepted full state (I,J) with |I|=k, |J|=r and

    lambda(I)+lambda(J)=q.                                (3)

An atomic packet event simultaneously changes selected full records between
two accepted arrivals. Its support is the number of changed slots, rather
than the number of individually legal one-slot transitions.

**Theorem.** Under the interface above, Omega(k,q) is connected by events
with support at most 2+2t. More generally, every constant-margin core path
with consecutive supports at most c lifts to a guarded packet path with
consecutive supports at most c+2t. The prescribed full helper records at
both endpoints can be retained, not just their charges.

**Proof.** A constant-margin core path I_0,...,I_m needs helper charge
q-lambda(I_l) at arrival l. For every interior arrival choose
J_l=sigma(q-lambda(I_l)); at the ends use the prescribed J_0,J_m. Equation
(2) proves (3) everywhere, and the compatibility hypothesis makes these
actual feasible records. Any two r-subsets of an h-set differ in at most
2 min(r,h-r)=2t positions. Add the core support to obtain the bound.

For arbitrary core endpoints, pair the elements of I_0\I_m with those of
I_m\I_0 and exchange each pair once. This is a Johnson-graph path of
m=|I_0 symmetric_difference I_m|/2 events, each of core support two. Every
interior set still has cardinality k. If m>=1 the preceding construction
uses exactly m events and the prescribed helper endpoints. If the core
endpoints coincide, change only the helper configuration, in one event of
support at most 2t (or zero events if the whole state is unchanged). QED.

In particular the number of events is at most max(1,min(k,n-k)), since
m<=min(k,n-k). The support bound
does not grow with n. A helper may change again at a later time, but it is
never assigned to two slots or counted twice at one arrival. Temporary
helper charge is explicitly present and returns to its prescribed value.

This is capacity-one reversible borrowing. If instead a released resource
may NEVER be occupied again at any later time, a helper required to return
to its original resource cannot change. The complete-section controller
then does not apply; primitive endpoint atoms can remain indivisible. The
incoming `session46_native_packet_decomposition_20261024/
GUARD_AND_REUSE_CONTRACT.md` proves this temporal-nonreturn boundary explicitly.

For actual clock costs, let d_i be the distance between the supplied records
of slot i. Each event costs the sum of d_i over its changed slots. Thus its
cost is bounded by the changed-core cost plus the sum of the largest 2t
helper distances (or all h if fewer). This counts slot-clock occurrences;
shared physical coordinates are not silently treated as independent ports.
Converse transport uses the literal reversed paths and reflected clock
indices. It does not obtain a new capacity from the same resource orbit.

## 3. Native specialization: forty new helpers, eight-source packets

The protected core is the complete 28-source binary menu in
`native_packet_decomposition_20261023/native/domains.json`. Its old endpoint,
indexed-arc and additional common-clock pins are unchanged. The targets are
the merged order-eight classes U=1965 and V=2234. The original paired current
group is

    Q = Z/128 x (Z/2)^4,
    lambda(s)=(X(s) mod128,Y1(s) mod2,...,Y4(s) mod2).

The supply proof in `../supply/README.md` authenticates 40 NEW source classes,
their chosen fixed common clocks and complete U/V domains. It gives a section
of every one of the 2,048 charges using exactly three V-helpers. All alternative
normalized resources are disjoint from one another and from the protected
core reservoir. Their literal black conductance is two, authenticated rather
than inferred from a common-neighbor label. The full source automorphism
kernels happen to be trivial here; the general theorem does not require that.

Both targets have X=146 and the same paired Y parity. The inherited physical
current equivalence therefore applies to these newly authenticated edges:
at fixed margin, (3) is precisely preservation of all original dyadic
aggregate guards. The chosen original per-slot X phase mod8 is also fixed.

**Corollary.** For every core margin k=0,...,28 and every q in Q, the complete
68-slot family with exactly three V-helpers has connected guarded fibres
under events changing at most EIGHT sources. Each of its states uses all
68 slots and 68 distinct normalized resources. The complete alternative
reservoir has 136 resources / 544 closed endpoint words; each arrival uses
272 such words. The reserved forty helpers are not three helper slots.
The two largest old-slot clock costs sum to eight, and the six largest
helper costs sum to 42, so every event costs at most 50 owner-clock occurrences.

This is a selected packet family at tournament order eight, not a tournament
on 68 vertices, nor a modification of the fully occupied canonical atlas.
Any unchanged background must avoid the entire reserved reservoir. New helper
pins are explicitly chosen; no protected core pin has been released. Extra
undeclared cross-source owner identifications are outside this corollary.

The finite code preserves EVERY original dyadic harmonic current: the
authenticated two-primary critical group has exponent 256, so its paired
test saturates. See `injective_native_blocks_20261023/electrical/
ELECTRICAL_CODES.md` Section 2. This does not mean the same printed integer
X values can be tested at an arbitrarily higher modulus without constructing
the required harmonic lift or scaling. Exact integer X current and
odd-primary guards are not claimed. Nor does aggregate current preservation
imply every ordinary per-port phase. The genuine
character 4e_U separates U and V at mod8. Those distinctions are retained in
the supply audit. Exact integer-current variants require a further condition,
not a reinterpretation of the same 2,048-charge section.

## 4. Mechanism and the genuine improvement

The earlier 264-source bound used monotone zero-sum blocks of an independent
difference list. Here a real, preallocated helper bank absorbs EVERY residual
charge after a two-source core exchange. It replaces a large zero-sum block
by a short sequence with recorded intermediate debts. This trades forty
helper slots, reserving eighty alternative normalized resources, for a
uniform support bound of eight.

Special circuits already admit six-source repairs with fewer helpers. The
new assertion is coverage of every margin and current signature of this
protected core, including the complete four-occupied-source fibres where
six-source connectivity fails. It is not a claim that eight or forty is
globally optimal. The bank is inclusion-minimal under single deletion; its
global minimum cardinality remains open.

Late incoming `sixslot570_20261024/FIBRES.md` proves the original28-source
fixed-current fibres connected at sharp support8 and22 clocks WITHOUT these
helpers. That is the stronger bound for that smaller state space. The present
bank instead allows the core charge to vary through every k-subset while
retaining the full combined charge; its coupled shared-owner extension also
retains an additional changing compatibility factor. Do not advertise the
helper bank as a cost improvement over the incoming fixed-inventory rule.

All 72 available new-source labels together fail to cover Q by pair sums.
Thus a two-active-helper complete-section template is impossible in this
alphabet. This is not a lower bound against every adaptive six-source
controller. Three active helpers require C(h,3)>=2048, hence h>=25; the
authenticated construction has h=40.

## 5. Shared compatibility is a retained factor

The global compatibility assumption can be weakened substantially. For each
internally admissible core state I let A(I) be a supplied family of r-subsets
of the SAME h helper slots. Require each (I,J), J in A(I), to be a genuinely
compatible full assignment, retaining the immutable remainder, actual owners,
kernel records, local guards and capacity. Assume

    {lambda(J) : J in A(I)} = G for every I.                (4)

Only the allowed family varies with I; the full record assigned to each
helper bit is fixed. Unrecorded changes of a kernel lift at an unchanged bit
would have an additional cost and are not permitted by this interface.

**State-dependent section theorem.** Every connected authenticated
constant-margin core-move graph with edge supports at most c lifts, at each
total charge q, to a connected guarded packet graph with support bound c+2t
(and helper-only support2t). No common charge section across core states is
required. Use sigma_I(q-lambda(I)) at each arrival. The proof of Section2
applies word for word: both endpoints of each simultaneous event are feasible,
and any two r-subsets still differ in at most2t slots. Prescribed endpoint
helper configurations remain unchanged. QED.

Condition(4) is exactly the projection-completeness criterion: every supplied
core I extends to a feasible state with J in the SUPPLIED A(I), for EVERY
total charge q, iff that family's helper-charge image is G. For necessity
among all feasible helper configurations, A(I) must also be a complete list.
This does not say projection completeness is
necessary for connectivity of one particular q-fibre. For a Johnson core
graph, the support bound is still2+2t. For a shared-pin core that already needs
c core slots to change together, that cost is retained as c, rather than hidden.

Charge coverage alone is insufficient when changing core owners constrain
helper records. The companion theorem
`../compatibility/NEUTRAL_FIBRE_LIFTING.md` compiles each actual helper
configuration into an exact interval of old implication upsets, retains all
charged or heavy components, and eliminates only light charge-neutral ones.
Its boundary distance is the sum of costs of coordinates forced on in one
interval and forbidden in the other. This gives an exact bounded-connectivity
test, including simultaneous neutral changes; it is not a statistical proxy.

Our bank supplies actual charged vertices to that compiler. A permitted
section must then be selected inside the appropriate compatibility intervals;
surjectivity of the unrestricted subset-sum table does not guarantee this.
The native shared-C3-owner family in that proof forces arbitrarily large
events despite constant charge and disjoint resources. Thus the unresolved
all-size target is an authenticated supply with sufficient COMPATIBLE charge
coverage, or a variable-budget theorem that pays its forced boundary costs.

Source: full core/helper records. Target: constant-margin core exchanges with
a charge section. Preserved: all declared guarded arrivals and exact normalized
capacity. Forgotten by the charge projection: owner and kernel choices.
Required sidecar: every literal record, actual pins, resource orbit, clock
distance and compatibility factor. Cheapest hostile test: hold charge fixed
while a shared actual owner forces an entire neutral component to change.

## 6. Reproduction

`check.py` independently enumerates the 9,880 helper triples, replays the
saved complete charge section, compiles routes across all core margins,
and evaluates physical Laplacian rows from actual resource edges. It retains
complete native arrivals for the smallest eight-source hostile, including
their actual paths, normalized resources, shared clocks and converse paths.
Run normally and with `python -O`; explicit exceptions keep checks active.
