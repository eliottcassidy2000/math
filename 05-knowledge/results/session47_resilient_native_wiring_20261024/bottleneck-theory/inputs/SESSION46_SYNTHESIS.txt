# Session 46 synthesis: charge atoms, charged buffers and native context

This session moves beyond the old sixteen-source example, but its more useful
result is a separation of different obstructions. A packet can be indivisible
inside its original pins, decomposable after adding actual charged resources,
and still resistant to a smaller move alphabet even when that alphabet spans
the entire integer charge lattice. These are different mathematical questions.

All new finite native certificates below retain the original marked order-eight
network. Their enlarged supply means additional actual source options from that
network, not an invented copy of its ports or an unproved transport to another
order. The conditional general theory and inherited higher-order comparison
are identified separately. Dyadic preservation means retention of the original
nodewise harmonic current banks by inclusion, not equality of kernels or
preservation of every ordinary individual-port harmonic phase.

## 1. The expanded supply is complete, not a selected binary fiction

The 28 inherited source options have 57 actual Hamiltonian paths after every
common clock/owner pin is fixed. Twenty-seven sources have exactly two paths.
Source 1923 alone has a third path, targeting 1559. When the original target
margins require only U=1965 and V=2234, that third path is impossible because
its nonnegative target count cannot cancel. The resulting complete margin-fixed
supply has exactly 56 options, one U path and one V path at each source.

The full raw supply has 57 distinct normalized black resources and 228 distinct
closed ports. All source automorphism groups are trivial, checked from the
original matrices. Actual source/path/owner maps and converse are retained;
no kernel lift is counted as extra capacity. Declared shared pins are preserved,
without falsely asserting that all source matrices form one common host.

[Complete pin domains](pinned-domains/THEOREM.md).

## 2. An eighteen-source atom, with a sharp 37-clock minimum in this supply

There is a complete pinned eighteen-source circuit with exactly two states
preserving the full original nodewise dyadic bank. The first certificate uses
43 changed owner clocks and has 792 accepted mod-8 states but only two at full
precision. Direct full-Laplacian integer-image witnesses authenticate the
accepted endpoint difference, and genuine harmonic separators exclude every
other balanced choice.

A stronger circuit has 37 changed owner clocks. Its positive source set is

    {1394, 1425, 1457, 1480, 2176, 2185, 2581, 2835, 2945},

and its negative source set is

    {1477, 1549, 1691, 1921, 1974, 2837, 2905, 3000, 3369}.

There are 262,144 complete binary assignments, 48,620 with the required target
margins, 754 full-mod-8 states and exactly two fully dyadic states. Each endpoint
uses 72 closed ports; their union has 144 distinct ports. The chosen integer X
current debt is exactly zero for this circuit, while the 43-clock example has
opposite debts of magnitude 512. Neither fact replaces the genuine harmonic
scaling beyond the original exponent 256.

Every changed permutation costs at least two owner clocks. A 36-clock
18-source circuit would therefore use only the twenty cost-two options, where
the previous exact maximum primitive support is 16. An independent new census
also excludes every such 36-clock primitive relation. At cost 37 there is exactly one
primitive relation up to reversing all signs, the one displayed above. Thus 37
is the sharp minimum for eighteen-source primitive circuits in this fixed
28-option, one-slot-per-source supply. It is not a global minimum for every
order-eight source menu or every target pair.

[First 18 certificate](larger-circuits/THEOREM.md),
[sharp 37-clock extension](larger-circuits/cost-optimal18/THEOREM.md),
[independent audit](audit-larger/cost-optimal18/AUDIT.md).

## 3. A real S2/S3 connection, with its exact obstruction

The minimizing circuit uses seventeen two-owner transpositions and one actual
three-owner cycle. At source 2835, clocks 4, 5, 6 change their owners from (6,4,0)
to (0,6,4), giving the relative permutation (0 4 6). This is an authenticated
permutation statement, not a connection inferred from the number 3.

Applying that same three-cycle a second time produces a permutation which is
NOT a Hamiltonian path: an actual required arc fails. Hence an ambient
order-three permutation does not give a three-state legal native orbit or an
S3 action on the permitted paths. The domain is still binary. Native legality
is a condition on partial transformations of presentations, not automatic
closure under every ambient group operation.

This also supplies exact generating functions. Let P_h(z,t) count accepted
relative flip sets at dyadic precision h, with z marking changed slots and t
marking changed owner clocks. Complementing the flip set gives

    P_h(z,t) = z^18 t^37 P_h(1/z,1/t).

At full precision P(z,t)=1+z^18 t^37. At every tested precision, exactly half
the accepted states have odd clock cost, because the only odd-cost slot is
2835 and complementation toggles it. In particular the 754 mod-8 states split
377/377. The univariate specialization P_h(1,t) has factor 1+t; that factor is
not asserted for the full bivariate polynomial. The coefficient tables and
literal failed iterate are independently checked.

## 4. Charged borrowing can decompose a primitive atom

The word 'primitive' is relative to the supplied pinned slots. It does not
forbid adding genuinely different charged slots and returning them afterward.

For the original sixteen-source circuit:

- One cost-two outside source can reduce the sharp largest event from 16 to 8.
- Source 3015 alone, with a cost-three alternative, permits six-slot events.
- Two cost-two buffers 1394 and 2905 also permit six-slot events.

The eighteen-source 43-clock circuit admits a six-slot decomposition using
buffer 1394. The sharp 37-clock circuit admits one using buffer 2071. Every
compensator returns to its exact initial record. Full source/path matrices,
all fixed owner coordinates, normalized resources, converse and actual current
rows are checked at every intermediate state.

Six is sharp for these endpoint pairs throughout the complete 28-source supply,
not merely an observed route size. All 27 possible four-slot trades preserve
b_1638+b_1974 as an INTEGER, whereas these endpoint pairs change it. No number
of such four-slot moves, and no initial arrangement of the available outside
slots, can overcome that invariant. Distinct labels exclude two-slot trades.
The six-slot constructions attain the lower bound.

Clock cost and event count are separate. In the complete 38-state menu for the
37-clock circuit with buffer 2071, four events and total work 49 are both optimal
at the sharp six-slot budget; one route attains both and peaks at 13 clocks.
The sharp peak-clock budget is instead 12, which requires at least five events
and at least 57 total clocks in this menu; a different route attains both.
These event and work minima are not claimed across all possible larger menus.
Thus reducing peak work need not reduce total work. The
[exact route trade-offs](compensated-decomposition/variants/THEOREM.md) and
[independent audit](audit-compensation/VARIANT_AUDIT.md) retain these distinctions
rather than calling a smaller simultaneous event an unqualified speedup.

These statements use capacity-one reversible borrowing. If a released port
may NEVER be occupied again, a changed buffer cannot return, and an original
binary slot cannot flip twice. Every permitted route is then conformal and
the original primitive 16/18 atom must move as a whole. Both no-reuse conventions
are addressed explicitly in the [guard and reuse contract](GUARD_AND_REUSE_CONTRACT.md).

## 5. Why full integer cancellation still need not make a route

In the same 28-source supply, the 27 four-slot trade vectors have rank 23. Adjoin
all 1,917 six-slot trade vectors. Their integer span is now the entire zero-margin
paired-charge kernel, with index 2048 in the ambient zero-margin lattice.
This is proved by direct integer combinations of 48 actual moves producing a
27-dimensional basis with determinant 2048, together with the complete original
group-image calculation. No additional additive
conservation law remains beyond target degree and the marked group charge.

Nevertheless six-slot routing fails in some complete binary fibers. At exactly
four V slots and charge (parity 7, cyclic 15), the only possible selected sets are

    {176, 1425, 1549, 2087}
    {1394, 1923, 2837, 3369}.

All 28 slots are allowed to vary, not just these eight. There are still only
two legal states. They differ in eight slots, so the sharp event budget is 8.
The integer lattice expresses their difference using small moves, but that
expression cannot be scheduled inside this capacity-one fiber.

The complete cardinality-four census has 20,475 states across 2,048 fibers.
Under events of at most six slots, 1,851 fibers are connected, 175 have two
components, 21 have three and one has four. Thus 197 fibers expose the issue.

This is the distinction between an algebraic generating set and a generating
set valid in every actual context. A binomial multiplied by a repeated port
need not be a legal native move. The correct state graph retains the occupancy,
shared-coordinate and guard predicates at both ends of every event.

[Integer lattice and native capacity boundary](move-lattices/THEOREM.md),
[independent audit](audit-move-lattices/AUDIT.md).

## 6. A conditional all-size decomposition rule

Take an authenticated binary supply whose two common targets have equal
marked paired charge, whose pins are jointly compatible, and whose alternatives
have pairwise disjoint normalized resources avoiding the fixed remainder.
Assume the only changing acceptance conditions are the exact target margin
and complete nodewise dyadic current equations. Every accepted endpoint
difference then partitions into primitive balanced equal-charge subtrades. They can be executed
in any order, each changed slot switching once, using only endpoint resources.
This is a conformal decomposition, not a promise that each primitive atom is small.

Pair the positive and negative slots of an atom. Their label differences form
a minimal zero-sum sequence in the finite paired-charge group G. Therefore an
atom has at most 2D(G) slots, where D is the Davenport constant. The classical
p-group formula, proved here by a self-contained augmentation-ideal argument,
gives

    D((C2)^4+C128)=132,

hence the conditional bound 264 for this exact order-eight common-target type.
That bound is not claimed sharp for native packets, is not the bound for an
arbitrary multi-target interface, and is not constant across tournament orders.

There is also a precise charged-buffer supply criterion. Match an atom into
pair differences d_j and take their partial sums s_j. If actual, jointly
compatible outside slots exist with labels h-s_j and distinct available
resources, one selected buffer can circulate through these slots while each
original pair changes. Each event has at most four slots and the buffer returns.
The missing native labels, pins or resource availability can prevent this
construction, as the current 28-source four-slot invariant does. The theorem
does not manufacture those resources from a formal group expression.

For general overlapping physical roles, shared owners and multiple targets,
retain all rows of the actual incidence change, applying its paired multiplicity
exactly once. A lossless finite one-hot encoding
can include actual record choices, resource capacities, target margins,
shared-owner values, full kernel lifts and complete residual factor tables.
Conformal Graver atoms of its exact equality/congruence lattice then give a
finite endpoint-union decomposition. The encoding may be enormous, and an
atom may involve every slot. This is an exact conditional rule, not an efficient
universal compiler or a classification of every tournament metagraph.

[Conditional decomposition theory](decomposition-theory/THEOREM.md),
[proof review and independent controls](audit-theory/AUDIT.md),
[final synthesis and incidence review](audit-synthesis/AUDIT.md).

## 7. Cohesive explanation of the obstructions

There are three different graphs in play:

1. The tournament supplies actual Hamiltonian paths and legal owner permutations.
2. The presentation metagraph supplies the original marked harmonic current guards.
3. The guarded assignment graph supplies legal simultaneous moves and remaining
   connectivity components once pins and capacity are enforced.

An automorphism, harmonic equality or lattice relation in one graph does not
by itself provide a path in another. Our examples locate the missing data:

- A second application of the permitted three-cycle leaves the Hamiltonian domain.
- Identical charge labels can represent incompatible or already occupied ports.
- Pairwise constraints or a full integer lattice can miss a binary state component.
- Neutral workspace can release vacancy constraints but cannot pay a projected
  current defect; charged workspace must carry and return that exact debt.
- A larger legal move can have fewer events yet higher peak clock cost.

The [inherited guarded-coverage construction](../guarded_coverage_20261023/PORTS.md)
gives a complementary all-size native warning: a target-fixed Gray-cycle
assignment can require 2^r simultaneous slots despite unchanged class-incidence
currents. Its alternative resources
are shared between slots, violating the globally disjoint-alternatives hypothesis
of the common-target decomposition theorem. Its authenticated four-layer repair
changes the actual resource inventory by enlarging the source and releasing
the stated auxiliary pins. It applies when injection is the only changing
cross-slot constraint, not to arbitrary shared equalities or original-network
electrical debt. This is a proof-level comparison with already inherited work,
not a claim of newly received commits or a fresh replay of all its native
certificates. No contradiction is hidden by calling both phenomena 'connectivity'.

There is also a precise limit to a matroid-exchange interpretation. In a fixed-cardinality,
fixed-charge fiber with distinct group-charge labels on the slots, a one-for-one exchange
would require the entering and leaving labels to be equal. None is available.
Therefore any non-singleton such fiber fails the basis-exchange axiom for a
matroid on these same labelled slots. This does not exclude recoded or generalized
matroid constructions; it explains why ordinary basis exchange cannot supply
the native routes claimed here. An event in this session means simultaneous
replacement by the specified actual paths, not an unverified sequence of
single-arc or single-owner changes.

## 8. What remains open

The complete all-size wiring classification remains open. The new work provides
larger certified atoms, an exact finite cost minimum, charged decompositions,
explicit lower bounds and a lossless conditional decomposition framework.
It does not prove native supply at every size, uniform small event width for
all guarded fibers, or globally injective replacement of a fully occupied atlas.

The next useful target is an authenticated supply theorem for compensated
interfaces, with actual shared coordinates, available normalized resources and
integer clocks. Rank and charge counts are effective filters, but the native
context is what decides whether the proposed mathematics can actually move.
