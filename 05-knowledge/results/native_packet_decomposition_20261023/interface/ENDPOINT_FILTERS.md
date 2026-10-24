# Endpoint mixtures are charged filters, with actual resources retained

**PROVED, scoped.** This is an exact representation of the interval between
two supplied feasible native assignments, not a construction of an unknown
second assignment or a classification of the unrestricted metagraph.

## 1. Input and the loss boundary

Fix slots, their actual source frames, and a fixed background. Supply two
feasible assignments A and B. A local record contains the literal path,
kernel/owner decorations, all shared-coordinate values, and the complete set
of normalized capacity-one resources it occupies. Converse and inverse-target
resources must be included. Distinct raw kernel representatives of the same
normalized resource are NOT extra capacity. Retain their raw records anyway.

The only changing nonadditive constraints allowed here are pairwise
compatibility predicates, such as equality of shared actual owners and
disjointness of normalized resource sets. Higher-arity predicates require an
additional retained factor; they do not become implications automatically.
Local pins and compatibility with the fixed background hold for both supplied
records. Unchanged slots can be absorbed into that background.

For each remaining slot i choose A_i or B_i, encoded by b_i=0 or1. Pairwise
compatibility accepts 00 and11 because the two full endpoints are feasible.
Its only possible prohibitions are therefore 10 and01. Forbid10 by the
implication i -> j; forbid01 by j -> i. Both prohibitions identify the bits.

**Endpoint-filter theorem.** Contract the strongly connected components of
this implication graph. Its condensation is a partial order, oriented so an
arrow points toward a required successor. Pairwise-compatible endpoint
mixtures are exactly upper sets of this order, with all bits of a contracted
component selected together.

Proof. A selected bit forces precisely its implication successors to be
selected. Mutual reachability forces equal values; after contraction the
remaining condition is precisely upward closure. Conversely every upper set
avoids every forbidden mixed pair. This proves both directions. No port or
owner quotient beyond the explicitly retained implications is being taken.

## 2. The electrical and margin factors remain higher-arity

Let the additive observable take values in an abelian group H, with exact
integer coordinates included as Z factors. Set

    q_i = observable(B_i) - observable(A_i).

For aggregate electrical guards, the observable is the ENTIRE physical-row
current residual, with the needed original character coordinates. Source and
target contributions meeting at the same owner must be combined there. A
global scalar checksum is insufficient in general. Exact degree/target
margins are also coordinates. Subtract the fixed background once only.

Both endpoints have the same required observable iff sum_i q_i=0. In this
case the legal endpoint mixtures are exactly the upper sets U satisfying

    sum_{i in U} q_i = 0 in H.                          (1)

The charge of a contracted component is the sum of its member charges. The
finite component order plus these charges is an exact decision representation
of this endpoint interval. For reconstruction keep every original full record,
normalized resource key and slot-clock coordinate as a sidecar. The partial
order alone discards the electrical factor; the charges alone discard the
port/shared-owner dependencies. Neither alone is sufficient.

This does NOT contradict Session45's 196-state mod8 family whose pairwise
projections are all full. Equation(1) is a genuinely joint factor, not a 2-SAT
encoding of that family.

## 3. A finite monotone packet decomposition always exists

Start from U=empty. If U is not the full order, the full complement is an
available increment: adjoining it reaches the feasible endpoint B and has
zero charge. Choose an inclusion-minimal nonempty increment D such that
U union D is upper and sum_D q_i=0. Perform all changes in D simultaneously.
Repeat. Each step reaches a feasible full-record assignment and changes at
least one previously untouched slot, so at most s steps suffice for s slots.
Each slot changes at most once. Resources may be held at successive times but
are never occupied twice at any arrival; obsolete and replacement records are
updated atomically, not superimposed as a fictitious extra inventory.

The irreducible object is therefore a **zero-charge upper-set increment**,
relative to the current residual order. This constructive finite statement
gives no small or size-independent support bound. The whole remaining set
can be the only possible increment. Nor is the decomposition generally unique.

For a finite accepted-set table, minimum monotone bottleneck cost is given by
the acyclic recurrence

    B(empty)=0,
    B(V)=min_{U proper subset V, U legal} max(B(U), cost(V\U)),

over legal upper sets V, with the usual unreachable value infinity. Cost is
nonnegative and can
count changed slots or sum the supplied actual changed-clock costs. This is
exact for monotone endpoint-mixture routes. An unrestricted route may leave
this interval and can have a smaller bottleneck; equality is not asserted.

## 4. Native examples and the new repair

- In the earlier eight-slot packet the two literal shared-owner equations
  force the component {0,4}. The other six bits are independent before the
  current test. The four accepted masks give the independent two-slot and
  six-slot packets, with no discarded indexed owner.
- In the earlier sixteen-source circuit the pairwise order is an antichain,
  but equation(1), including the full dyadic bank, accepts only empty and full.
  Its monotone endpoint packet is genuinely of support16.
- The present [native repair](../native/README.md) retains that same interval
  obstruction. Extra source slots start and finish unchanged; hence they are
  absent from the original endpoint interval. Letting them temporarily change
  creates nonzero compensating charges at the SAME target owners and yields a
  smaller-support route through the enlarged state space. This is precisely
  why the old atomicity statement and the new repair can both be true.
- A shared-owner cycle can contract all s slots into one component even when
  every q_i is zero. Conversely an independent supply can have no implication
  edge and still need a large minimal zero-sum packet. These are different
  mechanisms. The [quantitative charge lane](../algebra/README.md) bounds only
  the latter under its explicitly independent-supply hypotheses.

## 5. Connections and the next obligation

Source: full native records. Target: a charged finite partial order. Map:
select before/after records, compile forbidden mixed pairs, sum row charges.
Preserved predicate: exact feasibility within that endpoint interval. Lost
information: alternate local paths, excursions through extra slots, and the
record geometry. Sidecar: literal records, normalized keys, shared owners,
kernel lifts, converse/inverse partners and integer clocks. Cheapest decisive
test: enumerate accepted endpoint mixtures and compare with the compiled
filter table, retaining the electrical equation as a joint factor.

This gives a precise interface between shared-pin assembly, source-charge
circuits and vacancy workspace. A universal native wiring theorem still needs
a complete, effectively supplied family of local records and compensators
that crosses every surviving interval obstruction while preserving the full
ambient current bank. This note does not supply that family.

## 6. Provenance and executable checks

The implication mechanism specializes the earlier
`injective_native_blocks_20261023/assembly/BINARY_ASSEMBLY.md`; full-resource
and per-owner charge discipline is inherited from Session45's
`injective-completion/THEOREM.md` and `dyadic-constraints/THEOREM.md`.
The added content is the exact endpoint interval, its retained joint charge
factor, monotone decomposition and the compensator interpretation above.

`check_filters.py` exhausts all4096 four-bit pairwise compatibility systems
accepting both endpoints, crossed with all8 mod2 charge assignments of total
zero. It separately checks literal resource/owner records, exact integer
charges and the two inherited native accepted-set tables. It compares literal
constraints with a reachability-based compiler and a separate subset dynamic
program. These are finite controls of a proof, not evidence for unrestricted
all-size coverage.
