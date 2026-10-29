# Recovery, correction precedence, and the live concept board

**Source recovery and scoped synthesis.** Historical titles are not truth
claims; current corrected mechanisms have precedence.

## 1. March enumeration: retain the admissible symmetry universe

The recovered March16 experiment
`04-computation/fast_a000568_v2_s90di.py` groups equal cycle lengths and
precomputes gcd interactions. Those techniques remain useful, but both
functions enumerate **all** partitions. Their outputs at n=2,3 are2,4,
the graph counts, instead of the tournament counts1,2.

The following day's `04-computation/a000568_turbo.py` already supplies the
repair: only odd cycle types can fix a tournament. The current session
propagates a maintained warning to the older script and records the
lineage in `01-canon/MISTAKES.md`. It does not claim to newly discover
the odd-cycle filter. The checker executes AST-extracted original function
definitions, avoiding the old script's long top-level benchmark.

The mechanism matters for the new precision filter: a permitted **number
of cycles** does not imply prime or small **cycle lengths**. For example
the one-cycle term at odd composite n must remain. Repeated parts also
carry a factorial denominator whose binary part controls the first
correction. These are two separate admissibility guards.

## 2. Older tail peeling and current response recurrences

[THM-514, unit-tail core kernel and prime-tail ladder](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
and [THM-516, Burnside core-kernel phi reframe](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md)
already contain tail elimination and gcd interaction laws. The
[October20 future-response packet](../tournament_future_boundary_20261020/README.md)
compresses histories by responses to still-possible insertions and expires
coordinates beyond the remaining degree. The new agenda contracts the
identity transitions already implied by that theorem; it is not a new
future-equivalence relation.

[Session87](../session87_four_vertices_enumeration_20261028/FAST_COUNTS.md)
adds degree buckets, packed integer responses and the terminal fold. It
also honestly benchmarks March Python and C/GMP implementations. The
current matched comparator is its faster `degree_bucket_dp`, not a slower
historical implementation or a richer symbolic task.

[Session98's fixed-degree projective quotient](../session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md)
retains h_j-j h_1 and absorbs2^(r h_1) at a **fixed** remaining mass r.
Its exact individual-continuation minimality is already proved.
[Incoming Session105](../session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md)
makes its negative binary shifts exact in an integer scale. It improves
the Fraction implementation substantially but does not uniformly beat the
whole-prefix method. This changed the present direction: retain that
quotient and seek a different requested-output reduction, namely bounded
binary precision, instead of claiming another smaller exact-degree state.

## 3. Session89's first bit becomes a finite-precision recursion

The decisive inherited sidecar is
[Session89's low-bit theorem](../session89_consecutive_tournament_pairs_20261028/counting/LOW_BITS.md).
It isolates one-cycle terms at odd n and coprime two-cycle terms at even n.
P2 keeps its valuation bounds at arbitrary precision, where only bounded
cycle number survives. P3 then resolves the first correction sharply by
retaining the repeated-unit triple and its2! denominator.

The connection contract is explicit: source=positive Burnside summands;
target=normalized residue; map=odd-denominator inversion modulo2^b;
preserved predicate=the exact requested bits; loss=all discarded high
valuation contributions; sidecar=cycle multiplicities and the normalization;
hostile=the invalid odd-k factorial bound at k=4. This is a projection of
an arithmetic observer, not a classification of individual tournaments.

June's THM-516 gives a second accepted map: its nonnegative totient
collision energy becomes P2's finite precision budget. Only bounded
common divisors can occur in a retained tuple at fixed precision, even
though each individual odd cycle can be arbitrarily long. This is why
the cutoff must not become a fixed wheel that deletes large prime cycles:
the one-cycle contribution1/n is the leading term at every odd n.

The aligned-pair channel requires independent weights. Session89's
same-response cores(3,7) and(5,5) have different energy/centralizer data;
squaring an already-summed ordinary weight is invalid. P3 works directly
with F_n(2,1), preserving this correction.

## 4. The root-switch map supplies an actual metagraph consequence

[Session79's parent-only root switch](../session79_tournament_enumeration_20261027/ROOT_SWITCH.md)
and [the order-ten bank](../consecutive_flags_order10_20261028/README.md)
connect counts to actual pointed states. A root-switch block records all
vertex orbits of one child, with orbit-size weights, and has trace equal
to its order. P4 turns the sharp channel correction into two universal
congruences on that specific reset operator and its pointed-state count.

These checks can reject a faulty aggregate result cheaply, but cannot
certify a whole representative bank: changing a count by the tested
modulus evades the check. They also do not describe native path complement.
The literal full-bit universes through5 and inherited complete bank through10
are separately labelled in the receipt.

## 5. Previous precision, phase and realization lessons

The [12/22 clock packet](../twelve_twenty_two_trace_20261029/POWER_CLOCKS.md)
separates clock-generator recovery from readout recovery. The counter is
unimodularly observable, while arithmetic trace has a ramified7 kernel.
Here the observer is an explicitly requested residue of a global count;
precision truncation deliberately changes what is recoverable. It does
not produce a universal free attachment coordinate.

The [overnight periodic/terminal observers](../overnight_rank_return_20261029/TERMINAL_SELECTOR.md)
show that new terminal operations can expose amplitudes without changing
the clock polynomial. Likewise, the new channel residue cannot be reused
unchanged for a different r,s or a structural port query. The binary
precision parameter is not a Frobenius phase; there is no field-theoretic
identification merely because both use powers of2.

The81/27 controls are replayed with their actual groups and transporter,
not only their group orders. Their common cyclic-triple set cannot repair
the transporter3 to7, which has27 elements in one and none in the other.
This preserves the actual attachment boundary of all global counting
shortcuts in this packet.

Incoming `atomic_tail_reset_20261029` distinguishes a native return clock
from an unpaid source reset. Its precise relevance is the operation
boundary: a forced-zero insertion is a proved identity with a capacity
certificate; an unproved arithmetic or metagraph reset is not. No Collatz
claim follows from the new schedule or precision recursion.

## Concept-board update and stopping reasons

| Concept | New evidence | Decision |
|---|---|---|
| Future gcd responses | Exact forced-zero contraction | Keep the same response quotient; save only unnecessary visits |
| Fixed-degree gauge | Session105 already supplies integer realization | Do not rediscover it or claim prefix minimality from it |
| Binary precision | Uniform cycle-number cutoff and sharp triple correction | New anchor: exact partial-output recurrence |
| Root transport | New trace/pointed congruences pass literal and inherited controls | Add arithmetic checks, retain full root records |
| Joint realization |81/27 transporter mismatch survives all scalar summaries | General strong-atom closure stays open |

The used research cards are “Keep observer, recurrence class, and finite
head” and “Search the statement before the method.” No new universal
meta-pattern is promoted from one counting application. The scheduler's
small timing differences are a recorded stopping reason, not a failed
correctness theorem. The next useful arithmetic target is a fast exact
three/five-cycle convolution; the next structural target needs a native
joint-response valuation bound before it can inherit this technique.
