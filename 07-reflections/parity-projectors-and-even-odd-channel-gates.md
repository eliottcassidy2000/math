# Parity Projectors And Even/Odd Channel Gates

codex-2026-06-13

The prompt's slogan has a clean algebraic core:

```text
midpoint scalar gate -> odd channels survive
reversal tournament gate -> even channels survive
```

The mechanism is not mystical.  It is the two projectors for an involution:

```text
P_+(f)=(f+fJ)/2
P_-(f)=(f-fJ)/2
```

For scalar midpoint pairing, the involution is `z -> -z`.  Pair-differences
keep odd powers.  The Faulhaber balance has one extra fixed atom `c^p`, then
only odd moments `S_1,S_3,...`.

For tournaments, encode every arc by a sign.  Converse sends every sign to its
negative.  A Walsh monomial of degree `k` changes by `(-1)^k`, so
complement-invariant functions can only see even Walsh levels.  This is the
exact version of the old HYP-1387 fact that `H(T)=H(T^op)` kills odd levels.

## The Marked-Viewpoint Lesson

The best computational finding was not merely that `H` is even and writhe is
odd.  It was this:

```text
start0(T) = end0(T^op).
```

So a rooted perspective is not naturally even or odd.  Reversal transports it
to its dual perspective.  Only the paired channels are honest:

```text
start0 + end0 -> even
start0 - end0 -> odd
```

That feels directly relevant to the observer-blind versus observer-coupled
LRC split.  If an LRC statistic has a source/sink, left/right, start/end, or
owner/co-owner dual, we should split it before throwing it into a quotient.
Otherwise we average away the channel that might carry the proof.

## Why This Helps LRC14

The current Q27 program has many clocks.  They should not all be compressed in
the same way.

```text
even scalar:
  denominator shell, total reset period, observer distance, unmarked blocked q

odd marked:
  owner support, carry residue, deletion derivative, pair sum/difference cut

transported:
  source/sink, left/right endpoint, start/end observer perspective

compatibility:
  OCF alpha packets, Q27 owner/carry packets, code72 support incidence
```

This suggests a practical algorithm.  First use even scalar gates to shrink the
state space.  Then split transported viewpoints into sum and difference.  Then
use odd marked channels to prove pressure or find a deletion/opening.  Finally
attach OCF-style compatibility packets before trusting an atom inventory.

In short:

```text
even channels compress;
odd channels witness;
transported channels must be paired;
compatibility channels decide liftability.
```

## Broad Synthesis

This same pattern explains why so many repo threads kept saying "scalar shadow"
and "hidden lift."

In code72, the Type II weight enumerator is an even scalar gate; the hard part
is the support/design lift.  In unit distance, the metric equality graph is an
even scalar gate; the unit spine and endpoint ears are marked support data.  In
polynomial irreducibility, sign-cube and residue summaries are scalar gates;
convolution factor grids and Newton slopes are hidden lifts.  In OCF, odd
cycle counts are atom inventory; `alpha_k` compatibility packets make the
Hamiltonian-path count.

The new connection is that parity tells us which data are safe to quotient and
which data are proof-bearing.  The midpoint odd channel and the tournament even
channel are not rivals.  They are two faces of the same projector calculus.

## Next Object

The next useful LRC artifact should be a typed Q27 ledger.  Each field should
declare one of:

```text
even_scalar
odd_marked
transported
compatibility_packet
```

Then AP, `Vstar`, `2AP`, one-stranger rows, and multi-stranger rows can be
compared after the same protocol.  My guess is that the remaining hard rows
are exactly the rows where the even quotient looks closed but an odd marked
owner/carry channel has not yet been forced to choose a side.
