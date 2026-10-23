# Independent final audits of the concurrent root packets

October 22, 2026. These were read-only audits from the order-15 lane. The
reviewer changed only this note in the owned folder. The packet owners made
the interface repairs described below before final replay.

## Delayed pulse replacement — PASS

Reviewed [the complete proof](../lagged_tail_pulses_20261022/README.md),
production consumer, supplied native paths and all-future interpretation.

For an authenticated frontier-avoiding prefix from pending node $u$ to
pending node $v$, the two first-arrival times differ by the **native edge
length**, while the transported sign uses the **summed integer clock**.
If the unknown tail of $v$ never reaches the frontier, neither owner arrives.
Otherwise their separate old/new vectors cancel after the second arrival
and equal the early vector during the intervening interval. This is a
four-coordinate identity, so cross terms with every other owner are retained.

Allowing the resulting pulse to be arbitrarily present/absent at a single
cut is a sound relaxation. It is not a new synchronous native snapshot and
does not inherit the ordinary state-update monotonicity theorem. Disjoint
pending owners prohibit spending an original source fibre twice; intersecting
known paths do not themselves violate ownership. The input prefix, mod15/mod35
source registers and actual clocks remain bound to the inherited receipt.

The initial consumer accepted floating aliases for the final path target,
edge clock and summed clock. The owner repaired these with literal-integer
guards and added validation of the input state before removing a late node.
All three independently supplied regression examples now reject. This was
an interface correction, not a change to the pulse theorem or native rows.

Independent normal and optimized saved replays both pass **15,017** checks:

`c0be85b63ccb3c00d7263b7180c63a88f3c4fcc6d5c4225728578a189c90a4ef`.

The finite positive bounds are 521233, 1039528 and 3856740; the N7 hostile
remains 1153 above its 1071 budget. No uniform source-entry conclusion follows.

## Linked signed-cycle capacity — PASS

Reviewed [the scoped proof](../incidence_tail_bridge_20261022/README.md),
all-cycle certificate consumer, edge ownership and exact mixed controls.

A negative signed cycle of allocated weight $p$ loses at least $2p$ from
the unsigned edge maximum. If all its vertices are present at one terminal,
some edge must be frustrated. If absence or different destinations break the
cycle, at least two edges disappear, giving the same bound. Equivalently,
the proof's independent random signs on terminal groups and omitted vertices
recover the actual pair objective in expectation. Summing cycle allocations
is valid exactly when every edge's total allocation is within its absolute
capacity. Shared vertices are permitted; a shared edge cannot be overspent.
Previously consumed cluster edges must remain excluded.

The coefficient is two because the kernel already counts each unordered
pair once. The proof does not require the linked four-dimensional form to be
positive. It does not authorize subtracting the same saving from another
optimized angular certificate. The exact triangle-versus-longer-cycle
control and the all-positive-lambda mixed comparison are sound.

The initial reader accepted numeric aliases for cycle node IDs and the
claimed saving. The owner added literal-integer owner/saving guards and
regressions. Independent node/saving alias controls now reject. The corrected
relative theorem link was also reported to the owner.

Independent normal and optimized saved replays both pass **697,921** checks:

`5e73e0f3a4fec33cbea799e3c042c08f7ffed2ef59493ae92597934de0aa9744`.

The native cycle savings are 5320, 19526 and 99200. The last is a bounded
greedy certificate, not a packing optimum. The independent exact angular
bound remains smaller in that case.

## Other reciprocal reviews completed in this lane

- [Center incidence](../center_incidence_20261022/README.md): full proof/API
  and optimized saved replay **PASS**, 19,331 checks,
  `13effe30119b6397fb58f3e247aa9e53f1bc292cbc91bcb8c6c404f89ff56b1a`.
  Reviewed the anti-involution orbit pairing, singleton-necklace family
  count, split-run separation, proper-shore deletion, and two equal odd-orbit
  obstruction. Its seven literal sources are controls, not censuses.
- [Weighted tail](../weighted_tail_20261022/README.md): full proof/API and
  optimized saved replay **PASS**, 19,435 verifier checks,
  `2d17a1124daa2c4859c2a3bee039da8ad728e667cac918b5e93b167113fdac2d`.
  Reviewed exact unary allocation, optional pair choices, shared-terminal
  fields, convex endpoint certificates, the fixed-parameter LP dual and
  typed-pulse composition. The resulting bounds change a relaxation, not
  physical energy, original sources or already-owned pair credits.

These audits concern the displayed elementary bridges and exact consumers.
They do not independently re-prove cited solvability/Hamiltonicity inputs,
the supplied incidence papers, or the unrelated lattice/Dehn literature.

## Subsequent root boundary hardening (October 23)

After the reciprocal review, the pulse consumer was tightened to require
the supplied frontier to equal the state's literal terminal keys. A caller
cannot bypass an intermediate terminal by supplying an empty frontier:
the hostile path 7 -> 11 -> 17 is rejected when 11 is a terminal.
Nonliteral integer aliases are rejected too. The root's optimized replay
now passes 15,052 checks with the same mathematical receipt hash. This
additional API check was performed by the root after the independent audit;
it does not change the frozen pulse pairs, paths, clocks or numerical bounds.
