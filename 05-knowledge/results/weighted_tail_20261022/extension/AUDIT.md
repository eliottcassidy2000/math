# Bounded extension audit

Status: **FINITE-EXACT** consumer audit and **PROVED** finite transport
justification. October 23, 2026. This audit concerns only the extension
packet and its use of the inherited weighted gate, not the uniform stopping
conjecture.

The potential failure points were ownership duplication, past-frontier
transport, loss of an old channel during coalescence, and interpreting a
pulse relaxation as an ordinary synchronous native state. The consumer
starts only from the frozen ordinary prefix. It checks a complete bijection
between supplied owners and its pending-node rows; replays every edge; and
retains the full four-vector until a synchronous sum is zero in every
coordinate. Once an owner arrives, its contribution is rebuilt once at
each checked cut rather than repeatedly deposited. An independent aggregate
step agrees with this owner construction at every cut.

The final bound follows from the already proved weighted finite-state
theorem, applied to the exact final native state. The producer does not
decide that theorem's future assignments, assume arrival, or choose
unavailable negative credit. The optimizer recomputes a complete objective
and supplies exact rational convex-optimality data. No gain from a pulse,
angular gate, or earlier cutoff is added to the new gain.

The grounding conclusion uses the additional actual finite paths and the
independently checked strict frontier ranks. The energy conclusion instead
uses every intermediate physical cut and the final all-future gate. These
are not equivalent: N7 fully grounds but fails the energy target.

Both normal and optimized Python reproduce the same receipt and bounds.
The initial probe values were reproduced exactly, rather than imported as
expected constants into the mathematical calculation. The remaining
inherited assumption is the existing authenticated source-prefix receipt;
the extension does not rerun its full original multi-million-source census.
Every original source remains decodable through the frozen carrier, and
every newly transported owner has its appended full clock recorded.

## Independent audit of the optional-pulse boundary repair

The separate `lagged_tail_pulses_20261022/check.py` consumer now requires its
frontier argument to be a set or frozenset of literal integers, exactly equal
to the complete marked terminal universe in `state.ports`. This closes the
gap where an omitted terminal could make a path that crosses the stopping
boundary appear to avoid it. It also rejects floating-point and Boolean
aliases before set equality can identify them with integer terminals.

The hostile path 7 -> 11 -> 17 has two actual native edges of clock 1 each.
With terminal 11 present in the state, it cannot support the claimed
delayed-cancellation relation: the late source stops at 11 and cannot be
continued to 17. The repaired consumer rejects an empty supplied frontier,
the alias set `{11.0}`, and the literal set `{11}` (in the latter case by
rejecting the path crossing itself). Both ends of every supplied edge are
tested against the bound frontier. Thus the proof's common stopping rule
and the public consumer's stopping rule now agree.

An independent optimized-Python replay passes 15052 explicit checks and
preserves the original pulse receipt SHA256
`c0be85b63ccb3c00d7263b7180c63a88f3c4fcc6d5c4225728578a189c90a4ef`.
Audited repaired code LF SHA256:
`95132e89c79e7d3d84f4edf6fb3ceb38c225491b7276015f91622e27021a2f08`.
The repair changes admissibility of false boundary inputs, not the saved
native paths, optional-pulse theorem, or four existing numerical bounds.
