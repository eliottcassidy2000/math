# Chan–Pak: choose the carry from the target, and program only its residue

**Exact universal placement remains open.** The work below changes
which parts of the inverse are solved, and supplies all-order physical
counts as well as a deterministic solver for the remaining local task.

1. **[No all-ones alignment](ABSOLUTE_CARRY_DECODER.md).** The actual
   high sum selects one of the two binary carry branches. For
   M=2,3,4,5,6,7,8,9,10,12,16, every sufficiently large odd N
   has two strong tournament counts bracketing it, congruent to N
   mod 2M, with error O_M(N^(3/8)). The [second certificate bank](EXTENDED_MODULI.md)
   adds five moduli. The finite certificates are exact over all
   input states; the bracketing theorem has no target-size cutoff
   other than “sufficiently large.” The earlier general fixed-modulus
   O_M(N^(2/5)) result still applies at other moduli.
2. **[The carry is a weighted cut](CARRY_CUT_LAW.md).** It chooses
   the graph of free bit positions; assigning the digits chooses
   the cut. At the first output bit the only frozen components are
   isolated vertices and runs of three consecutive vertices. The
   six-column test gives a forced local choice for the target's
   first bit without changing the entire carry. At higher precision,
   [guarded XOR gates](XOR_CARRY_GATE.md) can force a prescribed bit
   without altering lower bits. Shared guards and absolute cost
   still matter. A separate linear-length cut-star construction
   realizes an arbitrary **relative** dyadic correction, but not
   an arbitrary absolute baseline.
3. **[Choose the nearest real phase with a 54-state carry scan](REVERSED_CARRY_SELECTOR.md).**
   The reversed carry gives the leading checksum exactly up to a
   certified real envelope. The [target-only terminal solver](RESIDUE_ROUTED_SELECTOR.md)
   runs the same scan over just the requested residue, with no
   all-syndrome alignment. An empty finishing set or model distance
   above the envelope rules out the entire interval. A physical
   congruent result with total certified error below 2M is exact.
   These are conditional certificates, not a proof that every
   target has a viable interval.
4. **[Uniform phase on quartic bands](QUARTIC_CARRY_SHELLS.md).**
   For each fixed modulus, explicit bands at height of order q^4
   give congruent physical brackets with error O_M(q) for order
   q^3 actual target positions per dyadic scale. The exact inverse
   there has at most one high-address candidate per fixed common
   phase. The bank has too few exact values to cover all targets
   in its own intervals. Another coarse selector is necessary.

[Rule 30](https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/)
suggests solving one input from an output and guarded neighbors.
[Quasi-Riemann's alternate proof](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex)
suggests duplicating the hard row through invariant copies before
using an average. The cut law produces literal modular copies and
exposes the exact constant that their real heights still read.
[Planar non-properness](https://arxiv.org/html/2101.05245) supplies
the complementary question: can locally compatible branches
escape every absolute bound? The present construction has a
specific example and an exact interval exclusion test, not an
application of the Jacobian conjecture. The prior repo sources
and their mathematical scopes are linked in the proofs.

The [Collatz graph](FUNCTIONAL_GRAPH_BRIDGE.md) has one outgoing
arrow per positive integer; removing the outgoing arrow of 1
would give a spanning tree toward 1 exactly if Collatz holds.
The target carry decoder also has one prescribed outgoing arrow
on its viable states, but it has a proven terminating rank:
the number of columns remaining. Finding a viable initial state
at every absolute target is the still-unproved step.

## Replay

Run the seven `replay_*.py` scripts with Python 3; the three
phase and selector replays use the adjacent predecessor packet
`chanpak_absolute_steering_20261115`. Receipts are included.
Use the instructions in [extended moduli](EXTENDED_MODULI.md)
for the independent C++ certificate.
