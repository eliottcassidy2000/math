# Twelve and twenty-two: power clocks, integral recovery, and mixed observations

**2026-10-29. PROVED scoped / FINITE-EXACT / CITED classical. Same-agent
verification; no independent audit or external novelty claim.**

This packet follows the user's request to recover earlier appearances of12
and22, starting from Session104's synthetic exterior model. The new object is
the information retained by powers of its actual counting matrix. It is not
a new identification of the six-state legality machine with an exterior
square, nor a solution of general tournament complement closure.

## Main results

Write

\[
C=\begin{pmatrix}1&0&1\\1&1&1\\1&1&0\end{pmatrix},\qquad D=C-I.
\]

- **P1--P3, integral reconstruction.** The identities
  \(C=C^4-5C^2+2I\) and \(D=3D^2-D^4-I\) show that replacing either
  primitive clock by its square loses no generator in its integer algebra.
  More generally, two unit clocks generate the clock at their gcd. In a
  rank-three integer algebra their separate power observations already
  span that generated algebra. For clocks12 and22, five prescribed readings
  therefore recover all primitive-clock observations.
- **P4, exact latency.** The safe-word counter has a unimodular ordinary
  observer. Four readings at times0,12,22,44 recover its entire state at
  every modulus. Among times12a+22b, the earliest complete horizon is44;
  modulo13 proves the lower bound. Three readings first suffice at44,46,48.
  This is a tradeoff between latency and number of linear observations.
- **P5, an all-size precision recurrence.** For clock12 times13^r, the
  13-primary Smith exponents are \((0,r+1,2r+2)\). Adding clock22 restores
  the full integer observer for every r. The same matrix's arithmetic trace
  observer has a separate, irreparable trace blind spot at7.
- **F1--F3, actual transport.** The dual roots are seventh-root traces;
  angle doubling is \(d\mapsto d^2-2\). At11, clock12 is a Frobenius
  conjugate of the exterior dual and clock22 a conjugate of the square.
  Reducing at2 gives a literal F8/Paley7 chart: Frobenius preserves its
  arcs, whereas changing the generator from C to D reverses every arc.
- **P6--P7, hostile boundaries.** The Pell unit behind Ramanujan's g22
  leaves an extra dyadic coordinate under these two clocks. In the
  biquadratic factors of g130, all separate-clock traces miss a mixed
  amplitude. One mixed reading detects it. The rank-three theorem cannot
  be generalized to rank four without more hypotheses.

Proofs and exact scopes are in [POWER_CLOCKS.md](POWER_CLOCKS.md) and
[FROBENIUS_DUAL.md](FROBENIUS_DUAL.md). The historical sweep and correction
lineage are in [ARCHAEOLOGY.md](ARCHAEOLOGY.md).

## Inheritance and concept board

**Closest proved mechanisms:** Session104's integer exterior conjugacy,
the October28 quadratic trace-clock Smith formula, and the common-source
observer kernels in the overnight terminal-selector packet. **Hostile
controls:** the actual legality word `ee` has rank2 and is not an exterior
square; Pell dyadic loss; a rank-four mixed trace; and the twelve-vertex
81/27 pair, whose equal masses do not identify joint port laws.
**Corrected near miss:** June's self-converse fixed-class count is a trace,
not the whole even eigenspace. **Least-used relevant sidecars:** March's
complete divisor-prime selector and August's folded C7/K4 common-drift
obstruction.

The final clean fetch also recovered the concurrent
[atomic-descent wedge packet](../atomic_descent_wedge_20261029/WEDGE.md),
which independently proves the same exterior conjugacy and the native
modulo-seven binomial count law. Both exterior certificates are credited
as one inherited result; P3 now compares its visible nilpotent moment
with the arithmetic trace's blind kernel.

The live board has five concepts: integer clock lattices; Frobenius phase;
joint versus separate observations; port realization; divisor-prime
selection. The anchor is clock recovery, the niche is historical selector
and phase corrections, and the wildcard is the quadratic/biquadratic
Ramanujan boundary. Their common question is which *operation* exposes a
missing coordinate. A shared integer alone does not supply that operation.

The archive search found9462 documents containing a digit-block12 or22,
of which3692 contain both, in three specified roots at cached HEAD
`2cc2ce454a`. This is a mechanical inventory, not9462 semantic readings.
The atlas records the selected source chains actually inspected and the
transfers accepted or rejected.

## Verification and scope

From the repository root:

```text
python -X utf8 -B 05-knowledge/results/twelve_twenty_two_trace_20261029/check_clocks.py
python -X utf8 -B -O 05-knowledge/results/twelve_twenty_two_trace_20261029/check_clocks.py
python -X utf8 -B 05-knowledge/results/twelve_twenty_two_trace_20261029/check_precision.py
python -X utf8 -B -O 05-knowledge/results/twelve_twenty_two_trace_20261029/check_precision.py
```

The scripts use SymPy and python-flint, explicit failure checks, and saved
receipts `clocks.json` and `precision.json`. Checks include a literal
five-letter grammar enumeration through length7; determinantal checks of
Smith factors; exhaustive three-reading banks through time48; trace
carriers modulo1 through24; all unramified prime controls through100;
all7! Paley permutations and Hamilton paths; and nine13-adic tower levels.
Those finite checks support the stated proofs; they are not the proofs of
the all-r or all-modulus claims. `search_archive.py --write` deliberately
refreshes only the archive inventory against the advancing checkout.

No new A000568 terms, arbitrary strong-atom classification, integral
twelve-vertex complement recovery, LRC theorem, or Collatz termination
claim is made. The 81/27 pair retains the inherited two-run ranks24/22 and
three-run ranks430/435 modulo3; its full complemented relations and
selected-path phases are still required. Linear reconstruction of a clock
observer is not permission to reset a tournament class or realize an
arbitrary signed tensor by a physical attachment.

## Next decisive questions

1. **OPEN:** For actual complemented strong-atom transfers, characterize
   when their integer response algebra has the rank-three sum property.
   Test the paired kernel against the full target relation, not just a
   scalar recurrence or characteristic polynomial.
2. **OPEN:** Classify minimal mixed-word probes when separate subalgebras
   generate a rank-four or larger carrier. The g130 tensor example gives
   a one-probe witness, not a universal bound.
3. **OPEN:** Find clock schedules minimizing both observation count and
   latency for integral recurrence modules, including their local Smith
   obstructions. Here the two optima separate: four readings by44 versus
   three by48.
