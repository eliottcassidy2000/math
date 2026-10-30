# Synthesis and remaining proof obligation

**Session 118, 30 October 2026.** Anchor: strengthen shared-B2 gate deletion.
Niche: exact global inventory constraints. Wildcard: incidence antichains and safe
Boolean factoring. LRC(14), general Ramsey records and unrestricted superlinear
explicit-NP circuit bounds remain open; no result here settles them.

## Inheritance and information retained

The pinned source is the provisional ternary-pair elimination manuscript reviewed
in [Session 117](../session117_ramsey_circuit_interfaces_20261030/SOURCE_MODELS.md).
Its original row-potential region was completely classified there. This session
changes one structural row, not merely its weights. We do not republish third-party
source text, execute submission code, or submit a challenge entry.

The closest proved mechanism was the one-bit bottleneck in the point-defect circuit
classification. The hostile was unrestricted sharing and the impossibility of
counting the same gate twice. The corrected near miss was the proposal that AP2
might force five deletions; the minimum-circuit four-deletion witness forbids that
under the local hypotheses. The least-used sidecar was output ancestry: keeping
an actual path to the live output, rather than only gate counts, forces the fourth
deletion. The concept board was sharing, safe-core dimension, usage inventory,
constant antichains and semantic resynthesis.

The maps and their losses are explicit:

- A shared DAG maps to deletion rows, preserving guaranteed distinct gate losses
  but forgetting ancestry, which is exactly the data used by the strengthened AP2
  lemma. Scalar equality schedules are not assumed realizable.
- A full trace maps to endpoint inventories. Pair births/deaths and usage changes
  survive, so AP1/AP2 frequency is bounded by terminal safe-core dimension. Gate
  coincidence information must already have been checked before taking this map.
- A legal restricted circuit is not necessarily a minimum circuit, even when the
  original circuit was minimum. The tight-prefix hostile exposes this missing
  hypothesis and its residuals admit an explicit safe factoring improvement.
- Constant-seed sets retain a usable deletion charge only with ancestry control.
  An antichain gives a sharp factor-two bound; a comparable chain refutes dropping
  that condition.

## Concurrency and attribution

The worktree began at cd6d98ba58a2418fb81054659bc37962676c8dd0, preserving concurrent
Ramsey capacities, counting-compiler and arithmetic/LRC results. The separate
[Ramsey/circuit interface packet](../ramsey_circuit_interfaces_20261030/README.md)
was read for its shared-circuit compiler and primary-input accounting caveat.
Those are credited source results, not re-proved in this session. A parsimonious
encoding preserving shared wires is useful for a future exact checker; its auxiliary
variable count is not an unrestricted gate lower bound or free input-length change.

## What would change the leading coefficient

The exact strengthened slack identity gives

    s0 = 4n-8h + sT +4(IT+J+U)+2qT+A+2Z+3H+D+K+E.

Notation and hypotheses are in the accounting proof. At terminal h=O(log n),
a positive linear lower bound on the remaining surplus would raise the leading
coefficient. We do not prove such a bound. The inventory theorem instead shows
that all but O(h) relevant events must be GP and P2/PE, focusing that obligation.

The genuine long-prefix example refutes a universal positive surplus claim for
arbitrary reduced circuits. Its linear h and nonminimum circuit presentation do
not refute a low-dimension/minimality-sensitive theorem. Conversely, finding a
rewrite for that example is not proof that every tight cycle admits the rewrite.

The next justified target is a dichotomy for the high-frequency GP/cleanup blocks:
either certify extra deletions or usage loss, or identify a safely replaceable
subcircuit with a quantified gate saving. The replacement must preserve every
external use of shared intermediate values and the restricted-domain semantics.
The factoring companion provides one exact sufficient pattern. On the hostile
family it saves k ordinary gates, then a legal M/D4 role exchange with the retained
map y=w+a saves k further charged equations. The resulting representation is
minimum in its new coordinates; restoring original coordinates costs k XOR gates.
This is not a complete dichotomy. That complete structural theorem remains OPEN.
