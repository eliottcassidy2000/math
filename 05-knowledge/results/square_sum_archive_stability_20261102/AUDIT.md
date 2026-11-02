# Audit scope and reproduction

Status: same-session mathematical review, with two different full exact
graph-counting paths. No independent agent was requested, used or claimed.

## Order-ten certificate

Universe: every unlabelled simple graph on ten vertices, generated without
filters by nauty 2.9.3 `geng -q 10`; total 12,005,168. There are no inherited
connectivity, degree, orientation or Hamiltonicity filters. Both kernels in
`graph10.cpp` evaluate every input, using matching recursion versus explicit
edge-pattern enumeration and two wedge counts. Output `graph10-dual.json`.
The complete-generation guarantee is inherited from `geng`, not inferred
from agreement of the two consumers or merely from the graph count.

`verify_reduction.py` independently generates the full cycle partitions and
distinguished-cycle lifts. Its exact rational matrix controls count H by
subset DP and compute determinants by elimination. They include transitive,
all-fair and five-fair-pair matrices, sixteen fixed-seed dyadic matrices and
six perturbations of the sharp witness. Normal and optimized Python outputs
are byte-identical (`reduction.json`, `reduction.optimized.json`). The actual
all-real-weight proof uses the interpolation reduction, not these 25 controls.

The exact maximum 9583/18432 and positive reserve 8849/18432 are checked as
Fractions. No floating-point inequality or removable Python assertion is
used. Integer scales and the fifth-matching removal are derived explicitly
in STABILITY_10.md. Pattern universes and graph6 length/alphabet/padding are
checked, as is the expected full universe count. Graph6 parsing and the
generator are shared dependencies of the two kernels.

## Prepublication correction

The first draft of the minus four-cycle group omitted `-2 M244`. Selecting
one four-cycle from a collection containing two counts the collection twice.
The symbolic residual checker rejected the draft with a negative coefficient
at partition (2,4,4). The repaired group retains `M44-2 M244>=M44/2`, since
the remaining two vertices have variance at most 1/4. The full symbolic and
matrix controls were rerun. No earlier published theorem or frozen receipt
is changed. The mechanism is recorded in the repository mistakes ledger.

## Native arithmetic contexts

The proof preserves an actual finite numeric source, its induced adjacency,
the selected forests, and the same named ports. It allows only nontrivial
path components with distinct endpoint vertices; zero-length incidence-slot
paths would require the stronger inherited slot model and are not silently
included. Fresh bridge vertices are degree two in the final induced graph.
This forces all connector edges in any Hamiltonian cycle and makes the
rejected component pairing a genuine nonextendability result.

The 65 exact controls in `verify_native_contexts.py` cover four square scales,
both accepted/rejected pairings, both connector orders and all connector
orientations. Every pair of labels is tested for inducedness, and both full
resulting 2-factors are checked for their components. A further eight-port control joins four source fragments into one16-cycle while the rejected pairing leaves cycles of sizes4 and12. Normal/-O outputs agree.
The theorem is all-size by proof; 65 is only the implementation test universe.

Incoming verification is separate: `check_incoming.py` checks every prime
label and every square edge of the50396 witness, independently reconstructs
the complete Q25 path-cube spectrum, and records the exact failed direct
order-eleven extension. It does not rerun the prime minimality exclusions.
Its normal/-O receipts agree; INCOMING.md retains full scope and provenance.

## Recovery controls and honest limits

`catalogue.py` records a pinned Git baseline, all matching file paths, Git
blob identities, excerpts and heuristic topic tags. Audit copies are not
counted as independent results, and historical messages are provenance only.
`HISTORY_DISCOVERY.json` additionally records the reachable-ref filename
history scan (2,641 paths; six absent at the baseline). Those six are other
algebraic uses of 'square', not missing square-sum Hamiltonian sessions.
This is not a claim to exhaust untracked conversations or every synonym.
RECOVERY.md separates reviewed proofs, inspected scoped reports and discovery.

## Commands

From this directory:

```bash
clang++ -O3 -std=c++17 graph10.cpp -o /tmp/square-stability-graph10
set -o pipefail
geng -q 10 | /tmp/square-stability-graph10 dual 12005168 > /tmp/graph10.replay.json
cmp graph10-dual.json /tmp/graph10.replay.json
python3 verify_reduction.py > /tmp/reduction.replay.json
python3 -O verify_reduction.py > /tmp/reduction.optimized.replay.json
cmp reduction.json /tmp/reduction.replay.json
cmp reduction.json /tmp/reduction.optimized.replay.json
python3 verify_native_contexts.py > /tmp/native-contexts.replay.json
python3 -O verify_native_contexts.py > /tmp/native-contexts.optimized.replay.json
cmp native_contexts.json /tmp/native-contexts.replay.json
cmp native_contexts.json /tmp/native-contexts.optimized.replay.json
python3 catalogue.py
python3 check_incoming.py > /tmp/incoming.replay.json
cmp incoming_checks.json /tmp/incoming.replay.json
python3 seal.py --check
```

The first command pair can take several minutes. Do not rerun a complete
unchanged universe merely for another routine progress check.
