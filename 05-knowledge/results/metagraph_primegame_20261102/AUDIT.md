# Audit, exact universe, and reproduction

The universal spectrum remains OPEN. There was no subagent or formal-proof
audit. The independent methods below were written and checked by one agent.

## Reproduction

From the repository root:

```bash
python3 05-knowledge/results/metagraph_primegame_20261102/replay.py
```

The script verifies frozen source hashes, restores the seven class banks
from public commit8286d710498a3b4e9b9f1ac683c6eb77716f3f4b, and compiles
with `c++ -std=c++17 -Wall -Wextra -Werror` under both `-O2` and `-O0`.
Python checks run normally and with `-O`; validity gates use explicit
exceptions rather than removable assertions. Receipts must agree byte for
byte. Temporary inputs and binaries live outside the repository. No nauty
installation or external network download is required for replay.

## Declared universe and independent paths

- `ports.cpp`: all unlabelled tournament classes at orders3..9, respectively
  2,4,12,56,456,6880,191536. Every vertex mark is tested. Source H is compared
  with the inherited separately implemented census on **every class**.
  Starting/ending subset DPs compute q,g and each port charge.
- The complete filter is **strong and H<125**:687 sources. Its output is
  reproduced as low_sources.txt; BINARY_SUMMARY.json is derived, not hand
  selected. SCC reduction plus the inherited f(10)=125 supply the all-order
  interpretation. Source H>=3 and raw H are retained throughout.
- `audit_ports.cpp`: a forward-extension DP counts the full child adjacency,
  independently of the q,g formula. It checks5210 children and327 literal
  permutation source counts. The two-birth and multi-port probes use actual
  labelled graphs, not abstract positive-coordinate states.
- `audit_two_birth.cpp`: independently enumerates each source's Hamilton
  paths, constructs all endpoint-moving injections, and explicitly checks
  legality and uniqueness of14228153 inserted paths for45216 ordered birth
  pairs. Full adjacency counts verify all lower bounds.
- `induced_branches.cpp` and `audit_induced.py`: subset DP versus literal
  enumeration of every directed simple path. All153914 nonempty proper
  subsets and713827 paths are checked. The latter supplies the explicit
  noncomposable47->73->117 witnesses. The all-order negative conclusions
  also require the prime-source reduction, not merely these numbers.
- `padded_ears.cpp`:438 source/padding choices and297236 complete cuts, on
  the nine declared one-port-dead H counts, with total padding<=2. Negative
  results from this probe remain bounded; they are not all-order ear bounds.
- `primegame_code.py`: integer/rational and finite-field ranks, the exact
  kernel, all128 square-sigma equivalences, simplex/Hamming duality, all14
  Fano circuits, and actual priority-positive controls11<->13 and the
  six-cycle starting4301. The complete priority-cycle classification is
  inherited from Session129, not reproved by those two finite controls.

## Hostiles, proof review, and stopping boundary

The review checked source/sink cases in the birth theorem; old/new and
repeated/distinct ports; uniqueness after deletion; oddness rounding;
nondecreasing induced path counts; prime SCC factorization; the strict
q<3p cutoff; p=7's unavailable source; and the difference between existence
of compatible embeddings and composability of their numerical images.

Positive controls include the two-port examples55->111,61->123,73->145.
Hostiles include31->63 (even arbitrary extension fails),47->73->117 (lost
middle representative), and instruction C (formal parity differs from
actual periodic execution). The Fano relation code is Hamming; the return
observation code is simplex. Neither code retains a word's priority guards.

No known false canon was discovered or altered in this packet; no MISTAKE
entry is needed for a failed research candidate that was never promoted.
The algebraic labels remain OPEN as an operator-identification question.
The two-birth theorem supersedes the finite depth limitation of
WORD_SEARCH.json; that search remains a reproducible independent probe.

## Completed replay and incoming integration

All six C++ programs and receipts passed under both O0 and O2. Python
normal/-O outputs also matched. After integrating the incoming packets, the
C++ sources and receipts were unchanged. The expanded PRIMEGAME checker was
rerun normally and under -O: the19-entry incidence factorization, rank13,
rank-three circulation quotient, seven exact probe classes, and the line
determinants(2,-2,0,0,0,0,2) passed. INTEGRATION_AUDIT.json records its hashes;
INCOMING_SOURCES.json pins the exact eight incoming documents read. The replay
command now checks both source-pin sets. Documentation and whitespace checks
are also part of Git close-out.
