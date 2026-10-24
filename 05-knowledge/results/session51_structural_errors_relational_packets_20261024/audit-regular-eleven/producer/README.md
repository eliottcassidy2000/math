# Targeted regular-eleven correction

The complete nine-parent action atlas yields13,333 almost-regular ten-parent
classes. Their parent-only reroot trace is13,453, hence1,223 regular-eleven
classes and correction133,210. Combined with the independently proved
nonregular packet: Es11=142,376 and E11=152,690. The full U11 and a11 remain open
in this work. No full ten-atlas or eleven-vertex canonicalization is used.

- `regular_trace.cpp`: eligible subset orbits, complete cyclic-cell full-matrix
  parent normal forms, full parent actions, reroot components, converse control.
- `eleven-O2-parents.txt`: normal-form parent rows and complete actions.
- `eleven-O2-states.tsv`: every exact reroot row and converse parent index.
- `eleven-O2-components.tsv`: child components, full-group orders, orbit sizes.
- `eleven-O2.json`: exact result and operation counts; O3 files are identical.
- `degree_mass.py`: independent exact labeled residual-degree dynamic program.
- `THEOREM.md`: universe, proofs, all distributions, scope, and limitations.
- `SOURCE_PINS.json`: authenticated input lineage, including structural packet.
- `combined-correction.json`: derived complete order-eleven defect distribution.
- `verify_packet.py`: portable manifest, source, arithmetic, and component checks.
- `hostile_controls.py`: six corrupted-input negative controls.
- `reproduce.py`: complete O2/O3 replay, with byte-for-byte artifact comparison.

Reproduce from a clean directory:

    python3 -B reproduce.py --output-directory /tmp/regular-eleven-replay
    python3 -B verify_packet.py
    python3 -B -O verify_packet.py

Each parent text file begins with its class count; each record is n row masks
followed by an action count, then that many old-label-to-new-label permutations.
State TSV destinations have format `state:multiplicity`, comma-separated.
Compiled binaries and caches are not part of the frozen certificate payload.
