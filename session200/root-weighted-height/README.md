# Root-weighted height and a short modular constructor target

`PROOF.md` derives exact binary-address laws for W=U−d, a graph-order bound
charged by root changes, and a logarithmic-variable finite-field formulation
with a linear prime threshold for the at-most-three-root-change family.
The all-target constructor remains open. The finite census refutes a tempting
two-length restriction and supplies an explicit target49,937 witness.

Reproduce the complete finite census:

    python3 probe.py > /tmp/root-weight-normal.json
    python3 -O probe.py > /tmp/root-weight-optimized.json
    cmp /tmp/root-weight-normal.json /tmp/root-weight-optimized.json

The finite universe is all source response/root-count/retained-letter-count
triples with balanced output at most50,001 and at most three J letters.
Cost coordinates are retained throughout. The census has641,253 states.

The two files under `dependencies/` freeze the inherited physical source
and C1 interpretation, the midpoint pairing, and the existing four-block
conjecture. The four-block conjecture is not assumed for the proved size or
finite-field statements. No current index or proof status outside this
packet is changed.
