# Reproduction and scope

From the repository root:

```sh
python3 05-knowledge/results/session288_twin_lengths_20261112/reproduce.py
```

The driver freshly runs ten engines normally and with Python optimization,
compares each with its stored receipt, and writes `REPLAY.json`. Python3,
NumPy and a C++17 compiler are required for the stable lane. Other main proof
engines use the standard library. The optional wildcard discovery scans need
NumPy/Numba as specified in `wildcard/REPRODUCE.md`; they are not assumptions
of the six-core universal identities and are not part of the aggregate replay.

The stable-cylinder compilation uses the inherited atlas's
`five_stress_connection_atlas_20261112/cyclic_port_receipt.json`, pinned by
SHA256 in the output. It does not rerun the full tournament census that
produced that input. It does recompute the signed all-length defect tensor
of every retained cylinder. Source/data fingerprints are recorded by the
independent review; no unobserved computation is claimed.

The ten comparisons cover the one-clock theorem's tests, exact valuation
controls, complete two-port skeleton and literal arbitrary-module witnesses,
stable tensor/canonical controls, all98 cylinder certificates, the finite
contraction DAG, a scoped deterministic12,000-extension screen, and an
independent permutation-plus-cut counter, the new three/four-clock tensors,
the exhaustive320-template next-clock obstruction, and independent literal
base counts. Finite controls are not substituted
for the all-order proofs in the notes.
