# Square-sum machines

**New explicit artifacts, not a new Hamiltonicity theorem.**

* A hand-written **14-fraction FRACTRAN** square recognizer, with exact terminal
  states and all-input step formulas.
* Actual deterministic **six-tape Turing machines** for square-sum Hamiltonian
  path/cycle certificates:25/36 nonhalting states. Full ordinary transition
  tables are mechanically expanded from compact disjoint clauses; frozen
  row counts and SHA256 values are included (no pseudocode/universality appeal).
* Exact cycle-verification clock
  `5n²+14n+11+3v_first+6Σ_edges sqrt(u+v)`, isolating a root-weight invariant.
* The named-vertex sidecar cannot be dropped: `(1,8)^4` passes local square
  tests, length, total sum, and root moments at n=8, but is not Hamiltonian.

Read [THEORY.md](THEORY.md) for contracts, proof, inherited boundaries, and
reproduction commands. Gerbicz's cycles for every n>=32 remain **CITED**;
finite witness verification is not a replacement all-order existence proof.

`build_tm.py` regenerates the JSON and fully expanded gzip transition tables.
The gzip tables are reproducible local derivatives, not committed binaries.
`check.py` materializes missing tables and checks their frozen expanded hashes.
`check.py` writes exact finite-scope verification to stdout. `simulate_tm.py`
executes only finite transition rules; `square_fractran.py` offers literal
integer and exponent-vector FRACTRAN execution.
