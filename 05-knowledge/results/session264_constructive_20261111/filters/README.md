# Positive component compiler and retained-endpoint consumer

[Proof and boundaries](PROOF.md).

New physical ingredients:

- Triangle wrap turns weighted rooted paths into ordinary Hamiltonian paths
  with the SAME old endpoint; gate position records the old path length.
- The full coloured cover response has an exact positive odd-component merging
  operator, with a finite skeleton compiler for arbitrary fixed uniform cores.
- The explicit11-vertex mod8 consumer has H7623, endpoint count3696, and462
  paths per suffix-length residue. This is conditioned readout, not a scalar
  projector onto one residue.
- The return-lane H29 modules become two8-vertex modules with H349 and equal
  complete ordinary/root-deleted profiles, but start counts52/50 and end41/43.

`controls.py` tests190 new-operation gates and produces `controls.out`.
`coalescence.py`, `interface.py`, and `explore.py` are reusable compilers.
`primitive_scan.cpp` and the retained input/output files document unsuccessful
primitive16 searches in precisely bounded one-, two- and three-wrap grammars;
`scanner_positive8*` gives the successful scanner control. No general primitive16
nonexistence conclusion is drawn. Sources and reproduction are in the proof.
