# Session145 magic/frequency lane

Read `PROOF.md` for full proofs and map/loss contracts. The main result is an
exact bridge, with an equally exact obstruction:

- The actual nine prior labels are 174,285,396,417,528,639,741,852,963.
- No permutation or nonconstant affine transform of them is fully magic,
  because their required center 555 is absent.
- Their Lo Shu rank placement has seven sums 1665 and one sum 1584.
- In orthogonal Latin coordinates it is
  `528+111U+(567/2)V+(81/2)V²`; the requested grid is `126+4U+18V`.
- Every real-valued map on the nine labels that is magic in this placement is
  exactly `c+pU+qV`. Strictly increasing maps satisfy `p>0, q>2p`.
- The order map and the trace map are distinct. Both retain the three fibers,
  but transported decimal rotation is not a numerical magic symmetry.
- There is a complete fixed-trace polynomial classification; the next mode
  after the exceptional cubic first occurs in degree seven.

`verify.py` recomputes the exact result and uses both center-fixed exhaustive
permutations and an independent constructive census. `test_controls.py` rejects
nine artifact mutations, checks eleven mathematical hostiles, and verifies
parameter criteria on 625 integer pairs. Both scripts pass with ordinary and
optimized Python, with byte-identical result files.

`SOURCE_PINS.json` pins three upstream blobs to revision
`7b76044ca13e6684d1746923285b83dbf83828d2`; `verify_source_pins.py` validates the
archived source bytes and Git blob hashes. The older atlas contains unrelated
preprint discussion: none of it is imported as a theorem here.

There is no health or physical-effect claim, no curve-level C3 action, no
canonical identification with another carrier of dimension 126, and no Ramsey
or knot theorem established by this lane.
