# Order-nine structural predictor

The next tournament count is predicted from a certified order-eight parent
atlas and its actual vertex automorphism actions, with no target atlas.

- `THEOREM.md`: exhaustive nonregular strong-nine error classification and
  the regular-nine reduction to eighty-five almost-regular-eight parents.
- `prepare_inputs.py`: reconstructs ordinary parents from the pinned merged
  source; exhaustive parent canonical forms and labelled orbit mass certify
  completeness. This is input preparation, not target prediction.
- `predictor.py`: standalone full predictor reading only the two pinned parent
  input files. Its bulk pass only extends and refines; only the regular sector
  performs parent-eight isomorphism tests.
- `regular_trace.py`: independent detailed regular-sector trace and component
  control. Child canonicalization is guarded out.
- `upper_trace.cpp`: independently structured compiled full-refinement control,
  with no deletion, isomorphism, or target automorphism implementation.
- `normal.json`, `optimized.json`, `upper-compiled.json`: full results including
  the per-parent state count, upper trace and root-cell histogram.
- `verify_packet.py`: manifest integrity, cross-controls and optional complete
  manifest-only relocation replay.

Run:

    python verify_packet.py
    python verify_packet.py --full --output /tmp/order-nine-replay.json

The full replay needs Python 3 and a C++17 compiler. It does not access a
repository checkout, sibling packet, network, external tournament database,
or hidden target list. All source inputs are archived under `archive/`.
`SOURCE_PINS.json` and `manifest.json` record bytes and SHA-256 digests.

Result: 1,716,608 pointed extensions; refinement upper trace 1,724,582;
structurally derived correction 758; exact predicted ordinary count 191,536.
The known target value was held out of executable inputs. This is not a claim
of prediction from a bare unmarked metagraph, nor an all-order closed formula.
