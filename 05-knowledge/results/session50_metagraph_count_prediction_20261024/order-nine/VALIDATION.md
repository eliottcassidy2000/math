# Validation

## Controls already passed

- Input preparation separately canonicalizes all reconstructed order-eight
  parents, exhausts their actual vertex actions, proves pairwise class
  distinctness, and verifies total labelled orbit mass exactly 2^28.
- The standalone predictor is run normally and with `python -O`. Their entire
  JSON outputs, including all 6,880 per-parent summaries, are identical.
- The separate regular-sector implementation is run normally and with `-O`.
  Those complete receipts are identical. It uses explicit root-transition
  components as well as the trace, and agrees with the predictor's regular
  sector on all four reported aggregate values.
- The C++ bulk control uses full color-vector refinement, unlike Python's
  explicit-block refinement with a root-singleton early exit. It matches the
  Python state total, upper trace, root-cell histogram, and every per-parent
  summary. Its action validation checks permutations, arc preservation,
  identity, closure, and complete attachment-orbit coverage.
- Python validation uses explicit exceptions, never `assert`, so no proof
  check disappears in optimized mode. C++ uses explicit runtime checks.
- The predictor installs an I/O audit hook before reading its inputs. It can
  read only its own source and the two pinned parent files, and write only
  the requested result path. It forbids network, subprocess and ctypes use.
- The only canonicalization routine in the predictor requires exactly eight
  vertices at runtime. Its only call receives an actual deleted parent from
  one of the eighty-five regular-sector extensions. The C++ bulk control has
  no canonicalization or deletion routine at all.
- The verifier checks the predictor's AST contains no integer literal equal
  to its derived target count, target upper trace, exact trace, full correction,
  regular-sector trace/error/count, or pointed-state count. These values are outcomes,
  not calibration constants.

## Reproduction

`python verify_packet.py` checks manifest hashes and the exact cross-controls.
`python verify_packet.py --full --output /tmp/order-nine-replay.json` first
copies only manifest members to a fresh temporary directory, then regenerates
the input, runs normal and -O full predictors, reruns the regular-sector
control, compiles and reruns the C++ control, and compares complete receipts.
No source file outside that temporary manifest-only copy is needed.

The frozen `replay-receipt.json` records the completed full relocation run.
Independent audit results are separately preserved by the Session50
order-nine audit packet; this packet does not substitute its own controls
for that independent review.

## Exact output and work budget

The regular sector has 85 pointed states and 765 parent-eight canonical calls.
The exhaustive candidate-label total is 2,130. Its trace is 135, giving fifteen
regular classes and error 630. Its regular classes have vertex-orbit shapes:
seven of shape 1^9, two of shape 1+1+1+3+3, three of shape 3+3+3, and three of
shape 9. These contribute errors 504+72+54+0=630.

The 1,716,608-state bulk pass gives upper trace 1,724,582. Python performs
1,497,092 nontrivial refinement rounds after safe early exits; the C++ full
refinement control performs 4,131,890 rounds. Both have zero target
canonicalizations and zero unrestricted target re-root passes. The 765
regular-sector deletions are explicitly counted, not hidden by the wording.

The correction is derived as 630+14+72+42=758, and the exact count is 191,536.
No closed formula at order ten or above is implied.
