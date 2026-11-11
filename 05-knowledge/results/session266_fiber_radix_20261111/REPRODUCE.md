# Reproduction and scope ledger

Run from this directory with Python3; no nonstandard packages are required.

    python3 -B fibers/build.py --n 30
    python3 -B fibers/build.py --h 63 --n 100 --output query_63_100.json
    python3 -B fibers/finite_census.py
    python3 -B radix/replay.py
    python3 -B collatz/probe.py
    python3 -B smooth/replay.py
    python3 -B geometry/replay.py

The query output is optional and is not part of the manifest. The first command
regenerates `fibers/report.json`. All-H-through8 uses `finite_input.json` rather
than performing any tournament enumeration. JSON receipts and CSVs use LF.

## New-constructor controls

- Fiber compiler:5704 active checks on the order30 report, including the weighted
  integer channel, cumulative summation, convexity and inverse reconstruction.
  H1 cumulative outputs exclude the empty object. Invalid scope queries reject.
- Finite atlas compiler:5312 checks, all odd H1..661 and sizes1..8; class totals
  1,1,2,4,12,56,456,6880 and cumulative7412. Larger H and all even H vanish in
  this finite range. This is inherited data processing, not a fresh census.
- Radix algebra:17154 controls of the stated elementary constructors.
- Collatz/puncture bridge: all64 labeled four-vertex tournaments; all8 fixed-path
  completions; five punctured profiles; four explicitly known cycle witnesses;
  500 signed table transports,1600 code products,499 raw-code steps and literal
  hole controls. These are not a complete Collatz cycle search.
- Smooth classification:1000 exponent triples0..9, comparing the explicit
  degree/leading formulas against an independent ordered-six-atom recurrence.
  The all-exponent result is proved in `smooth/PROOF.md`, not inferred by fitting.
- Electrical construction:26 checks on five atom networks and six literal
  padded composites, with rational resistance and independent path enumeration
  for the small atom inputs. General statements have direct energy/series proofs.

Every program was run normally and with `python -O`; corresponding deterministic
outputs agree. These bounded controls concern new constructors and compilers,
not a re-verification campaign on inherited mathematics.

## Publication

Checkpoint1:ba8d0d22ef3f443a451ec01114738f244407d7fd.
Checkpoint2:03576e3f93527803c6d06cc93e800d7293750b57.
Both were committed to main with expected-head checks, no force push, and each
published file was byte-compared with fetched main. The first checkpoint
normalized three CSV files from CRLF to LF; numerical content was unchanged.
The stale dirty local checkout was not reset, staged or committed. All unrelated
local changes and concurrent remote commits were preserved.

The final synthesis/navigation form a third checkpoint. `MANIFEST.json` hashes
the final packet plus its navigation page, excluding itself. The final fetched
main comparison is required before reporting completion. This session has no
Lean proof claims or software deployment. Repository-wide verification is not
claimed; the final record distinguishes scoped navigation/link checks.
