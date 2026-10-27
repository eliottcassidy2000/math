# C3 portable replay adapter

`c3_replay.py` replays the **expanded original logical files**, not gzip
transport paths. Run the checkpoint transport verifier/materializer first.
Python 3 (standard library only), a C++17 `g++`, and an ordinary writable
temporary directory suffice; no nauty installation or repository checkout is
needed for these replay paths.

```sh
python3 c3_replay.py --self-test
python3 -O c3_replay.py --self-test
python3 c3_replay.py --materialized /tmp/c3-logical --work /tmp/c3-normal --mode normal
python3 -O c3_replay.py --materialized /tmp/c3-logical --work /tmp/c3-optimized --mode optimized
python3 c3_replay.py --materialized /tmp/c3-logical --work /tmp/c3-ubsan --mode ubsan \
  --component tc --component tc-window --component boundary --component boundary-audit
python3 check_cpp_rejections.py --normal /tmp/c3-normal \
  --optimized /tmp/c3-optimized --ubsan /tmp/c3-ubsan
python3 -O check_cpp_rejections.py --normal /tmp/c3-normal \
  --optimized /tmp/c3-optimized --ubsan /tmp/c3-ubsan
```

Every `--work` directory must be new and outside `--materialized`. Compilation,
generated ledgers, modified adapters, and run logs stay there. The source packet
is not modified. No original `/tmp/session75-*` directory is needed.

For parallel replay, repeat `--component` with any of `boundary`,
`boundary-audit`, `twins`, `twins-audit`, `root-twins`, `joint`, `joint-audit`,
`tc`, or `tc-window`, using a different work directory per invocation. A
complete normal/optimized acceptance needs all nine components in each mode.
The native command timeout defaults to two hours **per command** and fails
closed; it never omits a slow source or scan job.

## Exact scope and comparisons

- **Boundary producer and independent audit:** authenticate every ordinary
  tournament class at orders 4–9 using C++, then rerun every frozen finite
  witness, family, survivor, and rejection-control receipt. All regenerated
  stdout is byte-compared with the frozen matching output. UBSan redoes the
  complete C++ universe authentication; it does not claim a new Python audit.
- **Twins:** the producer, separate auditor, and root checker each rerun only
  their declared two order-10 sources. The audit includes the complete named
  cut comparison and the distinct odd-cycle collection polynomials. This is
  not an additional census or equality of labelled triangle sets.
- **Joint producer:** regenerate all 19 complete native ledgers and every one
  of the 245,826 declared rectangles, the literal controls, pair-birth control,
  boundary invoice and mechanism reports. All regenerated files and the
  summary are byte-compared. This producer replay is not an independent
  numerical audit of every rectangle.
- **Joint independent audit:** reconstruct all 1,000,908 source records,
  structurally check all 245,826 receipt rectangles, independently numerically
  recompute 42,425 rectangles, and check 43 direct cutoff controls. The whole
  independent and supplementary JSON receipts must be byte-identical.
- **TC analytic/Cayley:** rerun the exact repaired analytic/group checks,
  all 526 Cayley connection sets, 16 degree-15 substitutions and the inherited
  degree-15 source program. The original multi-million order-12 scan is retained
  historical evidence and is **not rerun** by this adapter.
- **TC window:** regenerate the complete 83-job original window and 45-job
  alternate order-12 gate from scratch, reconstruct all orbitals, regenerate
  and compare the portable exact automorphism scanner source, and rerun every
  job. Counts are 185,642 and 46,968 strong generated tournaments respectively.
  Exact automorphism controls compare against all 33,867 labelled tournaments
  through order six plus cyclic/wreath positive controls. Only variable
  per-run elapsed times and raw JSONL hashes are removed from the summary
  comparison; every original-window scanner line's entire mathematical suffix
  is separately matched by the frozen summary checker. UBSan additionally
  checks the specified order-22 module and order-16 nonmodule jobs, not the
  full 128-job window under instrumentation.

## Two explicit compatibility adapters

1. The frozen independent joint audit authenticates four inputs by `git show`
   at commit `96e35642110104cc547516a63f37eebcf930f707`. Those exact four source
   snapshots are already embedded. An isolated adapter subprocess intercepts
   **only** the exact command shape, sentinel repository name, commit and four
   paths, then returns bytes whose hashes match the unchanged `SOURCE_PINS`.
   All other calls fail. This is an **offline snapshot integrity check**, not
   independent evidence of repository origin; the checkpoint's separate
   source-origin receipt establishes that binding. The audit source hash is
   checked before execution and no mathematical code is changed. Five hostile
   tests reject a wrong repository sentinel, commit, path, command shape, and
   corrupted snapshot, under both normal Python and `-O`.
2. The frozen TC `verify_receipts.py` names an original `/tmp` sibling. In a
   disposable copy, one and only one expected literal is changed to the
   supplied materialized sibling. The original source hash and the adapted
   hash are recorded. This auxiliary command checks **frozen receipt
   consistency**; fresh analytic/Cayley computations and fresh TC-window scans
   are separate checks above. It must not be described as independently
   rerunning the original large order-12 computation.

`check_cpp_rejections.py` separately mutates the authenticated ordinary-universe
input streams: outside order, wrong automorphism count, duplicated class,
deleted class, and deleted whole order. Both C++ authenticators must reject all
five with their declared failure message under normal, `-O3 -DNDEBUG`, and
UBSan builds: 30 fixed reject-path controls. A crash or sanitizer error does not
count as the expected controlled rejection. This is not a generalized parser
correctness claim. Its entire receipt must agree under Python normal and `-O` against the same
three compiled workspaces. Observed executable hashes may differ across
workspaces or compilers (especially UBSan debug paths); cross-workspace binary
identity is not required or claimed.

The C3 driver hashes its own source at startup and fails if that source changes
before completion. `REPLAY.json` records commands, hashes, exact scopes and the
five adapter rejection controls. Compiler output and sanitizer errors are
fatal. Normal C++ uses `-O2`; optimized C++ uses `-O3 -DNDEBUG`; UBSan uses
`-O1 -g -DNDEBUG -fsanitize=undefined -fno-sanitize-recover=undefined`.
