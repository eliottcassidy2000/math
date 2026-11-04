# Complete named marked-profile reduction

The proof gives an all-order sufficient finite-data reduction modulo eight for arbitrary fixed quotient substitutions, retaining both physical ports, shared names and the complete cubic/quartic response. It does not prove universal noncancellation or a collision.

Reproduce with `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The standard-library verifier compiles the local C++17 run-count program with the installed compiler; no packages or installs are needed. Every mathematical check is explicit and remains active with `-O`.

Controls include 78 actual child tournaments (all labeled orders1..4 plus TT5, TT6 and an ordered sum of two C3s), complete endpoint-prescription tables with at most four contact roles, and 96 full physical queries on six four-variable faces. The physical controls retain exact counts by first-child run number and all face coefficients, including cubic/quartic terms and shared physical names. Every exceptional-name cutoff has nonzero high-run controls; their cancellation is divisibility, not an absence of words.

Frozen allowlist: PROOF.md, README.md, verify.py, runs.cpp, normal.json, optimized.json, SHA256SUMS. Independent review is a separate artifact.
