# Reproduction and release scope

Run from the repository root. Python3 uses standard-library exact arithmetic.
The completion comparison also requires g++ and GMP development libraries,
which are used by the explicitly inherited C++ counters.

```bash
P=05-knowledge/results/session283_grammar_20261112
python3 "$P/enumeration/residual_compiler.py"
python3 "$P/enumeration/thresholds.py"
python3 "$P/enumeration/complete.py" --cases 64:11,100:19,150:27,150:35,200:49,250:59
python3 "$P/grammar/verify.py"
python3 "$P/metagraph/verify.py"
python3 "$P/paths/verify.py"
python3 "$P/paths/switches.py"
python3 "$P/audit/independent.py"
python3 "$P/audit/enumeration_independent.py"
python3 "$P/audit/switch_independent.py"
```

Repeat with `python3 -O`. Every mathematical gate uses an explicit exception,
not an assert removed by optimization. Deterministic receipts match exactly;
performance-bearing receipts match after removing only timing fields. Timing
is informational and is not a proof or a machine-independent cost guarantee.
The exact source snapshots and comparison checks are in the independent audit.

The complete pipeline constructs the interval and its independent low-bit
residue before computing/reading the inherited full answer. The order150,
200 and250 interval-only cases use no low-bit oracle at all. Exact-cache
agreement is an after-the-fact control. This is not a discovery claim for
these already known values.

## Dependencies and boundaries

- Session279 ground-moment consumer is an explicitly inherited proved kernel;
  the independent audit also reconstructs a separate consumer.
- Session277 modular counter and Session270 full counter/cache are inherited,
  hash-pinned in `enumeration/completion.json`.
- The new local metagraph proofs rely on the inherited homogeneous-pair path
  forest and unique weighted-core decomposition.
- Path-profile contraction and triangle-cap coefficients are inherited; the
  superincreasing single-probe decoder and scoped switch deductions are new
  applications.

The full repository documentation checker reports two unchanged pre-existing
budget failures: META-PATTERNS has401 maintained lines against400, and the
startup surface has180862 bytes against180000. This packet does not modify
those maintained surfaces. No full-repository test or dependency re-audit is
claimed.

A local isolated checkout initially showed33 tracked historical message files
as EOL-only dirty (the diff ignoring end-of-line whitespace was empty). They
were neither staged nor altered by this session. Publication stages only this
packet and its own navigation route; unrelated work is preserved.
