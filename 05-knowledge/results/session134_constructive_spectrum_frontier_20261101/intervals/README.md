# Session134: exact interval-coverage arc-budget gate

**PROVED necessary complexity bound, not a spectrum closure.**

[THEORY.md](THEORY.md) gives an exact core-order cap with at most b cyclic
SCCs, a finite binomial bound on values realizable with at most b backward
arcs, and the uniform consequence that density-one odd targets need
`(1-o(1)) log n/log log n` backward arcs in any realization.

For example, six backward arcs cannot cover an entire dyadic interval
already at `(2^63,2^64)`, by this counting certificate. The all-size
six-arc deficit-nine construction is retained and independently tested;
there is no contradiction and no globally forbidden new value.

## Reproduce

Python 3 standard library only; no repository, network, packages or data
files are required for the verifier. It also runs from another directory.

```sh
python3 /path/to/verify_arc_budget.py > normal.rerun.json
python3 -O /path/to/verify_arc_budget.py > optimized.rerun.json
cmp normal.rerun.json optimized.rerun.json
cmp normal.rerun.json /path/to/normal.json
python3 /path/to/verify_packet.py
```

The theory depends on the cited exact Moon–Busch strong minimum and the
pinned prior core work. Finite tests verify implementation and controls,
not the all-order input theorem. `SOURCE_PINS.json` identifies inherited
repository files and primary-source access; `MANIFEST.json` seals this
packet. Discovery-only `explore.py` is deliberately not a deliverable.
