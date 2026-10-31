# Session 125: hidden phase and successor carry costs

**Proved structural results with independently checked finite inputs.**
The graph has vertices 1 through n and an edge when two distinct labels sum
to a square. This packet studies the original ports 1, 3, 8 and Gerbicz's
whole-block radix-49 constructor; it does not claim new all-order cycle
existence, a Zeckendorf identification, or an obstruction for arbitrary cycle
families.

Read [the proofs](SUCCESSOR_PHASE_AND_CARRY.md) and
[the independent audit](audit/AUDIT.md).

| Result | Scope |
|---|---|
| A globally fixed position phase needs at least six balanced defects | Any original-port phase-coherent path family |
| No phase-coherent family spans orders 41–44 | Complete three-mask catalogue and short degree obstructions |
| Independent consecutive nice pairs have O(log n) deleted-edge cost | Valid original-port Q49 whole-block construction |
| Every Shorter/Longer bank selector has limsup edit ratio at least 4/87 | Any valid base pairs and valid whole-block words in this construction |
| Every bank selector has exact global limsup 22273/24042 | **Only the pinned original tables**, not arbitrary graph families |

A bank selector chooses at order n either the Shorter member of Pair(n) or
the Longer member of Pair(n−1). Closing a path adds its port edge. The key
carry recurrence is `H(49r+24)=49H(r)+K`, with nonnegative connector
correction; for the original tables K is 40 or 52 according to parity.

## Inputs and reproducibility

Python 3 standard library only. Two kinds of input are deliberately separate:

1. **Published sibling origin packet:**
   `../session125_origin32_numeration_20261031/origin_pairs/` contains the
   finite degree-subset checker and gzip/base64 certificate. The phase checker
   calls the published authenticated `read_reference()` helper, then recomputes
   the relevant catalogue. It does not need an unpacked certificate JSON.
   `--origin-dir` overrides the default sibling location.
2. **Separately obtained original C source:** required only for the original
   table audits. It is not in this packet. [SOURCE.json](SOURCE.json) gives
   public URLs, SHA256, attribution, and the prior reconstruction. No license
   permitting source redistribution was found. Treat the C file as inert data:
   **do not compile or execute it**. Both auditors reject any other whole-file
   digest before parsing literal integer initializers.

[DEPENDENCIES.json](DEPENDENCIES.json) pins the sibling inputs and source.
The universal proof does not depend on Gerbicz's particular base choices;
its finite three-mask input can also be independently recovered by the audit's
endpoint search.

From this packet's directory:

```sh
python3 check_fixed_phase.py --output phase.normal.json
python3 -O check_fixed_phase.py --output phase.optimized.json
cmp phase.normal.json phase.optimized.json

python3 q49_edit_audit.py /path/to/squares.c --output q49.normal.json
python3 -O q49_edit_audit.py /path/to/squares.c --output q49.optimized.json
cmp q49.normal.json q49.optimized.json

python3 audit/independent_audit.py /path/to/squares.c --output audit.normal.json
python3 -O audit/independent_audit.py /path/to/squares.c --output audit.optimized.json
cmp audit.normal.json audit.optimized.json

sha256sum -c SHA256SUMS
```

The included normal/optimized receipts match. [RELEASE_CHECK.json](RELEASE_CHECK.json)
records relocated runs using the default sibling layout, source-rejection
controls, the exact file whitelist, and cross-audit comparisons. The independent
audit imports no producer module and reconstructs the phase catalogue by a
different search. No external source, archive, arrays, executables, caches,
or unpublished origin-certificate JSON are included.
