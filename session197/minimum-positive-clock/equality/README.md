# Unrestricted amplitude-72 equality packet

**PROVED and independently VERIFIED.** There is one rooted strong equality core for the transitive-root-clone family H_m=72*3^m+1. Every unrestricted equality source is that rigid nine-vertex core with an arbitrary transitive prefix and suffix. At each n>=9 this gives n-8 unlabelled sources and (n-8)n! labelled sources/pairs. No self-converse hypothesis is imposed.

- [Exact frozen proof](PROOF.md)
- [Independent review](REVIEW.md)
- [Global minimum theorem](../PROOF.md)
- [Machine-readable source receipt](receipt-normal.json)
- [Exact file manifest](MANIFEST.json)

The proof and review retain their original scratch-era source names so their reviewed hashes remain unchanged. In this published layout the preceding minimum proof is `../PROOF.md`; the original name `../minimum-positive-clock/PROOF.md` in the frozen equality proof refers to that same theorem. Other scratch-era names identify provenance, not additional runtime dependencies.

## Reproduce the complete equality search

From this directory, run:

```bash
python3 verify.py > regenerated-normal.json
python3 -O verify.py > regenerated-optimized.json
cmp regenerated-normal.json regenerated-optimized.json
cmp regenerated-normal.json receipt-normal.json
```

The verifier requires Python's standard library and a C++17 compiler, regenerates its bank, and keeps temporary files in this directory. The checked source bank is included as `strong-bank.jsonl`; it is not needed to discover or validate the result.

## Reproduce the independent controls

`audits/positive-equality/check_catalog.py` compares the unique core against the included unchanged nine-vertex rooted catalog, and exhaustively checks its unmarked automorphisms. Run it with Python and compare its output to the adjacent `catalog-check.json`.

The independent `gap_check.cpp` takes a plain text stream: first the number of bank records, then one line per record containing its n, H and adjacency rows. Generate this stream from `strong-bank.jsonl`, compile the checker with a C++17 compiler, and compare its stdout to `audits/positive-equality/gap.json`. It reconstructs root-run polynomials from exterior path gaps, rather than counting explicit root clones.

All frozen proof, verifier, receipt, bank and audit bytes are bound by `MANIFEST.json`. The two transitive padding lengths are independent. Equal padding is a separate, narrower self-converse result.
