# Independent endpoint-586 closure audit

**FINITE-EXACT / independently verified relative to the serialized THM-4318
carrier. LRC(14) remains open.**

The endpoint-586 frontier was reconstructed from the frozen prior typed
partition, rather than accepted as an input list. It has twelve rows and FNV64
`a1b617faa2e7f63f`. An import-free C++ implementation rebuilt each pair wall,
recomputed activity of the serialized 3,925-mask THM-4318 carrier, and scanned
all `12*C(30,9)=171,685,800` labelled bodies under both O2 and O3.

The two optimization modes agree byte-for-byte. They find 4,090 failures, all
at `(50,586)`, with body FNV64 `14ce094f4ab4ba94`, aggregate failure FNV64
`ffb884b2b17e6ef4`, pair-ledger FNV64 `46a2c17caecc55df`, 525,048 exposed
bodies, and 12,125,735 hit incidences. The pair and failure CSVs are also
byte-identical to the scout packet.

A separate Python implementation rebuilt the 8,381 open wall cells and kept
all 5,802 pair-safe cells unaggregated while evaluating every failure body.
It finds 2,424 full failure classes and 2,371 rank-at-most-nine classes. Every
body has strictly positive truncated mass surplus

```text
63*L(B)-4D > 0.
```

The minimum is `136976666519138544` at body `011c6405`; truncated and full
masses agree for 3,602 bodies. Normal and `python -O` runs are byte-identical,
and the 4,090-row detail CSV is byte-identical to both scout implementations,
SHA256 `6dbfc6178c7da5c844abc95e980f2696cf02c3845a83bdfff95f881c12c9c4fd`.

Independently consuming the reconstructed frontier gives typed union 2,229
(FNV64 `035ebf12f02ecc62`), residual 20,418 (FNV64
`89b73e31224821c4`), and next endpoint 585 with twenty-three rows (FNV64
`8f1b7c8db8fd5e87`). These generated files are byte-identical to the scout
successor files.

This closes endpoint 586 only in the frozen finite pool/typed model. It does
not provide a physical owner/entry map, a terminating descent, or a proof of
LRC(14).
