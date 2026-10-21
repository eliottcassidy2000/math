# Independent frozen extension audit

This package audits extension manifest
`a6667f7a99d5fa186ba09b5b6df462c13d084599f1af2c0001aba9f9aff4c390`:

- the exact corridor through `delta<=2k-1`, for `k>=3`;
- the first `delta=2k` owner correction, for `k>=4`;
- separately scoped finite small exceptions and literal marked-port controls.

Read AUDIT.md for the proof review, exhaustive finite universes, inheritance
credit, and historical-source pin handling. No result beyond `delta=2k`
is included. This does not replace the prior main-packet audit.

Portable replay, using Python's standard library:

```
python3 verify.py --source /path/to/session18-zero-heavy-extension --base /path/to/session18-zero-heavy --repo /path/to/math
```

The two source packets must retain their exact frozen manifests and listed
payloads. The repository supplies their pinned dependency files, except
that this audit bundles the one historical base source subsequently edited
by incoming work. That capture is authenticated to the original hash and
commit; the incoming lineage retains its separate new hashes. No git
command, repository mutation, or network access is needed for replay.

The independent checker uses its own generic weighted canonicalizer and
bundled auditor-written support. It imports no producer checker. All
four/five-position cells and collision examples retain actual outside
ports. Included C/D-only results at k=4 are not mislabeled as a complete
mass-layer census.
