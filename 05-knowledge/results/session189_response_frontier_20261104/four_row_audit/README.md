# Four-row independent audit

See `REVIEW.md` for the exact approved theorem and its limits. The audit reads the frozen producer packet without importing or executing producer code and never modifies it.

Reproduce with:

    python3 audit.py > normal.json
    python3 -O audit.py > optimized.json
    cmp normal.json optimized.json

An optional first argument identifies the producer directory. Reviewed dependencies must be present as sibling directories. All code uses standard-library Python and the installed C++17 compiler.

The explicit publication allowlist excludes the generated `recount` executable and logs. Complete endpoint, anchor, run/cycle, literal physical, affine-certificate, and named-cube verification data are retained. The profile proof covers unbounded child sizes; the literal finite graphs are verification controls and interpolation inputs, not a replacement of arbitrary children by small ones.
