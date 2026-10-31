# Portable replay receipt

The curated packet was copied to a fresh temporary directory and its absolute
`verify_all.py` path was invoked from an unrelated `/tmp` working directory.
All eleven checks passed normally and with `-O`; each matched its bundled
receipt. The complete replay summary was byte-identical before and after
relocation and is saved as `receipts/full_replay.json`.

The packet has no workspace-specific paths, remote dependencies, floating LP
requirements, certificate-generation requirements, or imports of exploratory
research scripts. Six certificate/modular checks require only Python's
standard library; the complete transfer and independent-audit replay also
requires SymPy1.14.0. Verification is read-only unless `--write-receipts` is
explicitly supplied.

The source research history remains outside this curated whitelist. The
obsolete fixed-z experiment is not executed or shipped. Earlier open-length
statements are explicitly historical/superseded; the currently untested
longer-macro frontier in this packet is length≥5.
