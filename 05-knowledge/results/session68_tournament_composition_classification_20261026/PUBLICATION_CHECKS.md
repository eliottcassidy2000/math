# First-checkpoint publication gates

- Native classification and literal-line notes accepted by the independent proof auditor;
  exact reviewed hashes are in `session68-audit/CHECKPOINT1_AUDIT.md`.
- Native author files are preserved byte-for-byte in `NATIVE_SOURCE_MANIFEST.json`.
- Each producer checker and the independently written checker are replayed using
  normal and optimized isolated Python in a temporary packet copy. The envelope
  checks the frozen file universe before running and compares exact output bytes.
- Six integrity mutations (source, receipt, missing file, extra file, symlink and
  traversal path) must be rejected under normal Python and `-O`; unchanged controls
  must pass. `verify_mutations.py` operates only in temporary copies.
- Repository gates: `python3 agents/check_docs.py` and `git diff --check`.

The local computations are exact bounded controls, not formalization or all-order
exhaustive enumeration. The complete proof scopes and classical-source attribution
remain in the native notes. Source commit pins are evidence pins, not publication
parents. Online publication and exact-commit CI status are checked separately;
an empty check list is not a passed CI result.
