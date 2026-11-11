# Reproduce and use the results

Run from the repository root. Python needs its standard library except exact
rank verification, which uses SymPy. Compiled checks require g++ and GMP
headers/libraries. No compiled binaries are published.

    D=05-knowledge/results/session270_chanpak_enumeration_20261111
    python3 "$D/enumeration/cache.py" 250
    python3 "$D/chanpak/compile_target.py" 999999
    python3 "$D/reproduce.py" --full

The first query reads the checked table. The second emits a literal tournament
adjacency list encoded as row bitmasks, its defining word, and its scope. The
full replay validates the final sources, small direct graph controls, rank
certificate, every atlas witness and the independent atlas subset controls.

For fresh counting, exact cache reconstruction and all performance comparisons,
see [enumeration/README.md](enumeration/README.md). To regenerate the exhaustive
477-million-visit word census rather than replay its certificate, see
[chanpak/PROOF.md](chanpak/PROOF.md). Sector benchmark reruns are in
[terminal/PROOF.md](terminal/PROOF.md); they can take over a minute and are not
part of the quick replay.

## Scope of validation

- Full A000568 table: all 251 coefficients through 250 independently agree
  between two exact recurrence implementations. Direct partition checks cover
  all orders through 30, including the historical floating-error boundary14/15.
- Million-target atlas: all 499,998 words replay;626 additional independent
  subset-DP graph recounts. Every graph is a literal two-core/spine construction.
- Positive transfer:2,205 exhaustive small comparisons against subset DP.
- Rank: exact invariant closure and nonzero minor, not just sampled rank.
- Marked carrier: whole small labelled universes canonicalized under actual
  core permutations and independently path-counted.
- Three normalized low bits: direct full Burnside controls at every even
  order4–40; huge-index formulas use supplied factorizations.
- Sector kernel: exact original-gcd multiplicity controls and identical-state
  gauge ablations, with single-run timing status explicit.

Checksums protect artifacts; they are not mathematical proofs. No Lean
formalization or unrestricted full-repository mathematical verification is
claimed. Inherited theorems are credited rather than re-audited.

## Repository-wide documentation check

`python3 agents/check_docs.py` was run before and after this session. Both
runs report the same pre-existing global failure: startup surface180146 bytes
exceeds180000. This session did not alter existing startup surfaces. Its new
S270 navigation page and packet-relative links are checked separately, and
owned staged files pass `git diff --check`. The old global budget issue is not
reported as a pass and unrelated files are preserved.
