# Manifest

This compact scratch promotion packet freezes the detached hostile audit that
the canonical THM-4281 final carrier closes `(256,663)`.

## Files

- `src/reconstruct_and_literal_audit.cpp`: carrier reconstruction, direct
  literal partition, independent exposed-body derivation, vulnerable-response
  check, and exhaustive 14,307,150-body scan.
- `results/audit_{O0,O3,UBSan}.out`: byte-identical semantic transcripts.
- `results/audit_{O0,O3,UBSan}.err`: three empty diagnostic controls.
- `results/reconstructed_final8951_carrier_{O0,O3,UBSan}.txt`: byte-identical
  ordered carrier reconstructions.
- `results/reconstructed_vulnerable_bodies_{O0,O3,UBSan}.csv`:
  byte-identical independent 71-row derivations.
- `controls/canonical_dependencies.sha256`: direct canonical inputs and the
  full three-source detached-literal include chain.
- `controls/replay_band_files.sha256`: all 59 canonical replay transcripts.
- `controls/mode_and_dormant_comparison.txt`: compact mode identities, local
  dormant-input comparisons, and singleton edge ledger.
- `README.md`, `REPRODUCTION.md`, and this manifest.

`SHA256SUMS` is generated last and is intentionally not self-listed.  No
binary, dSYM directory, progress stream, full rank-eight active-mask list, or
full active-mask response-pattern list is included.

## External canonical dependencies

Replay requires the repository tree at commit
`6ff6ee322b40e806d8d512baa9fd1df5730ce8b9`, or a later tree in which both
control manifests still pass.  All paths are repository-relative.  The source
uses only C++20 and the standard library.

The dormant packet at
`/private/tmp/thm4281-sig663-active.TwcWXY/endpoint663-active-response-packet`
was used only after the independent run for byte comparisons.  It is not a
reproduction dependency and is not copied here.

