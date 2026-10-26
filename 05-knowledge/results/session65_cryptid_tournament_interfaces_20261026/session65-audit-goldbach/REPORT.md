# Session65 independent Goldbach and multiperfect audit

2026-10-26. **PASS within the scopes below, with one prose provenance correction applied.**

* [Goldbach report](goldbach/REPORT.md): exact labelled-control isomorphism, independent native prefix and comprehensive local bounce-padding checks; source-inspected and upstream-CI-observed reduction, **not** a local Lean rebuild or axiom audit.
* [Multiperfect report](multiperfect/REPORT.md): elementary exact arc-family transport and odd valuation ledger, independently checked by divisor census; **not** odd nonexistence or sufficiency of the necessary Euler form.

The producer `session65-twentyfive/REPORT.md` previously called the 21-state development's endpoint “Mathlib's GoldbachConjecture.” At the parent's explicit request this was changed to the project-defined `ConjectureMachines.GoldbachConjecture`, using Mathlib's `Nat.Prime`. The producer report entry in `MANIFEST.json` and affected `SHA256SUMS` entries were refreshed. This was prose-only: no producer checker, machine source, result JSON or execution receipt changed. Audit input pins reflect the corrected report and the root's additional r=10 example.

Both independent authored checkers use explicit exceptions, not Python assertions; normal and optimized outputs are byte-identical. Neither imports nor executes the audited external software or producer checkers. No repository file, Git ref, commit or push was changed by this audit.

## Evidence and reproduction

Run from this directory, passing the two source roots if they are not in the original sibling layout:

```sh
python3 -B goldbach/check.py ../session65-twentyfive > goldbach/normal.json
python3 -O -B goldbach/check.py ../session65-twentyfive > goldbach/optimized.json
cmp goldbach/normal.json goldbach/optimized.json
python3 -B multiperfect/check.py ../session65-root/multiperfect > multiperfect/normal.json
python3 -O -B multiperfect/check.py ../session65-root/multiperfect > multiperfect/optimized.json
cmp multiperfect/normal.json multiperfect/optimized.json
sha256sum -c SHA256SUMS
```

`INPUTS.json` records every sealed producer input. `REPO_PROVENANCE.json` records historical-input comparisons (13 byte-identical and four differing only in trailing newlines) to `eliottcassidy2000/math@6a8ad5d3bf01a232da4de41bcbff68390a9337ee`, independent of the moving read-only checkout. `SOURCE_RECEIPT.json` records fresh read-only upstream observations, the exact successful CI run, and the correction's primary definition. The source producer has 48 files; its SHA256 list and manifest each cover 45 payload files, with the three inventory/seal files excluded by its convention.

The Leng source includes an explicit Apache-2.0 `LICENSE` at the pinned commit, and the saved DeepMind source includes its license and copyright headers. Retain these when redistributing saved source. The producer source selection is deliberately partial: it contains the local Tm25 import cone but not the Tm31 files imported by the root or the toolchain/mathlib dependencies. It is an audit snapshot, not a complete standalone Lean build bundle.

## Remaining verification boundary

Neither `lean` nor `lake` is installed in this executor. Local Lean compilation and endpoint `#print axioms` were **NOT RUN**; no toolchain was installed and no downloaded software was executed. A read-only web rendering of the CI page initially returned a cache miss; the supported GitHub connector independently returned the exact run and successful build-step metadata. This resolved the observation gap, not the local-build gap. No permission-denial retry or unresolved access block occurred.
