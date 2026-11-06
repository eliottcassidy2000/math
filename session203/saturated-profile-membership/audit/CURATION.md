# Portable audit curation

This is a publication-ready subset of the completed independent audit.
The original local audit was preserved without editing any of its files.

## Retained premises and controls

The independent checker and original mathematical review are byte-identical
to their audited originals. The frozen candidate proof, both implementation
files, normal and optimized 137 receipts, and both resource-guard receipts
are also unchanged.

The frozen producer manifest is retained. Its three actual dependency files
were added from the producer's dependency snapshot:

- `dependencies/FINITE_INTERFACE_PROOF.md` is the published theorem premise.
- `dependencies/UNIT_START_INVENTORY997.json` is only an independent historical
  inventory control, not a theorem premise.
- `dependencies/PROVENANCE.json` records their exact commits, paths and hashes.

Every one of the ten entries in `frozen_source/FROZEN_SHA256SUMS` verifies.
The dependency proof hash agrees with the original review's cited proof hash.
The independently computed 41 unit-start bases through 137 also agree with
the historical inventory restricted to that bound.

## Omitted material

This curated package excludes the unrelated before/after shared-checkout
status snapshots, all four `frozen_source/PROBE*.json` exploratory diagnostics,
and the unused producer README draft. In particular, the old concatenated
`PROBE1.json` diagnostic is not a production receipt and is not needed for any
test or mathematical premise.

The independent checker never loaded these omitted files except to hash
top-level source files generically. No computational check was removed or
weakened. No candidate proof, implementation, production receipt, or reviewed
dependency was replaced by an exploratory file.

## Regenerated receipts

`INPUT_SHA256.json` describes the retained top-level frozen inputs, including
the producer manifest. That manifest in turn covers the three dependency
files. The package-wide `SHA256SUMS` covers every retained file and dependency.

Both normal and optimized independent audits were rerun in this directory:

```
python3 -B independent_audit.py > RESULT.json
python3 -B -O independent_audit.py > RESULT_OPTIMIZED.json
(cd frozen_source && sha256sum -c FROZEN_SHA256SUMS)
sha256sum -c SHA256SUMS
```

Each regenerated result agrees with its original result after removing only
`input_hashes` and `elapsed_seconds`. The regenerated normal and optimized
results agree after removing only `elapsed_seconds`. `CURATION_CHECK.json`
records these comparisons and confirms that the original audit files stayed
unchanged. This is a packaging change only. The two required release
clarifications in `REVIEW.md` and `AUDIT.json` remain in force until the exact
promoted proof delta is reviewed.
