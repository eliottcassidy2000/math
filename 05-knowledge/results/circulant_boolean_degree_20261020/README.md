# Boolean degree boundary for circulant type bijections

This packet contains the unchanged [mathematical proof](BOOLEAN_DEGREE_BOUNDARY.md),
its [validation receipt](VALIDATION.md), reproducible exact controls, and the
[independent proof review](review/REVIEW.md). The proof establishes dyadic-part
degree lower bounds and unbounded necessary degree across the applicable primes.
Its explicit witnesses attain minimum degree two at p=13 and three at p=17.
Finite checks corroborate the proofs; they are not substitutes for all-order
arguments. No sparsity-optimality or exhaustive cubic-census claim is made.

## Reproduce the publication layout

Run these commands from this packet directory, in a writable working copy.
Python 3.10 or later suffices for the Python checks; they use only the standard
library. The optional exhaustive programs require a C++17 compiler supporting
`__uint128_t`.

```sh
scratch=$(mktemp -d)
python3 verify_degree_boundary.py --output "$scratch/verification.json" > /dev/null
python3 -O verify_degree_boundary.py --output "$scratch/verification.optimized.json" > /dev/null
cmp verification.json "$scratch/verification.json"
cmp "$scratch/verification.json" "$scratch/verification.optimized.json"

python3 review/independent_checks.py --output "$scratch/independent_checks.json" > /dev/null
python3 -O review/independent_checks.py --output "$scratch/independent_checks.optimized.json" > /dev/null
cmp review/independent_checks.json "$scratch/independent_checks.json"
cmp "$scratch/independent_checks.json" "$scratch/independent_checks.optimized.json"

python3 audit_basis.py
python3 audit_small_enumeration.py
sha256sum -c SHA256SUMS
(cd review && sha256sum -c REVIEW_SHA256SUMS)
sha256sum -c PUBLICATION_SHA256SUMS
```

The independent checker regenerates `review/input8.independent.txt` and verifies
that it matches the producer input. These commands reproduce 268,536 producer
checks and 249,374 independent checks. Both checkers retain their failures with
Python optimization enabled.

Optional full quadratic enumerations, each visiting all 536,870,912 maps:

```sh
g++ -O3 -std=c++17 review/independent_exhaust8.cpp -o "$scratch/independent_exhaust8"
"$scratch/independent_exhaust8" review/input8.independent.txt > "$scratch/exhaust8.independent.json"
cmp exhaust8.json "$scratch/exhaust8.independent.json"

g++ -O3 -std=c++17 exhaust8.cpp -o "$scratch/exhaust8"
"$scratch/exhaust8" input8.txt > "$scratch/exhaust8.rerun.json"
cmp exhaust8.json "$scratch/exhaust8.rerun.json"
```

The review's original exhaustive outputs, including the three previously
unpackaged rerun JSON receipts, are included. The optional seeded cubic search
is not needed to validate either witness; its original reproduction commands
remain in the proof.

## Packaging and evidence provenance

Only the review checker's producer-directory locator changed for publication:
`PRODUCER = HERE.parent` replaces the original sibling-scratch-directory path.
Normal and optimized runs of this relocated script were completed again. Their
receipts are identical to one another and differ from the original independent
receipts only in `script_sha256`. The producer checker was also rerun in both
modes and reproduced its original receipts exactly. The mathematical proof is
unchanged, with SHA-256
`cbdf5e4566ddec0be6e337302df683d18bf7524db7eff8d5f518a8866b300b55`.

The original independent script, both original receipts, and original review
checksum inventory are preserved under [review/originals](review/originals).
The current `review/REVIEW_SHA256SUMS` uses paths relative to `review/` and hashes
the current relocated checker and regenerated receipts.

`review/SOURCE_SHA256SUMS` is an unchanged archival inventory of the
pre-publication audit snapshot, not a portable current-layout checksum file.
Its original scratch/repository paths are historical, and its VALIDATION.md
hash predates the update recording the completed external review. The proof
source hash did not change. Likewise, the scratch-layout commands in the
unchanged review are historical; use the relative commands above for this
published layout. `SHA256SUMS` checks the current producer payload, while
`PUBLICATION_SHA256SUMS` checks every file in this publication packet except
the checksum file itself.

See [PACKAGING_RECEIPT.json](PACKAGING_RECEIPT.json) for the exact adaptation,
original approved manifest pin, and replay comparisons. No mathematical claim
or original audit evidence was revised during packaging.
