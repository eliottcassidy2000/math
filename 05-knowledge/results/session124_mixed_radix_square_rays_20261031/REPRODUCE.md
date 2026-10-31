# Portable reproduction

Python 3 standard library only. No MILP solver, SciPy, network, Git database,
or discovery script is needed.

Place this directory beside the existing repository dependencies:

- `square_sum_B_base_interval_20261019/` (compressed old base bundle);
- `session124_square_prefix_threshold42_20261031/SMALL_BASES.json` (eleven small pairs);
- `session124_recursive_square_prefix_20261031/` (Q49 certificates, only for operator C).

The base bundle and Q49 tables are referenced, not copied. `DEPENDENCIES.json`
pins the small-base file and all referenced Q49 certificate hashes;
`BASE_DEPENDENCIES.json` pins the base manifest and compressed/decoded chunk hashes.
`CERTIFICATES.json` pins the 200 new Q25 certificates included here.

## Constructor and producer checks

```sh
python construct_cycle.py B
python construct_cycle.py ABC
python construct_cycle.py AB 26563
python verify_package.py > replay.normal.json
python -O verify_package.py > replay.optimized.json
```

The optional integer is the interval order N. If omitted, the exact certified
threshold for the requested word is used. `core.py` also exposes `prefix(word)`
and `threshold(word)`. The constructor verifies the entire cycle while streaming
its vertices; its receipt gives a SHA256 of the decimal-lines stream without
storing a huge cycle file. No construction outside the specified theorem domain
is claimed.

For a different checkout layout, all constructor/producer CLIs accept
`--base-dir PATH --small-base-file FILE --q49-dir PATH`. A/B-only words do not
read Q49 certificates; an explicit nonexistent Q49 directory is tested.
The empty schedule suffix uses an A-only weak-B supply, preserving this property.

## Independent universal audits

The unchanged independent scripts import no producer module. Run from this directory:

```sh
python A/audit_independent.py --pair-dir A/paired --cycle-dir A/cycles \
  --small-base-file ../session124_square_prefix_threshold42_20261031/SMALL_BASES.json \
  --base-dir ../square_sum_B_base_interval_20261019
python B/audit_independent.py --pair-dir B/paired --cycle-dir B/cycles \
  --small-base-file ../session124_square_prefix_threshold42_20261031/SMALL_BASES.json \
  --base-dir ../square_sum_B_base_interval_20261019
```

Repeat with `python -O` to reproduce optimized-mode controls. Do not pass
`--partial` to B's script when checking the complete theorem. Each full audit
checks all 50 paired and 50 terminal states, all 1,020 base pairs through 1061,
100 paired expansions, and 103 streamed cycles including N=16,601,562.
Their receipts are in each family's `independent.normal.json` and
`independent.optimized.json` and are byte-identical across modes.

The scripts are named `audit_independent.py` and the receipts
`independent.normal.json`/`independent.optimized.json` in this package. Script
and receipt bytes were copied unchanged; prose artifact references were
adapted to these portable names.

## Scope and dependency boundaries

A and B use base pairs 42–1061 and recurse from 1062. C uses bases 42–2081
and recurses from 2082. Returning arbitrary weak bases above these cutoffs
would invalidate the stronger prescribed-prefix induction. All operators
consume and produce the same weak B type; no historical uniform head/tail
shift property is required of the inputs.

The new ray prefixes differ from the earlier Q49 ray except where explicitly
shared. The package does not claim arbitrary prescribed-path extension,
simultaneous embedding of incompatible prefixes, or tournament-spectrum
infinite completeness. The latter remains OPEN.

## Independent ternary composition and frontier audit

```sh
python ternary/audit.py --a-paired A/paired --a-cycles A/cycles \
  --b-paired B/paired --b-cycles B/cycles \
  --c-paired ../session124_recursive_square_prefix_20261031/paired \
  --c-cycles ../session124_recursive_square_prefix_20261031/cycles \
  --small-base-file ../session124_square_prefix_threshold42_20261031/SMALL_BASES.json
```

Repeat with `python -O`. The unchanged independent auditor uses closed-form
prefixes rather than the producer's recursive routine, checks 30 mixed cycles,
exact decoding and antichain counts, and shares all 27 producer cycle hashes.
Its 31,877,061-check receipts are in `ternary/normal.json` and
`ternary/optimized.json`; [the independent note](ternary/AUDIT.md) also proves
the coded-family count's precise exponent and checks the portable constructor.
This exponent is not an upper bound on all Hamiltonian cycles.

Finally run `python verify_integrity.py` to check the published whitelist and
all non-self-referential file hashes.
