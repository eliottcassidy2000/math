# Two fixed ears: exact transfers and limits of scalar accumulators

This portable research package studies tournaments consisting of two marked
vertices x→y and an ordered transitive spine. It does **not** prove the full
Hamiltonian-path-count attainment conjecture or exclude any new integer count.

## Certified conclusions

* An exact nonnegative36-state transfer, independently checked against complete
  run partitions and direct adjacency/subset Hamiltonian-path counting.
* A common triangular filtration: every length-L word has characteristic
  polynomial `z^8(z−1)^16(z−2^L)^10(z−3^L)^2`.
* No complete stateless centered affine accumulator `H→BH+d` for **2≤B<243**
  on this carrier, even if unequal template lengths were initially allowed.
  All binary powers are excluded; the certified ternary cases are3,9,27,81.
* Exact occupancy relations and separating certificates for the stated finite
  template families. Complementary two-letter periods have a forced nonzero
  linear spectral term, established by an invariant modulo3.

“Complete stateless” means one fixed physical template for every digit
`−(B−1),−(B−3),…,B−1`, with its assigned action valid at every descendant of a
seed with H>1. The separate binary theorem also treats explicitly defined
finite-phase equal-length systems. Do not extend these conclusions to changing
ears, general tournament gadgets, partial digit domains, arbitrary changing
physical banks, or unrestricted count-dependent rules.

## Read in order

1. [Exact transfer and binary obstruction](TWO_EAR_TRANSFER_AND_COMPILER_OBSTRUCTION.md)
2. [Ternary certificates and common spectrum](TERNARY_AND_SPECTRAL_RESTRICTIONS.md)
3. [Radix27/81 certificates and combined small-radix theorem](RADIX27_REPEATED_TEMPLATE_OBSTRUCTION.md)
4. [Occupancy and modular lower modes](OCCUPANCY_RELATIONS_AND_MODULAR_MODE.md)
5. Independent reviews: [initial audit](audit/AUDIT.md), [extended audit](audit/EXTENDED_AUDIT.md)

The notes preserve correction lineage and clearly identify historical finite
searches. Exploratory LPs, seed censuses, obsolete affine-ray calculations, and
numerical solver dependencies are intentionally **not** shipped. Exact witness
verification—not the success status of a floating solver—is the proof gate.

## Reproduce

Python3 is required. The complete replay additionally uses SymPy1.14.0, listed
in `requirements.txt`; no network access or external data is needed at runtime.

```sh
python -B verify_all.py
```

This runs eleven checks normally and with `-O`, compares each output pair, and
checks both against the bundled exact receipts. It works from any current
working directory. Default verification does not modify package files.

For the standard-library-only certificate and modular checks:

```sh
python -B verify_all.py --stdlib-only
```

To deliberately regenerate bundled receipts after reviewing a code change:

```sh
python -B verify_all.py --write-receipts
```

The independent audit scripts reconstruct matrices and recurrence rows without
importing producer code. The lightweight `code/verify_*.py` certificate checkers
use integer/Fraction arithmetic and include hostile mutation checks where stated.

## Integrity and transport

`WHITELIST.txt` is the explicit publication inventory. `SHA256SUMS` hashes every
listed payload file except itself and the transport-only `BATCHES.json`:

```sh
sha256sum -c SHA256SUMS
```

`BATCHES.json` proposes ≤100,000-byte raw-file batches, keeping individual remote
publication requests below the140KB transport limit after ordinary encoding.
Verify actual serialized request size before sending. The largest witness file
is under84KB. Neither the manifest nor batching metadata is a theorem dependency.

No source worktree, global ledger, branch, or remote was modified while preparing
this package. Publication is owned by the parent research session.
