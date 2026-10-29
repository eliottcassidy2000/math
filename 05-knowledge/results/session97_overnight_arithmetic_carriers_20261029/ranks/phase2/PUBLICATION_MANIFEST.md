# Phase2 publication manifest

Publish the nineteen payload files below, plus this manifest and
`PUBLICATION_MANIFEST.json` if desired. Payload total: **745553 bytes**.
Do not publish the local HTML cache, raw database, PDF, `__pycache__`, or other
unlisted generated outputs. JSON certificates listed here are essential proof
inputs, not disposable execution receipts.

Runtime: Python3.12; **SymPy1.14.0 only for the six symbolic scripts named below**.
All other code uses the standard library. `named_family.py` is the only helper
module: exclude it from the runnable-check glob. Run the eleven explicit
checks in the JSON manifest, rather than globbing arbitrary Python files.

The first checkpoint is unmodified. This phase is self-contained together
with five ALREADY PUBLISHED files in the parent `ranks/` directory:
`check_characters.py`, `curve273.json`, `curve302.json`, `curve4239.json`,
`specialization_certificates.json`. Preserve that directory layout.

## Exact payload

|File|Bytes|
|---|---:|
|`README.md`|9807|
|`COMMON_A3_CARRIER.md`|2595|
|`check_common_a3_carrier.py`|929|
|`check_elkies_948.py`|1589|
|`check_fibration_geometry.py`|836|
|`check_named_sections.py`|910|
|`check_odd_primary_certificates.py`|2209|
|`check_secant_bit_quotient.py`|1829|
|`count_k3_reductions.py`|1155|
|`verify_highrank_snapshot.py`|1721|
|`verify_k3_counts.py`|953|
|`verify_odd_primary_certificates.py`|3176|
|`verify_all_rank_coverage.py`|1271|
|`named_family.py`|1066|
|`highrank_inputs.json`|502174|
|`all_highrank_character_certificates.json`|172881|
|`odd_primary_certificates.json`|32374|
|`elkies_948_mathdata.json`|7550|
|`k3_reduction_counts.json`|528|

## Dependency graph

- `check_common_a3_carrier.py` → `named_family.py`; SymPy.
- `check_fibration_geometry.py` → `named_family.py`; SymPy.
- `count_k3_reductions.py` → `named_family.py`; SymPy.
- `check_elkies_948.py` → `elkies_948_mathdata.json`; SymPy.
- `check_secant_bit_quotient.py` → `elkies_948_mathdata.json`; SymPy.
- `verify_k3_counts.py` → `named_family.py`, `k3_reduction_counts.json`; SymPy.
- `check_named_sections.py` → `named_family.py` and all five existing parent
  dependencies listed above; standard library.
- `check_odd_primary_certificates.py` → `highrank_inputs.json`; standard library.
- `verify_highrank_snapshot.py` → `highrank_inputs.json`,
  `all_highrank_character_certificates.json`; standard library.
- `verify_odd_primary_certificates.py` → `highrank_inputs.json`,
  `odd_primary_certificates.json`; standard library.
- `verify_all_rank_coverage.py` → both `verify_*primary/snapshot` programs and
  their three JSON inputs. It runs both independent verifiers and checks
  coverage: 64 binary certificates + 2 ternary repairs = all66 curves and
  2005 independent point directions; 140334 verifier checks. Standard library.

Every file read is resolved from `__file__`; ordinary imports find their
same-directory helper by Python script-directory semantics. Symbolic scripts
without data files still import only that same-directory helper. No script
requires the caller's working directory to be this folder.

## Fresh-copy test

Copied precisely these nineteen payload files and the five existing parent
dependencies into a new isolated directory. Ran all eleven checks with the
caller working directory `/tmp`, both normally and under `python -O`.
All eleven exited successfully and produced byte-identical stdout across the
two modes. Exact stdout SHA256 and file SHA256 values are in the JSON manifest.
No raw database, downloaded HTML, PDF, generated receipt, or ambient source
checkout was present in the isolated copy.

## Scope guards for synthesis

- The 32768 secant choices are bases of the rank17 **generated subgroup**;
  saturation in the full Mordell–Weil group is not independently proved.
- The common A3(2) sublattice is not claimed primitive.
- K3 inequivalence is overQ, not overQbar.
- The phase2 README's binary-certificate gap in§5 is explicitly closed by§6.
- These files prove rank lower bounds, not unconditionally exact ranks.
