# Independent physical and structural audit of the 27-name boundary

This audit independently counts one actual coarse separator for every one of the 2,348 endpoint-zero prime-nine converse-merged representatives with nine C3 children. It uses exact physical-subset orbits under the C3 rotations, not the producer's run-word recurrence. The companion census audit independently establishes exhaustive coverage. The review combines this finite certificate with the separately reviewed original-module structural restrictions to exclude nontrivial source-relative collisions on at most 27 old names.

## Reproduction

Use the installed C++17 compiler and Python standard library only. From this folder:

```
g++ -std=c++17 -O3 -Wall -Wextra -Werror orbit_physical.cpp -o orbit_physical
g++ -std=c++17 -O3 -Wall -Wextra -Werror literal_physical.cpp -o literal_physical
python3 controls.py > controls.normal.json
python3 -O controls.py > controls.optimized.json
cmp controls.normal.json controls.optimized.json
python3 verify.py --recount > normal.json
python3 -O verify.py > optimized.json
cmp normal.json optimized.json
```

The 2,348-pair large recount takes approximately five minutes in the audit environment. The complete literal controls are much smaller. `verify.py` without `--recount` validates retained large counts, every physical row and mask, dependency pins, and the separate census certificate. No check depends on Python assertions. The census sibling can use either its original local name or the publication name `c3_census_audit`.

## Retained evidence

- `REVIEW.md`: derivations, scope, and structural integration.
- `orbit_physical.cpp`: independent exact modulo-eight physical recurrence, 1,048,576 subset orbits per 29-vertex side.
- `all_recount.in/out`: all 2,348 source/target witness inputs and independently computed residues.
- `physical_records.json`: complete actual augmented source and target rows, physical child sets, port names, masks, and residues.
- `literal_physical.cpp`, `controls.py`, `small_physical_controls.json`: 1,160 exact literal physical graphs, including common relabeling and arbitrary child converse controls.
- `input/`: frozen residual candidate list, fixed battery, first-witness records, reviewed dependency proofs and reviews.
- `normal.json`, `optimized.json`: byte-identical final audit receipts.
- `PUBLICATION_FILES.json`, `SHA256SUMS`: exact publication allowlist and content checksums.

The large independently verified counts are residues, not claimed full integers. The small controls retain exact integer counts. Only one separator per residual representative is independently recounted; this audit does not claim to independently recompute the producer's entire twelve-query matrix. No new arbitrary-ambient or arbitrary-child prime-nine theorem is claimed.
