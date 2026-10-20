# Verification scope and freeze record

20 October 2026.

## Structural proof review

A separate reviewer checked and approved:

- the degree-eleven transitive-image proof from the four cyclic-tournament representatives, neighborhood rigidity, and distinct in-neighbor signatures;
- the complete subgroup assembly with exactly 8 ambient-conjugacy classes at degree eleven and 10 at degree twelve;
- the diagonal/product fixed-cube identity, exact normalizers, pronormal-tower argument, small-f multiplicities, and type-count factors;
- the all-f relative laws in README §5.5, explicitly without inferring universal equality;
- FAMILY_SCOPE.md's extension of equality to f≤p−4 via the inherited one-quotient-prime theorem, and its distinction between four selected degree-13/14 probes and a complete census.

The color-partition formulas provide an additional direct arithmetic control. The main finite counts do not depend on the separate, still-candidate double-triple descent argument, which is excluded from this publication inventory.

## Computations completed

1. `prime11_certificate.py`: all 32 skew connection sets; independent out-neighborhood permutation/signature reconstruction. PASS, with automorphism-order histogram 11:30 and 55:2.
2. `extension_census.py`: all 18 structurally classified degree-11/12 rows. PASS. Normal and optimized-Python JSON outputs agree byte for byte.
3. `pair_action_replay.py`: all 18 rows recomputed with uncompressed pair actions and the original unbucketed isomorphism classifier. PASS. It compares all scalar counts, all three multiplicity histograms, normalizer orders, pair-orbit counts, and fusion/escape defects. The largest row took about 761 seconds; it completed, rather than being omitted for cost.
4. `relative_burnside_control.py`: 12 rows' relative tournament-orbit counts recomputed solely by odd-partition Burnside sums, importing no fixed-cube or isomorphism engine. PASS. Normal and optimized outputs agree.
5. `family_scope_probes.py`: four explicitly selected product/diagonal embeddings at degrees thirteen and fourteen. PASS, including every all-f relative identity. Normal and optimized outputs agree. No full degree-13/14 subgroup census or separate uncompressed probe replay is claimed.

## Independence limits

The main and pair-action paths share the affine-normalizer generator, pair-orbit constructor, full-isomorphism primitive, and twin-quotient global-sign primitive. They differ in action coordinates, signed-action compression, enumeration order, and isomorphism bucket filtering. They are not wholly independent programs.

The degree-eleven neighborhood certificate and arithmetic-only Burnside control are separate implementations of their stated, narrower tests. Neither alone verifies global graph signs or the entire finite table.

## Engine provenance

The copied predecessor engine is unchanged. Both source and copy have SHA-256:

    46cdad3ae5c4008dc912e83122383bfd6abd7f3a761688bc934e4705412d1684

Source: `odd_mixed_prime_census_20261020/mixed_symmetry_census.py` in the intended published sibling packet. Copy: `fixed_cube_engine.py` here.

`SHA256SUMS.txt` pins the publication payload. Logs with wall-clock timings, temporary selected-cube files, bytecode caches, and the unapproved double-triple candidate are not part of that payload.
