# Review and verification scope

20 October 2026.

## Proof review

A separate reviewer within the session read the complete n=11,12 H-extendible-copy descent and the n=13 exact-defect proof. The review checked:

- restriction of the relative q=3 matching to the original synchronized and independent double-triple predicates;
- correction of a class-preserving automorphism into the C5 normalizer by sign-even internal independent-twin permutations;
- eligibility, not arbitrary C5 presence, in the 5+6 and 5+7 cases;
- the complete small weighted-quotient reduction, all four globally even copy-fusion types and both sign-escape types;
- the distinction between K-based restricted defects and the H-normalizer census defects;
- the source Sylow argument, using the regular-module bound instead of assuming a complete transitive degree-thirteen classification.

The reviewer also checked the regular-module argument and its equality geometry conceptually before its final formulation here. A further separately delegated adversarial check of the small quotient exhaustiveness was requested before publication; no claim that this additional review has completed is made in this freeze record.

The structural claims are marked PROVED on the basis of the displayed proofs and completed session review. Numerical counts have the separate status FINITE-EXACT.

## Completed computations

1. `weighted_copy_audit.py` enumerates every edge setting on weighted quotients (5,5,3), (6,6,1), (5,8), and (6,7): twenty presentations, collapsing to eleven candidate global graph types. It merges independent-twin blocks with equal open neighborhoods, directly computes the quotient edge-reversal character and eligible-class orbits, and returns exactly four fusion losses and two escaping relative-copy orbits. Its six explicit checks also preserve the hostile distinction between one and two triple cycles on K_(5,5) disjoint union I2. Normal and optimized-Python JSON agree byte for byte.
2. `double_triple_n13_probe.py` computes precisely the two specified degree-thirteen fixed cubes, with original permutation generators and all global type/multiplicity/sign sidecars retained in JSON. For H_d the cube has 2^16 points, normalizer order 1440, relative count 1024, and full counts 768 and 762. For H_i it has 2^14 points, normalizer order 2880, relative count 256, and full counts 256 and 250. Normal and optimized-Python JSON agree byte for byte.
3. The degree-eleven/twelve counts are inherited from the complete, independently replayed predecessor packet. They were not newly recomputed for this companion.

## Independence and limits

The weighted audit imports only the Python standard library. It does not import the fixed-cube enumeration, full-isomorphism classifier, affine normalizer, or twin-sign engine. Its finite universe is justified by the written support argument; the program alone is not a proof that this universe exhausts arbitrary degree-thirteen graphs.

The full-cube probe reuses the predecessor normalizer, signed-cube, full-isomorphism and graph-sign engine. Normal versus `-O` execution checks that optimized Python does not disable the controls; it is not an independent algorithm. No separate uncompressed pair-action replay of the degree-thirteen cube totals is claimed. The exact deficit-six proof does not depend on those numerical totals.

## Pinned dependencies

The probe imports the sibling predecessor extension engine, without modifying it:

- `mixed_prime_census_extension_20261020/extension_census.py`, SHA-256 `73a502ec5fc0bc299ac07e99706da9ea41a6882ae692431ee267fb627de36bc0`;
- `mixed_prime_census_extension_20261020/fixed_cube_engine.py`, SHA-256 `46cdad3ae5c4008dc912e83122383bfd6abd7f3a761688bc934e4705412d1684`.

The fallback import also supports the original staging-child layout. All main predecessor payload files stayed unchanged. The old separately staged candidate remains a historical precursor; this companion contains its reviewed successor `DOUBLE_TRIPLE_DESCENT.md`.

`SHA256SUMS.txt` pins this companion payload. The publication manifest is coordination metadata outside its own checksum set. Wall-clock logs and bytecode caches are excluded.
