# Boundary fibres and exact endpoint transport

**PROVED, independently audited.** This packet extends the actual two-return Hamiltonian-path decoder to arbitrary permitted interior deletion, endpoint-growth kernels, and the full unimodular endpoint matrix. It also retains a literal 625/576 ordered-join model and an explicit marked F5 path table. It does not assert generic deletion safety or canonical field adjacency.

## Read in order

1. [Main proof](BALANCED_DELETION_EXTENSIONS.md): rectangular fibres; all-depth black corrections; general finite product fibres; complete endpoint transport; determinant/permanent one; minimal twin orphan; ordered-join colours; marked F5 coordinates.
2. [Independent audit](AUDIT.md): all-order arguments, hostile controls, exact bounds, parsed printed-witness regression, and publication provenance.
3. [Exact marked paths](marked_f5_paths.json): all 25 source paths with their states, field vectors, and blue predicate.
4. [Cross-packet consistency receipt](cross_packet_receipt.json): equality of all 25 path/state/vector triples with the prime-field packet, not merely equality of cardinalities.

## Adjacent results and the marking boundary

- The [incoming pointed-product decoder](../zero_heavy_ports_20261020g/microcosm.md) supplies the normalized common-coordinate convention. The right coordinate here is the reflection of the raw right-shore subset, with infinity fixed.
- The [two-return source theorem](../strong_blue_games_20261020f/double_return_counts.md) gives the inherited count, rigidity and blue diagonal.
- The [marked F5 determinant/icosahedron bridge](../session17_prime_shapes_20261020/F5_PATH_PAIR_ICOSAHEDRON.md) uses exactly the same five-state labeling and all 25 actual paths. Its determinant graph is an additional **marked relation on sign-orbits**, not a canonical adjacency on the original seven-vertex paths. The explicit marking is required by both packets.
- The [general insertion-response discussion](../zero_heavy_ports_20261020g/profile_growth.md) explains why this family-specific one-component endpoint closure does not replace the distinct-two-component sidecar needed for arbitrary insertion.

The 625/576 join counts are respectively 25/600 and 24/552 blue/black. Removing the diagonal instead would leave 600, not 576. The distinct arithmetic identity 625+125+25+1=776 has not been supplied a four-stratum path interpretation.

## Reproduce

Run `python3 check_fibres.py` and `python3 -O check_fibres.py`; compare the emitted JSON. Run `python3 audit_balanced.py .` and `python3 -O audit_balanced.py .` to check this publication proof and supplied table. Run the audit without a directory argument to check the preserved original proof/table. See the audit for full commands and finite universes.

The source script and producer results are unchanged. `frozen/` preserves the original proof/table and six-file source manifest. [The exact wording patch](SUGGESTED_NOTATION_PATCH.diff) adds only the common right-coordinate convention and the finite-set scope. No theorem ID, global canon, navigation file, or previously published file is changed by this add-only packet.
