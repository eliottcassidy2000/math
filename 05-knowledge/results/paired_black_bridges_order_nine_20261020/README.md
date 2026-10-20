# Paired-black bridge classification through order nine

Status: independently verified finite-exact census and proved local nonbridge criteria. Full all-order nonleaf bridge classification remains open. Independent audit is complete; root review and publication are separate pending steps.

Main report: `BRIDGES_THROUGH_NINE_AND_LOCAL_CYCLE_LAWS.md`.

Results:

- Complete order-nine universe: 268,435,456 tilings, 191,536 ordinary classes, 97,144 merged nodes.
- Exactly seven paired-black quotient bridges, all non-self-converse rigid-word leaf edges.
- The separately audited D(z)=z^4 theorem closes the entire pendant-bridge series; the remaining nonleaf bridge series vanishes through order nine.
- Every other Q edge lies on a loop, digon, triangle, or quadrilateral; all 27,524 strictly quadrilateral witnesses are retained.
- The actual blue graph has 24 bridges, all incident to blue-degree-one nodes. These are not the Q bridges.
- A score-separated, arc-disjoint Hamilton-path criterion proves a genuine black quadrilateral without isomorphism testing.
- An explicit four-vertex insertion gives a black Q quadrilateral at every order n>=8.
- Reid's cited Hamilton-path reversal theorem gives all-order combined connectivity, hence every NS node has a black-only path to some mixed self-converse node. It does not prove that the active Q component is unique.

The literal doubled black graph K has no ordinary bridges. The seven bridges belong to its paired quotient Q, not K.

## Reproduction

Run from this packet directory with a C++17 compiler and Python3. The full map uses 1GiB; allow roughly 3GiB RAM for the independent full-adjacency checks.

    g++ -O3 -std=c++17 census_stream.cpp -o census_stream
    ./census_stream 9 census_n9
    g++ -O3 -std=c++17 sparse_bridge_certificate.cpp -o sparse_bridge_certificate
    ./sparse_bridge_certificate 9 census_n9
    g++ -O3 -std=c++17 short_cycle_audit.cpp -o short_cycle_audit
    ./short_cycle_audit 9 census_n9
    python3 construct_black_quadrilaterals.py

The generated `census_n9.map` is a little-endian int32 array indexed by tiling mask. It is omitted from the compact publication because it is exactly 1GiB; its reproducible hash is retained in `LARGE_ARTIFACT_SHA256SUMS`. The compressed node atlas and sparse certificate are included. The full map and atlas are reproducible rather than inferred from the compact result table.

For audit replay, first decompress `census_n9.json.gz`, `census_n9_sparse_certificate.tsv.gz`, and `audit/expected_nodes.tsv.gz` if the corresponding uncompressed files are not present. Then:

    g++ -O3 -std=c++17 audit/verify_map.cpp -o audit/verify_map
    ./audit/verify_map > audit/map_audit.json
    g++ -O3 -std=c++17 audit/verify_graph.cpp -o audit/verify_graph
    ./audit/verify_graph
    python3 audit/verify_certificates.py
    python3 audit/verify_score_construction.py

The core n9 and score-construction checks above are self-contained after regenerating the map. Supplemental replay scripts retain their original research-stage paths: `audit/verify_controls_and_squares.py` requires the sibling `session14-metagraph-census` atlas plus regenerated n7/n8 maps and node JSON; `audit/verify_pendant_consequence.py` also requires the sibling `session15-deficit-centers` packet and its frozen manifest. The complete `audit/run_audit.sh` expects that original sibling layout. The development `commuting_square_probe.py` runs from the parent staging directory and uses the same session14 atlas; its n4 through n7 outputs are retained. These supplemental dependencies are not needed to reproduce the full n9 map or direct graph verification. Published companion locations are recorded in the publication manifest.

Timing fields in generated reports are observational and need not reproduce byte-for-byte. The node map, mathematical records, bridge masks, and witnesses are deterministic.

## Honest limits

A first nonleaf Q bridge, if one exists, must have order at least ten. There is no proof here that all bridges are leaves at every order, nor that every nonleaf edge always lies on a cycle of length at most four. The new uniform construction supplies examples, not universal edge coverage. The separately audited center classification closes the leaf alphabet and pendant-bridge series. It does not exclude nonleaf bridges.
