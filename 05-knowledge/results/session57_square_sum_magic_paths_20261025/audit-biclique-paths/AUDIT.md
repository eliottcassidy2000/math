# Independent audit: decorated biclique bridge and sharp path certificates

Verdict: PASS, with the execution caveat below. No mathematical defect found in the reviewed claims. This audit does not assert an integer nine-distinct-square solution, arithmetic independence of nine square tests, or literature priority.

## Evidence

- `audit.py` independently builds all paths by vertex-level DFS, deduplicates cycle edge sets, and checks every edge and pair incidence. It does not import either producer verifier.
- It enumerates all 31,684,800 unordered triples directly using edge-mask union, independently of the producer's incidence-bitset intersection method. It obtains 675,264 edge covers.
- Counts: 576 unoriented paths, 72 cycles, 20,624 label-covering pairs, and 204 optimal nine-edge pairs. Union-size distribution: 9:204, 10:1852, 11:5600, 12:7720, 13:3760, 14:1488. There are 36 cycle-cover pairs.
- The displayed sharp paths each use all eight shared vertices, have five common edges, and their nine-edge union has all nine distinct formal labels.
- Independently checked reconstruction edge indices, all eight symbolic magic lines, label multiplicities, and the exact two-element fixed-formal-label automorphism group.
- Independent polynomial arithmetic verifies the finite-field witnesses, checks supplied extension moduli for irreducibility by exhaustive monic-factor division, and exhaustively checks every gauge. Good-gauge counts are 20, 50, 72, 2178 at q=29,59,81,2187.
- `normal.json` and `optimized.json` are identical outputs from normal and optimized Python. All audit checks use explicit exceptions. An altered expected triple count is rejected under both modes; see `audit_negative_controls.json`.
- Producer bridge verifier was executed only on an audit-owned copy. Normal and optimized modes both pass and both reject a tampered certificate; see `producer_bridge_runs.json`. All pinned local source hashes match.

## Proof review and boundaries

The normal-form derivation is valid in characteristic 3: counting the four center lines gives T=3C without dividing by 3. The S+S=(x+y)+(S-S) identity gives exactly nine forbidden gauges whenever the nine offsets are distinct in an odd field. Integer square-scale 4 and gauge 1 produce positive disjoint vertex classes with residues 3 and 1 modulo 4. Reconstruction is correct. Neither an arbitrary bare biclique nor arbitrary un-decorated path union implies an MSOS.

The two-path conclusion requires shared numerical decoration and nine distinct numerical labels, plus integer squarehood and positivity in the integer setting. Nine tests are minimal for direct label coverage, not necessarily logically independent arithmetic tests. Vertex values need not be squares; inducedness is not required; zero-sum edges must be permitted when a field witness uses zero, or otherwise excluded by a stronger hypothesis.

## Execution caveat

The path-count producer `verify.py` uses Python `assert` and hard-coded absolute output paths. It was not run in place, and its checks must not be advertised as validating under `python -O`. Independent audit verification does remain active under `-O`. This is a reproducibility caveat, not a failure of the enumerated counts. All audit writes are confined to this audit directory; no producer or remote modifications were made.
