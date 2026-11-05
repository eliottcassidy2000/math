# Independent review: sharp separating-vertex bound

**PASS, no correction.** The all-order argument proves that a bridge-free strong tournament with a separating vertex has at least nine Hamilton cycles and at least seven vertices. The order-seven example is unique up to isomorphism within that specified class.

## Frozen source reviewed

All source files are in session197-lowcycle and were left unchanged.

- SHARP_CUT_BOUND.md: 116460a800b344556637707f5c3de448483461e494597fa9c20739b66a7e39ec.
- verify_sharp_cut_bound.py: 88bacf3e3128a728a829ea12c1d30a200f408682a25a5d4b33ec3f07f8eb7ded.
- sharp_cut_bound_certificate.json: ce46ea18a2fc477049c180747490437c9b581eaeb8f79a6666cfa3dc21a3c49b.

The proof uses the independently reviewed start-two indegree lemma and C1 marked-pair lemma in BRIDGE_TWO.md, and the exact separating-vertex product. It is independent of START_THREE and its enumeration.

## Mathematical review

In the SCC order of T-r, the first component's entry-weight sum A is positive. A=1 forces a unique entry arc from outside that component, which is a strong bridge. A=2 has exactly the single start-two entry or two start-one entries stated in the proof. The former supplies an indegree-one vertex different from the entry; the latter implies the component is C1 and supplies an indegree-one vertex outside both entries. Since r is not an in-neighbor of this supplied vertex and every later component is dominated by it, its indegree stays one in the full tournament. Thus A cannot be one or two. The dual argument gives B>=3, so the exact cut product is at least nine.

Singleton endpoint components are excluded by A,B>=3. Each endpoint therefore has at least three vertices, and the separator makes seven. Equality in order forces exactly two directed triangles. Each triangle vertex has one starting and one ending Hamilton path, so the entry and exit sums force all three corresponding incidences. This proves uniqueness without assuming a finite census.

The explicit nine-cycle witness and all five types of alternative paths are correct. Replacing a deleted arc with an alternative directed path preserves strong connectivity, so this proves no strong bridge. The conclusion for at most eight cycles is precisely 2-strongness, not a claim about all nine-cycle tournaments.

## Reproduction and independent controls

A copied source verifier was run normally and with Python optimization. Both outputs match the frozen source certificate byte-for-byte: all 33,867 labelled tournaments through order six and 37,094 separating-vertex cut products pass.

The new independent.py uses adjacency sets, depth-first reachability, and bottom-up anchored subset Hamilton-cycle DP. It checks all 64 incidence patterns of a new separator over two ordered directed triangles. Exactly 49 are strong. Their cycle counts equal the product of the numbers of entries and exits; exactly one pattern is bridge-free, the stated witness. It also constructs an alternative simple directed path for each of that witness's 21 arcs. Normal and optimized outputs agree.

- independent.py SHA256: d25c89fb60513e9112f5d6ea65113d60eb944d1f62d62e19212d063bed2e07cc.
- INDEPENDENT.json and INDEPENDENT_OPTIMIZED.json SHA256: fb97a42d3868cda69f8feb886f72383f0f3e0dbc70527c29556f9502104fadcd.

These finite controls corroborate the reviewed all-order proof. No additional classification or priority claim is made.
