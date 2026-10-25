# Independent audit: character graphs and path separations

Status: mathematical claims PASS; verification-hardening caveat below.

The independent audit imports no producer code. It uses explicit exceptions, so its checks survive Python optimization. It exhaustively re-enumerates normalized grids over all 30 supplied fields, recovering 202 grids, checks all supplied representative edge counts, degree sequences and Hamiltonian-path counts, and checks the supplemental witnesses. This audit does not independently recalculate characteristic polynomials or group every enumerated grid by spectrum; those remain the producer check's remit.

## Confirmed claims

- Opposite-pair sums are 2c and magic totals are 3c. For nonzero square c their nonzero-square status is controlled by the characters of 2 and 3 respectively. Characteristic three and zero center require the explicitly stated zero exceptions.
- The characteristic-three parametrization covers exactly c + span(u,v). Nine distinct entries force dimension two. For distinct entries x,y, the member z = -x-y yields x+y = -z. When -1 is nonsquare, a square-entry plane cannot contain zero because that would force both w and -w to be squares. Thus the graph is empty for odd extension degree and complete for even extension degree with the zero-allowed convention.
- If the even-degree plane contains zero, excluding zero sums removes four disjoint opposite edges; otherwise it removes none. Disjointness of the two cosets A and -A follows because an intersection would force c into the direction space in characteristic three.
- F17 is Parker under the exhaustive normalized enumeration. Brute force over all 9! vertex orders gives 1,152 unoriented paths with zero sums and 112 without.
- The printed F59 matrix is a valid nine-distinct-square magic grid of total 3. Vertex 57 is isolated under both edge conventions, so neither graph has a Hamiltonian path.
- The F2187 witness is valid in the pinned polynomial field, has no zero entry and gives an empty graph. The field modulus was independently tested for irreducibility by all possible monic factors through half its degree. The F81 and F729 supplied planes give complete graphs under zero-allowed adjacency.
- On each fixed F81, F729 and F2187 plane, exhaustive 9! permutations yield exactly 432 magic arrangements; the zero-sum triple hypergraph has exactly 12 edges. This agrees with the 9 times 8 times 6 ordered affine-frame count. Graph adjacency does not recover these designated eight magic lines.
- The injective, non-induced hypergraph embedding formulation is an actual equivalence. Neither direction of a magic/HP equivalence survives the stated finite-field examples. No integer-existence inference is warranted.

## Source boundary and pin

Primary source directly inspected on 2026-10-25: O. M. Cain, arXiv:1908.03236v2, https://arxiv.org/html/1908.03236v2 . A privately inspected HTML snapshot had SHA256 3644671bee87716cf1a3300517e1f718207aab7856cef4dcd5d29f876d402a33.

Theorem 2.1 and Corollary 2.1 support the center and opposite-pair identities. Corollary 5.4 covers F17; Corollary 7.1 gives an F59 construction. Section 6 treats scaling normalization. Within the 2019 Cain paper, the broader prime and strict-prime-power Parker lists are expressly labeled Conjectures 7.1 and 7.2. This describes that source paper, not the current mathematical or repository status. The inherited Session 55 finite-field and residue-ring classification results, and its subsequent-source context, are outside this audit and are not downgraded or reopened here. The present packet needs only the explicit small-field claims and correctly confines its own census to q <= 101.

## Verification-hardening caveat

The producer check.py and its field_arithmetic.py depend on Python assert statements. Running those scripts under -O removes their substantive assertions while allowing PASS output. This is not a mathematical counterexample, and a normal producer replay may be valid, but an optimized producer replay is not evidence. The independent audit.py uses explicit RuntimeError checks and verifies frozen SOURCE_PINS.json, all 12 producer file hashes, before calculation. The primary-source receipt records a prior inspection; public replay does not fetch or require the article. Both normal and optimized independent runs must pass and produce identical result JSON.

No producer files or remote repositories were modified. This publication edition excludes the full article. The primary-source receipt records provenance, not a bundled source file.

## Public replay

Run `python audit.py /path/to/character-graphs` and `python -O audit.py /path/to/character-graphs`. The input directory must be the producer packet pinned by SOURCE_PINS.json. The article is neither bundled nor needed. Compare normal.json and optimized.json. SHA256SUMS covers this publication edition only; producer dependencies have their own frozen manifest.
