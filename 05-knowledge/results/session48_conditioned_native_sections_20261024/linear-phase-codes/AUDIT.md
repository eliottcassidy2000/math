# Exact verification and scope ledger

The all-size conclusions rest on AFFINE_PHASE_CODES.md and the vendored native carrier proofs. Finite checks do not stand in for all-size proof. The checker uses standard-library integer arithmetic and explicit exceptions, not Python assert or floating eigenvalue approximations.

The deterministic receipt contains107,676 checks before comparison with the frozen receipt. Each frozen normal/optimized run executes107,677 checks; outputs must agree exactly.

A final coordinate-precision repair explicitly distinguishes canonical phase k=owner(8r)-8r from raw owner residue o. The raw equation is A o=b+(8r)A1 modulo3, while the converse-head phase equation has the separate shift2A1. All prior computations already used canonical phases correctly. Added nonzero-translation controls include r=1,row(1,0),canonical right side0, where raw right side2 is required; r=2 distinguishes raw shift1 from converse-head shift2. Arithmetic is in the declared suffix numbering; other owner embeddings require a transported phase chart.

Finite coverage:

- All249 ternary linear codes on0 through4 labelled slots, generated once by reduced row-echelon basis form. Their2,663 affine cosets are compared with direct predicates on actual indexed suffix owners, and independently solved by Gaussian elimination.
- All1,195 full budget graphs on those codes, for budgets0 through the number of slots. Exact graph components are checked against low-weight-span ranks and cosets, not just their counts. Character spectra, zero multiplicities, integral divisibility by3, and trace identities are checked independently of graph search.
- Constructive routes use a basis selected from the original bounded-support kernel vectors, never possibly heavy elimination rows. Every code difference is tested for route membership, accepted arrivals, endpoint correctness, rank event bound, budget bound, and owner-clock-work bound.
- The no-return stage-DAG algorithm is compared with an independent literal visited-record-mask breadth-first search on811 endpoint/budget cases on up to3 slots.
- The five-resource, nine-state native counterexample has ordinary connectivity threshold3 and a specified endpoint's no-return threshold4. Its original budget3 route actually revisits the fifth full record while continuously holding every normalized resource.
- Actual source matrices, complete prefix-pinned menus, distinct four-word normalized reservoirs, source rotations, target partner maps, and converse reflection are reconstructed for1 through4 prefix blocks, of source orders11,19,27,35. The rigid seed's automorphisms are checked over all8! permutations. An unauthenticated same-source terminal reflection is explicitly rejected.
- All three helper phases are retained in the pair/helper code for L=1 through5. Its complete weight enumerator, sharp first mode bridge, budget2 components, and selected/full-budget graph distinction are verified. The27-state L=2 example is additionally enumerated from literal native suffix-owner predicates and all729 ordered record pairs have their exact clock costs checked.
- Every character of the selected graphs for L=1 through5 is evaluated by exact counts of exponents0,1,2. Its closed spectral formula, gap, multiplicities, and total dimension are checked. The actual27-by-27 L=2 Laplacian is independently checked through an annihilating split polynomial and exact traces of powers0 through3. No numerical eigensolver is used.
- The sharp normalized-gap bound is checked with rational arithmetic for L=1 through64; its all-size proof is algebraic. General affine-monomial equation/cost transformations are tested in1,296 coordinate cases. Their mathematical validity is not confused with native owner symmetry.

Every slot retains its own actual full records and one fixed distinct normalized resource. Full retained guards must factor through the fixed incidence data or be separately authenticated and retained. Arbitrary affine equations here are explicitly declared higher-arity predicates; the work does not promote them to pairwise owner equalities, one global observer, native affine reflection holonomy, arbitrary electrical banks, or microscopic metagraph edges.

The selected spectral graph is not the full budget-(L+1) graph once L>=3. Its gap estimates cost growing packet budget and exponential degree. None of its Laplacian statements refers to the original tournament metagraph.

Reproduction closes over the manifest only. portable_check.py first verifies the original package, copies precisely the listed files plus MANIFEST.json into a fresh /tmp directory, and runs normal and optimized checks in Python isolated mode from /tmp. This disallows reliance on a repository checkout, sibling package, working-directory import, or network input. Every shipped file is hashed with LF normalization; unlisted files and symlinks are rejected. Normal, optimized, and portability output logs are retained.
