# Independent review: primitive unit ternary source witnesses

**PASS, no correction.** The two explicit original sources are strong, rigid and module-indecomposable. They have the exact all-order ternary clone families and finite arithmetic witnesses stated. These are counterexamples to a universal module/twin extrapolation, not a minimality or coverage theorem.

## Frozen source binding

Directory: `session197-prime-atoms/primitive-unit-witnesses/`.

- PROOF.md: `f73f24bf6f3089dbeb61bd437d75700871bf5496b844ca1d4cde539c7a3d8425`
- verify.py: `bbaf1e7c4b20f54978b3515050cb07a6891c2696f8b22cf9775ea6e34132b19f`
- Both source receipts: `1b4f8efbd02a6119d74e0722de3658ef9a1c642fb25d33614f0754e6d742b78b`

The original source files remain unchanged. Copied normal and Python -O replays reproduce the frozen receipt exactly.

## Mathematical and source review

The displayed rows agree with the explicit cross-word construction and its i=k=1 specialization. The root's crossing edges have a two-vertex cover and contain two disjoint edges. This gives exactly rank two, and the positive clone sequence therefore uses only the three positive run coefficients. The reported vectors cancel the binary term and have constants -1 and +1. The separately counted deletion values establish the m=0 statements; no positive-clone interpolation is extended across zero without evidence.

The module enumeration visits every subset of size2 through n-1 and checks the defining uniform outside incidences. The automorphism routine exhausts every degree-preserving permutation; all automorphisms must be among them. It finds only the identity. These are complete finite certificates for the two literal source graphs. The proof explicitly does not extend module indecomposability to clone multiplicities m>=2, whose clone blocks are modules.

The four prime claims have exactly the stated integer-square-root trial bounds. The 7p factors and source factorization 779=19*41 are correct. There is no unbounded primality, least-order or least-amplitude assertion. The rank-two examples do not conflict with the separate rank-one unit-boundary/twin theorem.

## Independent checks

A separate standalone C++ checker starts from the two literal adjacency arrays and builds larger graphs by replacing the last vertex with actual transitive clones. It does not import the source's construction or counting code.

- Hamilton counts for all five displayed source/member graphs are computed by inclusion-exclusion over allowed vertex sets of unrestricted walks of length n-1, with signed 128-bit intermediate arithmetic. They are 779,1561,2339,7019,189539.
- Full root-gap polynomials are computed by literal permutations of the exterior vertices and multiplication of every admissible gap factor. All coefficients beyond degree three vanish, and the full vectors are (259,779,780,260) and (521,1561,1560,520). This is separate from the source's subset gap recurrence.
- Every module subset is tested using the outside vertex's neighborhood intersection with the whole candidate set. No nontrivial module exists in either original source.
- A separate incremental adjacency-preserving backtracking algorithm exhausts all possible automorphisms, with outdegree filtering but without enumerating the source's degree-block Cartesian products. Each source is rigid.
- A complete sieve independently verifies the four primes2339,7019,223,27077 and the exact 7p products.

Optimized and unoptimized library-debug runs agree byte for byte.

- Independent checker: `0ebdddf2608e229d786b2a07cf8dea8dfe6f75ed2fc5b6a7a91f0198d502cdad`
- Independent receipt: `1d79e35bf7d51114c739fba30b77fe6ab4b820a7406b57d7dd6b36d62f8a943c`

The full work artifacts are under `session197-hybrid-constructions/audits/primitive-unit-witnesses/`. No unresolved issue remains.
