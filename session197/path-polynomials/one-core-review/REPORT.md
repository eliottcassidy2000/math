# Independent review: one stationary core between two binary histories

**PASS on the corrected frozen snapshot.** One wording correction was requested and applied before this signoff: the bounded prime corollary now says odd primes and seven-times-odd-primes, so it does not accidentally include2 or14. No mathematical formula, verifier or receipt changed. The source's separately corrected prefix/suffix padding wording was also checked in this final snapshot.

## Exact binding

Source files in session197-odd-cycle-path-polynomials/:

- ONE_CORE_BINARY_INVERSE.md: 701d44aeb13e4f8c26ff44c3d25c3166350371232e7e00e3019e70b5aeab0321
- verify_one_core_binary_inverse.py: 74b92f51e4e4386880e1416af9dce144096d79c292eb4f5f6e9952ab19dc942c
- one_core_binary_normal.json and one_core_binary_optimized.json: 8d3883aa53fb454388a2f35c254ba5fd1e457ef3e4e63f22fdd1f812e9b80400

Named dependencies, read and checked against the receipt:

- RETAINED_ROOT_SEWING.md: c5a7a18611363985b49d4ab455ad42160150ab012157c7d4c01aa6f838332655
- RETAINED_SEWING_SOURCE_SCOPE.md: e301d1085e9d66104ecfac6f514f4f2f2c1efeb9df40c5426aa237f82bc874c6
- verify_retained_root_sewing.py: 62d824eb95905c6c408b38079eb7f56a3c7e8226539f9bc7588f3d1c3855beed
- verify.py: 10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb

## All-order mathematical review

The prefix recurrence follows from the exact O/Z four-state transfers, with O=1 and Z=0 read most significant bit first. Leading Z steps preserve all four coordinates only at the singleton seed. They may be removed before a later stationary core because that continuation uses precisely these coordinates. The proof does not remove internal padding.

The weighted inversion formula counts every odd cycle exactly once: its first chain vertex is an earlier O, its last a later Z, and an even-sized subset of the intervening chain gives the remaining cycle vertices. The gap-one boundary has weight1. All cycles use the retained root, so the odd-cycle collection formula contains at most one such cycle. Reverse-complementing a fixed-width word is an actual converse operation on the root-plus-chain graph, yielding the scalar identity without claiming equality of endpoint states.

The suffix has the opposite bit convention and orientation. Its R coefficient is t, its initial-e coefficient is minus the low bit of t, and the singleton evaluation fixes the d coefficient. Linearity proves the formula on every actual input. Only chronological terminal O steps can be removed at this stage, since the construction ends there. A post-core Z step is not neutral. The corrected transitive-block discussion keeps these cases separate.

Substitution of the exact stationary sewing rows gives the stated two-integer formula. The canonical graph retains every core vertex and all prescribed binary births. The root-deleted graph is the ordered join of the prefix chain, core and suffix chain, so its deletion response is H(core), as used.

The finite inverse bounds are valid for every nontransitive core. In particular C>=(n+1)H and Q<=C/2 give a=C-Q-H>0 for n>=3. The initial core count is at least L(q+1), and A_q=C(q+1)-Q-H is at least a(q+1). Thus the q cutoff and the t cutoff are necessary. Every enumerated pair is physically realizable, proving sufficiency. The harmonic-series work estimate and logarithmic canonical order bound handle q=0 or t=0 separately. They are numerical-target bounds, not bit-polynomial complexity claims.

The no-core inverse correctly strips final O steps and then bounds the even index2v by2M-2. It is essential to do this separately because arbitrary histories have infinitely many neutral representations. Deduplicating actual core signatures is exact for the two-integer formula, without claiming graph uniqueness.

The all-order inverse is for the specified one-core language. The million-bound support is a finite-exact consequence of that inverse, not an unrestricted spectrum or infinite prime-coverage theorem. The two-core H462661 example is a physical escape from this language, not a contradiction to its complete inverse.

## Independent reconstruction

The standalone C++ checker imports none of the submitted code.

1. It inventories every labelled tournament through order5 by literal Hamilton vertex permutations and all labelled path cuts. It recovers all1099 graphs,5405 rooted candidates,65 distinct signatures and50 nontransitive signatures. It counts the designated start-bin by the prefix convention; swapping the two labels identifies this with the source's equivalent suffix convention. The full signature sets match exactly.
2. It constructs2,000,001 prefix states by applying the physical O/Z rules directly to their bit words, rather than using the scalar F recurrence. Their endpoint, cover and deletion coordinates satisfy the proposed formula.
3. Its million-bound census uses the broader product cutoff t<=M/(a(q+1)), not the source's sharper (M-h_q)/A_q loop. It tests21,248,736 core pairs and obtains exactly499,997 odd outputs. The missing values are precisely7,21,462661. The no-core support has249,077 values. The entire sorted support, not merely its cardinality, has SHA256 cfd4ec3cb57c2bae029d42f5bd780728b3c4847c3bd1f713c95eeedb20b87213, exactly matching the submitted certificate.
4. For every nontransitive signature and q,t in0,...,3 it builds the actual adjacency graph and counts all Hamilton words by literal DFS. All800 counts agree with applying the physical four-state transfers.
5. It independently reconstructs both supplied stationary cores, their signatures and the full19-vertex two-core graph. A backwards subset recurrence and a separate literal Hamilton-word DFS both give462661. All19 adjacency rows agree with the source receipt. The middle state is(415,0,4128,15). The elementary factorization127*3643 and nondivisibility by7 are verified.

The source was copied before replay. Its normal and -O runs both reproduce the frozen receipt byte for byte. The independent optimized and debug builds also agree. No source file was changed during this audit.

## Reproduction and artifact hashes

Compile check.cpp with C++17. Run it as `./check support.bin > independent-normal.json`; the optional filename receives the full support bitset used by compare.py. A debug build with -O0 -D_GLIBCXX_ASSERTIONS gives the same JSON. The binary support file is a generated intermediate, not a required published input.

- check.cpp: 1f088d7abad1d9b7bf168e76451525581b0de60b58de4f2f90bbf8f87c076a72
- independent-normal.json and independent-debug.json: f4a40220aa8c1c72eac805947edc9de8bc2450d3cdaf93416d58d99f35f10684
- compare.py: f5f6c8c80fe13b26ed7084028fdfc45b64bf186178631f5df3b3f280688f8d32
- comparison-normal.json and comparison-optimized.json: a06eabcf1cb0909d1ea2e6e6b7a2dde727273f85c24acaa975ff45617c15f380
- source-normal.json and source-optimized.json: 8d3883aa53fb454388a2f35c254ba5fd1e457ef3e4e63f22fdd1f812e9b80400

compare.py reads the adjacent generated support.bin, independent-normal.json and source-normal.json. It compares every signature, every escape adjacency row and the full support digest. Its normal and -O outputs agree, with explicit exception checks.

No issue remains in the corrected snapshot. The inherited stationary scalar/state laws and the inherited binary scalar recurrence remain attributed to their named sources; this review makes no priority claim for them.
