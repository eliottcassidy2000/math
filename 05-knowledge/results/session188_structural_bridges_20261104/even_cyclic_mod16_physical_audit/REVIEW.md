# Independent normalization and physical review of the even-child modulo-sixteen theorem

**PASS.** Frozen local proof SHA256 `d759e62510d0469dd9e864a46117f2bccc02e8256892162664a021821f07a78c`, checker `b64c97ad9784657af1ed8478a957deb71b129418d5507cb16d526a6bae8ebee0`. This audit covers domain completeness, literal physical/profile controls, complete named coefficients and the local root-family consequence. The independent 512-by-64 formal run-word table reconstruction is a separate audit by the other reviewer; it is not counted twice here.

## Complete normalized upper domain

H is odd, so its inverse modulo sixteen exists. Reduction of `q2/H` and `q3/H` modulo eight restricts their modulo-sixteen lifts to `{6,14}` and `{4,12}`. No information about the odd scalar is lost, since it is retained as the common multiplicative factor.

For any k, dividing q_k by k! counts unordered path covers. Its parity is the number of set partitions into k nonempty parts: every induced tournament has an odd Hamilton count, so each partition contributes one modulo two. Thus `q4/24=S(N,4) modulo2`. When N is even and at least four, the Stirling recurrence reduces this parity to `S(N-1,3)=1`. It follows that q4 is eight modulo sixteen, unchanged by division by odd H. Divisibility by 120 leaves q5 either zero or eight, and `16 | k!` for every k>=6 forces the entire remaining tail to vanish.

This proves the stated eight-element upper domain for every actual even-order child of order at least four, including the N=4 boundary where q5 is zero. The proof does not require every formal vector to be realized. Enlarging the actual domain to all eight vectors is safe for a universal non-annihilation certificate.

The full Hamilton run functional contains one complete-profile factor for every child. Its normalization factor is exactly the odd product of their Hamilton counts. All profile coordinates above five vanish modulo sixteen, so truncating the denominator series above shift five is exact. This retains all visits and all physical two-port terms, including the cubic and quartic parts.

## Physical counts and named data

The audit imports no producer implementation. Complete profiles for all 120 supplied actual tournaments are independently obtained by predecessor Hamilton DP and unordered set partitions, followed by slot factorials. Every normalized vector and every tail agrees. The actual unequal-child systems, both quotient orientations, all block-mask expansions and the complete child profiles are reconstructed independently.

Every retained physical Hamilton count is recomputed with a separate predecessor DP. In particular the TT4-cubed fixed query `(1,3)` is exactly `12625089` versus `13471001`, a nonzero residue-eight contrast. Three full coarse actual banks and the two larger unequal-child controls are checked.

The complete TT4-cubed named bank contains exactly all 6,241 pairs of masks having at most two a names and two c names. The audit recomputes every literal value, every Boolean finite-difference coefficient including zero and shared-index terms, the complete degree histogram and all 32 full-mask reconstructions. These checks retain the full bidegree-(2,2) observer; no low-degree approximation is used.

## Root-family direction and scope

The fixed core query pair is `(1,3)` and its converse-complement `(4,6)`. Their baseline integer differences are exact negatives. The first has negative source-minus-baseline difference. Coordinatewise growth of target child profiles can only increase the target count at a child-uniform query, so the same first query strictly separates every allowed target from the source.

The finite records enumerate all 64 original-module families inside the first TT4 child, with the root also reversed. All selected intervals, physical target labels, complete target profiles, forty strict increases and common negative-query inequalities are independently verified. This is explicitly a one-child finite subfamily, not a census of every original-module family of the entire substitution. The universal root-family conclusion uses the separately proved hereditary monotonicity and dual-query identities.

The frozen packet's new conclusion is local to the cyclic core. Nonzero modulo sixteen alone is not treated as an arbitrary-ambient transfer theorem. Any later special full-cover lift is a separate result with its own hypotheses and certificate. The local conclusion, combined with the earlier five-type theorem, covers arbitrary nonempty cyclic children and profile-preserving target replacements, but does not reconstruct arbitrary prime quotients or identify unrelated physical partitions.

## Reproduction

Normal and optimized audit receipts agree byte-for-byte and pin the final frozen source. The script compiles its independent count engine into a temporary directory and uses only the standard Python library and available C++17 compiler. Literal Hamilton values fit in unsigned 64-bit arithmetic at the checked size bound. The separate full formal-table audit is still required alongside this physical/domain audit for publication.
