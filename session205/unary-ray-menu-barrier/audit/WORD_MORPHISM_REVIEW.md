# Independent audit of the word-valued morphism extension

Status: PASS for arbitrary finite-word-valued morphisms with fixed end contexts.

This supplements rather than retrospectively expands the earlier letter-or-erasure review. Reviewed proof SHA256: `7d81d1f713f58f56f8e959de8ecc0b6f28714df11a674ebc8ecec19e9a203859`. Verifier SHA256: `73bd847a6d25c49c2859a983f07de5456708a4ee0f6b514d9b778e1a4f2d3389`. Receipt SHA256: `4f66ea01625fe7bdf2e11d40e6ad5deb6771a998f1ae869e1d009ff5142c6999`.

Both normal and optimized verifier runs exactly reproduced the receipt. I also separately checked every off-diagonal nonzero of each regenerated matrix against the proposed common permutation: all four matrices are lower triangular in that ordering, and all four diagonals coincide. The respective numbers of off-diagonal nonzero entries are 35,46,46,45.

For any word v of length k, chronological multiplication uses the reverse letter order in the matrix expression, as required by column states. Common triangularity and shared diagonals make the product-spectrum conclusion valid in either multiplication order: its diagonal entries are the k-th powers of the common entries. Thus its eigenvalues are precisely 0,1,2^k,3^k with the advertised multiplicities.

For k>=2, an arbitrary repeated-word sequence has, after the finite zero-eigenvalue transient, an exponential-polynomial representation using bases 1,2^k,3^k. None equals 2. Subtracting the desired target therefore leaves the nonzero coefficient -2 on base 2. Positive distinct bases ensure a nonzero expression has an eventual fixed sign by its largest nonzero base and polynomial degree, so there are only finitely many successful exponents. This argument does not extrapolate past a nilpotent transient and needs no repeated-word rank-stabilization assertion.

The k=1 case retains the separately audited rank-stabilization and forbidden-seven proof, including its valid n=1 endpoint. The k=0 case is constant. Finite unions and finite compositions preserve the conclusion. In a composition, the new prefix and suffix remain fixed words and the composed map remains a word-valued morphism.

No mathematical changes are required. This rules out the full eligible-language selector-free menu from the Session204 criterion, including word substitutions, but not a restricted witness language, stateful transducer, variable context, free insertion position, or unbounded normalization procedure. It proves no obstruction to the broader eventual odd spectrum conjecture.
