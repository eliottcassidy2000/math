# Independent review of the unit-clock reduction and quaternary obstruction

**PASS, no correction.** Frozen source: `session197-orbit-wildcard/quaternary-obstruction/PROOF.md`, SHA256 `f0bf7a6bfb64a13b6b606bb346b7b28b64395b5d18ae56c4022b36adce59289c`; `verify.py`, `3174441a728009f0abd4b68f9c875281fa8a173d9b8c2c08470c1e16203abe2b`; `clone_counts.cpp`, `6e92cc8d2e6461f985aab4b310c28c0abe1e70f80dafa88404732f3b562ea0aa`; `gap_counts.cpp`, `74e4f4a2b190c614a173d4cc1d30fb15bb3ceff846db82a4c74d69fece96a066`; receipt `7162385bbf4f323e1f721d79e996a8a68a444aa220cd069a3d331634e275c720`.

## Mathematical scope and proof

The run expansion has a positive highest coefficient a_K. Thus an eventual positive real-base pure clock A b^m+epsilon must have b=K and A=a_K; lower bases cannot cancel the highest exponential rate. Distinct-base Vandermonde independence extends the identity to every positive multiplicity, without asserting anything at zero.

Permuting the j-1 nonempty internal root-return fragments is a free action on skeletons. Every fragment begins in the root's out-neighborhood and ends in its in-neighborhood, so every permutation remains legal, while labelled disjoint nonempty words make the action free. Therefore (j-1)! divides a_j, giving integer A divisible by (b-1)!.

The common divisor g of all positive-m values divides b-1 by subtracting successive values. The factorial divisibility also gives b-1|A, so every H_m is epsilon modulo b-1. Hence g=1, including b=2. This is sufficient to force all SCCs outside the nontrivial root component to be singletons. A singleton root SCC would give a constant sequence and is excluded. The proof therefore correctly removes arbitrarily many neutral prefix/suffix vertices before applying a strong-order bound.

For b=4, only A=6,12,...,54 need exclusion below60, with strong-core H1<=217. The exact arithmetic 5^11>217^3 rules out order12 and higher, and Busch's sharp order11 minimum225 separately excludes order11. The weaker uniform bound alone would not exclude order11. I checked the cited published equality proof and its adjacent independent review at commit a7d7c95bb3532507b63025b0f74d702ced0e692c: their complete bank contains every strong H<=217 graph through order10, before any root/ternary filter. The copied bank hash agrees exactly. Bank completeness and the external Busch theorem are inherited reviewed inputs, not reproved by the present finite counters.

The finite clone equations at multiplicities2 and3 are necessary for any such clock. Their failure on every eligible root suffices for exclusion, without requiring a rank filter or an m4 calculation. The gap-polynomial target (4A+epsilon,6A,4A,A) is also correct by the Stirling expansion of4^m. Saturation at1024 is safe for every compared target coefficient; each 12-bit lane receives a pre-cap sum at most2048, and discarding terms of degree greater than4 cannot affect lower nonnegative coefficients. The matching-rank3 filter is justified by the separately reviewed maximum-run theorem, but the literal clone counter does not depend on that filter.

The conclusion is only the necessary bound A>=60 for fixed-root transitive quaternary unit clocks. It does not assert existence at60, minimal attainability, a bound for other graph-changing operations, or a full scalar-spectrum result.

## Independent reproduction and a third counter

All work is in `session197-lowcycle/quaternary-audit/`; frozen source was copied and left unchanged. Normal and Python -O source runs reproduce the receipt exactly.

`root_revisit.cpp` is a third counting method. It walks through subsets of old vertices, with an explicit root-placeholder state that may be revisited but never consecutively. Old vertices are visited once. The number of root visits records the clone-run coefficient. Unlike the source gap DP it does not multiply exterior/interior gap polynomials, and unlike the source clone DP it never constructs a clone graph. Every internal return needs two distinct old vertices, so the finite visit bound floor((n-1)/2)+1 is valid.

On all410 eligible graph/sign choices and all3426 roots, this checker verifies a1=H1, checks15255 factorial-divisibility instances, and computes H2=a1+2a2 and H3=a1+6a2+6a3. It finds exactly the same two H2 survivors: order8, H1=145, A=36, epsilon=+1, roots3 and4 of the same bank representative. Both have run coefficients

    (a1,a2,a3,a4)=(145,216,78,0).

Thus H2=577 matches, but H3=1909 differs from the required2305. No candidate remains. Optimized and assertion-enabled unoptimized C++ runs agree byte-for-byte.

There is no unresolved issue. The all-order reduction and the independently reproduced complete finite exclusion support the claimed unrestricted A>=60 obstruction.
