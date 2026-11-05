# Independent review: quaternary unit-clone amplitude at least 96

**PASS, no correction required.** This review concerns the frozen cap361 packet, separately from the earlier amplitude60 result. It proves a necessary amplitude bound for eventual physical transitive-root-clone identities with either unit offset. It does not assert existence at amplitude96 or exclude arbitrary graph-changing constructions.

## Frozen source

Directory: `session197-orbit-wildcard/quaternary-96/`.

- PROOF.md: `91e7d5b395c1a60a4bb4679df68afcfe7f4abf938298b0affa47a741c34c2165`.
- verify.py: `0484f6c0c7e2b435e0b329b08db532be3cad72ba91de4d877d853800c41f87d8`.
- generate.cpp: `716a5494f1ed241cd7b64425490869e5833a3eff6ff68c073b486bbf8b98a4c5`.
- clone_counts.cpp: `e4b241fd36e9b026191e95c0b533cad924ee2a3a467607f13df7f858adc4be0f`.
- gap_counts.cpp: `74e4f4a2b190c614a173d4cc1d30fb15bb3ceff846db82a4c74d69fece96a066`.
- strong-bank.jsonl: `fd502c3c82507fe72991f2db9bcd56f3c399fc4d39b005df557cb627562a7a25`.
- Both source receipts and both copied-source reruns: `5253e608960aaed87dc90182b83bff66b61ae87798d22dffb116ccfd9efa1dac`.

The original files were left unchanged. Normal and Python `-O` reruns each regenerated the bank from the triangle, compared its bytes, and executed both supplied counters.

## Mathematical scope and completeness

1. The finite root-run expansion converts an eventual fixed-base identity into the exact positive-clone identity. Its leading amplitude is an integer and the free action permuting the three nonempty internal fragments gives 6 dividing A. This is a free action on disjoint, nonempty labelled vertex sets, not an assumption of graph symmetry.
2. For A divisible by6 and epsilon in {-1,1}, any exterior SCC factor divides both 4A+epsilon and 16A+epsilon. Their gcd divides3, and 4A+epsilon is not divisible by3, so the gcd is1. Thus all nonroot SCCs have H=1 and are transitive singleton SCCs. The root SCC is nontrivial and strong because the sequence is not constant. Arbitrarily long prefix/suffix padding therefore causes no missed orders.
3. A<96 leaves A=6,12,...,90 and H<=361. The exact comparison 361^3=47,045,881<48,828,125=5^11 with Busch's strong bound excludes order12 and higher, while retaining order11. The theorem and its scope were checked against the previously reviewed named Busch source; no inaccurate order6 minimum is used.
4. The strong-deletion proof is valid, including the case of exactly two vertices outside a maximum proper strong set: after deleting a vertex of that set, the remaining three nonempty cyclic blocks are strong even when the first block is internally weak. Thus every strong graph in the cap has a strong one-vertex predecessor until the triangle. Mixed extensions enumerate all possibilities.
5. The deterministic insertion injection proves induced-H monotonicity. Rejecting a partial graph when an induced-subset Hamilton count exceeds361 cannot remove a retained graph. Forward DP remains exact on retained inputs; additions are bounded before overflow. Degree-group canonicalization enumerates every within-group permutation and minimizes the entire adjacency tuple. Hash-table collisions are resolved by full key equality. There is no heuristic isomorphism or fixed embedded seed assumption.
6. The source clone counter tests every eligible root, without matching restrictions. Saturation at target+1 is safe because all transitions are nonnegative. The second packed gap counter has disjoint 12-bit fields and saturates safely at1024; intermediate field sums cannot carry. Its targets lie strictly below the cap. Its matching-rank restriction is unnecessary for the first counter's complete exclusion.

## Independent controls

The independent `check.cpp` imports no source implementation. It uses separate, uncompressed four-coordinate root-gap arrays, a subset-state bipartite matching computation, and literal strong-connectivity and tournament checks. It recounts H on all8,740 graphs and tests every12,263 eligible root. It does not use matching rank to filter those roots. There are19,605 root-gap calls in total. Truncating coefficients above degree3 cannot affect these lower nonnegative coefficients, and saturation at1024 is above every tested target.

Results in both optimized and assertion-enabled unoptimized builds:

- Complete strong inventory by orders3 through11: 1,1,6,35,353,3343,3153,1517,331.
- 1,398 eligible graph/sign choices;12,263 eligible roots.
- 13,456 roots of matching rank3.
- Six roots satisfy a2=6A, equivalently the second-clone equation.
- None also satisfies a3=4A, equivalently the third-clone equation after the first two are fixed.

The six surviving tuples (order,H,root,A,epsilon,a0,a2,a3) are recorded in the receipt. Their final a3 values are78,244,244,78,198,198; the required values are144,288,288,144,264,264. Thus the independent lower-coefficient check already excludes every candidate, regardless of higher matching rank or higher run coefficients.

The cap217 subset of the new bank agrees record for record with the independently reviewed2,456-class equality bank (SHA `3d4ddbe02807c7c6bf962147c582e0c6707eb4d221134a715cdd8a84bf730951`). This is a separate consistency control, not the proof of completeness for the new orders/counts.

Independent checker SHA: `378b2a395de24ce1ff474d2b1bb03f36cf5f9594e035ff649cd04e5708f75165`.
Independent normal/debug receipt SHA: `27aa9bd8092c23413cf3241bb36af4c24401cd85f5db7f09b5449ce08711af9b`.

## Reproduction

The local `source/` directory is an unchanged copy of the six source inputs. Run `python source/verify.py` and `python -O source/verify.py` there. Run `python check_bank.py` to check the named previous-bank dependency and regenerate `bank.csv`. Compile `check.cpp` with C++17, then run the executable with `bank.csv` as its only argument. Optimized and unoptimized outputs are identical. The optional previous-bank comparison can be omitted without affecting the main checker, which only requires the supplied cap361 bank converted to CSV.

Conclusion: the all-order reduction, exact capped census, and independent exclusion all pass. The necessary lower bound is A>=96 for both unit signs, with no symmetry hypothesis and no existence claim at the bound.
