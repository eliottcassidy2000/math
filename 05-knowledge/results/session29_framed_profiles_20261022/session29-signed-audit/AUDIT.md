# Independent signed-controller audit and pair-mask compression

Audit date: 2026-10-22. No repository or Git state was changed.

## Verdict

The signed wedge/fold formula, mean m!, star-phase reduction, universal-source binary recurrence, all-positive-odd blue-count existence theorem, and finite fixed-path maxima through five pairs pass this audit. The even/odd blue maxima193/195 at orders10/11 pass independently. The unsigned Springer relaxation is not falsely promoted to a realizable maximum.

The Gluck-dependent general odd-order spectrum transfer has a correct mathematical application. Primary-source verification was completed independently by root after this auditor's DOI fetch failed; see the source-closure addendum below.

## Independent implementation and results

`audit.cpp` constructs actual signed-vertex adjacency directly from controller orientation and edge labels. It counts half-paths by the pair-mask endpoint DP below. It does not import or use the producer's wedge evaluator. It recovers phase coefficients by baseline and single-star-flip evaluations, then optimizes phases for every local bipartition of every controller through m5.

Fixed-path maxima:1,3,9,39,195.

Max-per-controller histograms:
- m3:9 on8 controllers.
- m4:37 on32;39 on32.
- m5:185 on24;189 on480;191 on240;193 on240;195 on40.

At m5 this used237568 baseline/single-flip DP calls, plus independent evaluation of every phase-optimal witness. Each of the40 maximizing bipartitions has no zero phase coefficient on a nonempty star, so its maximizing phase assignment is unique. Every resulting maximum witness has signed-controller automorphism group of order3, checked over all120 pair permutations. Thus all f195 gauges are covered, not merely40 selected examples.

`factor_audit.py` additionally compares direct adjacency counts to a separately written wedge formula for every controller/sign gauge through m4:4165 exact comparisons. It checks mean m! separately for each controller, all-minus equality with H(Q), positivity/oddness, and strict failure to attain U at m4.

`literal.py` reads producer adjacency witnesses as input data, without importing producer executable code, and independently counts fixed paths, ordinary total Hamiltonian paths, and fold-preserving pair automorphisms. Every supplied centralizer agrees. In particular:
- n6: f=blue=7, centralizer1, ordinary total H=37.
- n10: f=blue=21, centralizer1, ordinary total H=913.
- n10: f=blue=193, centralizer1, ordinary total H=15397.
- n10: f195, centralizer3, blue65, ordinary total H15417.
- n11 after one mate-bit recoloring: f=blue=195, centralizer1, ordinary total H84517.
- Literal repunit witnesses31,63,127 and their odd-order transfers also pass.

No conclusion about a forbidden value of ordinary total H follows from blue7 or blue21. These are distinct observables.

## Pair-mask endpoint theorem: O(m²2^m), not O(m²3^m)

Fix an anti-involution with pairs i+,i−, i=1,...,m, and optionally one fixed center. Let D[S,x] count directed paths using exactly one signed vertex from each pair in S and ending at x, where the pair containing x belongs to S. Initialize D[{i},i+]=D[{i},i−]=1. For every unused pair j and sign t, update

    D[S union {j},jt] += D[S,x] [x -> jt].

The exact fixed-Hamiltonian-path count is

    sum_x D[all pairs,x] [x -> sigma(x)]

in even order, and

    sum_x D[all pairs,x] [x -> center]

in odd order. In the normalized fold gauge these both reduce to summing positive endpoints.

Proof: each directed half-path chooses one vertex from each pair, and its reversed sigma-image gives the second half uniquely. The fold criterion is exactly the remaining middle adjacency (or middle two adjacencies in odd order, the second following from anti-invariance). Every partial half-path has a unique previous endpoint and previous pair set, so the recurrence is disjoint and exhaustive. Previously chosen signs affect neither future pair availability nor the next adjacency once the current endpoint is known. Summing them out is therefore exact, not a heuristic. There are O(m2^m) states and O(m) transitions per state, giving O(m²2^m) arithmetic operations and O(m2^m) integer storage. Bit complexity additionally depends on count bit lengths.

This is an addendum improving state compression. It does not silently modify or claim an error in the earlier valid ternary implementation.

## Proof audit details

For a half-ordering, each adjacency pins the sign of exactly its controller. A repeated pin is either an internal vertex controlling both neighbors or a final vertex controlling its predecessor and the fold. There are m pin occurrences on m sign variables. Each consistent doubled pin releases one free variable, giving its factor2. Opposite required pins give0. This proves the wedge factors1-b_vu b_vw and fold1+b_vu. No edge variable appears twice in a monomial, and all nonconstant terms have zero sign-average, proving mean m!.

Changing all labels on one outgoing star changes no internal wedge. There is at most one final fold per permutation, so the phase polynomial is affine A+sum epsilon_v B_v with no mixed phase terms. Maximization is exactly A+sum |B_v|.

To attain the unsigned relaxation for every ordering, putting an edge at the final fold forces its label+. If one controller has two outneighbors, an ordering placing it between them then has wedge0, contradicting relaxation saturation. For m>=4, a tournament necessarily has a vertex with outdegree>=2. This is a genuine global sign-frustration obstruction.

For a new universal source-controller u with uniform label epsilon, interior placements vanish. Initial placement preserves old f; terminal placement contributes(1+epsilon)A. For the wedge-only sum, both endpoint placements contribute A and all internal placements vanish. Hence f'=f+(1+epsilon)A and A'=2A. Starting at f=A=1 gives all odd1,...,2^m-1, with unique binary row-sign choices.

The transitive controller has trivial automorphism group. An automorphism commuting with sigma preserves its pairs and the positive fold representative; thus no pair sign can flip, and no pair permutation remains. Its centralizer is1. Therefore the constructed blue count is exactly f.

For h>=3, choose the minimal m=ceil(log2(h+1)); the highest binary digit forces the last added source label+. In even order, u− dominates all old signed vertices, all old signed vertices dominate u+, and u+ dominates u−. This makes the tournament strongly connected, irrespective of old connectivity. In odd order, the route closes through u+ -> center -> u−, with the same domination relations to old paired vertices. Thus every positive odd h>=3 has strong realizations at both minimal orders. For h1, C3 supplies a strong odd-order realization by choosing its otherwise unused mate arc to close the cycle. No strong even-order h1 claim is made here.

## Finite blue maxima

Since every blue count is f/c for a positive odd centralizer order c, the verified fixed maximum195 bounds blue. All f195 even witnesses have c3 and yield65. Any smaller odd f is at most193, and the explicit f193,c1 witness realizes193. Hence max blue at n10 is193. At n11, mate-bit recoloring breaks the cyclic centralizer while preserving f195, proving maxblue195. At smaller certified orders the maxima have a c1 witness, yielding even-order1,3,9,39 and odd-order1,3,9,39.

These finite enumerations are not a formula for all m, and do not identify the user's unresolved M_k unless its definition is separately recovered.

## Reproduction

Compile audit.cpp with a standard C++17 compiler and run in this directory. Run python literal.py and python factor_audit.py. literal.py defaults to the sibling session29-coherent-halfpaths/blue_checks.json and accepts --witnesses PATH for another packaging layout. Compiled binaries and __pycache__ are not deliverable evidence. results.txt, literal_results.json, factor_results.json, and top_m5.txt are exact audit outputs. The auditor uses the opposite q-bit convention from the producer and direct required-sign b labels, so raw controller/sign codes must not be compared without conversion.

## Source closure (2026-10-22 07:34 UTC)

Root independently read the primary Cambridge PDF, not merely an abstract or search snippet. Definition(*) on printed page60 is the existence of a subset with trivial setwise stabilizer. Corollary1 on printed page66 (PDF page7, extracted lines341–344) applies to any odd-order permutation group. This closes the external-source qualification above and verifies the general odd-order mate-bit regular-subset application.

Primary PDF: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/9667B85D9DDD0FE911C8382287B86167/S0008414X00022550a.pdf/trivial-set-stabilizers-in-finite-permutation-groups.pdf

Source provenance: direct primary-source reading reported by root to this auditor at07:33 UTC. This audit does not misattribute that independent reading to its own failed DOI retrieval.

## Final frozen-source closure (2026-10-22 07:44 UTC)

The producer's final eleven-artifact SHA256SUMS verifies without discrepancy. The check_model.py import change is confined to adding os, selecting a sibling-relative DP dependency with optional HALFPATH_DP_PATH override, and supplying that selected path to the existing import loader. Compared with the earlier code inspected during this audit, the counting/formula logic is unchanged.

Section13's strong-realization proof agrees with the independently proved wrap argument above and correctly limits the h1 strong construction to odd C3. The claim about the faster pair-mask algorithm is explicitly credited to this separate audit rather than to the unchanged ternary verifier.

The new check_strong.py was inspected and executed in a separate copied directory, using the documented HALFPATH_DP_PATH override. It checks reachability and fixed-half-path counts for127 constructed graphs. The regenerated strong_checks.json is byte-identical to the frozen producer output. check_model.py also reran there, regenerating independent_checks.json byte-identically. No producer file or repository file was mutated by this recheck.

REVIEWED_SOURCE_SHA256SUMS records the exact final producer sources/results and its DP dependency reviewed here. These are provenance hashes, not a claim that every producer file was independently rederived. The independent adjacency enumeration and factor controls described above remain the mathematical cross-check.
