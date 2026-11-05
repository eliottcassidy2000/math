# A stronger quaternary unit-clone obstruction: amplitude at least 96

Session 197. **PROVED by an all-order reduction and a complete capped strong-tournament census:** if a rooted finite tournament has an eventual identity H_m=A*4^m+epsilon, with A>0 and epsilon in {-1,1}, then A>=96. There is no symmetry assumption. This is a separate strengthening of the previously reviewed A>=60 packet, not an edit to that frozen checkpoint. No source at amplitude96 is asserted.

## 1. Exact all-order reduction

The general unit-clock theorem in the published quaternary-obstruction packet proves that an eventual real-base identity is exact for every m>=1, its base and amplitude are positive integers, (b-1)! divides the amplitude, and every SCC outside the nontrivial strong root component is a singleton. Its proof uses the free permutation action on nonempty internal root-return fragments and the resulting gcd1 identity. Thus arbitrary weak padding is removed before any order bound is applied.

For b=4, 6 divides A. If A<96 then A is one of6,12,...,90. The strong root core has

    H_1=4A+epsilon<=361.

Busch's strong-tournament bound gives 5^(n-1)<=H^3. The strict integer comparison

    361^3=47,045,881 < 48,828,125=5^11

therefore excludes n>=12. Every candidate core has order at most11. Order11 is included in the census; it is not incorrectly discarded by the earlier cap217 cutoff.

## 2. Complete generation under the cap

The generator uses the same exact procedure previously reviewed for the cap217 equality bank, with cap361 and maximum order11. Its elementary completeness argument is restated to distinguish the expanded bank from a mere filtered list.

Hamilton-path counts are monotone under induced extension: insert a new vertex at a deterministic valid position in an old Hamilton word. If neither end permits insertion, a change in its incidences along the word supplies an internal position. Deleting the new vertex recovers the old word, so the map is injective.

Every strong tournament of order at least4 has a strong induced subtournament one vertex smaller. Choose a maximum proper strong subset S; a cyclic triangle exists. If at least two exterior vertices remain, no one exterior vertex can be mixed toward S, because that would enlarge the proper strong subset. Partition the exterior into A dominated by S and B dominating S. Strongness forces both classes and some a->b from A to B. Then S union{a,b} is strong. With at least three exterior vertices this contradicts maximality. With exactly two, deleting any vertex of S leaves a nonempty S-part together with a,b; this whole set is strong through the cyclic block orientation S-part->a->b->S-part, even if the S-part is internally weak. Hence an order n-1 strong subset always exists.

Repeated deletion yields a chain from the directed triangle. Induced-count monotonicity keeps every predecessor under the cap. Conversely every added vertex with at least one incoming and one outgoing incidence to a strong predecessor makes a strong graph. The generator therefore enumerates every mixed orientation at each order and loses no strong H<=361 graph.

The isomorphism key groups vertices by outgoing degree and minimizes the complete relabelled adjacency tuple over **all** within-group permutations. Every isomorphism preserves these degree groups, so this is an exact key, not a heuristic graph hash. All new orientations are considered after deduplication, so no marked seed embedding needs to be preserved.

The path counter may reject as soon as one induced-subset total exceeds361. This is safe by the same injection. Retained totals are exact; there is no floating-point comparison. The complete sorted bank has SHA256

    fd502c3c82507fe72991f2db9bcd56f3c399fc4d39b005df557cb627562a7a25.

Its strong class counts at orders3 through11 are

    1,1,6,35,353,3343,3153,1517,331,

for8,740 classes. The11-vertex level is nonempty and essential. The full bank is included, and the verifier regenerates it from the triangle and compares every byte.

## 3. Exhaustive quaternary exclusion with two counters

The actual-clone counter independently recounts H for all8,740 retained graphs. There are1,398 graph/sign choices with H_1=4A+epsilon for a listed A, and12,263 roots in those graphs. Every such root is tested, without a matching filter. Any unit-quaternary identity must satisfy

    H_2=16A+epsilon, H_3=64A+epsilon.

Exactly six roots pass the first equation, and none passes the second. This is already a complete exclusion. No sufficiency claim is made from just these two equations, and the unused optional fourth-multiplicity test is not needed.

A second counter evaluates the exterior-gap root-run polynomial, independently of constructing clone graphs. Each old word receives a factor1 for a valid unmarked gap and a factorz for a valid root insertion; both choices are kept when both apply, including the two exterior gaps. At matching rank3, a unit-quaternary identity requires the run coefficients

    (a1,a2,a3,a4)=(4A+epsilon,6A,4A,A).

All13,456 matching-rank3 roots in the bank are tested and none has that profile. Nonnegative coefficients are saturated at1024, which is strictly above every target coefficient for A<=90. Five12-bit fields preserve additions and shifts without inter-field carries. The positive leading-base/matching theorem rules out other matching ranks; the first counter does not need that theorem or filter at all.

Combining Section1's all-order reduction with the complete finite exclusion proves A>=96 for both unit offsets. The conclusion does not concern other graph-changing operations, arbitrary nonunit offsets, other bases, or the unrestricted Hamilton-path spectrum.

## 4. Sources and reproducibility

The already reviewed general unit-clock reduction is published at:

https://github.com/eliottcassidy2000/math/blob/f633072db124b398e2563b9e3ba7b12e4fc2c7ec/session197/quaternary-obstruction/PROOF.md

Its SHA256 is `f0bf7a6bfb64a13b6b606bb346b7b28b64395b5d18ae56c4022b36adce59289c`. The earlier cap217 bank procedure and its independent review are the named equality packet at commit `a7d7c95bb3532507b63025b0f74d702ced0e692c`, path `session197/minimum-positive-clock/equality/`. The present cap361 bank is newly regenerated and checked; completeness is not inherited merely because its smaller sub-bank was reviewed.

The external order theorem is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006)N3, Theorem1, DOI10.37236/1141, already checked in those named source packets. There is no literature-priority claim.

Run `python3 verify.py` and `python3 -O verify.py`. They compile all three included C++17 programs in a temporary directory, regenerate and hash-check the full bank, execute both different root counters, and verify the strict arithmetic cutoff. Outputs are deterministic. No network access, heuristic isomorphism, or probabilistic count is used.
