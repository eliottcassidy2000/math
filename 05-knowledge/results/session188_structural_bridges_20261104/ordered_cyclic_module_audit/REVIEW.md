# Independent review: ordered cyclic-module quartic

Result: PASS. No mathematical correction was required. This review pins the author's frozen proof `d528ad7a713b16bb88971466bf1d7ba1f4e2260a1efe7abeafb790ddbcea01cf` and all sources/data recorded in the paired audit receipts. It does not include any later strengthening.

## Analytic checks

I independently checked all eight exterior incidence interfaces, the effective two-run transition matrix, all eliminated endpoint terms, and the exact direction/sign convention in the inherited physical C3 kernel. Clearing the common denominator gives an identity over ten independent integer indeterminates. Grouping the numerator into alpha beta, alpha b, a beta and ab gives exactly the four groups displayed in the proof. No special polynomial relation among the exterior interfaces is needed.

The exterior variables are squarefree but the two run counters are not. A product of exterior excursions enforces disjoint physical support. Every nonconstant denominator term consumes positive exterior support, so its formal expansion terminates on each full-support extraction. The complete marked cover profiles are used. The marked endpoint occupies one prescribed run slot, and permutation of cover slots introduces no extra factor equal to the number of runs. Endpoint balance is used only where required to exchange the start/end role of that physical mark.

All three final series have nonnegative integer coefficients. Their xy coefficient is 2 alpha beta + alpha b + a beta. Splitting a Hamilton path at the incidence singleton identifies the full-support coefficient of alpha beta with the asserted positive Hamilton count, including the empty-exterior case. This proves the strict lower bound and not merely nonnegativity.

I checked restoration after whole quotient converse against the inherited exact disjoint quartic and marked-profile statements. The scalar-equality root-family branch uses hereditary equality of every child q profile, while vertex transitivity ensures preservation of physical marked profiles under all original internal module families. The untouched triangle branch uses the actual exposed-anchor battery. The scalar query is already among the sixteen quartic evaluations. The stated m+20 upper bound is valid. The all-C3 orbit argument applies from each state with the physical witness roles exchanged when necessary; it does not identify responses up to relabeling.

The local statement permits arbitrary exterior E satisfying its two module hypotheses. The quotient-converse and root-family statements do not claim an arbitrary strong outer lift. The transitive interface X->A->Y is excluded. Critical-family classification is an explicitly cited dependency, not inferred from the finite tests.

## Independent executable checks

`independent.py` uses a different exact integer polynomial representation, explicit enumeration of all exterior paths, and a predecessor endpoint recurrence followed by ordered-cover convolution. Its physical Hamilton counter is predecessor-based, rather than the producer's forward recurrence.

- Complete ten-variable numerator equality, retaining all 113 monomials.
- All 549 exterior cases through three exterior names, every exterior tournament and every independent incidence pattern toward the triangle and ordered module.
- 48 additional balanced-mark controls, including a non-vertex-transitive W5 with a balanced hub mark, and a transitive-attachment boundary whose quartic is exactly zero.
- All 19,136 physical integer values for those 598 cases are retained and reproduced in both Python modes.

`audit.py` binds the frozen source and dependencies, reproduces that independent packet, reconstructs physical names and targets in all 21 producer local/quotient controls, checks every source small-case record against the independent result, and checks all 45 critical structure batches. It independently enumerates the complete proper module lists of the four all-C3 root-family examples, reconstructs the exact XOR targets, finite-difference masks, anchor masks, and battery unions, and compares the full integer response vectors in all 104 families. A separate predecessor-DP program recomputes all 3,466 distinct retained augmented graphs. It also freshly replays the producer in the corresponding Python mode and requires byte equality with its complete retained receipt.

All counts fit unsigned 64-bit integers because every queried augmented graph has at most 17 names and each Hamilton count is at most 17!. Explicit exception checks, not Python assertions, enforce every audit condition under both normal Python and `-O`.

## Replay and allowlist

Run `python audit.py > normal.json`, then `python -O audit.py > optimized.json`, and compare the files. No producer code is imported by the independent mathematical/physical routes. An existing standard C++ compiler is required; no installation is performed.

Audit allowlist: REVIEW.md, independent.py, count.cpp, independent.normal.json, independent.optimized.json, audit.py, replay.cpp, normal.json, optimized.json, SHA256SUMS.
