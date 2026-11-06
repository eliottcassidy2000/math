# Independent frozen audit: PASS

Subject: session203-twoport-digit-obstruction/PROOF.md, verify.py, literal_forest.py, CERTIFICATE.json, and accompanying receipts. SOURCE_SHA256SUMS is the copied source manifest; all entries were verified before and after review. Candidate bytes were preserved.

Proof SHA256: 7a5349e2ad6763380e44745cc4def28e0df29369b17fb08aab7ffe66aa658a1e.
Verifier SHA256: b6abefeca64d752e78b281c583fe07a203efa272568a0761b8c3b192c5b306c1.
Certificate SHA256: 2e08fcf99832d07085ed6f6ddf56bf5f3e9317d0acf914a2dbd3df7a5ad2ea2a.

## Mathematical review

The theorem is valid in its stated one-state, one-update-per-digit architecture, with two retained ports and peak width three. It does not claim true peak width two, general transducer impossibility, or a counterexample to F.

### Universal finite forest semantics

A retained forest has at most two bulk-start components and at most two active-start components, hence at most four components. The endpoint incidence filter in patterns() exhausts the possible active endpoint types, including active singleton components and active vertices internal to paths. Formal unreachable patterns are harmless: actual forests occupy their correct coordinates and need not realize arbitrary nonnegative vectors.

Introduction chooses a predecessor endpoint beating the new vertex and a successor start beaten by it. Bulk endpoints can be predecessors, but bulk starts cannot be successors. Selecting the same component on both sides would form a cycle and is correctly prohibited. Component multiplicity is counted by choosing indices before passing back to an unordered multiset. Deleting the new vertex recovers the old forest and the selected attachments, proving exact transfer rather than merely an upper bound.

Forgetting an endpoint changes its label to bulk. The number of bulk-start components cannot subsequently decrease: no future introduction can enter one, and a merge with an active-start component preserves any bulk start in its predecessor component. Thus discarded states with more than two bulk starts cannot contribute to the final one- or two-component readouts. The truncation is safe at arbitrary order and through arbitrary histories.

The closing readouts are exact. One-component forests give h, root-starting ones give s, and an unordered one- or two-component forest ending at the root has a unique designated component for L. There is no omitted factor of two. Direct marking and introduction of a final root are both covered.

### Necessary physical constraints and bootstrap

When the bulk is nonempty, the four active-singleton-only patterns are impossible. Coordinate 41 is exactly H(bulk), since the active vertices are singleton components and all remaining vertices form one path. Every nonempty tournament has such a path, so x41>=1 independently of seed order, strongness, and root-start counts.

An empty-bulk seed necessarily has exactly two vertices. Either one-digit child has three vertices with two retained ports and therefore nonempty bulk. All descendant identities and the same kappa are inherited. This reduction is valid even if the newly introduced vertex is immediately forgotten.

### Spectral necessity

For the two one-digit child states, the difference of closing readouts is exactly 6. Repeatedly applying either fixed digit doubles that difference because its common affine term cancels. Therefore ell M_i^n(x1−x0)=6·2^n for all n>=0. Cayley–Hamilton then forces the characteristic polynomial at 2 to vanish. This is an all-seed argument; it does not presume a stationary or eigenvector seed. The exact determinant computation excludes precisely the ten advertised matrices.

### Necessary rows and exact contradiction

The quantities r0 x and delta x are respectively constant −kappa and 6 on all prefixes. Therefore every listed word row is a necessary homogeneous constraint. The 62 words of lengths one through five for each of two rows, plus four impossible coordinates, give 128 rows. Their order in the checker agrees with the supplied multipliers.

G=6e41−delta satisfies Gx>=0 for every candidate physical state. For each retained ordered operation pair and closing row, exact rational verification establishes

    H^T y − zG^T − delta^T >= 0,    z>=0.

Multiplication by x>=0 gives 0<=−zGx−6<=−6, a strict contradiction. The multiplier y is correctly unrestricted in sign. Neither integrality of the whole state nor a bound on s or kappa is used. The complete case set is exactly 6·14·14=1,176, including equal digit operations. The discarded ten operation choices are excluded by the independently verified spectral argument.

No LP result, floating-point tolerance, rank estimate, sampled seed range, or search log is required by the proof. The raw rational CERTIFICATE.json data are a necessary part of the mathematical proof and are frozen by hash.

## Code review and replay

Reviewed the full checker and literal_forest.py. All author acceptance gates use explicit exceptions rather than assertions, so optimization does not suppress them. Normal and python3 -O reruns agree byte-for-byte with the author NORMAL.json and OPTIMIZED.json. The checker regenerates 46 patterns, 24 matrices, six closures, determinants, constraints, and all exact rational slacks. It verifies complete case coverage and runs 216 literal transition and 54 closure regression checks.

## Independent controls

1. independent.py independently regenerates all 24 endpoint-update matrices, compares their entries, computes all determinants by Fraction Gaussian elimination rather than Bareiss, explicitly enumerates each constraint word, and verifies every rational slack in all 1,176 cases. All PASS. Maximum observed multiplier denominator is 1633051061078656200; no rounding is used.
2. closure_rows.py independently enumerates local attachments and final readouts for all 276 formal closure columns. All match the author's closure rows.
3. literal_permutations.py uses complete vertex permutations and all cut sets, deduplicating actual unordered path forests, rather than the supplied subset-partition DP. On every labeled seed of orders two through four, it checks all 24 literal updates and six literal closures: 74 seeds, 1,776 transitions, and 444 closing readouts. It includes empty-bulk seeds and checks the necessary physical coordinate conditions on nonempty bulk. All PASS.

These finite controls support the exact code and semantics. Universality follows from the forest bijections, safe truncation, spectral argument, and complete rational certificates, not extrapolation from bounded seed tests.

No mathematical or implementation gap found. Reviewer performed no publication or repository edits.
