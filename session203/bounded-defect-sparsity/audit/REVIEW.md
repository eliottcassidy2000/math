# Independent audit: bounded absolute Z-defect sparsity

**PASS in the stated scope.** The frozen proof establishes the bound for all actual nontrivial positive-start sources with physical Z defect at most a fixed positive integer K. This is separate from the equality-face theorem and does not prove or refute F(8,137).

## Inputs

- `reviewed-inputs/PROOF.md`: SHA256 `51fe374ede02e44fb98f8e36d51efbf04f195da67361e775bb4ac3bd975bef6d`.
- Producer `verify.cpp`: SHA256 `917026ce4b131f5a574be1531a87717e335c8accf1f72286c19561f9234c490c`.
- The three mathematical dependencies are the reviewed cover-equality proof (`680ee2d8e97f700d33b008e2a42f9e61a9287da03d25668e4b7d430e8db72ea7`), corrected sharp-defect proof (`582687a3d966a928efb7355cca6d9f1181d76d09d4065da17a0815d1343f457d`), and published rooted-cover/physical-edge-extension proof (`8e78951f9dd15decff0d23cfb3ebe64f358b3c9d4561b269b527a1e9aa902168`). They are frozen beside the producer. No later or unreviewed strengthening is used.

The cover-equality dependency has now independently passed its own audit. The main proof requires only that dependency's exact bad-cover split and two-terminal lower bound, not its density addendum or complete normal form.

## All-order verification

### 1. Bounded residual Hamilton count and marked exceptions

For root outdegree k>=2, the nonnegative defect decomposition gives D>=s+(k-1)d, hence d<=D<=K and k<=D<=K. The retained marked vertices are precisely the at most K outgoing neighbors of the root.

For k=1, write a for its unique successor and T=U-{v,a}. The rooted directed edge is handled separately: (B,D)=(1,3). For nonempty T, let I be the vertices beating a. Then start_v(U)=H(T)-sum_I start_x(T). The exact L excess and unique-successor identities give

    D = H(T) + 2 sum_I end_x(T) + K_T(I).

The two-cover term has distinct specified end/start vertices and is at least binomial(|I|,2). Thus H(T)<=K and |I|<=K, including |I|=0 or 1. This step is an exact identity on physical graphs, not a formal-state inference. It also handles weak T and nontrivial weak U.

### 2. Uniformly bounded core and ordered runs

The nontrivial SCC counts multiply to H(T)<=K and each is at least 3. The sum of these counts is at most their product when there are at least two; the zero- and one-component cases are immediate. Each such SCC has order at most its Hamilton count, since each vertex supports a starting Hamilton path. Their total vertex count is at most K.

Retaining these vertices plus all exceptional marked singleton SCCs adds at most K further vertices. Unretained singleton SCCs form at most 2K+1 ordered transitive runs. Every retained vertex, the root, and the optional second exceptional vertex sees a run uniformly. In the physical edge extension E(U), the new vertex dominates all old nonroot vertices, so it also sees each run uniformly. Therefore the fixed core has at most 2K+3 vertices and all variable runs remain linearly ordered. There are only finitely many such labeled patterns for fixed K. Allowing zero lengths or patterns that do not preserve the original defect restriction merely enlarges the set used for an upper bound.

### 3. Exact positive Hamilton-path expansion

A skeleton records the order of the fixed core vertices and the strictly increasing list of visited run types in each core gap. No gap can revisit a run after moving to a later run, because every inter-run arc points forward. Repeated visits across distinct gaps are allowed and are counted as separate ordered occurrences.

For run i, assigning its n_i linearly ordered vertices onto q_i nonempty occurrences specifies exactly one path: sort within each occurrence. Every assignment is allowed because core incidence is uniform, and every Hamilton path recovers its unique skeleton and assignment. Thus the coefficient is the onto-map count Surj(n_i,q_i), with distinct skeletons retained at their full positive multiplicity. Both q_i=0 and n_i=0 boundaries are correctly included; no inclusion-exclusion sign enters the skeleton sum.

### 4. Normalizing invisible long runs

If a skeleton is active at a length vector and q_i>=2, its positive integer contribution is at least 2^(n_i-q_i). Since q_i<=c+1, any n_i>floor(log2 X)+c+1 would alone force H>X. Thus all active skeletons at an oversized positive run have q_i=1.

Shortening that run to N=floor(log2 X)+c+1 preserves every factor's zero/nonzero status: both lengths are positive and at least all possible q_i. Inactive skeletons cannot become active, and every active factor for that run stays one. The full count is unchanged. Repeating for all runs is legitimate. This normalization does **not** assert that the shortened source retains D<=K, its original s, or any profile other than the counted B=H(E(U)). Such preservation is unnecessary for the upper bound, and the proof does not use it.

Consequently each fixed pattern contributes at most (floor(log2 X)+c+2)^r distinct counts. With c<=2K+3 and r<=2K+1, finitely many patterns give exactly the displayed exponent and additive constant. The extra 1 safely includes the rooted-edge B value even when K<3. The argument also handles real X>=1.

## Scope of the conclusion

For fixed K the support is O_K(log^(2K+1) X), hence zero-density among odd positive B. Fixed finite shifts remain zero-density. A bounded-start staggered construction covering an odd tail must therefore allow absolute defect to grow without a fixed cap. Neither this theorem nor its controls supply missing profiles, establish a target gap, bound defect when K grows with X, or refute F(8,137). The singleton source is excluded throughout the physical-defect statement, avoiding the corrected singleton boundary from the prior packet.

## Independent controls

- `check.cpp` uses literal directed-path-prefix enumeration, not the producer's subset endpoint recurrence. All 164,030 positive marked tournaments of orders 2 through 6 pass: 141,374 have multiple successors, 22,654 have a unique successor and nonempty T, and two labeled rooted edges hit the explicit boundary.
- For all 164,028 nonedge inputs, the checker independently extracts residual SCCs, marked vertices, fixed cores and 43,222 homogeneous runs. It checks the asserted count/size bounds, every run's transitivity and external homogeneity, the run order, and the physical E(U) count H(E(U))=B. This tests the connection from physical source to abstract pattern, not just the two numerical inequalities.
- `check_skeletons.cpp` builds skeletons by a depth-first walk on unused core labels and the current gap's last run type. This differs from the producer's permutation-and-gap-subset enumeration. It retains multiplicities by occurrence vector and evaluates onto counts by exact integer inclusion-exclusion.
- Across all 636 core/run patterns with c<=3 and r<=2, it compares the positive sum with 9,096 literal graph counts (run lengths 0 through 3, order at most 9). It independently enumerates Hamilton paths and checks 18,078 occurrence-histogram entries, so the test verifies more than the final total.
- Across 93,768 length-vector cases through length 12, exact simultaneous normalization preserves counts; 17,096 cases actually shorten at least one run. Empty runs, no core, no runs, inactive skeletons and invisible runs are included. Declared integer ranges fit signed 64-bit arithmetic.
- Both independent C++ verifiers have identical optimized/debug receipts. The frozen producer optimized/debug replays reproduce its saved receipts exactly.

## Reproduction

From this audit directory:

```sh
g++ -std=c++20 -O3 check.cpp -o /tmp/bounded-defect-audit && /tmp/bounded-defect-audit
g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check.cpp -o /tmp/bounded-defect-audit-debug && /tmp/bounded-defect-audit-debug
g++ -std=c++20 -O3 check_skeletons.cpp -o /tmp/bounded-defect-skeletons && /tmp/bounded-defect-skeletons
g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check_skeletons.cpp -o /tmp/bounded-defect-skeletons-debug && /tmp/bounded-defect-skeletons-debug
g++ -std=c++20 -O3 reviewed-inputs/verify.cpp -o /tmp/bounded-defect-producer && /tmp/bounded-defect-producer
sha256sum -c SHA256SUMS
```

No mathematical correction to the frozen proof was needed. Status-only promotion may proceed with preserved original bytes and a separate packaging review.
