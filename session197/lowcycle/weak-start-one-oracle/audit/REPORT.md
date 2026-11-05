# Independent weak start-one oracle audit

Verdict: PASS, with the final complexity wording clarification incorporated. The original mathematical statements needed no correction.

## Frozen inputs

- WEAK_START_ONE_ORACLE.md: 55048003e63b68ef21a51d424117145a4afaeb4bdfb43b5946d19524b6006ec0
- verify_weak_start_one.py: ff38a88467e86c58a65922891230486bc42b5e68e6e81e6f5a5991330c50df37
- weak_start_one_certificate.json: aefeb206c02d96fb7c0c9dfb6bb543e17f2c88e78d0a2401867e6b85fce43041
- BACKWARD_PROFILE_INVERSE.md: 8a07f1bb59f7ee3af9c172a1cfec8551c87125747eeb95dddef7335df33df1fb
- verify_backward_profiles.py: 3f11b79a5470736e01ce1596a5d7b5fe9b75be9a78489cc285a9cfc2e22cf4b4
- ROOTED.md: ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e

The six inputs are copied under source/. The reviewed strong inverse is an explicit dependency, not re-proved by this addendum.

## Proof review

The SCC product forces the marked vertex into the first component and every later component to be a singleton. A nontrivial strong later component has at least three paths by Camion and cycle cutting. Thus every weak start-one source is exactly a strong start-one core (or singleton) followed by a nonempty transitive tail. Appending a sink preserves H and the deletion count and doubles the reserved-start two-cover count, including the empty first path. The ending count becomes zero.

The half-open normalization intervals are disjoint. The strict upper core bound U_core<2D is justified by its positive endpoint count, so no dyadic endpoint ambiguity remains. Exact divisibility and k>=1 are essential and are enforced. The endpoint erased by the tail must be existentially recovered; it is not determined by the weak triple. The accepted strong inverse witnesses provide both directions of the membership proof and every possible deletion value. E>0 correctly excludes every weak source. Singleton, two-vertex transient, empty-deletion and even-D boundaries are handled consistently.

The final statement distinguishes O(log U) exact halvings from their bit cost and makes no polynomial-in-log D claim. This is a numerical finite oracle, not a polynomial input-bit complexity claim or a scalar spectrum completeness theorem.

## Independent controls

Copied source runs under normal Python and python -O reproduce the frozen certificate byte-for-byte. The new independent.py instead builds the bounded source image by forward four-coordinate recurrences, retaining the deletion coordinate, then compares against the inverse oracle. It checks 324 forward quadruples through H=101, 323 strong/singleton core tests, and all 82,416 weak integer candidates D=1..101, U=1..16D. These bounds suffice for the independent forward bank because a weak source has the same H as its core.

A separate literal permutation-and-cut counter checks all 1,695 start-one rooted occurrences among every labelled tournament through order five. It directly counts H, ending and starting paths, reserved covers and deletion paths, rather than using the inherited subset DP. The independent tests also use tail exponents 20,80,200, reject nondivisible neighboring responses, and test the singleton boundary. Independent normal and optimized outputs are identical.

The collision (9,0,18) has exactly deletion values {1,3}. Its two source graphs and resulting adjacent first-clone profiles (61,17,89) and (57,15,87) agree with literal counting. This is a real loss of physical information and prevents treating the weak triple as a complete interface.

## Reproduction

From this directory run python independent.py and python -O independent.py. The copied source verifier writes its certificate beside itself; replay only the copy. The source and independent receipts are included. No source file outside this audit directory was modified.
