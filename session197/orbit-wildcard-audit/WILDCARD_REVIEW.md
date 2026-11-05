# Independent review: parity-context divisors

Review result: **PASS** on 2026-11-05, with no correction requested. This review owns the all-order parity-context theorem and its arithmetic/physical scope. A separate reviewer is independently rebuilding the 14-vertex graph certificate.

## Frozen source

Path: `session197-orbit-wildcard/`.

- PROOF.md: `69497d610f4ffdb3530ed260cec3618af3d973aadaadbbddd23aaff343424f4a`
- verify.py: `9f354c23dfe8616f73f58cc15f56b6054a5f197bdb51fc944c6550b04d6ceba9`
- receipt.json: `50f7e00c1375e818bdc01f419d0c13b2340d59caee197cd2e2aaeae3e5c31aef`
- witness.json: `b4358dec35824a42a5756bd73a2cd180a3b7fd242870628fb1e54ae5430d36d7`

Both normal and optimized verifier runs completed, their emitted JSON agreed byte for byte, and all four hashes remained unchanged. The verifier rewrites its two certificates deterministically; their bytes were verified unchanged after replay.

## Proof checks

1. The labelled-bin OCF expansion is exact. A fixed disjoint odd-cycle collection contracts every selected cycle to one bin-assignment object and leaves uncovered vertices as singleton objects. Its exponent is positive for nonempty input and has the same parity as input order. This proves a polynomial identity, not merely an evaluation at one bin count.
2. The whole-module context is a fixed integer linear functional on falling-factorial coefficients. Maximal module runs are genuinely nonempty; the factor k! orders an unordered path cover among the k occurrences of the module marker. The coefficient a_1 equals H(Q). Applying this functional to the positive transitive-profile expansion is therefore valid.
3. The two infinite gcds agree in both directions: the transitive-parity gcd divides every arbitrary-input count by the expansion, and each transitive graph is itself an allowed input. The input class ranges over all positive orders of one parity, not one fixed order.
4. The stride-two recurrence is monic and integral, with characteristic polynomial product_(b=1)^D (x-b^2). Hence every later value is an integer combination of the first D values. The claimed finite gcd test follows without a positivity or realizability converse.
5. The even sequence correctly starts at clone size two, never at a formal size-zero term. In the prime-modulus criterion, bases divisible by the modulus disappear because input lengths are positive. For the remaining distinct nonzero square classes, the Vandermonde criterion has weights beta_b*b for odd lengths and beta_b*b^2 for positive even lengths. This is valid also when square classes merge.
6. The numerical odd-context atom exclusions distinguish divisibility from coverage. The C3/H21 boundary checks orders one and three separately and uses the positive leading h_n term thereafter. The C5 and C7 minimum-size bounds rule out their small exceptional multiples. The even-input examples appropriately defeat an unqualified parity-free statement.
7. The context extremal gap follows from the empty cycle collection and one triangle. Equality forces exactly one nonempty odd-cycle collection, hence exactly one odd cycle, necessarily a triangle. Extra cycles or collections have strictly positive contributions in every context.
8. The carousel/Paley switch is valid in every whole-module context because the complete cover vectors differ only at f_1, by 14. This is not a claim about nonmodule automorphism orbits.

## Finite evidence

The supplied verifier exhausts labelled input profiles through order five, complete legal run words for C3/C5/C7, actual transitive clones, arbitrary small substituted inputs, the seven-orbit cover vectors, and physical switch examples. Its finite checks support the displayed coefficients; the all-order statement follows from the proof above rather than these tests.

The 14-vertex section's conceptual strongness, module argument, and group-order deduction also read correctly. In particular, equal total degrees force any module of a regular seven-tournament to have odd order; the only possible proper nontrivial sizes three and five are excluded. A cross-pair in the 14-vertex graph cannot be a module because its vertex degrees differ by three, whereas a two-vertex module has degree difference one. Independent numerical certification of that graph is assigned separately and is not represented as new evidence from this review.
