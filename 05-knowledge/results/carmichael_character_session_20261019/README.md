# Carmichael numbers, character transfer, and the structures that survive a quotient

Research session begun 19 October 2026 and completed after midnight UTC. This is the readable map for the user's request to read the attached Carmichael paper carefully, connect each equation with prior repository work, and develop several independent cross-field themes. No major open conjecture is claimed solved by this session, and historical novelty has not been established.

## Start with the paper, then follow the strongest transfers

The attached paper is Howe–Huffman–Kohlmeyer–Kruer, *The exponent of growth of the Carmichael counting function*, dated October 8, 2026. Its finite-group transfer takes a supply of smooth shifted primes and an ambient character estimate, prunes poorly supported directions, uses a positive Fejer kernel to separate phases, then counts distinct identity-product subsets and applies Korselt's criterion.

[Equation-by-equation reading](EQUATION_BY_EQUATION.md) covers every numbered equation, including the analytic inputs, proof-budget inequalities, Bernoulli transform, size collar and passage to all large scales. It distinguishes the supplied version's formal-verification account from later public reports. The paper's analytic inputs (S),(A) are named assumptions in our new asymptotic arithmetic corollaries; they were not independently proved here.

## 1. Exact-cardinality transfer: a genuinely stronger finite-group tool

[Proof and arithmetic corollaries](../carmichael_affine_transfer_20261019/FIXED_CARDINALITY_TRANSFER.md).

Conditioning the paper's independent selections on exactly k labels changes the relevant invariant. Identity-centered energy must be replaced by phase dispersion

    Delta_chi = n − |sum_i chi(g_i)|.

A coefficient-extraction proof gives an explicit Fourier bound for uniform k-subsets. Affine pruning retains a coset aH, not merely a subgroup. The same ambient cancellation bounds its index and therefore the period of the compatible subset sizes. For compatible k in the mixing range, every target in a^kH has at least binom(n,k)/(2|H|) preimages.

This is an elementary finite-group theorem. Under the paper's analytic inputs, it produces Carmichael constructions with exactly k prime factors, including k divisible by any prescribed fixed integer. It removes the Bernoulli size collar and keeps the exact binomial entropy.

[Prime-power budget refinement](PRIME_POWER_BUDGET.md): retaining logarithms gives the sufficient condition eta>2delta+2sigma instead of eta>2delta+4sigma. The conditional count exponent becomes 1−sigma−eta. This sharpens the transfer, not the analytic supply theorem or the finer Pomerance asymptotic.

## 2. Binary automata and tournament cuts are exact code instances

[Applications and exact codes](../carmichael_affine_transfer_20261019/APPLICATIONS_AND_EXACT_CODES.md).

- The biased binary strip has a Pascal-code weight enumerator. Its complete-period return law and total-variation bound follow exactly from the paper's product-character formula.
- Lucas's binary monomial description gives complementary duality C_F^P⊥=C_(P−F)^P, hence a MacWilliams transform. This is an explicit connection between automata, Hasse/Pascal digits and classical code duality.
- Fixing the number of input ones forces a terminal parity equation. An independent-input mixing theorem cannot simply be conditioned on weight without that affine sidecar.
- Reversing tournament edges changes the endpoint score parities. Characters are graph cuts, and exact-k selection gives Krawtchouk polynomials. For a connected reversal host, bipartiteness is exactly the additional affine parity obstruction.
- Binary-digit arithmetic lifts are a subset-sum problem only after the first matrix lifting obstruction vanishes. Pruning changes the number of available coordinates; the count retains that number explicitly.

These formulas concern the stated finite operations. They do not assert a random prime-gap model or recover Hamiltonian statistics from score parity.

## 3. Square-zero labels separate subset products from closed walks

[Squarefree cycle inventory](../squarefree_cycle_inventory_20261019/SQUAREFREE_CYCLE_BRIDGES.md).

Commuting variables with x_i²=0 retain precisely the objects that must not repeat. Marking prime labels makes regular-representation traces count distinct identity-product subsets. Marking graph vertices instead makes a log determinant count simple directed cycles by their vertex supports.

The map is literal, but the predicates differ: a walk with distinct generator labels may revisit a group vertex. Minimal product-one subsets are characterized by every ordering giving a simple Cayley cycle, not merely one ordering. A naive exterior-algebra replacement cancels the desired terms; a naive logarithm of the return polynomial has negative coefficients because minimal block decompositions are not unique.

This bridge exposed a concrete inherited error. THM-1870's ambient Hamiltonian-length spectral boundary is false: co-spectral six-vertex tournaments with different c6 retain that difference after adjoining a universal source, so c6 remains nonspectral at every ambient order at least six. The positive survivor is exact: in every loopless digraph without digons, c3,c4,c5 are single traces tr(A^k)/k. Historical bodies were preserved with current-truth corrections. No complexity-class claim follows from a co-spectral witness.

## 4. Carmichael power maps have an exact arithmetic/jet clock

[Power retractions and hidden sheets](../carmichael_frobenius_20261019/CARMICHAEL_POWER_RETRACTIONS.md).

For a Carmichael number n, the power map x↦x^n is the identity on Z/n. On Z/n^a it instead converges to a unique multiplicative Teichmuller section; the uniform a−1-step bound is sharp. It is generally not additive.

Retaining polynomial infinitesimals gives a sharper structure. On

    (Z/n^a)[T]/(T^b),   b≥2,

the exact worst-case iteration count is

    a + floor(log_(p_min)(b−1)),

where p_min is the least prime divisor of n. More generally, for a monomial coefficient-precision ring, an explicit maximum over surviving p-power axis monomials determines the clock. Mixed monomial relations do not affect this worst case once the axis precisions are fixed.

The first surviving local Hasse jet has degree p at a unit. The power map can hide both inseparable collapse and separable extension-field sheets while acting identically on base-field points. Its fixed locus is nevertheless etale. These distinctions explain exactly why pointwise equality cannot be transported into a polynomial Jacobian statement.

The additive carry of the composite lift has an explicit CRT-normalized local formula. Integer-valued divisibility is distinguished from coefficientwise polynomial divisibility by a concrete binomial-coefficient hostile.

## 5. Carmichael completion and the correct Collatz carry interface

[Inheritance, completion and carry barriers](../carmichael_collatz_carry_20261019/INHERITANCE_AND_BRIDGES.md).

The earlier one-prime extension theorem is generalized to every squarefree core N. A new prime q completes it exactly when q−1 divides N−1 and Nq is one modulo the old predecessor lcm. This is a complete finite classifier. Examples include 33→561 and the two completions of 91.

The all-unit Fermat-test kernel has an exact index, separating repeated-prime loss from missing field-period divisibility. The marked local factors distinguish causes that a single scalar index cannot.

For the finite Collatz affine group, every one-dimensional group character kills every translation. At the cycle clock, all fixed-shape words have multiplier one, so those characters erase precisely the carry that needs to vanish. A fixed-composition additive-Fourier recurrence restores it.

Integer occupancy gives a second obstruction: when there are far fewer words than residues, their full carry distribution is far from uniform in total variation. A positive deterministic relative asymptotic tending to zero cannot count an integer number of cycles. The repaired heuristic therefore needs an explicit ensemble, clock-gap control and primitive-period bookkeeping. Burnside/Möbius corrections and the sign guard are retained.

## 6. A correct finite-clock spectral comparison, and a new obstruction

[Commutator mixing and smooth Fourier rectangles](../carmichael_collatz_carry_20261019/clock/COMMUTATOR_CLOCK_MIXING.md).

The actual affine permutations E(x)=x/2 and O(x)=(3x+1)/2 have commutator translation by −1/4. This gives a proved lazy symmetric gap at least sin²(pi/D)/16 and a separate forward eigenvalue bound. The symmetric chain includes inverse and holding steps; the directed operator is nonnormal. Neither is a fixed-composition carry theorem.

The converse obstruction is also explicit: sine-weighted Fourier modes on a rectangle of frequencies 2^i3^j give a gap upper bound O((log D)^−2). Therefore there is no uniformly positive lazy symmetric gap across these clocks. The incoming 3-smooth-number theme supplies the actual almost-invariant functions, rather than merely matching the numbers two and three.

The lower and upper bounds do not match in order. They do identify why a changing-level radius near one-half cannot be inserted as a uniform fixed-clock estimate.

## 7. An Apollonius gap closed by a discriminant code

[Generic squareclass proof](../apollonius_squareclass_code_20261019/APOLLONIUS_SQUARECLASS_CODE.md).

Incoming work had the four-discriminant parity relation but inferred the full generic group from a finite census. Direct symbolic elimination now gives each discriminant as a product of three oriented tangency defects divided by an area square. The six defect factors are independent prime divisors in the rational function field. Their four incidence rows have rank three.

That matrix is exactly the binary incidence matrix of K4. It proves the full generic parity-cube Galois group for labelled circles with radii in the base field. Uniform group elements have fixed-circle counts 0,4,8 in proportions 1/8,6/8,1/8. A geometric specialization error term remains a separate problem; finite sample frequencies are not used to prove the rank.

More generally, a connected loopless graph with independent edge squareclasses gives one component parity relation and a degree 2^(v−1) multiquadratic field. This is a reusable graph-to-Galois construction.

## Connection ledger: independent themes introduced or extended

1. Affine cosets and exact-size entropy, replacing unconditioned subgroup mixing.
2. Pascal monomial codes, complementary duality and Krawtchouk forcing laws.
3. Commuting square-zero inventories, distinct-label exclusion and cycle-cover information loss.
4. Arithmetic/jet retraction clocks determined by p-power axis precisions.
5. Finite affine commutators and 3-smooth Fourier almost-invariance.
6. Discriminant incidence codes and generic multiquadratic Galois groups.
7. Quantized occupancy as an obstruction to sparse-reservoir equidistribution.
8. Exact prime-power lcm envelopes, retaining support and height separately.

Every connection retains a specified map, predicate and failure boundary. Numerical coincidences motivated some routes, but do not carry the proofs.

## Next proof targets

- Determine when a prescribed Carmichael core's inverse target lies in a large bank's actual affine product coset. Identity-product abundance alone does not settle this.
- Classify which nonmonomial mixed-characteristic quotients retain the axis-clock formula and where relation-induced cancellation changes it.
- Find a target-specific fixed-composition carry estimate. Full scalar affine characters and sparse global equidistribution are now excluded as substitutes.
- Sharpen the finite-clock symmetric gap between the proved D^(−2) lower scale and logarithmic upper scale, without importing it into the wrong directed or conditioned operator.
- Classify specialized Apollonius squareclass-rank drops and separately prove quantitative arithmetic distribution.

## Publication and evidence

All intentional repository changes are ordinary commits directly on main. Incoming Apollonius and 2–3/Burnside work was preserved and incorporated. The interruption during publication was reconciled before continuing; no unconfirmed operation was blindly duplicated.

[Publication ledger](PUBLICATION_LEDGER.md) records verified checkpoints. Each result packet contains readable proofs, explicit scope labels, replay code and receipts, with independent review where stated. Integer controls and numerical spectral controls are labelled separately. No new Lean build, independent proof of the paper's analytic inputs, CI success without a status, or worldwide novelty is claimed.
