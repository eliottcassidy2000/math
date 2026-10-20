# Independent transfer and budget audit

Reviewed October 19–20, 2026. **PASS after one corrected application-scope issue.**

## Scope and verdict

Reviewed `FIXED_CARDINALITY_TRANSFER.md`, `APPLICATIONS_AND_EXACT_CODES.md`, and `supporting/PRIME_POWER_BUDGET.md` against the supplied paper, especially Sections 2–6. The elementary finite-group arguments, exact application identities, and stated conditional exponent refinements are sound. This is an independent mathematical reading and fresh finite-control implementation, not an independent proof of analytic inputs (S),(A), a Lean rebuild, or a historical novelty assessment. The source paper's claimed author-side formal verification remains distinct from this review.

The only substantive issue found was an ambiguous use of the original number of correction coordinates after affine pruning. It has been corrected in the reviewed application text. No outstanding substantive objection was found.

## Main transfer checks

1. **Fixed-size coefficient extraction (F1)–(F4): valid.** The kth coefficient of the independent Bernoulli generating product is exactly `b_(n,k) E chi(Y)`. Bounding its circle integral by the uniform modulus bound gives the reciprocal binomial-mass factor. The common phase must be optimized; its minimum energy is `n−|sum chi(g_i)|`. Character inversion supplies the stated factor `(|H|−1)/|H|`, and the exact binomial mode inequality supplies `b_(n,k)≥1/(n+1)`. No independence is assumed after conditioning.

2. **Affine pruning and pair dispersion: valid.** A dominant character fiber is a coset of the character kernel. Renormalizing by any retained label therefore permits the next step. At most `D` strict subgroup reductions lose fewer than `DM=N/2` labels. The terminal fiber cap implies at least `nM/2` unequal unordered pairs. Together with the pair identity and `n+|sum z_i|≤2n`, this gives the stated low-order constant `2N/(DR²)`.

3. **Translated Fejer control: valid.** A character extension has order at least the restricted character's order. Thus all powers `1,…,R−1` remain nontrivial in the high-order branch. Translation only multiplies each Fourier sum by a unit phase. The ambient total is at most `Q+(R−1)A0`, the central kernel is at least `R/4`, and the resulting `N/8` arc cap leaves at least `3N/8` retained labels outside. Normalization by the affine representative rotates the entire sum and does not change its modulus. The second constant `3N/(4π²R²)` is correct. The crossover at `D≥27` is correct.

4. **Index and cardinality-period repayment: valid.** Orthogonality on the quotient gives `dn≤Q+(d−1)A0`. Rearrangement requires, and explicitly assumes, `n>A0`. The class of the affine representative has order dividing `d`; an identity-product k-subset requires its kth power to be trivial. This is a necessary affine compatibility condition, not sufficient small-k attainability without mixing.

5. **Conditional arithmetic compiler: valid.** Keeping `U0=y^δ log(5y)` gives `D,log|G|=O(U0)` and `R=O(y^σ)`. The mixing-to-group-size ratio is bounded below at the scale `y^(η−2δ−2σ)/(log y)²`, so the strict condition `η>2δ+2σ` pays both group-size factors and the binomial coefficient-extraction loss. Ambient cancellation still needs `β+2δ+2σ<1`. The quotient period is `O(y^σ)`, hence rounding changes `k=y^η` by a relative `o(1)`. Unique factorization, squarefreeness, compositeness for large k, and Korselt's criterion justify the injection. The exponent `1−σ−η` and the passage through the continuous unbounded scale `X(y)` are valid. Fixed auxiliary divisibility and the stated bounded varying-modulus extension preserve the same reasoning.

6. **Abstract scale law: valid with its compiler interpretation.** Under `N≥y^(α−o(1))`, `Q≤y^(q+o(1))`, and `log|G|≤y^(d+o(1))`, the sufficient mixing threshold is `η>2d+2(q−α)`. The revised statement explicitly requires that a k-subset compiles to the product of its k bounded atoms. Its exponent is then `(α−η)/κ`, approaching `(3α−2q−2d)/κ` from below when positive. The revised statement correctly uses the upper bound `R≤y^(q−α+o(1))`. Neither strict boundary is settled by this proof.

## Applications and corrected issue

- The independent biased-input transform, complete-period return law, Parseval/total-variation bound, Pascal recursion, exact non-dyadic distance two, and fixed-weight Krawtchouk transform are correct. The last-coordinate parity is the complete affine-linear obstruction when `t≥F`; it is not a universal reachability assertion for each k.
- The complementary duality `(C_F^P)^perp=C_(P−F)^P` is correct. Pascal rows are Boolean monomials; two rows pair oddly only when their bitwise union is `P−1`. For the complementary initial row ranges, `h OR k≤h+k≤P−2`, and dimensions finish the proof. The displayed MacWilliams transform has the correct normalization `2^(F−P)`.
- Tournament reversals really change both endpoint degree parities. Cut characters give the stated independent and exact-k transforms. Uniform parity fibers count labelled orientations, not isomorphism classes. For a connected reversal host, an extra affine equation exists exactly in the bipartite case, giving dimensions `v−2` and `v−1` as stated.
- Binary correction digits require the actual syndrome equation after the first-lift obstruction vanishes. **Corrected issue:** pruning to n columns licenses `binom(n,k)/(2|H|)`, never replacement of n by the original a. Nine zero columns and one unit column over F₂, target zero, and k=8 give exactly nine corrections, below the erroneous unpruned bound `45/4`. The revised text separates the full-reservoir and pruned-reservoir hypotheses and includes this counterexample.
- Exponent-two reservoirs with independent probabilities at most one-half have nonnegative Fourier coefficients, so the identity has probability at least uniform and is modal. Equal-probability inverse pairs give the same conclusion in any finite abelian group. The distinction from exact-k selection is essential and is preserved.

## Prime-power envelope and logarithmic budget

The envelope is exactly the lcm of all Y-smooth positive integers at most X. For each prime q its exponent is the largest j with `q^j≤X`, provided `q≤Y`. Interchanging these finite prime-power incidences proves the theta-layer formula exactly, including real cutoffs. At `Y=X` this is the ordinary lcm envelope and its logarithm is ψ(X). This records maximal valuations, not only prime support or a sum of predecessor sizes.

The Bernoulli bookkeeping refinement is also valid: retaining logarithms changes the sufficient threshold from `η>2δ+4σ` to `η>2δ+2σ`, and retaining the constant one-half in the size bound preserves `1−σ−η` instead of weakening it to `1−2σ−η`. The Bernoulli collar still contributes `(1−ν)/(1+ν)`; its removal comes from the separate exact-size theorem. These changes do not independently establish any analytic input or improve the theorem to the finer Pomerance asymptotic.

## Independent replay and evidence labels

Run `python independent_controls.py` in this directory, or run the script by its full path. It imports and executes none of the reviewed verifier. `controls.json` contains exact counts, numerical tolerances, source hashes, and the script hash; normal and `python -O` outputs are retained separately and agree.

The fresh controls include 1,780 exhaustive labelled lists in cyclic and noncyclic groups; 500 affine-pruning cases; 462 exact index checks with a genuinely nonuniform ambient reservoir; 45 effective fixed-k bounds at 80-digit precision; 24 translated Fejer controls of orders 127, 257, and 509; 36 Pascal enumerators and complementary dual/MacWilliams checks; 24 exact rational biased-return/TV checks; 771 connected-graph affine-rank checks; and 4,950 exact envelope identities. Integer and Fraction checks are **FINITE-EXACT for those instances**. Complex and Decimal inequalities are **numerical controls**, not certified general proofs. In particular, an effective case with n=240 and k=120 has actual Fourier magnitude about `1.06e−36` and F1 upper bound about `1.70e−25`; these controls do not merely test bounds larger than one.

## Portable publication scope

The supplied PDF and extracted paper text are not redistributed or required at runtime. Their pinned SHA-256 values are recorded in `SOURCE_METADATA.json`. The independent replay checks the included proof files against the relative-path hashes in that metadata, then records its own script hash. The reported mathematical review used the supplied source; the portable replay only reproduces finite controls.
