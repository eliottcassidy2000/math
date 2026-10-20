# Reading the Carmichael counting paper through the repository

2026-10-19. This is an equation-by-equation research ledger for the supplied nine-page paper, *The exponent of growth of the Carmichael counting function*, by Will Howe, Grant Huffman, Jensen Kohlmeyer and Liam Kruer, dated October 8, 2026. The PDF was read in full; pages 3–4 were also rendered to check the character bounds, ceilings and constants in the extracted equations.

Source PDF SHA-256: `b1f34a345f63c952a0f1731682817be16ee86a70c8396ada27fa3c9884014f28`.

## Scope and version

The paper asserts C(X)=X^(1−o(1)) using two analytic inputs and an elementary finite-group transfer. Its attached version records an October 8 author-side Lean audit and explicitly says that a clean portable build had not then been independently reproduced. The [companion repository](https://github.com/JENW1N/erdos1057-lean), read October 19, describes later October 10 checks of two revisions and says prize review remains pending. These are different version-specific reports. This session does not rerun Lean or independently establish the analytic inputs.

The new finite-group, graph, ring and algebraic statements below have their own elementary proofs. Arithmetic asymptotic corollaries invoking (S),(A) are labelled conditional on those precise analytic statements. The paper's exponent-one conclusion is not silently promoted to an independently established result of this session.

## Inheritance and the live board

The direct predecessor is [Carmichael successor blocks](../carmichael_successor_blocks_20261014b.md): the lossless interface (N,lambda,e), complete one/two-prime extension classifiers, and the product-one compiler with its exact CRT guard. The hostile L=6,k=2, primes 7 and 13 gives 91, showing why independent congruences modulo L and k cannot be merged without coprimality. The new paper avoids that specific gap by choosing one lcm L covering every predecessor.

Other inherited objects are the integral tournament/Hasse lattice, finite Gilbreath control matrices, the prime-power norm-torsion carrier, the Collatz affine group and carry coordinate, the cycle-spectrum boundary, and the incoming Apollonius discriminant parity relation. The working board is: labelled subsets; cardinality/affine cosets; character phase dispersion; prime-power precision; square-zero cycle inventories; and discriminant codes.

## Section 1: the target and what a logarithm discards

### (1.1): log C(X)/log X → 1

Role: records the exponent of growth, not density or a sharp asymptotic.
Connection: logarithmic normalization is an information-losing representation, as in the repo's distinction between an invariant's order and its integral structure. Both X/log X and X exp(−sqrt(log X)) have exponent one but substantially different relative counts.
New question: which finite-group loss terms can be retained before this normalization? The exact-cardinality result keeps binomial entropy and the size-period sidecar. It still supplies no rate for the analytic-input thresholds.

### (1.2): Pomerance's finer predicted exponential correction

Role: exhibits information absent from (1.1).
Connection: the lower-bound proof spends entropy on pruning, Fourier mixing and size restriction. Write these losses separately before asking whether any can produce a second-order term. A matching functional shape is not enough: the paper fixes delta,sigma,eta before y tends to infinity, while a quantitative second-order analysis would require controlled dependence as those parameters shrink.
Cheapest decisive test: ask for a uniform threshold function in (S),(A). It is not supplied by the elementary transfer. This is a genuine unresolved analytic input, not an algebraic constant optimization.

## Section 2: analytic supply and its exact quantifiers

### (2.1): many primes in (y/4,5y] with all factors of p−1 at most y^delta

Role: gives many distinct labels whose predecessor obligations live on a small prime support.
Connection: this is the multiplicative analogue of selecting objects with a common finite certificate interface. It links to the repo's successor-block lambda and to smooth-number exponent lattices. The labels are primes, not merely their residues: distinct labels are needed for the injective subset-to-integer map.
Transfer: square-zero variables retain these labels even when several primes have the same group residue. The resulting exact trace formulas are in [the squarefree inventory packet](../squarefree_cycle_inventory_20261019/).

### (2.2): log L ≤ y^delta log(5y)

Role: permits exponentially large moduli, uniformly.
Connection: log L is a sum of maximum prime valuations. Its exact universal envelope is a truncated prime-power Chebyshev sum, proved in [the prime-power budget note](PRIME_POWER_BUDGET.md). This connects the paper to the earlier zeta/norm-lattice carrier through an actual formula, rather than a shared small prime.
Transfer: retain the logarithm instead of absorbing it twice into y^sigma; the sufficient mixing threshold improves to eta>2delta+2sigma.

### (2.3): uniform ambient character sum bound

Role: controls all nonprincipal characters of every allowed modulus over the whole prime interval.
Connection: this is an anti-concentration certificate for the ambient family, not a distribution theorem for the selected smooth-prime subset. Positivity of the later Fejer kernel is what makes restriction to a subset safe.
Transfer: the same bound controls arbitrary translated phase arcs and bounds the index of the affine subgroup retained by fixed-size pruning. Both are proved in [the affine transfer packet](../carmichael_affine_transfer_20261019/).
Failure boundary: no such transfer is valid for arbitrary signed weights or for nonlinear observables whose Fourier modes are not covered by the bound.

### (2.4): nonvanishing for Re(s)>23/24

Role: a stated analytic input from the selected formal development.
Connection: the boundary 23/24 is a version-specific dependency, not the 7/8 appearing in another manuscript's title. Neither number is forced by a tournament order or a Platonic count. The finite transfer only needs some character-saving exponent beta<1.
Research consequence: parameterize the transfer by beta and preserve its uniformity requirements. No new nonvanishing claim is made here.

### (2.5): |psi(x,chi)| ≤ C x^(95/96) log(qx)

Role: converts analytic nonvanishing into a prime-power weighted character estimate, then into the ambient prime sum.
Connection: psi counts prime powers with logarithmic weights; its removal of higher powers must retain a controlled error. This is the same prime-versus-prime-power distinction underlying the previous session's exact norm-torsion Dirichlet series, but the latter supplies no cancellation estimate for psi.
New question: can a sparse family of character powers suffice? The Fejer proof only consumes finitely many powers of each tested character, leading to the weighted/power-specific version in the affine packet.

## Section 3: constructing the correct finite group

### (3.1): L=lcm(p−1), G=(Z/LZ)*, separate ambient and selected banks

Role: makes product-one sufficient for every Korselt predecessor condition.
Connection: the group is a certificate carrier; its identity means something because of the lcm choice. Compare the Collatz affine group: abelianization forgets translation, so identity in its abelian character quotient cannot certify zero carry.
Transfer: the carry packet proves the precise blindness of one-dimensional affine-group characters and retains the additive carry transform instead.

### (3.2): log L=sum a_q log q ≤ pi(y^delta)log(5y)

Role: bounds the group size through the maximum predecessor exponents.
Connection: coordinatewise max is a tropical join. It differs from both the radical and the product of all predecessors. In the exact envelope E(X,Y), prime support and power height are independent cutoffs.
Test: compare actual bank lcms with E(5y,y^delta); equality need not occur. Do not identify the integer lcm with the elementary p-primary torsion groups that merely share its order.

### (3.3): D,U,R,A and their different jobs

Role: D budgets pruning depth; U bounds logarithmic group size; R resolves phase arcs; A bounds ambient cancellation.
Connection: these are four coordinates of a proof, not one generic size parameter. The repo's repeated quotient failures suggest keeping each until the final estimate.
Transfer: affine pruning introduces a fifth coordinate, the order h of the retained coset in G/H. Its size is repaid using ambient cancellation, so it cannot silently become an uncontrolled group exponent.

### (3.4): group, depth and phase-resolution inequalities

Role: makes all losses polynomial in y with small exponents.
Connection: asymptotic scales are compositional, but replacing log y by y^sigma twice needlessly enlarges the power budget. The sharpened note retains the original logarithms.
Boundary: improved constants in R or D do not by themselves change the exponent-one result or supply Pomerance's second-order asymptotic.

### (3.5): beta+2delta+2sigma<1

Role: ensures the analytic error is smaller than the selected reservoir.
Connection: this is a feasibility inequality among supply, conductor growth and cancellation. It is analogous to keeping a parameter cone rather than evaluating a single numerical example.
Test: if beta≥1, shrinking delta,sigma cannot satisfy it. Thus the finite transfer cannot manufacture analytic saving out of a non-saving estimate.

### (3.6): (R−1)A≤Q

Role: makes the total nontrivial Fourier contribution to the Fejer average no larger than the principal term.
Connection: a positive polynomial with many Fourier modes needs their combined error budget, not just one small coefficient. The weighted generalization replaces (R−1)A by the explicitly weighted sum of mode errors.

### (3.7): every nontrivial character moves at least M labels

Role: eliminates nearly invisible subgroup directions at controlled label loss.
Connection: this is a finite abelian analogue of retaining rank witnesses. For fixed cardinality, however, the bad condition is near-constancy at any phase, not just phase one.
Transfer: affine pruning removes a nearly constant character fiber and retains a coset. This corrects the naive attempt to condition the original Bernoulli argument on exact size.

## Section 4: the positive kernel and the missing phase

### (4.1): character energy E_chi=sum(1−Re chi(p))

Role: quantitatively separates every nontrivial mode from the identity phase.
Connection: Bernoulli selection naturally uses this energy; fixed-size selection instead uses Delta_chi=n−|sum chi(p)|. These are not interchangeable. If all labels equal one nonidentity group element, E_chi may be large while Delta_chi=0 and every fixed-size product is deterministic.
Transfer: Delta has the pairwise formula n²−|sum z_i|²=sum_(i<j)|z_i−z_j|². Thus pair ratios generate the correct fixed-size subgroup.

### (4.2): the Fejer polynomial

Role: a nonnegative localized detector whose Fourier support is explicitly bounded.
Connection: it links to existing Fejer localization of shifted combs, but its use here is ambient-to-subset anti-concentration. The same polynomial translated by a phase controls every short arc.
Transfer: high-order characters in the exact-size theorem can be centered at the argument of their retained sum. This produces dispersion rather than identity-centered energy.
Failure boundary: replacing a nonnegative kernel by a signed truncation loses the monotone restriction argument.

### (4.3): total Fejer mass ≤2Q

Role: limits how many ambient labels can occupy a short phase arc.
Connection: this is an occupancy bound. In the opposite sparse regime, the repo's integer-occupancy inequality forces Fourier energy to remain large: a tiny set cannot approximate uniform mass on a much larger residue space.
Transfer: the carry packet applies that exact obstruction to fixed-shape Collatz words. A hoped-for exponentially small expected cycle count is not the same as full carry equidistribution.

## Section 5: from independent choices to distinct integers

### (5.1): eta between the mixing and reservoir thresholds

Role: chooses a subset size large enough to mix but small enough to remain sparse.
Connection: it is a scale window, not a special numerical coincidence. Keeping logs lowers its sufficient left boundary from 2delta+4sigma to 2delta+2sigma. The right boundary remains eta<1−sigma.

### (5.2): lower bound K≥e^(−U)(n/t)^k/2

Role: converts probability mass into a count of distinct subsets using a bound on each subset's probability.
Connection: entropy is paid per label. Identifying equal group residues would erase multiplicity and destroy the count.
Transfer: fixed-size selection uses binom(n,k)/(2|H|) directly, while square-zero labels encode the same combinatorics algebraically. Neither substitutes repeated walks for subsets.

### (5.3): factorization of the Bernoulli character transform

Role: independence makes the transform a product of linear factors.
Connection: the group-algebra product product_i(1−rho+rho[g_i]) is exact. In commuting square-zero variables, exp(sum x_i[g_i]) becomes product_i(1+x_i[g_i]), retaining every label once.
Transfer: the same nilpotent trace-log mechanism counts simple directed cycles when the variables mark vertices rather than prime labels. The choice of what is marked is the essential sidecar.

### (5.4): character orthogonality gives identity probability

Role: isolates the identity coefficient in a finite abelian group.
Connection: Apollonius discriminant squareclasses form a literal binary incidence code. Its generic Galois sign group is the even-weight code, and the same orthogonality counts fixed-pair patterns. A new valuation-rank proof establishes that group without extrapolating a census.
Boundary: group uniformity and geometric Frobenius distribution are separate statements; the exact code law does not certify a finite-field sampling error term.

### (5.5): concentration of subset size

Role: the Bernoulli size collar keeps selected products below the claimed X(y).
Connection: conditioning on size changes the Fourier invariant and introduces an affine-coset obstruction. It is not lawful simply to remove the collar while keeping (5.3).
Transfer: coefficient extraction on a common auxiliary phase proves the fixed-size bound and pays an explicit binomial mode probability. This repairs the dependence rather than ignoring it.

### (5.6): X(y)=(5y)^((1+nu)y^eta)

Role: transfers a cardinality cap into a numerical size cap.
Connection: products become sums after logs, but their injectivity still comes from unique prime factorization. The square-zero inventory deliberately retains labels until this final map.
Transfer: an exact-k construction with k≤y^eta uses X(y)=(5y)^(y^eta), removing the collar factor while preserving an increasing continuous scale.

### (5.7): logarithmic Carmichael count

Role: subtracts the group-size cost from subset entropy.
Connection: every proof loss is now visible: retained-bank size, exact or approximate cardinality, and the logarithmic group order. The new budget note and affine theorem improve bookkeeping without adding an analytic theorem.
Next target: a uniform quantitative dependence of the source input thresholds would be needed before turning this into a second-order asymptotic.

## Section 6: parameter choice and all-scale coverage

### (6.1): explicit choice of e,delta,eta,nu,sigma

Role: produces a nonempty feasible parameter cone for every epsilon.
Connection: the displayed denominators 8,16,32 are sufficient slack choices. The general transfer depends on strict inequalities, not on those integers as intrinsic invariants.
Transfer: retain beta symbolically; the exact-cardinality version needs no nu and uses the sharpened eta inequality.

### (6.2): the final exponent exceeds 1−epsilon

Role: verifies that accumulated losses can be made arbitrarily small.
Connection: this is a certificate of feasibility, not cancellation of independent errors. Its optimized exact-size counterpart is 1−sigma−eta.
Boundary: fixed small parameters prove an exponent bound. A finer asymptotic requires their dependence on X to be justified rather than chosen after the limit.

### (6.3): liminf at X(y), then all sufficiently large X

Role: closes the often-missed passage from a construction scale to a full counting theorem.
Connection: X(y) is continuous, strictly increasing and unbounded; every sufficiently large X is covered. This differs from a sparse subsequence result.
Hostile: an estimate along an arbitrary sequence X_j with uncontrolled gaps need not extend with the same exponent. The monotone scale is retained in all new conditional counting corollaries.

## Section 7 and Appendix A: the dependency map is part of the mathematics

The exact formal target is the affirmative counting-exponent statement, and the attached account names its nonvanishing, character-saving, smooth-prime and selected-residue dependencies. These entries are not new finite-group lemmas. This session uses the analytic statements as explicitly named inputs and develops independently proved transfers around them.

The author's preserved audit, later public checking reports, and a fresh independent build are distinct kinds of evidence. No comparison of matching constants can replace matching the actual theorem type and its quantifiers.

## Six additional themes opened by this reading

1. **Affine rather than linear pruning:** exact size exposes a quotient cardinality period; ambient cancellation bounds it.
2. **Square-zero labels and cycle inventories:** trace and determinant identities become lawful subset counters only after repetitions are explicitly suppressed.
3. **Carmichael power maps on thickenings:** identity on residue points hides tangent collapse, higher Hasse jets and separable sheets; the thickened-ring theorem computes exact iteration time.
4. **Discriminant graph codes:** independent divisor factors turn graph incidence rank into a multiquadratic Galois degree and a concrete finite character law.
5. **Sparse occupancy as an obstruction to carry mixing:** the same Fourier language that proves a populous reservoir mixes also proves an exponentially sparse residue image cannot be uniform.
6. **Prime-power envelopes and logarithmic budgets:** the predecessor lcm is a two-cutoff Chebyshev sum, connecting the new paper to the previous integral torsion/zeta carrier while preserving the difference between multiplicity and exponent.

Each theme has an explicit carrier and a failure boundary. None supplies an automatic proof of Collatz, the planar Jacobian conjecture, Wilson-prime completeness or a new theorem about zeta zeros.
