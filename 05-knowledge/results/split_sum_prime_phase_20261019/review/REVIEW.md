# Independent review: split-sum branching, binary balance, and kernel recovery

Date: 19 October 2026. Verdict: **PASS for the stated mathematical claims and finite certificates.** No blocking correction was found. Neither packet proves greedy prime completeness, and neither is being endorsed as a globally new result.

## 1. Scope and disposition

Reviewed independently:

- `../PRIME_BRANCHING_PHASE_AND_FERMAT_BALANCE.md`
- `../verify_prime_branching.py`, its receipt, README, and manifest
- `../KERNEL_PHASE_RECOVERY.md`
- `../verify_kernel_phase.py` and both supplied outputs
- The inherited repository theorem THM-4644, specifically its prime-order circulant-tournament statement and proof.

Reviewed source hashes and presentation overlays are recorded in `SOURCE_PROVENANCE.json` in the parent directory. The supplied prime verifier writes its receipt beside itself.

The verdict covers mathematical correctness, the stated source hypotheses, finite examples, and reproduction. It is not a literature-wide novelty determination or a proof of the cited deep smooth-residue theorem from scratch.

## 2. Source and attribution audit

The following primary sources were opened directly on 19 October 2026:

1. [Booker–Pomerance, *Squarefree smooth numbers and Euclidean prime generators*, v4](https://arxiv.org/html/1607.01557v4): Theorem 1 applies to primes other than 5 and 7 and gives positive squarefree p-smooth representatives. For nonzero residues the factor p cannot occur, so all prime factors are strictly below p. Section 4 identifies residue 4 as the sole exception at 5 and 7. Corollary 8 already uses the eligibility test and nonresidue phase flip for freely selected branches; Section 6 removes the exceptional cases when signs are allowed and explicitly asks about greedy completeness. The packet credits these ingredients correctly.
2. [Booker, *A variant of the Euclid–Mullin sequence containing every prime*](https://arxiv.org/html/1605.08929): the introduction identifies Chua's sequence with choosing the least prime factor across all split sums. Its existence theorem allows choosing a branch and does not establish completeness of that greedy sequence. The packet preserves this distinction.
3. [Booker–Simon, *A generalisation of the Euclid–Mullin sequences*, v1](https://arxiv.org/html/2601.21901v1): the introduction separates least-factor and largest-factor N+1 sequences and states that completeness of the former remains open. Its results about the latter do not resolve the former.
4. [Akin, *Rock, Paper, Scissors, Etc*, v4, Theorems 3.13–3.14](https://arxiv.org/pdf/1806.11241): the cited multiplier-stabilizer classification and odd-order subgroup construction match THM-4644. The packet accurately treats full automorphism rigidity as the inherited Akin/Burnside/Alspach assembly, not a new theorem proved here. THM-4644's prime-degree Burnside argument supplies the additional full-automorphism step.

In particular, the obsolete 2016 status of squarefree residue coverage has not been imported into the packet as an unresolved conjecture. No unsupported source hypothesis was found.

## 3. Split-sum reduction

### Exact candidate count

For q prime, q not dividing squarefree N, multiplication by the invertible divisor d gives

    q | d + N/d  iff  d² = -N (mod q).

The two roots, when present, form exactly one signed class. Consequently the count is the signed divisor multiplicity at that class, rather than twice that multiplicity. The stated formula has the correct normalization.

For N > 1, complementing an actual divisor has no fixed point because a positive squarefree integer greater than 1 is not a square. The ordered witness count is therefore even, and division by two counts unordered splits. The empty source N = 1 is correctly excluded from this parity assertion.

The full-support premise cannot be dropped: the displayed N = 2, q = 11 example has a positive quadratic phase and no divisor in the required root class. Both the phase formula and this hostile pass independent literal enumeration.

### Smallest missing prime

Once q is the smallest missing odd prime, the cited support theorem gives every signed divisor class. A candidate below q cannot divide any split sum because it is already a factor of N. Thus eligibility of q is equivalent to its greedy selection. At a blocked state, multiplying N by the selected r changes the phase by the Legendre symbol of r. This proves the next-step trigger and the eventual-residue trapping statement exactly as written.

The equivalence between completeness and the absence of a permanently trapped quadratic-residue tail is sound. Distinct outputs, the eventual inclusion of every prime below a least omitted prime, and infinitude of prime nonresidues are sufficient. The elementary CRT proof of the latter is valid: prescribe a nonresidue modulo q and 1 modulo the finite product of alleged nonresidue primes; a prime factor of the resulting positive integer supplies a contradiction.

The Paley orientation is correct precisely for q = 3 modulo 4. At q = 1 modulo 4 it is an undirected graph. There is no hidden Hamiltonian-path premise.

The original-Mullin countermodel is correctly limited to residue dynamics. Alternating the nonresidues 3 and 5 modulo 7 gives cumulative products 3 and 1, avoiding -1. It demonstrates the failed implication for the N+1 rule without purporting to construct a Mullin trajectory.

## 4. Binary balance and the Fermat threshold

The cyclic binary-balance theorem is correct, including h = 1, empty lists, repetitions, and zero labels.

- Uniform integer counts require h to divide 2^k, hence h is a power of two.
- At h = 2^s, a character of index j with 2-adic valuation a is killed by a binary factor exactly when the corresponding label has valuation s-a-1.
- A product of complex numbers can vanish only through a vanishing factor. Thus the required layer coverage is both necessary and sufficient; no cancellation argument has been substituted for it.
- Changing the primitive root multiplies logarithmic coordinates by an odd unit, preserving all relevant layers.

For a Fermat prime q > 3, -1 is a square, and the product phase is the parity of the number of odd logarithms. One representative of every layer gives exactly one odd label and the blocked phase. A balanced generating source needs at least two odd labels, so at least s+1 prime factors. At equality there are precisely two layer-0 factors and one from each higher layer; every signed class has multiplicity 2, giving one unordered qualifying split.

The q = 3 boundary is treated separately and correctly. Existence of an abstract balanced residue list must remain distinct from its realization in a supplied set of primes. The note already makes this distinction.

## 5. Independent certificate for the 327941 optimum

The new controls do not use the draft's sieve or discrete-log table. They use trial division for primality/factorization and repeated squaring to calculate each prime's order in the signed quotient.

For q = 65537, the layer counts below q reproduce exactly:

    3260, 1631, 817, 444, 203, 93, 49, 26,
    9, 4, 2, 1, 1, 2, 0.

The only missing layer is 14. Its full residues are exactly 256 and 65281, the two roots of -1. There are only eleven positive integers in those residue classes at or below 327941:

    256 = 2^8
    65281 = 97·673
    65793 = 3·7·13·241
    130818 = 2·3·21803
    131330 = 2·5·23·571
    196355 = 5·173·227
    196867 = 11²·1627
    261892 = 2²·233·281
    262404 = 2²·3²·37·197
    327429 = 3³·67·181
    327941 is prime.

This is a complete finite necessity certificate for the largest-prime lower bound. The displayed 15-prime source supplies sufficiency, and its extension by 5 supplies the generating version with the same largest prime.

Literal enumeration of all divisor subsets independently verifies:

- The 15-factor basis has one divisor per signed class and no qualifying split.
- The 16-factor extension has two divisors per signed class.
- Its product is 36054720350329695914767176044436363972270.
- Its only ordered witnesses are 135625235 and 265840795411928288380603919649882.
- Their sum is 65537 times 4056346726458768152045471341.

There are 6542 primes below 65537. After appending the repair prime, the uniform multiplicity is therefore 2^(6543-15) = 2^6528. The proof obtains this multiplicity without enumerating those primorial divisors. Nothing in this certificate asserts that the source is reached by a greedy trajectory, or that 65537 is its least candidate.

## 6. Kernel recovery and its boundaries

### Stabilizer theorem

The proof is valid for every finite cyclic group. Complementing a subset gives a reflection about total T. Composing any other reflection with it yields a translation, and composing a period with the known reflection yields another reflection. Hence the entire center set is precisely T + Per(D).

Every period orbit has size |Per(D)|, so summing integer multiplicities gives |Per(D)| dividing 2^k. Together with Lagrange's theorem this forces a 2-subgroup. If h = 2^s m with m odd, the subgroup of order 2^r is generated by 2^(s-r)m. It lies in the even residues exactly when r < s. This proves the parity-center assertion. The draft correctly avoids claiming that every reflection center is realized by another binary list.

### Fermat phase recovery

For h = 2^s, a nonuniform kernel has a proper period subgroup, and every period is even. All reflection centers therefore have the same parity. The phase is recoverable from the labeled kernel.

For uniform multiplicity one, k = s and exactly one label occupies layer 0, so the total is odd. For uniform even multiplicity, the two explicit constructions have the same length and identical uniform arrays but opposite total parities. Thus the abstract-list ambiguity classification is exact.

The Fermat restriction matters. For example, at h = 6 the lists (0,3) and (3,3) both give the nonuniform kernel (2,0,0,2,0,0), with opposite total parities. This is covered by the independent exhaustive controls. Also, the q = 7 hostile in the companion packet is not a contradiction: h = 3 is odd and the quadratic character does not descend through the signed quotient.

### Root-class hostiles

Independent literal checks reproduce both:

- q = 5: N = 11 and 19 both have kernel (2,0) and positive phase for -N, but ordered witness counts 0 and 2.
- q = 17: N = 86 and 118 both have kernel (2,1,0,0,0,0,0,1) and positive phase, but ordered witness counts 0 and 2.

The relevant root classes differ: at q = 5 they are [2] and [1]; at q = 17 they are [4] and [1]. This verifies that recovering phase need not recover the relevant root class or eligibility without support. These sources do not satisfy the smallest-missing-prime premise, as the note explicitly says.

## 7. Reproduction and independent controls

The copied prime verifier passes in both normal and optimized Python with 79,609 always-active checks. Its JSON receipt exactly reproduces SHA256:

    596864b2b827961bd4cc37a531b81928d240f90ee8cfa674114d659874cadda6

The copied kernel verifier also passes in both modes, reproducing 74,612 cyclic lists, 54 ambiguity pairs, and both root hostiles.

The independently written `independent_controls.py` imports neither supplied verifier. It checks:

- All 118,754 cyclic residue multisets with 1 <= h <= 24 and 0 <= k <= 4, using literal binary addresses.
- The balancing equivalence, reflection-center coset, period-order restriction, parity criterion, and same-length phase-collision boundary.
- The complete q = 65537 prime layer census by signed multiplicative orders.
- Every possible repair candidate through 327941 by trial division.
- The entire 15- and 16-factor divisor certificates.
- All four actual arithmetic hostile pairs and the unsaturated q = 11 control.
- The first ten greedy outputs by factoring every actual unordered split sum, independent of any residue-support theorem.

It passes **553,351 always-active checks** in both normal and optimized Python, with identical output and JSON receipt. The independent receipt SHA256 is:

    ea311646bdac83621b54b732c42dc9d563964274b2671dd7b66a9c2d60c2b59d

Reproduce from this directory:

    python3 -B independent_controls.py
    python3 -B -O independent_controls.py
    python3 -B ../verify_prime_branching.py
    python3 -B -O ../verify_prime_branching.py
    python3 -B ../verify_kernel_phase.py
    python3 -B -O ../verify_kernel_phase.py

The finite checks corroborate proofs and certificates. They do not certify unbounded completeness, recurrence, waiting-time bounds, or literature priority.

## 8. Nonblocking presentation recommendations

1. Cross-link the kernel-recovery companion beside the first packet's discussion of the retained phase marker. The original statement is a valid global failure-of-reconstruction claim, while the companion identifies the regimes where that marker is already recoverable.
2. Preserve the qualifier that reflection calculations occur in multiplicative quotient/logarithmic coordinates. The display order 1,...,(q-1)/2 of signed representatives is not itself the additive cyclic coordinate order.
3. Keep the 327941 statement tied to squarefree prime-factor sources with q not dividing N, as in the surrounding setup; do not reuse it for unconstrained integer factors.
4. Retain all existing distinctions between chosen branches, greedy split-sum selection, and original Mullin. No evidence in this review upgrades either open completeness question.

No mathematical edit is required for acceptance of the reviewed claims.

## 9. Supplemental audit: appended odd-order Section 4

The author expanded `KERNEL_PHASE_RECOVERY.md` immediately after the initial review. The new Section 4 has now also been independently reviewed. Its current SHA256 is:

    83f2f7a3c0fbaa904c69d4e7527a8cb656de25665d4042d59e1e1bc62668061f

**Verdict remains PASS, with no required correction.** The initial and final source-hash records are preserved in `SOURCE_PROVENANCE.json` in the parent directory.

The extension to any finite abelian group of odd order is valid. In that setting, inversion and translation still compose to give the reflection-center coset T + Per(D). The period subgroup order divides both the odd group order and the power of two 2^k, so it is 1. Hence the center, and therefore total T, is uniquely determined. Cyclicity is unnecessary for this odd-order conclusion.

At q = 3 modulo 4 the group of signed residues has odd order, so the full signed divisor kernel does recover [N]. But the Legendre character takes value -1 on -1 and therefore cannot be a function on that quotient. This is a separate information loss from ambiguity in the reflection center. In the actual q = 7 example, 330 is 1 modulo 7 and 510 is -1 modulo 7; both have signed product [1] and kernel (6,5,5), while their phases for -N are -1 and +1. The example is fully consistent with unique signed-product recovery.

The added `odd_order_controls.py` directly enumerates 14,092 multisets of length at most three in eleven odd-order groups, including noncyclic groups C3 x C3, C5 x C5, and C3 x C3 x C3. It checks trivial periods and unique centers, then computes the q = 7 example directly in the multiplicative signed quotient. All **28,196 always-active checks** pass in normal and optimized Python; outputs and JSON receipts agree exactly.

    python3 -B odd_order_controls.py
    python3 -B -O odd_order_controls.py

The supplemental receipt SHA256 is:

    c34a3ce7c8dfea1e8a5db39495f1d3f82e3e32fe86dce4d41a627b44cdc06c5b
