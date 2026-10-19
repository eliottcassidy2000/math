# Independent review: sharp induced-book height

2026-10-19. **PASS.** Reviewed `BOOK_HEIGHT_ASYMPTOTIC.md`, SHA256 `08650c90529ab34069465cc25755dc34344eb652788f7ea9aeeee99b5c930628`. The complete positive-integer theorem, the exact construction for every m, the leading asymptotic constant, and the equivalent all-N capacity asymptotic are valid. No mathematical correction is required. No literature-priority claim is certified.

## Approved statements

For B_m=K_2 joined to m independent vertices and H(m) the least positive integer label ceiling of an induced realization,

    log H(m) ~ (log m)(log log m)/log 2.

For b(N) the largest induced-book page count within labels 1,...,N,

    log b(N) ~ (log 2)(log N)/(log log N).

Consequently a uniform polynomial height bound for all 2-degenerate graphs is impossible, even for non-induced realizations. The book family itself is 2-degenerate and induced-realizable for every m. This does not give a worst-case height asymptotic for the entire 2-degenerate class.

## Proof audit

1. **Exact common-neighbor map.** For A<B and a positive distinct common neighbor X, positive square roots r<s give (s-r)(s+r)=B-A. Strict positive factor pairs of equal parity recover r,s and then X uniquely. Positivity, the ceiling, and exclusion of the anchors are correctly applied after the factorization. The unrestricted capacities are floor(tau(d)/2) for odd d, zero for d=2 modulo 4, and floor(tau(d/4)/2) for 4 dividing d. The floor removes the diagonal factorization. Therefore every book, with or without inducedness, has at most D(N-1)/2 pages.
2. **Elementary uniform divisor bound.** With L=log N and y=L/(log L)^3, every exponent is at most L/log 2. At most y small primes contribute O(L/(log L)^2). The large-prime exponents sum to at most L/log y; a+1<=2^a gives their contribution at most (log 2)L/log y. These estimates are uniform over all positive t<=N, including t=1, and establish the required upper direction of Wigert's estimate without a prime-distribution theorem.
3. **Lower-bound inversion.** Distinct positive labels imply H(m)>=m+2. With X=log H(m), Y=log m, one has X>=Y and Y<=(log 2+o(1))X/log X. Hence X>=Y log Y/(log 2+o(1)). This justifies the leading constant; no unproved comparability of log X and log Y is being assumed.
4. **Anchor construction.** Q is odd, squarefree, divisible by 5 and coprime to 6. The chosen S is odd and not divisible by 3, with S>3 sqrt Q. Both anchors are positive integers, V-U=9Q, and U+V=S^2. The stated loose S-gap bound of six is valid. For Q>=5005, 18 sqrt Q+18<=Q, so V<=10Q.
5. **Pages, positivity, and collisions.** Retained divisors satisfy a<sqrt Q/2 and c=Q/a>2 sqrt Q. Thus r_a=3(c-a)/2>(9/4)sqrt Q and r_a^2>(81/16)Q>U. Both roots are integers. Distinct a give strictly decreasing positive r_a, so page labels are distinct. The square identities follow from s_a^2-r_a^2=9ac=9Q.
6. **Genuine inducedness.** Both anchors are 2 modulo 3, while every page is 1 modulo 3. Thus pages cannot equal anchors and every page-page sum is 2 modulo 3, a nonsquare. The center edge and both incidences of every page are squares. This verifies the entire induced graph simultaneously, rather than only its required edges.
7. **Height.** The exact second cutoff a^2 P^(2h)>=Q implies c<=sqrt Q P^h. Consequently X_a<(9/4)QP^(2h); the anchors are at most 10Q. The integer ceiling N_k=10QP^(2h) covers all labels.
8. **Quarter of divisors below the central interval.** Pair d with 5d for every d dividing Q/5. An interval with endpoint ratio four cannot contain both members of a ratio-five pair. Hence at least half the divisors lie outside [sqrt Q/2,2sqrt Q]. Reciprocal pairing gives equal lower and upper outside counts, leaving at least 2^(k-2) divisors with 4a^2<Q. Boundary conventions are consistent.
9. **Exact retained supply.** Uniform divisors correspond to independent Bernoulli prime selections. The logarithmic variance is (1/4)sum(log p_i)^2<=(k/4)(log P)^2. Failure of the second cutoff means a downward deviation greater than h log P. Chebyshev bounds its fraction by k/(4h^2)<=1/8 for every k>=4, since h=ceil(k^(3/4)). Thus at least 2^(k-3) divisors survive both tests. Choosing k=max(4,ceil(log_2 m)+3) supplies m pages for every positive integer m, with no hidden large-m threshold.
10. **Elementary prime growth.** The central-binomial coefficient is at least 4^t/(2t+1). Its p-adic valuation is a sum of zeros and ones, with at most floor(log(2t)/log p) nonzero positions. Therefore its logarithm is at most pi(2t)log(2t), yielding pi(x)>=c x/log x and hence the needed p_j=O(j log j). Since k!<=Q<=P^k and P is the (k+2)nd prime, log Q~k log k. Moreover 2h log P=o(k log k), so log N_k~k log k. Substitution k=log m/log 2+O(1) gives the matching upper constant.
11. **All-N capacity inversion.** For fixed 0<epsilon<log 2, the proposed m_N has log m_N~(log 2-epsilon)log N/log log N and log log m_N~log log N. The proved all-m height asymptotic puts H(m_N)<N eventually. This establishes the lower bound for every sufficiently large N, not merely along the explicit ceilings N_k. Letting epsilon tend to zero is a valid liminf argument.

## Primary-source checks and attribution

The classical maximal-divisor order is independently recorded in Elsholtz–Technau–Technau, *The maximal order of iterated multiplicative functions*, arXiv:1709.04799v2, abstract and Section 1.1, https://arxiv.org/pdf/1709.04799 . The original Rosser–Schoenfeld paper was also inspected: *Approximate formulas for some functions of prime numbers*, Illinois J. Math. 6 (1962), 64–94, printed page 69, Theorem 3 and its corollary, https://doi.org/10.1215/ijm/1255631807 ; scan https://denisevellachemla.eu/Rosser-Schoenfeld-1962.pdf . It supplies a stronger nth-prime estimate than required. Both arithmetic estimates needed by this packet are proved within the note, so these are contextual citations.

The closest two-anchor square-sum mechanism is explicitly present in Dujella–Elsholtz, *Sumsets being squares*, Acta Math. Hungar. 141 (2013), 353–357, first page, https://www.math.tugraz.at/~elsholtz/WWW/papers/papers44.pdf . The final note properly credits that background and makes no external-priority claim for its graph conclusion.

## Reproducible exact controls

The independent `independent_controls.py` imports no producer code. Normal and optimized runs give identical JSON and **105,750 exact checks**:

- Full divisor enumeration and all constructed labels for the first-prime products with k=4,...,14.
- Toggle-five and reciprocal counting; exact integer h and cutoffs; the finite one-eighth supply; positivity, square identities, modular separation, height, and injectivity.
- 7,562 direct pairwise page nonedge checks.
- 1,055 direct bounded common-neighbor factorizations, including injectivity.
- The exact unrestricted parity-capacity formula for every difference d=1,...,2000.
- Ceiling selection for all m=1,...,999 and immediately around powers of two through 2^100.

The producer checker was read in full and rerun normally and with optimization. Both outputs reproduce its saved receipt byte-for-byte: **4,081,492 checks**. Its direct common-neighbor/factor-pair comparison, exact parity capacities, full divisor enumerations through k=18, larger sampled parameters, and complete induced adjacency tests are correctly scoped. No check uses removable assertions. These finite tests are implementation controls; the unbounded theorem is supplied by the audited proof.

No files inside the repository checkout were edited, and no Git write was performed.
