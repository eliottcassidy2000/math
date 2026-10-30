# Independent audit: complete Collatz entry classes and Lucas phases

**Session 115. Verdict: PASS for both scoped packets.** This is independent algebraic proof review and a separately written exact checker, not formal proof-assistant certification. No arbitrary-controller, finite-word universal bound, or Collatz convergence theorem is inferred.

Audited: `collatz/ONE_RESET_CLASSIFICATION.md` and `collatz/CONSTANT_CONTROLLERS.md`.

## 1. Inherited ordinary means and word faithfulness

The Session114 marked-periodic theorem supplies the ordinary, untruncated mean and its O(log N) error; it is not merely a limit of fixed-depth means. Its common-root/LTE treatment of infinite-depth phases remains applicable to these positive words. All entries and exits are fixed while the controller repeats.

The faithful positive-word representation needed for equality classification can be checked directly. The affine slope recovers the length d and total clock A from its powers of3 and2. For d>1, the exact valuation of B−3^(d−1) at2 recovers the first letter, since the remaining bracket starts with an odd term and all later terms have positive additional powers of2. Stripping that letter recovers the suffix recursively. For d=1 the clock already determines the letter. Consequently commuting affine maps of positive words give commuting words, whose elementary centralizer consists of powers of a common primitive word.

## 2. One-reset family

For w=1^L2, the displayed coefficients, residue order4 and lifting depth

    r_L=v5(9^(L+1)+16·4^L)

are correct. The all-one subgroup has order2·5^(e−1). The commutator with its j-th power is an odd-unit multiple of3^j−2^j. Therefore membership in the controller subgroup forces that all-one power to be the identity modulo5^e. This proves the trivial subgroup intersection at every depth, including stalls.

A fixed entry selects at most one all-one residue class modulo that order. Its occupancy among0,...,L is at most1+floor(L/h_e), and the empty entry achieves this cap at every layer. This proves the all-layer maximum formula without exchanging limits or optimizing different entries at different layers.

Every layer deficit is nonnegative. An equality entry must have an infinite-depth phase, because otherwise it misses the strictly positive entire high-layer tail. Since w contains exactly one letter2, it is primitive. Word faithfulness therefore gives exactly the completion equation1^j u=w^b. Conversely such an entry selects the residue j in every layer, so it is maximal exactly when every occupancy is maximal. This proves completeness, not merely a list of sufficient entries.

The residue inequalities are genuinely equivalent to the Lucas digit condition. Writing L=2K+delta and j=2a+epsilon first separates epsilon≤delta. If any base5 digit of a exceeds the corresponding digit of K, the residue through that digit is too large: all lower-digit deficits together have magnitude less than the value of that one higher digit. Conversely digitwise comparison proves every residue inequality. The characteristic5 identity for(1+X)^K supplies the nonzero-binomial criterion and the product count.

The strict clock bound below1/2 is valid for all L: T_e≥4 gives the internal sum≤5L/32, and5^r≤25·9^L with9^8<5^11 bounds r≤2+11L/8. The resulting bound is deliberately not sharp. The observed9/32 maximum is not promoted to a universal theorem.

The unique nested stall residue follows from9/4 generating1+5Z5; the target−16/9 lies in that same principal-unit group. This is a finite-precision classification, not an unsupported ordinary average over L.

## 3. Constant-letter family

The multiplier orders2,4,1 in the stated residue classes are correct. In the residue1 multiplier case a≡3 mod4, the translation is a unit, so the **affine** lifting depth is0 regardless of multiplier stalls. The maximum is the full geometric subgroup-order sum. Strict increase of that sum in its cutoff, followed by word faithfulness, proves that precisely powers of the same letter attain it.

The entry commutator numerator is exactly2^a E_a(u); the supplied H=F_(a u) coefficients agree with the stated word-order convention. Thus the inherited cutoff rule gives every nonmaximal entry's exact rate.

For a≥2 in the t=2 or4 cases,2^(at)>3^t and

    5^r≤2^(at)−3^t<2^(at)<5^(at/2).

Since at/2 is an integer, r≤at/2−1. This proves the strict normalized bound below1/2. The remaining residue class has rate1/(4a), and a=1 alone gives5/8, with an all-one entry. Repetition of the same letter changes period and clock by the same factor.

The parameter-Haar3/8 statement is also correctly typed. Within each applicable residue class modulo4, the exponential step is a unit power of the generator2^4 of1+5Z5. Exactly one residue at precision5^(k−1) gives valuation at least k. Monotone summation gives the claimed integrable Haar expectations. This does not assert an ordinary initial-integer average or a clock-normalized parameter mean, and the draft explicitly separates those notions.

## 4. Independent finite checks

Run from this directory:

    python check_classes_independent.py
    python -O check_classes_independent.py

The checker imports no producer implementation. It builds affine coefficients from their explicit finite sums and generates literal cyclic subgroups. It verifies:

- One-reset L=0..20, all entries of clock at most6, and literal depths1..3.
- Trivial subgroup intersections, exact orders, phase membership and layerwise caps.
- Exact maximal-entry equivalence for the declared finite universe.
- Lucas residue equivalence for every phase at L=0..200.
- Constant letters a=1..100 with all entries of clock at most6; literal subgroup tests for a≤20 through depth3.
- Complete finite parameter residue classes modulo4·5^3, including exact tail proportions through level4.

**50,894 explicit gates pass. Normal and optimized outputs are byte-identical.** `independent_receipt.json` records the counts and universes. These are independent controls for the written all-parameter arguments, not a replacement for them.
