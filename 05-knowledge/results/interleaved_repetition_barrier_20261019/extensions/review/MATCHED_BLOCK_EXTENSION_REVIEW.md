# Independent review: matched factors and the repetition-poor greedy family

Date: 19 October 2026.

## Verdict

**Both frozen extension notes pass this independent proof audit, conditional on the already-audited primary height, carry, first-crossing-length, and short-range lemmas.** No additional external number-theoretic input is introduced by these extensions. The arbitrary-bit-start generalization, the exact-return alternative, the factor-complexity obstruction, and the explicit mechanical coding are valid.

The cycle alternative is essential and is preserved. Low factor complexity proves an internal periodic orbit for a hypothetical non-descending source; it does not prove the original source periodic, strict descent, existence of a non-descending source, or membership in a fixed positive-displacement set. The greedy note consistently separates its infinitely many lawful cylinder lifts from the unresolved diagonal endpoint equation.

One numerical scope point was identified during review and is resolved in the frozen text: the pure-expansion power exclusion has cutoff `m>900` for the displayed V,F parse, but the robust sufficient cutoff after arbitrary exact K^4F/V rewrites is `m>915`. The corrected proof and an explicit countercontrol for the stronger parse-independent reading are included below.

## Frozen material and audit boundary

- `MATCHED_BLOCK_EXTENSION.md`: `ca2838cdff301d3df75c3ce1054e74ad6e71173e5bcc448a57cb624c8774f799`.
- `REPETITION_POOR_GREEDY_BOUNDARY.md`: `1a363f9b3541535b90259034949cfd5021e33c3aabaf9c73da3fda5c382a9f18`.
- `verify_boundary_extensions.py`: `fbe51b6cf77eb0afd5adf60d8efa16cd3588ba66d6317e2cd23ab1b27d6e7516`.
- Submitted verification receipt: `d505b838a3fe1234ede8ca53027e1856a93d849fb13eecbec3713dce6e15f7fa`.

I read the frozen primary manuscript and its existing independent audit. This review accepts those audited lemmas, including their explicit expansion-first scope and the cited logarithmic-form theorem. It does not claim an independent proof of that cited theorem.

## 1. Matched cylinders give the required real gap

An r-bit chronological parity word has a single lawful source residue modulo 2^r. Therefore two unequal actual sources of that word differ by a nonzero multiple of 2^r, and their absolute difference is at least 2^r. Overlap of the occurrences changes none of this: only their individual lawful prefixes are used.

At macroblock boundaries preceding the final crossing, the primary carry identity and multiplier envelope give

`y_P <= lambda^(L/30) (x+N/6)`.

Both states are positive, so their nonzero difference is strictly smaller than their maximum. Hence equal factors at unequal states force

`x+N/6 > 2^r/lambda^(L/30)`.

With `100r>=L`, `log 2>1/2`, and `log lambda<1/12`, this is exactly `log(x+N/6)>L/450`. The primary height ceiling, short-range closure, and overlap contradict this under non-descent. The weaker quantitative matching budget and its prefix-specific expansion-length improvement are correct as well.

## 2. Arbitrary bit starts: the 2^122 adjustment is sufficient

Let q be the number of bits since the last completed bank letter. It is at most 122; using 122 instead of the sharper possible 121 only costs slack. For every positive integer z,

`T(z)+1 <= (3/2)(z+1) < 2(z+1)`.

Iteration gives `y_t < 2^122(y_P+1)`, including q=0. The last completed macroblock prefix is still before the first crossing, even if the bit start lies within the final reset. Thus its carry bound remains in scope. Since `lambda^(L/30)>=1`,

`y_t < 2^122 lambda^(L/30) (x+N/6+1)`.

The augmented height `Y=x+N/6+1` is smaller than the primary augmented height `Z=x+N/6+C`, since `C>1`. Consequently the existing logarithmic height ceiling applies immediately; alternatively the stated source-floor argument proves it directly.

An unequal matched pair therefore forces

`log Y > r log 2 - (L/30)log lambda - 122 log 2 > L/450-122`.

The short range is unchanged. At `L0=2^48`, the primary upper ceiling is strictly below `53B+52`, while the new lower bound exceeds it because `L0/450>53B+174`. Subtracting the constant 122 does not change the derivative; `L0>450(B+1)` makes the gap increasing for every larger L. Thus the arbitrary-start generalization is proved with the displayed constants.

## 3. Return alternative and aligned mismatch

If two actual states agree, determinism gives identical entire future trajectories. Therefore even a single observed difference at a later aligned parity position inside both occurrences rules out equality. It is not necessary to assume disjoint factors or a mismatch immediately after the matched word.

For a nonempty bank separator with multiplier at least one, positive affine carry makes its output strictly larger than its input. In particular this applies when the earlier macroblock start is the original source, because the later boundary precedes the first crossing. These are valid sufficient ways of eliminating the return alternative.

An exact positive return instead satisfies `(1-s)y=e` with positive carry, so its separator has multiplier strictly less than one. Contracting returns really can occur: `F=10` has affine map `(3y+1)/4` and fixes y=1. This reset-only example lies outside the global expansion-first theorem, but correctly demonstrates why equal matched factors alone cannot be converted into unequal states. The extension does not make that conversion.

## 4. Factor-complexity threshold and cycle parameters

Let `P=p_W(r)` and suppose `P<N-r+1`. Since these quantities are integers, this is equivalent to `P+r<=N`. Therefore all starts `0,...,P` support r-bit factors. Pigeonholing these P+1 starts produces equal factors at `i<j`, with `i<=P-1` and `j-i<=P`.

Under hypothetical non-descent, the matched-block theorem forces the actual states at i and j to agree. Thus the orbit is periodic from time i onward, with a positive return period j-i and hence a least period no greater than P. The source may have a preperiod. The strict inequality matters: at equality all available factors can be different.

The numerical sufficient condition follows because

`P+r <= (C+1)r+D < N+1`,

and `P+r` is integral. For fixed `C<99`, the usual nonnegative complexity coefficient satisfies `r<=L/100+1`, so the left side grows strictly more slowly than `L<=N`. Nonpositive coefficients cause no difficulty and are either easier or eventually vacuous.

The stronger continuation test is valid: every equal r-factor pair in a hypothetical non-descending word must have identical available aligned future bits. One witnessed continuation branch is enough to rule out that hypothesis. A small factor count by itself does not supply a branch.

## 5. Exact greedy construction and arithmetic scope

Independent reconstruction gives `v=3^77/2^122` and `v^9<4/3<v^10`. If the iteration begins with `1<=S<4/3`, one V followed by the prescribed optional F restores that same interval. A single optional reset suffices because `v<4/3`. Every optional reset finishes at or above one; the extra final F finishes strictly below one. The literal expansion prefixes of K and V also have multiplier greater than one, as directly checked.

Putting `theta=log(v)/log(4/3)`, the reset count is uniquely `floor(m theta)`, since taking logarithms of the multiplier invariant gives `0<=m theta-c<1`. Hence the nth reset indicator is exactly `floor(n theta)-floor((n-1)theta)`. The prime-factor argument proves theta irrational. From immediately after a reset, the state multiplier lies in `[1,v)`, so eight more V letters cannot trigger the next reset and ten must; complete runs have length nine or ten, initially ten.

For a contracting full affine branch, the positive source equation is exactly

`x=(b-Md)/(M-a)`.

Positivity and divisibility are necessary and sufficient. Full affine integrality authenticates the chronological parity cylinder by the odd-coefficient prefix induction. Since `gcd(M,M-a)=1`, the allowed nonnegative displacements form the stated residue progression `d=delta+j(M-a)` below `b/M`. This does not assert that the progression has any allowed member. The finite controls find none in their tested examples; unbounded fixed-displacement realizability remains open in the note.

## 6. Rotation partition and literal factor complexity

For a phase z in `[0,1)`, a length-k mechanical block has coordinates

`floor(z+j theta)-floor(z+(j-1)theta)`, for j=1,...,k.

The only interior discontinuities occur at the k distinct points `{-j theta}`, j=1,...,k. Along with the boundary at zero they divide the interval into k+1 pieces. The floor functions are right-continuous, so breakpoint values belong to the piece on their right. This proves at most k+1 factors without an appeal to a classification theorem.

The variable-length coding sends 0 to V and 1 to VF, of literal bit lengths 122 and 124. Every r-bit factor can be assigned a containing first image, an offset at most 123, and the next `k=ceil(r/122)+2` reset symbols. Those k images have at least `122k>=r+244` bits, more than enough to include that offset and the whole factor. There are at most `124(k+1)` such descriptions, so

`p_G(r) <= 124(ceil(r/122)+3) <= (62/61)r+496`.

Appending the last two bits introduces at most two additional factor starts. If the prefix is shorter than r, the final word has at most two factors in total. Thus `p_Wm(r)<=(62/61)r+498` holds for every positive r.

At `r=ceil(122m/100)`, the required sufficient inequality is `(123/61)r+498<N+1`. Bounding `r<=61m/50+1` and `N>=122m` proves this already at m=5, with a strictly improving linear margin thereafter. No description-length or automaton-size bound is being substituted for the literal factor count.

## 7. Resolved parse-threshold issue

In the displayed V,F parse, a pure expansion run has at most ten V letters, so a power factor has duplicated length at most `9*122=1098`. Therefore `100r<L=122m` for m>900.

An exact split replaces V by K^4F and reduces the expansion-bit accounting by two. Every resulting expansion run is V^qK^4, V^q, or K^4. A power containing both K and V cannot fit in the first shape because its second copy would place a V after a K. The same duplicated-length cap therefore survives, while `L>=120m`. Thus m>915 is a valid parse-robust sufficient threshold.

The distinction is real: in W_901, keep the first ten V letters and split all later V letters. The first run still contains V^10 with duplicated length 1098, but total expansion length becomes `10*122+891*120=108140`, and `100*1098>=108140`. The frozen note explicitly distinguishes the two cutoffs, resolving the issue.

## 8. Reproduction and independent controls

The submitted checker runs identically under normal Python and `python3 -O`, matching the frozen receipt with 49,321 explicit checks. It appropriately labels finite controls separately from universal proofs.

The separate `independent_extension_controls.py` imports no submitted code. It reconstructs carries by the closed sum over odd positions and constructs lawful residues by one-bit parity lifting. Its controls cover:

- all bank partial prefixes and the 2^122 envelope, with direct positive-integer trajectory checks;
- all overlap constants using exact rationals;
- an independent 96-term rational atanh bracket for theta, exact greedy/mechanical agreement through 2,048 iterations, and certified rotation-partition endpoint conventions;
- exhaustive dyadic-cylinder and contracting endpoint-residue tests for all binary words of lengths 1 through 9;
- exhaustive strict factor-pigeonhole tests for all binary words of lengths 1 through 9;
- matched-state divisibility and full available continuation equality on actual positive trajectories;
- reset-only exact return, nonzero preperiod, equal-slope/different-word, and mixed-parse threshold hostiles;
- greedy words through m=2,048, multiple literal factor lengths, carry inequalities, and the finite endpoint filter.

The independent run made **178,380 checks**: **1,022** exhaustive dyadic-cylinder pairs, **787** contracting endpoint words, **3,368** strict-pigeonhole cases, **860** certified rotation-partition intervals, and **29,388** matched trajectory pairs. The m=2,048 greedy word has N=250,272 and 2,614 distinct factors at r=2,499. The original independent run, shift-optimized rerun, and `python3 -O` rerun agree byte-for-byte. Exact fingerprints are recorded in the accompanying receipt. These finite controls support reproducibility and catch scope regressions; the arguments above, using the audited primary lemmas, supply the all-length quantifiers.
