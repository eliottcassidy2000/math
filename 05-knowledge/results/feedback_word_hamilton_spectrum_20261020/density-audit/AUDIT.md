# Independent audit: binary-product density extension

Status: ACCEPT. The all-cutoff word asymptotic and the consequent upper bound for the distinct intercyclic spectrum are proved with the stated normalizations and error order. No mathematical correction was required. Two minor source-reference/empty-word wording issues were corrected by the producer before freezing.

Audited proof: `session16-hamilton-spectrum/BINARY_PRODUCT_DENSITY.md`.

Final proof SHA-256: `e3d1cec97f2964ac3719d5bb41836cb2a1ee919eb38da2b1362be8f200515cc1`.

This is a separate audit. The accepted three-packet audit and its manifest are unchanged. Their original `SHA256SUMS` hash remains `3340469165eb26906d615a3b38bf4e2b750b2e3dfabd41f44815125d2b79115a`, and every entry was reverified successfully.

## 1. Accepted conclusions and strict scope

Let W(X) count reduced one-apex binary words of all lengths whose ordinary Hamilton-path count is at most X. Then

W(X) = X log 2 + O(sqrt X),

uniformly at arbitrary real cutoffs X tending to infinity. This is an asymptotic for presentations, not numerical values, tournament isomorphism classes, or automorphism orbits.

Complement-reversal preserves H and yields

|S_oneapex intersect [1,X]| <= (log 2 / 2) X + O(sqrt X).

The separately audited exact intercyclic spectrum adds only two dyadic rays, giving the same upper bound for S_intercyclic. Consequently its upper relative asymptotic density among positive odd integers is at most log 2, and the lower relative asymptotic density of its complement is at least 1-log 2, approximately 30.6853%.

These are one-sided bounds for distinct values. They do not assert that either distinct-value density exists, that the upper bound is sharp, or that any of the omitted intercyclic values is globally forbidden. An attainable value outside the intercyclic spectrum necessarily requires two vertex-disjoint directed triangles. The global odd-spectrum conjecture remains OPEN.

## 2. Exact binary-product identity

The scalings are correct. With N the usual m-bit integer and R the bit-reversed complement integer, the product NR/2^m has coefficient 2^(j-i-1) on b_i(1-b_j). The diagonal vanishes. At forward distances at least two this is exactly the coefficient in H-1. An adjacent inversion contributes one in the product but two in H-1. All backward terms form B. Thus

H - NR/2^m = 1 + A - B

is exact for every binary word, including the empty word under the natural conventions.

Adjacent 10 occurrences cannot overlap, so A<=floor(m/2). Summing the geometric backward weights gives B<=m/2. The repaired weak inequality includes m=0. Therefore the uniform error bound |H-NR/2^m|<=1+m/2 follows without a sign assumption on 1+A-B. Dividing by 2^m gives exactly the stated normalized error.

The reference to the original weighted-inversion formula has been corrected from `proofs.md` Section 1 to Section 2. Neither correction changes the reduced-word theorem.

## 3. Independent-prefix coupling and uniform CDF control

For a uniformly sampled reduced word, the top binary digit of N/2^m and R/2^m is one. With s=floor((m-2)/2), the next s digits come from two disjoint subsets of the independent interior bits. The selected prefix pair is therefore independent and uniform. This argument correctly does not declare the full finite variables independent.

Appending independent unbiased tails produces genuinely independent U,V uniform on [1/2,1]. Sharing s+1 leading digits gives coordinate error at most 2^(-s-1). Since all coordinates are bounded by one, the product error is at most 2^(-s). Adding the exact binary-product error yields the stated delta_m.

The product CDF is correctly normalized to density four on [1/2,1]^2. Both branches integrate to the displayed formulas, agree at 1/2, and meet the support endpoints 1/4 and 1. Their derivatives are 4 log(4z) and -4 log z, so the global Lipschitz constant is 4 log 2. The deterministic coupling therefore implies the two-sided CDF squeeze and a uniform Kolmogorov bound, including outside the support and at jump locations of the finite distribution.

Uniformity in z is established, rather than inferred from pointwise convergence. This is precisely what permits the later moving dyadic phase and the O(sqrt X) remainder.

## 4. Exactly two partially occupied word lengths

Every reduced word of length m satisfies 2^(m-2)<H<2^m. The outer inversion proves the lower bound, with m=2 checked separately. For the upper bound, a Hamilton path is determined by the set of chain vertices before the apex. At least the singleton set containing the first chain vertex is invalid, so H<=2^m-1. This argument includes all paths and does not use an insertion matching.

For X>=4, write X=t*2^L with L=floor(log2 X) and 1<=t<2. All lengths m<=L are full, no m>=L+3 contributes, and only L+1,L+2 are partial. The exact count is therefore

W(X)=2^(L-1)-1 + 2^(L-1) C_(L+1)(t/2) + 2^L C_(L+2)(t/4).

The full-bank constant -1 is correct. No word length or endpoint case is omitted at X=4, at exact powers of two, or just below the next power of two.

Multiplying each delta_m by its bank size 2^(m-2) gives

2^floor((m-1)/2) + (m+2)/8.

Thus the explicit aggregate error is exactly bounded by

4 log 2 * [2^floor(L/2) + 2^floor((L+1)/2) + (2L+7)/8],

which is O(sqrt X) uniformly in t. In particular, the linear-in-L term does not spoil the claimed remainder.

## 5. Exact phase cancellation

On 1<=t<2, t/2 and t/4 lie in the upper and lower nonconstant branches of the product CDF. Their substitutions give

F(t/2)=2t-3-2t log t+2t log 2,

2F(t/4)=2t log t-2t+2.

Adding the full-bank contribution one cancels both the constant and the phase-dependent terms, leaving exactly 2t log 2. Therefore

W(X)=X log 2-1+E(X)

with the previously stated uniform error bound. This is an exact algebraic cancellation followed by a controlled approximation, not evidence from a finite collection of phases or a subsequence argument.

## 6. Passing to distinct values

Deleting source/sink padding gives a reduced representative for every nontransitive one-apex value. Complement-reversal acts within the entire H-sublevel set. For even m it has 2^(m/2-1) fixed reduced words; for odd m it has none.

All relevant lengths are at most L+2, so summing the fixed-word counts gives J(X)<=2^floor((L+2)/2)-1. This is O(sqrt X) without assuming every fixed word in that bank qualifies.

The sublevel set has exactly (W+J)/2 involution orbits. Since H is constant on an orbit, the number of numerical values is at most that number, plus one for the transitive value. Different orbits can collide, so only an upper bound is justified and only an upper bound is stated. Substituting the error bound verifies every factor of two and the remaining fixed-point term in explicit inequality (13).

Adding O(log X) intercyclic ray values is absorbed by O(sqrt X). Dividing by X/2+O(1), the count of positive odd integers, yields the claimed limsup and complementary liminf at all real cutoffs.

## 7. Moment checks

For l=m-2>=1, the weighted distance sum sum_(d=1)^(l-1) (l-d)w(d)=2^(l-1)-1 follows by induction. Separating outer, endpoint/interior, and interior pairs therefore proves the exact mean

E H = 9*2^(m-4)+1/2, for m>=3.

The m=2 exception is correctly excluded. Bounded uniform coupling also gives first- and second-moment convergence. E U=3/4 and E U^2=7/12 imply limiting normalized mean 9/16 and variance 55/2304. These are valid checks, not replacements for the uniform CDF estimate.

## 8. Reproduction and genuinely independent controls

Producer normal and optimized replays byte-match both frozen receipts. Its floating-point logarithmic checks are honestly separated from exact combinatorial checks.

The separate `independent_density.py` imports no producer or repository mathematics. It uses a direct inversion sum, apex-first literal path enumeration, and arbitrary-precision rational arithmetic. For analytic controls it encloses logarithms by the positive atanh series with an explicit rational remainder. The enclosure of log 2 has width 1/39774420110079400828226875596. No floating-point comparison decides a check; floating point appears only in labelled approximate output summaries.

Normal and optimized independent outputs byte-match. The verified universes are:

- 2,047 arbitrary binary words through length ten, including the empty word, for the exact scaled identity and product-error ingredients.
- 8,191 reduced words through length fourteen, for duality, cutoff, coupling cells, exact disjoint-prefix frequencies, fixed-word counts, and exact moments.
- 31 literal Hamilton-path controls through reduced chain length six.
- 6,862 CDF left/right jump checks using rigorous rational logarithm enclosures, certifying the uniform CDF bounds for all tested lengths.
- Exact coefficient cancellation of the dyadic phase in the basis 1,t,t log t,t log 2.
- 8,353 integer and non-dyadic rational cutoffs below 8,192, for the exact two-length formula, rigorous analytic error bound, and distinct-value upper bound.

For these cutoffs, the length-fourteen bank is complete by the proved lower cutoff. No omitted later word can affect the comparisons. Normal and optimized commands are:

`python3 session16-feedback-density-audit/independent_density.py`

`python3 -O session16-feedback-density-audit/independent_density.py`

The prior mathematical dependencies were already checked against primary sources in the unchanged first audit: Bang-Jensen--Bessy for feedback two, Irving--Omar/Grinberg--Stanley for the odd-cycle formula, and Busch for strong minima. This new extension uses them only through the accepted intercyclic classification; its product asymptotic is elementary and self-contained.

## 9. Publication handoff

Approved new producer files are the proof, verifier, both receipts, and `BINARY_PRODUCT_SHA256SUMS`. The separate audit manifest also seals the imported `verify.py` and the two mathematical dependency files. Preserve the previous three-packet manifest unchanged, and use this separate audit as the approval record for the stronger all-cutoff density statement. No repository or remote writes were performed by the auditor.
