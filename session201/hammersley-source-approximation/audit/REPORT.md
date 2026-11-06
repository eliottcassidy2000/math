# Independent audit: constructive actual-source slope approximation

**PASS, without mathematical correction.** The exact index map, run-error identity, one-sided slope bounds, explicit approximate constructor, corrected auxiliary rectangle estimates, and exact J residual are valid in the stated domains. These are constructive approximation results for actual indices, not exact target occupancy.

## 1. Frozen input and independent controls

`INPUTS.json` freezes the seven source/manifest files. Principal hashes:

- `PROOF.md`: `e8e169f3a3af83ba87b4cc24387710e668353f287ff5d1876457daaff1c63454`.
- `verify.py`: `49b32ffefd693c41ba89137126a738e73b388343642f0c74389719d2914a074d`.
- Original manifest: `ccccf2ff4b87ba5a2a194ff96f22afb3cce6e691bf473d92090b5b2a6061a721`.

No source file was edited. The author normal and optimized runs match their frozen receipts. The independent checker imports no author code. It uses full physical source profiles, a separately proved weighted-digit evaluation of g, closed prefix formulas for run errors, exact strict-cutoff CDF limits, and sorted-slice rectangle counts instead of the author's two-dimensional prefix-sum algorithm. Its normal and optimized receipts match.

Independent controls cover 32,767 words through length 14, 114,688 complete Z-runs, 98,316 explicit constructors, 130,807 dyadic-block cutoffs, exact rectangle discrepancies through length 9, and 16,129 J pairs. Every overlapping exact error/discrepancy value agrees with the author output, recorded in `CROSSCHECK.json`. The all-order conclusions rest on the proofs below, not these finite ranges.

## 2. Physical index map and empty-word boundary

The chronological convention is essential: epsilon=1 means Z, and epsilon=0 means O. Iterating u->2u-epsilon gives

    u=2^(k+1)-sum epsilon_i 2^(k-1-i).

Iterating t->(t+epsilon)/2 from 1/2 gives t=(A+1/2)/2^k with A=sum epsilon_i 2^i. Thus the source index uses the fixed-length bit reversal of A, not A itself. Leading zero bits are retained.

Reversal is a permutation, so every index in 2^k<u<=2^(k+1) appears exactly once. All are genuine nontrivial sources. The independent checker reconstructs the entire profile (g(u),u mod 2,u,1), not merely its first two coordinates. At k=0 the unique empty word gives u=2, t=1/2, and zero error. The scalar recurrence base g(1) is not used to create a singleton source.

## 3. Exact error per run and all boundaries

Direct substitution gives O-invariance of e=g-ut and the claimed Z increment. At the beginning of each maximal Z-run, E=0 and the starting midpoint satisfies 0<s<=1/2. The seed attains 1/2; a preceding O halves a midpoint strictly below one.

The first Z increment is (s+1)/2. At the j-th subsequent position, the exact decrement is -(1-s)/2^j. Summing gives s+(1-s)/2^r, which is positive and at most (1+s)/2<=3/4 for every r>=1. The argument includes a one-letter run and a run ending the entire word. No terminal O is required, and individual negative increments are not incorrectly treated as positive.

At most floor((k+1)/2) nonempty Z-runs fit into k letters, so summation proves the stated Delta_k bound. At k=0 it is exactly zero. The all-O source has g=1 and zero error. For k>=1 the all-Z source is u=2^k+1 with g(u)=u and error (2^k+1)/2^(k+1). Both endpoints agree with the formulas.

## 4. One-sided star discrepancy and actual approximate indices

The displacement t_A<=s_A<=t_A+epsilon_k follows from the exact error and the minimum index u>=2^k+1. It is a pointwise matching of labeled lists; no monotonicity of the physical slopes in A or in u is assumed.

For strict-cutoff CDFs, t_A<z-epsilon_k implies s_A<z, and s_A<z implies t_A<z. This proves the CDF sandwich with the correct strictness. Applying the midpoint-grid discrepancy gives the one-sided bounds. When z<epsilon_k, the negative error follows directly from F_s>=0, as stated. Slopes equal to one do not create an exception: at z=1 they are excluded by the strict cutoff and contribute to the negative discrepancy. The independent CDF computation explicitly treats that endpoint and does not take a fictitious right limit beyond one.

For every real theta in [0,1], the chosen A lies in 0,...,B-1 and places theta within 1/(2B) of its midpoint. The cap at B-1 is necessary at theta=1. Combining with the one-sided displacement gives both stated approximation inequalities for the actual index u=2B-rev_k(A). The independent constructor controls use thirds of grid cells, in addition to the analytic endpoint proof, rather than copying the author's half-cell controls.

At k=0 the physical slope is 1/2, whose star discrepancy is 1/2, and the constructor works for every theta with error at most 1/2. As k tends to infinity, the displayed discrepancy bound tends to zero, proving uniform distribution of these equally weighted actual-index slope lists. It does not imply any exactly prescribed rational slope is attained.

## 5. Auxiliary Hammersley rectangles and the retained correction

For an aligned block of 2^ell indices, bit reversal splits into the reversed block label plus a fixed step times every ell-bit reversal. Thus its reversed ordinates are exactly an arithmetic progression, with offset ratio eta in [0,1).

The strict cutoff count is ceil(M*y-eta). Because y is in [0,1], this ceiling already lies in [0,M]; no missing endpoint truncation is needed. Its error lies in [-eta,1-eta), giving absolute error strictly below one. This holds at y=0, y=1, zero offset, and the singleton block.

A grid-aligned first cutoff a/B decomposes into at most k aligned blocks when 0<a<B. The a=B case uses one block, and a=0 is exact. For k>=1 this yields the stated k/B bound even when the second cutoff is an arbitrary real number.

For real x, the exact selected indices are A<ceil(Bx), and the area change is less than 1/B. This is precisely the missing term in an unrounded real-cutoff argument. The corrected (k+1)/B bound is valid. At k=0 the auxiliary grid is the single point (0,0), with real-rectangle star discrepancy one. This is a different object from the physical midpoint slope 1/2 and is correctly treated separately.

The stronger all-real k/B statement is false. At k=1 and exact cutoffs x=y=9/16, both points lie in the rectangle and the discrepancy is 175/256>1/2. The actual supremum is 3/4. The source proof retains that correction explicitly. The independent sorted-slice computation reproduces the exact real, first-grid-aligned, and both-grid-aligned discrepancies through k=9, including all singleton cases. No optimal asymptotic joint constant is claimed.

`SOURCE_CHECK.md` verifies the standard Hammersley naming at zero shift and coordinate exchange. No L2 estimate from the cited paper is used as a star-discrepancy premise.

## 6. Exact J residual and scope

Substituting g(u)=u*t_u+e_u and g(v)=v*t_v+e_v into

    Psi_J(u,v)=g(u)g(v)+uv-pq

gives exactly the displayed residual. Every term is nonnegative because both source errors are nonnegative. Replacing midpoints by the upper bound one and errors by their Delta bounds proves the stated upper estimate. Empty words, all-O sources, zero errors, and both index orders are included.

The residual generally grows with the source sizes. It is not known to be smaller than one and supplies no integer-rounding guarantee for a prescribed target. Approximate slope distribution, auxiliary rectangle counts, exact graph recognition, and exact target support remain distinct predicates. This positive approximation theorem is compatible with the separately proved exact-graph obstruction.

This audit does not include the separate union-witness or bounded-error companion proofs. Their reports and any status promotion are separate.

## Reproduction

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    cmp check-normal.json check-optimized.json
    sha256sum -c MANIFEST.sha256

Both the author and independent checkers are standalone standard-library programs. For publication, scratch-specific proof links and reproduction commands should be adjusted only in an explicitly reviewed promotion diff. Original source hashes and receipts should remain archived.

Exact target occupancy and the twelve-template all-target lemma remain **OPEN**.
