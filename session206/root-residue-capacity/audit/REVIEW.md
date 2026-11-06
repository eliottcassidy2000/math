# Independent root residue-capacity review

Verdict: PASS for ROOT_RESIDUE_CAPACITY.md in session206-root-residue-capacity.

## Finite-state proof

The state (m mod 2p,Q(m) mod p) carries both the residue and parity needed for the exact digit update. For every target state the target first-coordinate parity fixes the appended bit, and the doubling congruence has exactly two predecessor first coordinates. For each predecessor, the old Q residue is uniquely determined. Thus the transition matrix is doubly stochastic and uniform is stationary.

The proposed path to first-coordinate zero exists once 2^k>=2p. At first-coordinate zero the block 2(2^b-1), with 2^b=1 mod p, preserves that coordinate and adds -1 to Q modulo p. The block identity establishes this for every prior Q coordinate because the updates are additive in Q. The subsequent fixed suffix has a Q increment depending only on the old first coordinate, so any target Q coordinate can be chosen beforehand. This proves strong connectivity for every odd p>1, including composite p. The zero-state self-loop implies aperiodicity. Standard convergence of finite irreducible aperiodic chains therefore applies.

After a fixed suffix, Q is still uniform conditional on the resulting first coordinate: conditional on every old first coordinate the old Q is uniform, and suffix updates are translations. Mixtures over predecessor first coordinates preserve this uniformity. Thus H is uniform modulo p on each fixed endpoint cylinder as its middle length tends to infinity.

## Scalar-cutoff approximation

The established digit identity gives |H(m)-m(2-u_L(m))|=O(X2^-L+log(2+X)) uniformly on the relevant dyadic blocks, all contained in m<=2X. Within a fixed endpoint cylinder, low-bit u_L is constant and the leading-bit condition bounds the variation of m by O(X2^-L). Any cylinder straddling H=X therefore lies in the stated thin approximate boundary.

For each low-bit residue there is a single fixed positive slope, so the boundary consists of an interval of length O(X2^-L+log X), intersected with one progression of step 2^L. Summing the progression counts gives O(X2^-L+log X+2^L). This directly handles correlation between the residue condition and cutoff, without assuming independence. Cutting at epsilon X costs at most epsilon X+1; its boundary cylinders lie in an interval of length O(X2^-L). For fixed epsilon and L the number of dyadic blocks and endpoint cylinders is bounded independently of X, and all middle lengths tend to infinity. Their mixing discrepancies therefore sum to o(X). Taking X to infinity first, then epsilon to zero and L to infinity proves the asserted residue asymptotic. No uniform assertion for growing p is obtained.

## Distinct images and mixed bank

The root involution preserves H exactly, hence preserves each residue-filtered admissible set. Its total fixed-point count is O(sqrt X), also an upper bound in any residue class. The orbit bound therefore halves the residue-specific parameter coefficient. The physical JJ state is (5,6,2,3), whose outputs are 0 or 2 modulo 3, so it contributes nothing to odd class 1 modulo 6. The root-plus-JJ bank misses at least ((1-log2)/6)X-o(X) values. The comparison with the unweighted criterion is exact; its coefficient exceeds 1/2 by the stated rational integral lower bound on log2.

## Regression replay

The author verifier uses explicit RuntimeError checks. Independent optimized replay matches RESULT.json byte-for-byte: 8,988 Markov states over 14 odd moduli, and exact parameter counts through X=1,000,000. These finite counts need not be close to their limits at this scale, especially for slowly mixing moduli; they are diagnostics, not a bound on the asymptotic remainder. The proof establishes only fixed-p o(X), with no effective threshold claimed.
