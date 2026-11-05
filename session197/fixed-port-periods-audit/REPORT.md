# Independent review: continuously fixed-port periodic births

**Result: PASS**, for the final snapshot listed below. The only requested correction was to specify a nonempty initial tournament, matching the interface's omission of the empty cover. The author made exactly that scope clarification. No substantive mathematical correction was needed.

## Bound source snapshot

Files are from `session197-hybrid-constructions/fixed-port-periods/`:

- `PROOF.md`: `d3fd5253d89a842456ca2d5cca1622390740edf1ff53e4ccb513eb91a783bef4`
- `interface.py`: `733155121ac7a5936a38a0ec1f3a602759f868c1e4e547fd4b6f89bced25f9cb`
- `verify.py`: `5f6e76dd3dee5bf2fdfd30340bfb3e7a12fd3cb6173871a6f135805024bd9fda`
- `VALIDATION.json`: `b0c762412d4578faf47e988745fc2637d8d43187413a4ff0123d1c8a6691cea9`
- `VALIDATION_OPTIMIZED.json`: `b0c762412d4578faf47e988745fc2637d8d43187413a4ff0123d1c8a6691cea9`

Copies, replays, and independent checkers are in this review directory. The original source was not modified by the reviewer.

## All-order proof review

1. **Physical cover closure.** Deleting a rightward newborn cannot introduce a non-port start, because every old outneighbor of that newborn is one of the fixed ports. Multiple non-port ends may be represented separately by their components' distinct starts. Internal ports need not be assigned to components in the exposed shape: they cannot receive a new path-cover incidence, so their hidden placement does not affect future transitions. The reverse insertion/deletion correspondence is exact.

2. **Triangular potential.** Existing path-cover edges are never removed on insertion. Used port incidences are nondecreasing, and an existing non-port start cannot be consumed. Equality of the potential excludes attachment to any port endpoint and excludes creation of a new non-port start. The only remaining operation appends the newborn to one existing non-port end, leaving the identical shape. Its multiplicity is exactly the number of those ends and is independent of the mask. This proves simultaneous triangularity for arbitrary finite port count, including zero.

3. **Periodic spectrum and eventual recurrences.** A length-L product has diagonal d^L with d in {0,...,k+1}. Cayley-Hamilton works without diagonalizability. Restricting an annihilated sequence to a tail does not lose the annihilator. In an eventual affine tail with multiplier b>1 and nonzero exponential coefficient, linear independence of b^t and the constant sequence forces b to be a root. Hence b=3 forces L=1. An eventually constant tail is not covered by the stated nonconstancy premise.

4. **One-letter clone interpretation.** With one repeated mask, all newborn vertices have identical relations to the original graph and form an ordered transitive module. Treating the first newborn as the replaced vertex is legitimate on the nonempty tail; no claim is needed about the m=0 skeleton formula. Compressing maximal clone runs gives the exact ordered-set-partition expansion with coefficients j! S(m,j).

5. **The +4 obstruction.** Permuting the j-1 nonempty internal fragments is free because their ordinary labels are disjoint. It preserves both the fragment-internal edges and the root boundary conditions, proving (j-1)! divisibility of the skeleton coefficient. A nonconstant eventual multiplier-three clock has largest run index three, by the positive top coefficient and exponential independence. The +4 constant then gives H_m=a_3 3^m-2; a_3 is even, contradicting the standard odd Hamiltonian-path count of a tournament. This reasoning is unaffected by initial transients.

6. **Odd prime corollary.** For an odd prime p the same spectral argument forces L=1. The top coefficient a_p is even and the constant exponential coefficient is integral. Oddness of H_m makes that constant coefficient odd, giving c congruent to p-1 modulo 2(p-1). This is a necessary condition only, exactly as stated.

7. **Scope.** The theorem uses permanently fixed physical ports, rightward-only births, and a fixed periodic mask word. Edge reversal gives the leftward-only counterpart. It does not address mixed-sided words, changing ports, rebuilding, or state-dependent choices. No bounded-backward-span condition is assumed or concluded.

## Replays and independent controls

- Source normal and optimized runs agree byte-for-byte with both frozen receipts: 40,728 formal transition entries through four ports, 2,764 physical profiles, and 16,597 physical birth checks.
- `independent.py` reconstructs complete path-cover profiles using endpoint subset dynamic programming plus canonical set partitions. It does not use the source permutation/cut enumerator. All 2,764 profiles and 16,597 actual births agree with the source transfer. Its physical-data digest is `560080a521dbd9a38e1ceb28ba358828d3871f75351690032dc1819a440413a4`.
- The independent checker also verifies the potential/diagonal conditions at five ports: 2,930 shapes and 487,872 transition entries.
- `clone_checks.py` independently enumerates repeated-root skeletons for every labelled rooted tournament through order five: 5,405 cases, including 840 with a nonzero three-run coefficient. Every factorial divisibility test passes. It also compares the skeleton/Stirling formula with literal graph path counts for 1,725 clones.

These finite checks test implementation and conventions. The symbolic potential and fragment-permutation arguments establish the all-order claims.

## Additional checked consequence, not needed by the source

There is a convenient explicit annihilator bound. For a length-L product P, put

`q_L(x)=product_(d=0)^(k+1) (x-d^L)`.

The matrix q_L(P) has zero diagonal and can have nonzero entries only when the potential strictly increases. Since the potential takes values in {0,...,2k+1},

`q_L(P)^(2k+2)=0`.

Thus one may take an annihilator of degree at most `2(k+1)(k+2)`, independently of the much larger interface dimension. This is an immediate refinement of the reviewed proof, not an additional realization claim.
