# Independent audit of quantitative transitive-chart escape

**PASS, with the stated operation-language restriction.** This audit independently reviewed `../induction/PROOF.md`; its source digest is in `receipt.json`. No global spectrum obstruction or native construction is certified.

## Symbolic review

1. The physical directed joint formula is correct. A two-slot cover of a transitive tournament has forced increasing orders in its slots. For first endpoint k and second start i with i<k, the k−i−1 middle vertices may choose either slot. For first endpoint i and second start k, a cover exists exactly when k=i+1. This yields the asserted powers and adjacency term without discarding directed information.
2. Multiplying u by its reversed binary fraction separates the diagonal w, forward-index term 2(J−a), and backward-index term B. Hence J−uR/2=a−(w+B)/2. The bounds 0≤B<w and 0≤a≤w−1 imply the stated absolute error at most m, including one-bit masks.
3. The four low-three-bit cylinders have left endpoints 1/2,3/4,5/8,7/8 and width strictly less than1/8. Multiplication by3 interchanges residue classes1↔3 and5↔7. Each paired interval has separation at least1/8. There is no unnoticed issue with terminating fractions: finite tails are smaller than the stated upper endpoints.
4. Multiplying the two normalized error estimates by3u gives exactly3m+m', with m'≤m+2. Thus the lower bound3u/8−4m−2 is valid. It may be negative for small u, which causes no logical difficulty.
5. For bitlength m≥8, the bound at the start of its dyadic block is increasing with m, since the increment is3·2^(m−1)/8−4>0. At m=8 the bound is14. Therefore |d+2|≤6 really forces u<128; the finite classification does not hide a larger unexplored interval.
6. Comparing coefficients at two distinct legal clocks is sufficient to force2^k f'=3f and the constant relation. Arbitrary negative clock shifts cause no exception: legal clocks still form an infinite tail. Removing powers of2 gives u'=3u for every step. A bounded digit alphabet bounds u effectively, while an infinite lineage forces u_n=3^n u_0. The contradiction is exact, not asymptotic evidence.
7. The optional convergence argument is also valid: bounded d yields summable normalized increments, whereas the reversed-bit cylinders have a persistent separation and the normalization errors vanish.

## Independent executable controls

`verify.py` does not import the source checker. It enumerates actual two-slot assignments to recover every directed coefficient through suffix order8; independently counts full rooted tournaments with TT1 prefix through suffix order6; checks exact rational identities on small and up-to256-bit masks; checks the quantitative inequality on100,000 odd masks; and checks legal positive/negative clock shifts, cutoff monotonicity and the64-mask classification. Full attachments are independently checked for strongness.

All checks use explicit exception guards. Run `python verify.py` and `python -O verify.py`; their outputs byte-match. Detailed counts and the source proof digest are recorded in `receipt.json`.

## Boundaries and attribution

The result concerns transitive suffixes and affine relations at two input clocks (equivalently all legal clocks), not one-clock adjacency returns. It does not apply to a nontransitive suffix output, including the actual S183 missing44 witness. It does not forbid arbitrary same-count representations, nonuniform maps, or other growing physical grammars. No claim about native guards or inherited observer-phase transport follows. The inherited exact chart and earlier one-digit facts are explicitly attributed; the proof adds the effective quantitative bound and bounded-alphabet impossibility.
