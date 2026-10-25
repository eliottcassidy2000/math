# Independent audit: score forests and the recovered Cayley mirror bridge

Status: PASS on the mathematical arguments reviewed; finite controls separately scoped below. No new universal theorem is inferred from the finite computations. No remote writes or existing repository edits were performed by this auditor.

## Proof-level audit, not just executable replay

1. Graphic lattice rounding is sound. Every fractional-support vertex has degree zero or at least two: its incidence equation is integral and a sole fractional signed summand cannot be integral. A nonempty finite such multigraph has a cycle, including a two-edge cycle for parallel edges. Its signed circulation permits a nonzero perturbation until a coordinate hits an endpoint, strictly reducing fractional support. A loop is a constant integral translation in the endpoint convention and a zero generator in the signed convention; it is harmless after translation.

2. The potentially delicate projection recurrence is sound. For e=uv, the map from Z^V to Z^(V/e) which adds coordinates u and v has kernel Z(e_u-e_v) and is surjective. Thus its quotient lattice is exactly the contracted graph's lattice, not an index-multiple lattice. Fix any quotient lattice point of the projected deletion zonotope and any integral representative. The line fibre in that coordinate is a nonempty closed interval [a,b]; a or b need not be integers. Its lattice count is floor(b)-ceil(a)+1, which is nonnegative for every nonempty interval. Extending its upper endpoint by integral t increases this count by exactly t, including a fibre initially containing no integer. Projection does not gain any new quotient fibres. Summing therefore proves E_G(t)=E_(G-e)(t)+tE_(G/e)(t). Loops delete; parallel edges remain distinct. Forest deletion/contraction has precisely the same rules. This proves the all-size identity independently of finite tests or a literature black box.

3. Path inclusion-exclusion correctly contracts each selected path-edge set into consecutive blocks. Additional internal edges or two parallel contracted edges create cycles and are excluded. Weighted Prufer counting gives M^(h-1) times the product of component-block sizes for a component on h+1 blocks. Equal-sized blocks are still distinct: the binomial multiplicities in the multiset recurrence and k!/product a_i! in the outer compression are necessary and present. The component-root subset recurrence counts each forest once by its unique component containing the least vertex. For its tree weights, sI-J-L_H has nonconstant eigenvalues s-lambda_i(H); det(sI-L_H)/s^2 is thus the tree cofactor. This remains zero in disconnected cases and equals 1 for a singleton. Path-run continuants use the correct endpoint and internal diagonals.

4. Mirror branches are correctly separated. h(x)=log(x)/x increases up to e and decreases thereafter. T on the principal branch lies in [0,1], F=exp(T), and h(F)=T exp(-T)=z. Thus for n>e, F(log(n)/n) is the distinct small base; for 1<n<e it is n itself. The exact pair 2,4 guards this distinction. At e the two roots coalesce, and for n<=1 no distinct positive root exists.

5. Cayley coefficients are rooted forests, not unrooted forests. Removing vertex 0 from a labeled tree on m+1 vertices supplies a unique chosen root in each remaining component, and the inverse adds exactly those root edges. Hence m! a_m=(m+1)^(m-1). The determinant specializes to the same rooted count, whereas the score-fibre theorem counts unrooted forests. Triangle totals 16 and 7 are genuinely different objects, not a normalization error.

6. The coefficient ratio is below e; Stirling gives a_m asymptotic to e^(m+1)/(sqrt(2 pi) m^(3/2)). Consequently absolute convergence includes the whole boundary circle. The geometric tail holds inside the circle, while the displayed O(M^(-1/2)) explicit majorant is uniform on the closed disk. Rational log2 enclosures and the rational e upper bound in the producer verifier are valid. Floating-point conversions in its certificate remain illustrations; the separately named outward_lower_30dp and outward_upper_30dp fields are correctly rounded outward using exact integer division.

## Independent finite controls

Run `python audit.py` and `python -O audit.py`. Both give identical PASS JSON, with explicit RuntimeError checks and an exercised failing sentinel. The independent program checks:

- 54 multigraphs on three labeled vertices: each nonloop multiplicity 0..2 and a possible loop, at dilations 0..3. This produces 648 edge-specific deletion/contraction checks, plus equality between integer endpoint-allocation counts and forest polynomial evaluations. The fact these allocations exhaust lattice points is the proved rounding lemma, not something the finite test assumes is empirically established.
- Direct edge-subset forest coefficient enumeration through n=7, using merging explicit vertex sets rather than the producer's union-find implementation.
- Independent constant, linear, quadratic and cubic coefficient structure through n=14, including triangle subtraction in degree three.
- Degree-32 exact rational fixed-point iteration B -> exp(zB), distinct from the producer's differentiated exponential recurrence.
- Rational nonintegral fibre endpoints including intervals with no integer points.

The producer's full checks are additionally replayed normally and under -O. Their polynomial equality through n=14 is FINITE-EXACT; the run-polynomial and compressed-composition identities are proved for all sizes. Normal and optimized producer stdout match. The independent auditor initially encountered a trailing-zero representation mismatch at n=2, corrected by padding the producer's trimmed polynomial; this was an audit-harness issue, not a producer mathematical error.

## Priority and limitations

The tree-series coefficients, mirror input, asymptotics and stronger analytic bounds are already in the supplied Session31 TREE_MOMENT_BRIDGE.md; Session31 SYNTHESIS.md also explicitly states the complete rooted-forest/mirror connection. They must be labeled RECOVERED, not a first discovery. Rooted-forest combinatorics and graphical-zonotope/Tutte counting are classical; this audit accepts the self-contained proofs rather than claiming to establish bibliographic novelty. No A000568 identification, CM implication, MRB defect theorem, or transcendence consequence follows. External primary-source bibliographic claims have not independently been browsed by this audit; source pins here certify the exact local inputs reviewed, not historical priority beyond those inputs.
