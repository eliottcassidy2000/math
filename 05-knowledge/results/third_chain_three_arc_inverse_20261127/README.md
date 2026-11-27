# The third-chain inverse: 189 breaks three arcs, and a five-gap decoder

2026-11-27. This builds on the [fixed-budget path compiler](../constructive_path_compiler_backward_arc_budget_20261127/README.md), the [ternary chain](../chanpak_third_ternary_chain_shadow_20261127/README.md), and the follow-up background calculation. Let b(h) be the least number of backward arcs in any ordered tournament with exact Hamiltonian-path count h. An omission at fixed budget is not an unrestricted Chan–Pak omission.

## 1. Exact answer at 189, and the three-arc range through exponent 23

The chain has b(63)=3, **b(189)=4**, and b(567)=3. Thus a uniform budget of three fails, and its failure is not monotone along multiplication by three.

| h | Order | Reversed pairs in zero-based order | b(h) |
| ---: | ---: | --- | ---: |
| 63 | 8 | (0,3), (0,4), (0,7) | 3 |
| 189 | 8 | (0,4), (0,7), (1,4), (1,7) | 4 |
| 567 | 9 | (0,4), (1,7), (2,8) | 3 |
| 1701 | 11 | (0,6), (1,7), (1,10) | 3 |
| 15309 | 13 | (0,12), (1,4), (2,8) | 3 |

These are all strong tournaments; the replay checks the original arcs, strong connectivity, and the exact counts by independent subset-and-last-vertex DP. The lower bounds use the preceding all-exponent two-arc exclusion and the complete finite search below.

**Complete budget-three inverse search.** Label the endpoints of the three reversals in relative order, letting m be their number. There are 1,16,30,15 patterns on m=3,4,5,6 that use every endpoint. Compress uncrossed ordinary gaps to zero: the omitted sites form ordinal chains, which do not change the count. A reversed pair between consecutive exceptional endpoints must have a nonempty gap, since otherwise swapping these two adjacent vertices improves the order. This is valid when searching for a minimum three-arc realizer, and the previous note rules out budgets ≤2 for this chain.

The cap applies **to an entire interval, not separately to each gap**. For reversed (i,j), induce on its two ends and all ordinary sites between them, omitting all other exceptional endpoints. This has a single reversed extreme arc and therefore 1+2^(sum of the ordinary gaps from i through j−1) Hamiltonian paths. First-fit insertion injects its paths into those of the full tournament. Thus for target ≤X each reversed interval must have total ordinary length ≤floor(log₂(X−1)).

For target 189 the cap is 7. The four endpoint strata give 21, 2227, 13902, 19885 candidates, totaling **36035**. The Python replay evaluates every one with exact integers and none has count 189. It also recovers the surrounding counts 181,183,185,187,191,193,195,197 (a local diagnostic, not a claim about all cardinalities). The eight-vertex witness above proves the matching upper bound.

For the larger window X=7·3^23, the cap is 39. The C++ sweep evaluates 43,272,827 encodings. It compares residues at 65521 and 65519 before exact evaluation of survivors. To justify the exact stage, any nonzero exceptional permutation contributes at most 2^(3·39) assignments: a gap crossed by d reversed intervals has at most 1+d≤2^d eligible slots, then apply all three interval caps. Multiplying by at most 6! still fits within 128 bits. Signed monomials calculated modulo 2^128 therefore recover the actual nonnegative count without truncation. Only exponents **2,4,5,7** occur, with respectively 4,12,2,8 ordered pattern-gap encodings. These are not isomorphism counts. In particular b(5103)=b(45927)=4, by adjoining one directed triangle to the exponent-5 or exponent-7 witness. For the other exponents from 9 through 23 this establishes a lower bound of four only. Nothing here extrapolates the cutoff to arbitrary t.

## 2. Five independent parameters from three reversed arcs

Take the six-endpoint pattern (0,2),(1,4),(3,5) and five positive ordinary gaps (a,b,c,d,e). Permuting the exceptional subsequence and enforcing the nonempty slots by inclusion–exclusion, as in the path compiler, gives

    H−1 = 2^(b+c+d+2) + 2^(a+b+1) + 2^(a+c+d+2)·3^b
          + 2^(d+e+1) + 2^(b+c+e+2)·3^d
          + 2^(a+b+d+e+2) + 2^(a+c+e+2)·3^(b+d).

For a compact independently replayable coefficient certificate: discard every zero-base term because the gaps are positive; the constant coefficient is one. The other base vectors (in the order above) are 12221, 22111, 23221, 11122, 12232, 22122, 23232 with coefficients 4,2,4,2,4,4,4.

For m≥2 put each gap independently in [km, km+floor(m/10)] with slopes k=(8,1,4,1,16). In this box, the seven displayed terms have strictly increasing 2-adic orders: their leading coefficients in m are 6,9,13,17,21,26,28, and the gap perturbations contribute at most 4m/10. The constant offsets are either 1 or 2 and cannot swap the orders. This gives a decoder:

1. Take the 2-adic order of H−1, subtract 2 to get S=b+c+d, and strip that pure power of two. Repeat, subtracting 1 from the next valuation to get Q=a+b.
2. The next valuation minus 2 is R=a+c+d. Thus b=(S+Q−R)/2, a=Q−b, and c+d=S−b. Strip 2^(R+2)·3^b exactly.
3. Strip the next pure power, whose valuation minus 1 is d+e. The following valuation minus b+2 is c+e. Solve for c,d,e from the three pairwise sums.
4. Validate that all gaps lie in the union of the boxes and satisfy the full formula. For a coordinate x with slope k the admissible m interval is [ceil(10x/(10k+1)), floor(x/k)]. Intersect the five intervals with m≥2. This validates the chart without any enumeration, even if m was not supplied.

The decoder returns the actual gaps, not an approximation or just a congruence class. Distinct gap tuples in a fixed box give distinct counts. The formula gives H<2^(40m) in the box, so choosing m=floor(log₂ X/40) proves Ω((log X)^5) distinct values through X. The previous upper bound matches, giving **Θ((log X)^5)** for the three-reversal spectrum. The exponent counts independent lengths in one specific chart; it does not say that five unconstrained geometry parameters realize a dense set of counts.

## 3. Three exact reformulations of today's inverse problems

1. **Finite sensors versus phase constraints.** With the endpoint pattern fixed, the count is a linear combination of products of the functions 0^l,1^l,2^l,3^l,4^l (0^0=1). Evaluations at l=0,…,4 form a Vandermonde matrix, hence a tensor grid of these evaluations recovers the entire coefficient table uniquely. This is like [Neményi's stress sensors](../five_stress_connection_atlas_20261112/STRESS.md): one scalar count is not a reconstruction measurement, while a full-rank family of probes is. The five-gap valuation box is stronger in a different direction: a *single* scalar recovers its parameters because the five measurements can be peeled off at different scales. The [393/347 same-Moon-profile example](../nemenyi_moon_blue_cuts_20261127/README.md) illustrates what is lost if the observable does not distinguish such scales.
2. **A finite certificate for a proposed infinite family.** On a proposed affine gap ray L_i(s)=a_i s+b_i with a_i,b_i≥0, the count for integers s≥1 is a sum of at most 27 real exponential terms with positive bases (zero-base terms are either dropped or constant). The bound 27 is the maximum raw compiler support for the 62 patterns, replayed exactly. Comparing with 7·3^(ds+k) adds at most one base. A nonzero combination of r distinct positive-base exponentials has at most r−1 distinct real zeros: divide by the smallest base to obtain a constant plus r−1 exponentials and use Rolle's theorem and induction on the derivative. So **28 exact equalities on one fixed affine ray certify an identity**; otherwise it cannot produce infinitely many exact hits on that prescribed progression. This gives a precise replacement for guessing a stable pattern from a few residues. It does not assert that a successful ray exists or that all infinite families must be affine.
3. **Marking is not completion.** The [Moore 63-cycle scaffold](../moore_57_63_cycle_spanning_scaffold_20261127/README.md) also has a deterministic reconstruction rule *conditional on actual incidence data*: common-neighbor closure labels all vertices, but does not place the required nonlocal edges. Here an admissible endpoint pattern is the scaffold; a polynomial response and exact tests supply the remaining constraints. The modulus checks in the census are only filters; the signed sum must still hit the exact integer. This is the same distinction as the two channels in the [Cullen–Woodall trace ladder](../cullen_woodall_split_complex_trace_ladder_20261127/README.md): a channel may verify a proposed solution, but cannot produce a successful parameter simply because the channels can be separated.

Replay: **python replay.py** for the 189 certificate, five literal witnesses, formula, and decoder. With GCC or Clang supporting unsigned 128-bit integers, run **g++ -O3 sweep.cpp -o sweep; python replay.py --terms | ./sweep 23** for the extended census. The Python compiler alone uses no third-party libraries.
