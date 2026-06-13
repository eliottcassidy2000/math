COMPREHENSIVE CATALOG OF PER-TOURNAMENT-CLASS INVARIANTS

This catalog lists ALL invariants that have been computed or could be computed
for each tournament isomorphism class (at a given n). These are candidates for
node properties, edge properties, or visualization attributes in a metagraph viewer.

================================================================================
SECTION 1: CORE HAMILTONIAN PATH INVARIANTS
================================================================================

1. H(T) — Hamiltonian Path Count
   Measures: Total number of Hamiltonian paths in the tournament
   Data type: Positive odd integer
   Available for: n ≥ 3, computed exhaustively for n ≤ 8, sampled for n ≥ 9
   Source scripts: hamiltonian-paths.md, Gn_general_s167.py
   Notes: Always odd (Redei's theorem). Range: [1, max(n)]
          Fundamental invariant connecting to W-polynomial and OCF

2. degree (in metagraph G_n)
   Measures: Number of OTHER isomorphism classes reachable by one arc reversal
   Data type: Non-negative integer
   Available for: n = 3,4,5 (computed exactly), n ≥ 6 (computed for n=6,7)
   Source scripts: Gn_general_s167.py, tournament-tiling-explorer.html
   Notes: max degree = C(n,2) - 1 when all arc reversals change class
          min degree = 0 only for transitive tournament

3. gridSym (in tiling representation) / flip symmetry
   Measures: Whether tournament can be represented with grid-symmetric tiling
   Data type: Boolean
   Available for: n = 3,4,5,6 via tiling representation
   Source scripts: tournament-tiling-explorer.html, class_color_type_s20dc.py
   Notes: Only defined for tournaments expressible via tiling (base path fixed)

4. parentClass (at n-1 system)
   Measures: Which class at n-1 this class "comes from" by adding vertices
   Data type: Integer (class index at n-1), or -1 if n ≤ 3
   Available for: n ≥ 4, hierarchical only
   Source scripts: tournament-tiling-explorer.html
   Notes: Multiple classes at n can map to same parent at n-1

5. parentVertexOrbits
   Measures: Number of vertex orbits in the automorphism group representative
   Data type: Positive integer, ≤ n
   Available for: n = 3,4,5,6 (computed via permutation enumeration)
   Source scripts: tournament-tiling-explorer.html, aut_H_deep_s182.py
   Notes: Counts "perspectives" this class contributes via removal

================================================================================
SECTION 2: CYCLE-RELATED INVARIANTS
================================================================================

6. t_3 — Directed 3-cycle count
   Measures: Number of directed 3-cycles (triangles) in tournament T
   Data type: Non-negative integer
   Available for: all n ≥ 3
   Source scripts: cycle-counts.md, cycle_cover_deep_investigation.py
   Formula: t_3 = C(n,3) - Σ_i C(score_i, 2)
   Range: [0, C(n,3)]
   Key property: t_3 determines D_2, D_3 Fourier components of W

7. t_5 — Directed 5-cycle count
   Measures: Number of distinct directed 5-cycles in tournament T
   Data type: Non-negative integer
   Available for: n ≥ 5, computed exhaustively for n ≤ 7
   Source scripts: cycle-counts.md, alpha1_3_structure_n7.py
   Key property: Appears in Fourier coefficient formulas at n ≥ 7

8. t_7, t_9, ... (t_k for higher odd k)
   Measures: Number of directed k-cycles
   Data type: Non-negative integer
   Available for: n ≥ k, exhaustive for small n
   Notes: All k-cycles (k odd) participate in OCF decomposition

9. bc33 — Disjoint 3-cycle pair count
   Measures: Number of disjoint pairs of directed 3-cycles
   Data type: Non-negative integer
   Available for: n ≥ 6, computed for n = 6,7
   Source scripts: fourier-coefficients.md, cycle_cover_deep_investigation.py
   Key equation: H(T) = 1 + 2*(t3 + t5 + t7) + 4*(bc33 + ...) at n=7

10. C_k — k-cycle covers in permanent
    Measures: Count of cycle covers with exactly k cycles
    Data type: Dictionary/array: k → count
    Available for: n ≥ 3, computed for n ≤ 6
    Source scripts: barvinok_cycle_polynomial.py
    Notes: C_1 = ham(A) = Hamiltonian circuits
           sum_k C_k = per(A) = permanent
           At n ≤ 5: per(A) = C_1 only (no multi-cycle covers)

11. per (or per_A) — Cycle cover polynomial sum
    Measures: Permanent of adjacency matrix; total cycle covers
    Data type: Non-negative integer
    Available for: all n, computed for n ≤ 6
    Source scripts: barvinok_cycle_polynomial.py, cycle_cover_deep_investigation.py
    Formula: per(A) = Σ_k C_k
    Notes: Always ≥ ham(A)

12. ham (or ham_A) — Hamiltonian circuits count
    Measures: C_1 value; number of single-cycle covers (Hamiltonian circuits)
    Data type: Non-negative integer
    Available for: all n, computed for n ≤ 6
    Source scripts: barvinok_cycle_polynomial.py
    Notes: Different from H (Hamiltonian PATHS); ham ≤ per

13. per_odd — Cycle covers using only odd-length cycles
    Measures: Permanent restricted to cycles of odd length
    Data type: Non-negative integer
    Available for: n ≥ 3, computed for n ≤ 6
    Source scripts: cycle_cover_deep_investigation.py
    Decomposition: per(A) = per_odd + per_even

14. per_even — Cycle covers using at least one even-length cycle
    Measures: Permanent restricted to covers with ≥1 even cycle
    Data type: Non-negative integer
    Available for: n ≥ 3, computed for n ≤ 6
    Source scripts: cycle_cover_deep_investigation.py
    Notes: Rare at small n; becomes more common with n

15. cover_types — Cycle cover decomposition by partition type
    Measures: Maps cycle-length partition to count
    Data type: Dictionary: tuple(partition) → count
    Available for: n ≤ 6
    Source scripts: cycle_cover_deep_investigation.py
    Example: {(n,): 1, (3,n-3): c, (3,3,n-6): d, ...}

================================================================================
SECTION 3: ODD-CYCLE CONFLICT GRAPH (OMEGA) INVARIANTS
================================================================================

16. Ω(T) — Odd-cycle conflict graph
    Measures: Graph where vertices = odd directed cycles, edges = shared vertex
    Data type: Graph (adjacency list / matrix)
    Available for: n ≥ 3, computed for n ≤ 7
    Source scripts: omega-graph.md, cycle_cover_deep_investigation.py
    Key property: OCF theorem: H(T) = I(Ω(T), 2)

17. alpha_k (or α_k) — Independent set counts in Ω(T)
    Measures: Number of ways to choose k vertex-disjoint odd cycles
    Data type: Dictionary: k → count
    Available for: n ≤ 7 (computed exhaustively), n > 7 sampled
    Source scripts: fourier-coefficients.md, cycle_cover_deep_investigation.py
    Key equation: H(T) = Σ_k α_k · 2^k = I(Ω(T), 2)
    Notes: α_0 = 1 always (empty set); α_1 = total odd cycles

    Specific variants:
    - alpha[0] = 1 (empty set, always)
    - alpha[1] = number of odd cycles = t_3 + t_5 + t_7 + ...
    - alpha[2], alpha[3], ... = higher independent sets

18. n_odd_cycles — Total number of distinct odd-directed cycles
    Measures: |V(Ω(T))|, the vertex count of the conflict graph
    Data type: Non-negative integer
    Available for: n ≥ 3, computed for n ≤ 7
    Source scripts: cycle_cover_deep_investigation.py, omega-graph.md
    Notes: Grows exponentially with n

19. Q(T) (also Q_T or Q_indep) — Odd-cycle covering sum
    Measures: Σ over vertex-covering disjoint odd-cycle partitions of 2^k
    Data type: Non-negative integer
    Available for: n ≤ 7
    Source scripts: cycle_cover_deep_investigation.py
    Formula: Q(T) = Σ over COVERING independent sets of 2^|set|
    Decomposition: H(T) = Q(T) + R(T)

20. R(T) — Partial odd-cycle packing remainder
    Measures: H(T) - Q(T); contribution from non-covering packings
    Data type: Non-negative integer
    Available for: n ≤ 7
    Source scripts: cycle_cover_deep_investigation.py
    Notes: Represents "surplus" from incomplete cycle collections

================================================================================
SECTION 4: FOURIER / W-POLYNOMIAL INVARIANTS
================================================================================

21. c_k — Fourier coefficients of W-polynomial
    Measures: Coefficients in W(T,r) = Σ_k c_k · r^(n-1-2k) [odd n]
    Data type: Integer or rational polynomial coefficients
    Available for: n ≤ 9 exhaustively
    Source scripts: fourier-coefficients.md, W-polynomial.md
    Key equations:
      - c_0 = constant term (varies, at even n always 0)
      - c_1 (D_1) = universal for all T at given n
      - D_2 depends on t_3 only
      - D_k (k ≥ 3) depends on higher cycle invariants

22. c_0, c_1, c_2, ... (or D_0, D_1, D_2, ...)
    Measures: Individual Fourier components
    Data type: Integer or rational
    Available for: n ≤ 9
    Source scripts: fourier-coefficients.md
    Special properties:
      - D_0, D_1: universal (same for all T at fixed n)
      - D_2, D_3: depend on t_3 linearly
      - D_4+: depend on t_5, bc33, bc35_w, higher invariants

23. W-polynomial (full, not just coefficients)
    Measures: W(T,r) = Σ_P Π_{e∈P} (r + s_e), evaluated or symbolic
    Data type: Polynomial in r
    Available for: n ≤ 8 explicitly
    Source scripts: W-polynomial.md, fourier-coefficients.md
    Evaluations:
      - W(T, 1/2) = H(T) / 2^(n-1)
      - W(T, 0) encodes c_0
      - W(T, 1) = sum of F_k

24. F-polynomial (forward-edge polynomial)
    Measures: F(T,x) = Σ_k F_k x^k where F_k = #{HPs with k ascents}
    Data type: Polynomial in x with nonnegative integer coefficients
    Available for: n ≤ 8 explicitly
    Source scripts: F-polynomial.md
    Key properties:
      - F(T,1) = H(T)
      - F_k(T) = F_(n-1-k)(T^op) [complement duality]
      - Worpitzky expansion: F = Σ w_k C(x+k, n-1)

25. F_k — Forward-edge (ascent) counts
    Measures: Number of Hamiltonian paths with exactly k ascents
    Data type: Dictionary/array: k → count
    Available for: n ≤ 8
    Source scripts: F-polynomial.md, worpitzky_*.py
    Notes: Each F_k sums to a different aggregate across all T

================================================================================
SECTION 5: STRUCTURAL/ALGEBRAIC INVARIANTS
================================================================================

26. scores (or score sequence)
    Measures: Sorted out-degree sequence [score_1, score_2, ..., score_n]
    Data type: Tuple of n non-negative integers, sum = C(n,2)
    Available for: all n
    Source scripts: tournament-tiling-explorer.html, aut_H_deep_s182.py
    Notes: Coarser than isomorphism class; constrains adjacency in G_n

27. |Aut(T)| (or aut, aut_size)
    Measures: Automorphism group size
    Data type: Positive integer dividing n!
    Available for: n ≤ 7 exhaustively, n ≥ 8 sampled
    Source scripts: aut_H_deep_s182.py, class_color_type_s20dc.py
    Relationship: |fiber| = |class| = n! / |Aut(T)|
    Notes: Determines compression from individual tournaments to class level

28. isSC (or sc) — Self-complementarity status
    Measures: Whether tournament is isomorphic to its complement T^op
    Data type: Boolean
    Available for: all n
    Source scripts: Gn_general_s167.py, tournament-tiling-explorer.html
    Notes: Partitions all classes into SC and NSC pairs
           NSC classes come in pairs: (T, T^op)

29. complementClass (or compCI) — Complement class index
    Measures: Index of the isomorphism class of T^op
    Data type: Integer (class index), or self if SC
    Available for: all n
    Source scripts: Gn_general_s167.py, metagraph-explorer.html
    Notes: Maps each class to its complement's class

30. size (or |class|) — Class fiber size
    Measures: Number of labeled tournaments in the iso class
    Data type: Positive integer, divides 2^C(n,2)
    Available for: all n
    Source scripts: all computation files
    Formula: size = n! / |Aut(T)|

31. blueCount / gs (grid-symmetric tiling count)
    Measures: How many tilings in the class are grid-symmetric
    Data type: Non-negative integer
    Available for: n = 3,4,5,6 (tiling-based)
    Source scripts: tournament-tiling-explorer.html, class_color_type_s20dc.py
    Notes: Only defined when tournament expressible via tiling

32. blueFraction (or gridSymFraction)
    Measures: Fraction of tilings that are grid-symmetric
    Data type: Float in [0,1], or rational
    Available for: n = 3,4,5,6 (tiling-based)
    Source scripts: class_color_type_s20dc.py
    Formula: blueFraction = gs / |class|
    Interpretation: Measure of "symmetry" from tiling perspective

================================================================================
SECTION 6: METAGRAPH AND HIGHER-LEVEL INVARIANTS
================================================================================

33. eccentricity (in metagraph G_n)
    Measures: Maximum distance to any other class via arc-reversal
    Data type: Non-negative integer
    Available for: n = 3,4,5 computed exactly; n ≥ 6 bounded
    Source scripts: abstract_profiles_s304.py, Gn_general_s167.py
    Conjectured range: diameter = n-2

34. distance_profile (from metagraph)
    Measures: Profile p(d) = #{classes at distance d from this class}
    Data type: Array/tuple: d → count
    Available for: n ≤ 7
    Source scripts: abstract_profiles_s304.py
    Generating function: P(x) = Σ_d p(d) · x^d

35. H_rank (H-value rank within all classes)
    Measures: Percentile or rank of H value in the sorted distribution
    Data type: Integer rank (1 to |V(G_n)|)
    Available for: all n
    Source scripts: Inherent from H
    Notes: Low rank = low H (transitive ~0); high rank = high H (Paley-like)

36. spectralGap (of metagraph)
    Measures: 1 - λ_2 where λ_2 = second eigenvalue of G_n's adjacency/Laplacian
    Data type: Float in (0,1]
    Available for: n = 3,4,5 computed exactly
    Source scripts: Gn_general_s167.py
    Conjecture: gap ~ 2/C(n,2) asymptotically

37. meanReversionCoeff (in arc-reversal Markov chain)
    Measures: Slope b where E[H(T_e)|class i] ≈ a + b*H_i
    Data type: Float in (0,1]
    Available for: n = 3,4,5 computed exactly
    Source scripts: Gn_general_s167.py
    Formula candidate: b ~ 1 - 2/C(n,2)

================================================================================
SECTION 7: SPECIAL INVARIANTS (LESS COMMON)
================================================================================

38. neutralArcCount
    Measures: Number of arc reversals that leave the class unchanged
    Data type: Non-negative integer
    Available for: n = 3,4,5 computed exactly
    Source scripts: Gn_general_s167.py
    Notes: F[i][i] in metagraph; related to class self-loop weight

39. pureBlueClass / classColorType
    Measures: Categorical: PURE_BLUE, PURE_BLACK, or MIXED
    Data type: String/enum
    Available for: n = 4,5,6 (tiling-defined)
    Source scripts: class_color_type_s20dc.py
    Definition:
      - PURE_BLUE: all tilings grid-symmetric
      - PURE_BLACK: no grid-symmetric tilings
      - MIXED: both types present

40. M(T) — Transfer matrix (polynomial matrix form)
    Measures: n×n symmetric matrix with HP-theoretic structure
    Data type: Integer matrix (at r=1/2), polynomial matrix in general
    Available for: n ≤ 7 computed explicitly
    Source scripts: transfer-matrix.md
    Key: tr(M) = H for odd n
    Property: All off-diagonal entries symmetric (THM-030)

41. M[a,b] — Individual transfer matrix entries
    Measures: Sum over subsets encoding path structure
    Data type: Integer
    Available for: n ≤ 7 computed explicitly
    Source scripts: transfer-matrix.md
    Formula involves: E_a(S) and B_b(R) path counts

================================================================================
SUMMARY TABLE: AVAILABILITY AND DIMENSIONALITY
================================================================================

Invariant Type                 Dimension   n≤5   n≤7   n≤9   Notes
─────────────────────────────────────────────────────────────────────
H, scores, degree             scalar        X     X     X    Always available
t_3, t_5, t_7, ...            scalar        X     X     X    Odd cycles only
C_k, alpha_k, c_k             dict/array    X     X     (X)  Dimension grows
per, ham, Q, R                scalar        X     X          Cycle covers
F-polynomial, W-polynomial    full poly     (X)   X     X    Symbolic computation
|Aut|, isSC, size             scalar        X     X          Enum to n=8
cover_types                   dict          X     X          Variable dimension
Ω(T), distance_profile        graph/array   X     X          Complex structure
gridSym indicators, colorType boolean       (X)   (X)        Tiling-dependent
metagraph degree, spectral    scalar        X     X          Only at G_n level

================================================================================
DATA AVAILABILITY BY n
================================================================================

n=3: Trivial. All 2 classes have H ∈ {1,3}
     Full computation of all invariants

n=4: 4 classes, H ∈ {1,3,5,9}
     Full computation of all invariants

n=5: 12 classes, H ∈ {1,3,5,7,9,11,13,15}
     Full computation of all invariants

n=6: 56 classes, H ∈ [1, ...], computation intensive
     Most invariants computed (H, scores, t_3, deg, |Aut|, isSC)
     Some partial (C_k, alpha_k)

n=7: 456 classes, computation challenging
     H, scores, t_3, t_5, degree, isSC computed exhaustively
     alpha_k, c_k computed via cycle enumeration
     Metagraph and some structural invariants computed

n≥8: 6880+ classes (n=8), exponential growth
     Sampling-based for most invariants
     Exhaustive at specific families (regular, transitive, Paley)
     Metagraph properties estimated or computed for subsets

================================================================================
NOTES ON VISUALIZATION AND SELECTABILITY
================================================================================

For a metagraph visualizer, the following are most suitable as:

NODE PROPERTIES (selectable for coloring/sizing):
- H(T): primary, continuous range
- degree: secondary, discrete
- scores: categorical (hash to colors)
- |Aut|: discrete, related to fiber size
- isSC: boolean, binary coloring
- distance from transitive: calculated
- eccentricity: pre-computed for small n
- t_3 (or other cycle counts): discrete

EDGE PROPERTIES:
- weight (in G_n): count of arc reversals between classes
- edge_type: spine (SC-SC), rib (SC-NS), sea (NS-NS)
- transition_probability: weight / (size_i * C(n,2))

ADVANCED PROPERTIES (for analysis layers):
- alpha_k distribution: sparkline or histogram per node
- F-polynomial shape: unimodality indicator
- Q/H ratio: partial packing efficiency
- blue fraction: grid symmetry indicator (if tiling-defined)
- spectral properties: for n≤5 only

================================================================================
