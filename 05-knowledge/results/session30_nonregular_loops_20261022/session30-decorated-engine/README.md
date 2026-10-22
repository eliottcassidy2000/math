# A finite decorated-core engine for first and second almost-regular order

Session30, 2026-10-22. Status: exact core classification and formal local
run-decomposition construction; independently checked finite coefficients and
implementations. Uniform asymptotic remainder remains OPEN. Numerical predictions
are parameter-free but are not certified error bounds.

## 1. Object and grading

Let S=A-A^T, M=S/n and d=S1. The forest polynomial is
Φ(z)=Σ_F product(S_e)(z/n)^|F|, and κ=2^(n-1)H/n!.
Its logarithm is the sum of connected run-position quotient graphs with their
partition-lattice Möbius weights. A graph has excess e=E-V and ell terminal
vertices of degree1. Suppress vertices of degree2 to obtain its reduced core.
The exception with no remaining core is a cycle, whose sum is X0=logρ(z).

For almost-regular even tournaments, |d_i|=1. Each terminal vertex can be
summed first, producing a score-defect vector of size O(1/n), rather than a
free vertex sum of size n. Consequently the correct near-regular grading is

    w=e+ell.

This also applies to families with any fixed bound |d_i|≤B, with constants
depending on B. It does NOT give this grading uniformly for arbitrary nonregular
tournaments with score deviations growing with n.

At internal reduced-core vertices degree is at least3. Since
3V_internal+ell ≤ 2E and E=V_internal+ell+e,

    V_internal≤2e+ell,     V_total≤2(e+ell)=2w.

Thus every required core through weight2 has at most4 vertices. Enumerating all
connected multigraphs, allowing loops and repeated edges, with these degree
conditions is a finite complete classification, not a heuristic graph list.

The exact inventory is:

- Weight1: 5 cores. Three vacuum two-loop cores (e,ell)=(1,0), one lollipop
  (0,1), and one path (-1,2).
- Weight2: 26 cores. Fifteen vacuum three-loop cores (2,0), seven bicyclic
  one-terminal cores (1,1), three unicyclic two-terminal cores (0,2), and one
  three-terminal tree (-1,3).

core-census.json contains every edge list and its full automorphism factor,
including parallel-edge permutations and loop flips. The census was independently
reconstructed by the theory lane with the same5+26 result.

## 2. Local construction and inherited source

loop_engine.py and loop_engine_fast.py are unchanged copies from the repository's
04-computation/experiments/sierpinski_ocf_20261022, at source commit
1b09e49efac420c4e78181a5777b00bd449dd44a. Their previous publication was not itself
an independent audit of the later three-loop claims. This packet independently
checks their extension to degree1 vertices; do not inherit the old asymptotic
claims as theorems.

At each core vertex choose a partial matching of incident half-edges. A paired
pair means that a run passes through the vertex. If the degree is k and there
are r pairs, the quotient block has size k-r and receives Möbius weight
(-1)^(k-r-1)(k-r-1)!. A degree1 vertex has exactly the one unpaired option of
weight1, so the construction extends without an extra rule.

Each subdivided core edge consists of one segment, or at least two segments.
One-segment paths use R_s=s zM(I-s zM)^(-1). Multiple-segment paths with end
orientations s,t use

    Q_st=-R_s R_t (I-z²M²)(I+z²M²)^(-1).

Internal free run orientations are summed, closed runs are excluded, and the
remaining macro-runs receive both orientations. Divide the resulting homomorphism
contraction by the full core automorphism factor. This is the same exact
Möbius expansion as the direct subdivided-edge matcher, with degree2 chains
summed first. These rational propagators converge for |z| in a neighborhood of1:
for every skew tournament M, paired eigenvalues give ||M||²≤(n-1)/(2n)<1/2.

Each propagator has O(1/n) entries uniformly. Its terminal row/column sum is
O(1/n) for |d_i|≤B, by writing the resolvent as I plus an O(1/n)-entry matrix.
Except for the isolated two-terminal path core, each leaf has an internal
neighbor. Summing those ell terminals and then the V_internal remaining labels
bounds each such fixed core by O(n^(V_internal-E))=O(n^-w). The path core is
handled separately by its exact value -δ^T(I+M²)^(-1)δ=O(1/n); its two terminal
sums must not be counted as independent factors on the same edge. This proves
the stated scale for each
of the finitely many included cores. It does not bound the omitted infinite sum.

## 3. Parameter-free inverse-factorial transform

Write A1(z) and A2(z) for the sums of the5 and26 cores respectively. Put
D=z d/dz, evaluated atz=1. Let c_j=D^j X0, and

    b1=c1,
    b2=c2+c1²,
    b3=c3+3c1c2+c1³,
    b4=c4+4c1c3+3c2²+6c1²c2+c1⁴,
    p1=(b2-b1)/2,
    p2=(3b4-2b3-3b2+2b1)/24.

The exact finite-forest identity uses weights n^L/(n)_L. Expanding that factor
formally through second order gives the following two predictions:

    Y1=A1+p1/n,
    Y2=A2+A1²/2
       +[A1 p1 +(D²A1 +(2c1-1)DA1)/2]/n+p2/n²,
    κ1=ρ(1+Y1),      κ2=ρ(1+Y1+Y2),
    y1=-nY1,         y2=-n(Y1+Y2).

There are no fitted coefficients. predict.py implements these expressions.
Its derivatives are computed by a16-point Cauchy sum in log(z), radius0.025;
its X0 derivatives use the absolutely convergent odd-power spectral series.
Changing derivative resolution/radius is a numerical stability check, not a
proof of the omitted-loop remainder. First-predictions.json records initial
values only; held-out exact-H comparisons live in the separate testing packet.

## 4. Independent checks and exact scope

1. check_kernels.py compares resummed descriptor coefficients to an independent
   expanded-edge local matcher (exact_kernel.py, developed by the theory lane).
   All31 cores and all1618 positive edge-length assignments with total at most8
   pass exactly as integers. The earlier cutoff6 comprises344 checks.
   This is independent of the macro-run resummation but shares the general
   Möbius combinatorial derivation, as it should.
2. The independently implemented theory engine imports no code from this packet.
   All124 per-core comparisons (31 cores on transitive4, transitive6, regular7,
   and a vertex-deleted regular7) agree within5.33e-15. That is a floating-point
   implementation cross-check, not an exact identity proof by itself.
3. The audit lane independently enumerates connected partitions of run positions
   at total lengths2,4,6 and groups them by(excess,terminal count). Exact rational
   references for a weighted3-by3 skew matrix and transitive4 agree with Cauchy
   coefficients from this engine to below2e-8. See the separate audit packet.

The kernel check remains active under python -O because mismatches raise explicit
exceptions. Reproduction:

    OPENBLAS_NUM_THREADS=1 python3 decorated.py
    OPENBLAS_NUM_THREADS=1 python3 check_kernels.py 8
    OPENBLAS_NUM_THREADS=1 python3 fixtures.py
    OPENBLAS_NUM_THREADS=1 python3 predict.py

## 5. Boundaries

A complete finite set of formal first/second-order terms is not a proof that
κ/ρ has a uniform asymptotic expansion. That requires controlling all higher
cores, the inverse-factorial transform, and the large-L tail simultaneously.
Neither the old regular conjecture nor its almost-regular extension is promoted
here. In particular no universal scalar second-order constant is claimed:
score-flow and rooted structure remain in the contractions.

## 6. Four, not five, propagator coordinates

The inherited five-matrix implementation has an exact linear redundancy:

    Q_++ + Q_-- + R_+ + R_-=0.

Indeed R_++R_-=2X²/(1-X²), while the two equal-sign Q terms
sum to -2X²/(1-X²), X=zM. Eliminating Q_-- gives an exact four-matrix
basis R_+,R_-,Q_++,Q_+-. The tensor size becomes4^E instead of5^E,
in particular4096 rather than15625 for a six-edge core. No spectral or
combinatorial information is discarded. evaluate() uses this reduced basis
by default; reduced=False retains the original five-coordinate calculation.
All124 fixture/core comparisons agree within1.78e-15.

export_kernels.py exports the exact integer coefficient tables in this
four-coordinate basis, with a separate automorphism denominator. This is a
stronger reproducibility object than floating-point evaluations at a few S.

The fresh independent theory implementation also matches every coefficient of
all31 reduced four-basis symbolic tensors exactly as rational numbers, before
substitution of any matrix. See its symbolic_comparison.json and
compare_symbolic.py. This removes dependence on selected numerical fixtures
for the local formal mechanism. The omitted infinite-core remainder remains
unchanged and open. check_prediction_stability.py additionally compares Cauchy
derivative radii0.025/0.05 and16/32 points on deleted Paley orders6 and10;
the differences in predicted quantities are recorded in prediction-stability.json.
