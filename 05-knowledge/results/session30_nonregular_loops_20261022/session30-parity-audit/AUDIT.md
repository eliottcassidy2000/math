# Independent audit: parity-adapted reduced-core kernels

**Status: PASS**, with the timing and finite-order qualifications below.

## Scope and provenance

This audit reads `session30-decorated-engine/symbolic-kernels.json` and candidate data files in `session30-parity-basis`. It does not import the producer's transformation, quotient, or evaluator code. The source engine was separately audited before this task; this audit checks the new algebra and contraction reductions, not a new derivation of the entire Hamiltonian expansion.

The audited inventory consists of the 5 weight-one and 26 weight-two reduced cores. No assertion about a uniform asymptotic remainder, correctness of arbitrary higher-weight tables, or a physical electrical-network equivalence is made.

## 1. Rational identities and domain

Put

- G = −M²(I+M²)⁻¹,
- H = M²(I−M²)⁻¹,
- O = M(I−M²)⁻¹,
- P = M(I+M²)⁻¹,
- R_s = sM(I−sM)⁻¹,
- Q_st = −R_s R_t (I−M²)(I+M²)⁻¹.

Then exactly

R+ = H+O; R− = H−O;
Q++ = −H+P−O; Q-- = −H−P+O; Q+- = −G.

One proof is to replace M by a scalar x, multiply through by the common denominator (1−x²)(1+x²), and compare numerators. `audit.py` independently verifies all five resulting rational identities with exact symbolic arithmetic. Matrix substitution is legitimate because all factors are rational functions of the same matrix, hence commute. Equivalently these are formal matrix power-series identities; for numerical evaluation every displayed inverse must exist. Replace M by zM without changing the proof.

For real skew M, G,H are symmetric and O,P are skew-symmetric. This transpose statement also holds for complex z without complex conjugation. Under M ↦ −M, G,H are even and O,P odd. The change from the original four formal symbols is invertible over the rationals: G=−Q+−, H=(R++R−)/2, O=(R+−R−)/2, P=Q+++R+.

## 2. Exact coefficient and parity certificate

For each original term, the auditor expands the above substitutions by a Cartesian product of sparse linear combinations and accumulates Python integers. No floating-point arithmetic enters the transformation comparison. Every resulting coefficient agrees with the candidate JSON, together with all core metadata and denominators.

Independently, the auditor checks that every source tensor polynomial is invariant under the original-basis involution

(R+,R−,Q++,Q+−) ↦ (R−,R+,−R+−R−−Q++,Q+−).

After the invertible basis change, this involution is diag(1,1,−1,−1). Thus a monomial containing an odd total number of O/P factors acquires a minus sign. Polynomial invariance in characteristic zero forces its coefficient to be zero. This proves the parity claim for all 31 audited tensors, rather than merely observing numerical cancellation after contraction.

Exact nonzero coefficient counts:

- Original four-basis tensors: 2364.
- Parity-adapted tensors: 553.
- Surviving odd-total-O/P monomials: 0.
- Surviving O/P loop entries: 0.

The last property is verified in the formal coefficient tables; independently, an O/P matrix loop always has zero diagonal.

## 3. Signed contraction quotient

`audit_quotient_palettes.py` independently enumerates the vertex permutations preserving each core's edge multiset. For each term it relabels vertices, introduces a minus sign when a skew edge reverses orientation, sorts propagators among parallel edges, and chooses a canonical index word. Terms with both canonical signs cancel. This signed orbit reduction exactly reproduces the candidate quotient after independent canonical normalization.

The 553 tensor terms reduce to **263 contraction monomials**. This is a contraction identity, not an assertion that the unsymmetrized tensor itself has only 263 coefficients. There is no extra automorphism division: equivalent terms are summed, while the original denominator is retained.

As a stronger numerical-free check, 93 contractions (31 cores at matrix orders 2,3,4) were tested using arbitrary small integer symmetric G,H and skew O,P matrices. Each parity-tensor numerator and quotient numerator agreed exactly. These tests do not impose rational-function relationships among the four matrices.

## 4. Palette support

Deleting an edge and testing graph connectivity gives an independent multigraph bridge test; loops and duplicate parallel edges are handled explicitly. In the 31 exact parity tensors:

- Loop-edge support is contained in {G,H}.
- Bridge-edge support is contained in {G,P}.
- Weight-one bridge support is contained in {G}.

Detailed supports for every edge are recorded in `quotient-palette-receipt.json`. These support constraints are certified for this finite inventory. They are not asserted here as an all-weight theorem.

## 5. Numerical engine comparisons

`audit.py` compares direct sparse contractions of the old and new JSON tables on four fixtures: a genuinely weighted skew order-3 matrix, a transitive order-4 tournament, an almost-regular order-6 tournament with scores 2,3,2,3,2,3, and a regular order-7 tournament. At z=1 and z=0.3+0.1i this gives 248 core comparisons. Maximum absolute difference: 1.3323 × 10⁻¹⁵.

`check_original_engine.py` separately imports the original five-basis optimized engine, with `reduced=False`. Deleted cyclic tournaments of orders 4 and 6, at z=1 and z=0.4+0.2i, give 124 further core comparisons. Maximum absolute difference: 4.4409 × 10⁻¹⁶.

## 6. Timing: do not confuse sparsity with overall speed

The exact term-count reduction is 2364/553 ≈ 4.275. This is not a measured predictor speedup.

Using the same generic sparse einsum evaluator on all 31 cores of one order-7 regular fixture, three repeats produced a median time ratio approximately 3.47 in favor of the parity representation. These are local empirical timings, recorded in `receipt.json`, including per-core path construction and excluding matrix construction. They show an advantage only against that identical, unoptimized evaluation strategy.

The producer's separate measurements in `session30-parity-basis/verification.json` show its naive sparse 263-term quotient slower than the optimized original four-basis engine at n=10,20,40. That is consistent with this audit: fewer symbolic terms do not necessarily improve the optimized contraction schedule. Those producer timings have not been independently reproduced here. A new compressed-palette implementation requires its own accuracy and performance checks; this audit does not certify it.

## Reproduction

Run from the workspace root:

    python3 session30-parity-audit/audit.py
    python3 session30-parity-audit/check_original_engine.py
    python3 session30-parity-audit/audit_quotient_palettes.py

Dependencies are Python, NumPy, and SymPy. JSON receipts are saved beside the scripts. `INPUT_SHA256SUMS` pins the data files read; `SHA256SUMS` covers this audit's deliverables. No Git changes or frozen-packet edits were performed.
