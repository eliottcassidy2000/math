# Independent review of TT3 degree tomography

**PASS, including the source-relative root-family integration in Section4. No mathematical correction is required.** This review binds the frozen proof and verifier pins in the paired receipts. It does not claim arbitrary-prime completeness or unrestricted reconstruction.

## 1. Exact three-inside-one-outside quartic

The eight-term quartic normal form specializes correctly. After eliminating the exterior, every H correction has at least four internal A variables and therefore vanishes when the full three-name support is extracted. The coefficient of C_p M_qr+C_q M_pr is the sum of the two arcs into r plus two copies of tau, namely d_A(r)+2tau. The opposite bracket gives 2-d_A(r)+2tau. Each triple product has coefficient two. Thus the displayed identity holds for all eight internal orientations, and changing the triple orientation changes the quartic by its indegree difference times K_i.

The virtual-singleton interpretation is exact. A first path ending at z contributes alpha before the distinguished endpoint, whereas a second path starting at z contributes beta after it. The squarefree full-support extraction enforces disjointness. Hence K_i=C_zi-C_iz, with ordered nonempty path slots and no product-of-marginals error.

## 2. Binomial positivity and equal-profile skew

The coefficient of x^j y^(N-j) in -(x-y)^2(x+y)^(N-2) is exactly binom(N,j)[N-(2j-N)^2]/[N(N-1)]. Here N>=2. For an interval-supported log-concave sequence, q_(c+j)q_(c+N-j) is symmetric and increases toward the center, including its zero boundary region. The same ordering holds for N-(2j-N)^2. Their covariance under binomial weights is nonnegative, and the latter factor has mean zero. This proves the claimed functional inequality. The r=0 value is exactly2q1^2.

Subtracting the two Schur formulas leaves antisymmetric multiples of x-y and (x-y)h plus the symmetric expression in equation(7). Equal complete profiles make the antisymmetric terms vanish. The remaining negative square is controlled term by term by the binomial inequality; all other terms are coefficientwise nonnegative. The r=0 square contribution is2q1^2[E](ab). The identity [E](ab)=H(E+z_M) follows by uniquely splitting a Hamilton path at z_M, with empty prefix or suffix allowed. It is strictly positive. The equal-profile and log-concavity hypotheses are genuinely used and are not removed.

The transitive-profile log-concavity argument is valid: adjacent convolution preserves interval-supported log-concavity by the stated expansion, and multiplication by k preserves it. The directed-triangle profile also satisfies the inequalities. No universal log-concavity assertion for arbitrary tournament profiles is made.

## 3. Independent resolution of the proposed family integration

Here m means the number of quotient names, |V(Q)|. The family is restricted to original modules inside the fixed core, with its prime root selected oddly.

If the core Hamilton count increases, the scalar physical query has all core bits zero and retains the root incidences outside. Its augmented core is x->y->core. Every ordered-cover coordinate of that ordered sum is a nonnegative combination of core cover coordinates, and the first coordinate increases strictly. The fixed ambient has positive first occupied-run coefficient, so the lifted scalar query separates. It need not be the globally all-zero port mask.

Otherwise hereditary equality preserves all child profiles. The target distinguished triple has Hamilton count one and is therefore transitive. Whole-converse invariance of every saturated quartic cover coordinate, followed by uniform restoration of the other child profiles, reduces the comparison to the source A versus the converse of the target A, with a common exterior. The grouped observer needs two distinct original a contacts and two distinct original c contacts. Bidegree at most(2,2) therefore excludes any additional alias route or contamination from fixed exterior bits.

The arbitrary ambient can be included in E: every added exterior vertex sees X union Y uniformly. The equal-profile positivity result gives K_X-K_Y>0, so these two numbers cannot both vanish. If all six quartic contrasts vanished, equation(2) would force all three named indegrees to agree between A and the converse of target A. A transitive triple is uniquely determined by its named indegrees. Therefore the target A is exactly A^op.

At that point, any query nonuniform only inside A and uniform on every other child restores the target exactly to the whole-core converse comparison. Partition the other m-1 children into their complete module blocks and leave the three A names singleton. There are m+2 blocks, including an uncovered singleton anchor. The inherited compressed anchor battery has at most m+4 queries and separates by2 modulo4. Its fixed-root ambient lift is valid. This does not apply a laminar-family theorem to crossing descendants: the prior exact profile restoration has already reduced the comparison to whole converse.

Consequently the stated conservative fixed battery bound96+1+(m+4)=m+101 is valid. This is source-relative detection, not pairwise reconstruction of arbitrary higher-profile outputs. The critical T/U/W all-TT3 family consequence follows from the already established cyclic ordered-pair witnesses and transitive-profile log-concavity. No arbitrary-prime closure follows.

These two clarifications, m=|Q| and the scalar query's fixed exterior bits, make explicit what the frozen proof uses; they do not change its hypotheses or conclusion.

## 4. Independent numerical controls

The auditor imports no producer implementation.

- Every one of the5,040 quartics is independently extracted by sixteen actual Hamilton queries, rather than the producer's path-matrix quartic formula. Its exterior expression is separately computed from literal simple-path words.
- Both complete1,024-query physical responses are independently recounted. Every one of their256 named coefficients, including all cubic/quartic terms and zeros, is recovered by direct finite differences and reconstructs every full query.
- All4,941 equal-profile skew records are independently recovered by a direct ordered-two-cover run automaton. It fixes the appropriate first-path endpoint and second-path start, rather than using Schur elimination. Its run bound n_E+3 is exact: each return to X or Y within a path requires an exterior or z occurrence, since X->Y cannot alone provide a return, and there are at most two path starts. There are549 complete exterior/incidence configurations.
- All3,456 binomial controls are independently reproduced by literal polynomial multiplication and exact profile contraction. Transitive profiles are constructed by the explicit surjection sum, not the producer's recurrence.
- Eighty complete root-containing original-module families are additionally checked, for TT3 with equal singleton or pair sides, both locally and in a cyclic fixed ambient. These include every internal transitive-triple order and crossing interval reversals. The fixed union of scalar, six quartic and compressed anchor queries separates every target. Actual query lists and full source/target values are retained.

In total49,787 distinct actual augmented graphs are independently counted. Both normal and optimized runs agree byte for byte. The physical engine uses a predecessor recurrence, distinct from the producer's forward recurrence. The other engine counts cut-constrained color words with ordered slots. Both use exact bounded integers; no floating arithmetic or installation is involved.

The finite controls support the all-order arguments and do not replace them. The unequal-profile boundary examples in the producer remain explicitly outside the positivity hypothesis; they are not treated as collisions.

Frozen audit allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `run_covers.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.
