# A000568 prediction across enormous scales

Status: **PROVED finite identities and symbolic model bounds; numerical forecasts
are arbitrary-precision evaluations, not certified interval arithmetic.** This
lane deliberately does not read or compare exact A000568 caches. No priority
claim, all-integer-digit certification, or Chan–Pak conclusion is made.

## 0. Inheritance and research board

The inherited input is the odd-cycle Burnside formula, documented in
`05-knowledge/results/session270_chanpak_enumeration_20261111/enumeration/PROOF.md`,
and its fixed-mass response quotient. The **defect-polynomial expansion,
degree3k, unique all-3-cycle leading coefficient, and fixed-order asymptotic
hierarchy are inherited**, not new in this session, from THM-4845:
`01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md`,
Sections D and proof D. That source explicitly credits Session98
`05-knowledge/results/session98_tournament_growth_20261029/limits/ALL_POLES_CORE_EXPANSION.md`.
The Session98 path is recorded as cited by THM-4845; it is not present in this
worktree, so we do not claim to have independently read that file. Our
falling-factorial basis is equivalent to the binomial basis in THM-4845.
This lane reorganizes those inherited results for executable giant-index
forecasting, adds the explicit support-endpoint envelope and adaptive budget,
and records the quantitative comparison to a naive Poisson gas. No external
priority is asserted for these practical derivations or comparisons. Here we change the
question from exact enumeration to prediction with quantified relative error.
The inherited hostile warning (Session199, repeated in Session270 synthesis)
is that a tiny relative tail does **not** justify rounding to an exact integer.
We preserve it explicitly. This lane did not locate the original Session199
file and attributes the warning through Session270 rather than inventing a link.
The least-used sidecar is the moved support of a permutation, which controls
exponent loss even without solving the response inventory.

Anchor: accurate giant-index predictions. Niche: finite defect polynomials and
interacting rare symmetries. Wildcard: structural rarity versus numerical
census accuracy. The live board has six entries:

1. Identity entropy, with factorial normalization retained.
2. Odd-cycle support, with exact defect retained.
3. Finite defect polynomials, with gcd interactions retained.
4. Stirling expansion, with its first omitted term retained.
5. Precision/scale separation, with model and numeric errors distinct.
6. Structural observers, with labelled/class sampling retained.

The meta-pattern used is “Correct the object before sharpening the technique”:
replace exact integer recovery by logarithmic/relative prediction. Its
counterindication is a request for low bits, exact rounding, or full digits.

## 1. Exact positive correction expansion

Let a(n) count ordinary tournament isomorphism classes, N=binom(n,2), and
b(n)=2^N/n!. Burnside implies a(n)=b(n)(1+δ_n), δ_n>=0. A permutation with any
even cycle fixes no tournament, because an appropriate power reverses an edge.
For a multiset μ of nonunit odd cycles, write

- m=sum of its cycle lengths (moved support), r=number of cycles;
- z_μ=product over odd j>=3 of j^(c_j)c_j!;
- e_μ=sum over cycles (j−1)/2 + sum over unordered pairs of cycles gcd(i,j).

There are (n)_m/z_μ permutations of this type, with n−m fixed vertices.
Their pair-orbit count is binom(n−m,2)+(n−m)r+e_μ. Thus their normalized
contribution is exactly

    w_μ(n) = (n)_m/z_μ · 2^(−D_μ(n)),
    D_μ(n) = (n−m)(m−r)+binom(m,2)−e_μ,
    δ_n = sum over μ with m<=n of w_μ(n).

This is not a new Burnside formula. It is a useful normalization and support
organization of the inherited formula. All coefficients are positive, so
retained terms give a lower bound without cancellation.

Examples, with (n)_m=n(n−1)…(n−m+1):

    δ_3 = (n)_3/3 · 2^(−2n+4),
    δ_5 = (n)_5/5 · 2^(−4n+12),
    δ_33 = (n)_6/18 · 2^(−4n+14).

Consequently δ_33/δ_5=10(n−5)/9 for n>=6: arranging correction types merely
by maximum cycle length misses the eventually larger two-3-cycle term.

## 2. Rigorous constant-cost relative envelope

For any nonidentity odd-cycle permutation of support m, r<=m/3. Each fixed
vertex contributes m−r>=2m/3 to the deficit. Among pairs of moved vertices,
no pair is fixed: fixed endpoints would contradict movedness, and a swapped
pair would require a 2-cycle. Every pair orbit has odd length, hence at least3.
Thus at most binom(m,2)/3 moved-pair orbits remain and

    D >= (2m/3)(n−m)+(2/3)binom(m,2)
      = m(2n−m−1)/3.

There are at most binom(n,m)m! <= n^m permutations with this support. Therefore
for any integer M with 3<=M<=n, the entire support tail is at most

    R_M(n) <= sum_(m=M)^n exp(f_n(m)),
    f_n(m)=m log(n)−log(2)m(2n−m−1)/3.

Since f_n is convex in m, its maximum on [M,n] is at an endpoint. This yields

    R_M(n) <= U_M(n)
      := (n−M+1) max{ n^M 2^(−M(2n−M−1)/3),
                       n^n 2^(−n(n−1)/3) }.

For M>n the tail is zero. The bound is valid, even if loose or greater than1,
for every indicated n. No monotonicity assumption is hidden in the endpoint
step. Integer ceilings would sharpen it but are unnecessary.

Let C_s sum every w_μ with m<=s. Then

    b(n)(1+C_s) <= a(n) <= b(n)(1+C_s+U_(s+1)).

If s<3 use M=3. In particular the relative error of b(n)(1+C_s), measured
against the actual a(n), is <=U_(s+1)/(1+C_s); the implementation conservatively
reports U_(s+1). The relative error of b(n) is <=C_s+U_(s+1). For logarithms,

    0 <= log10 a(n)−log10 b(n)
       <= (C_s+U_(s+1))/log(10).

The bound is symbolic and rigorous. The decimal strings evaluating the bound
are numerical approximations, not outward-rounded certified endpoints.
For strict machine certification, evaluate the same expression in a directed
interval package with sufficient precision.

### Useful asymptotic consequences

For fixed M, the M endpoint dominates eventually. Hence U_M is
O(n^(M+1)2^(−2Mn/3)). Retaining supports3,5,6 gives

    δ_n = δ_3+δ_5+δ_33+O(n^8 2^(−14n/3)).

The finite leading correction is exponentially more accurate than the identity
term, and the second layer is correctly organized by defect, not cycle size.
All asymptotic statements here are fixed truncation order as n tends to infinity.
Uniform growing-order statements require the displayed finite envelope instead.

## 3. Inherited defect transseries and a non-Poisson gas comparison

Every d=m−r is even. Put d=2k. Since every nonunit cycle has j−1>=2,
r<=k and m<=3k; also m>=2k+1. Grouping the exact sum gives

    δ_n = sum_(k=1)^floor((n−1)/2) P_k(n) 2^(−2kn),
    P_k(n) = sum_(μ: m−r=2k) [(n)_m/z_μ]
               2^[m(m−r)−binom(m,2)+e_μ].

Here a falling factorial with m>n vanishes for nonnegative integer n. Each
P_k is a polynomial in the falling-factorial basis with positive rational
coefficients, of degree3k. The finite inventory is generated by partitions of
k: a part h corresponds to a cycle of length2h+1. Its sole degree3k term is
all 3-cycles, with coefficient

    [ (n)_(3k) ] P_k = 2^(3k²+k)/(3^k k!).

For example

    P_1=(16/3)(n)_3,
    P_2=(4096/5)(n)_5+(8192/9)(n)_6,
    P_3=(16777216/7)(n)_7+(16777216/15)(n)_8
         +(536870912/81)(n)_9.

`defect_coefficients.json` supplies k=1,…,6 (defects2k=2,…,12) and each contributing cycle type.
This classifies every finite defect layer, not every tournament structure.
These polynomials are exact reorganizations, not a convergent infinite power
series at arbitrary real n. For integer n the original sum is finite.

**Inherited fixed-order hierarchy, with a self-contained weaker envelope proof.**
THM-4845 already supplies the sharper remainder order
O_K(n^(3K+3)2^(−2(K+1)n)); the support argument below gives a simpler,
slightly weaker estimate and does not claim a new hierarchy. For every fixed
K>=1,

    δ_n = sum_(k=1)^K P_k(n)2^(−2kn)
            +O_K(n^(3K+4)2^(−2(K+1)n)).

To prove this, divide omitted types at support M=3K+3. The finitely many types
below M have defect d>=2K+2, hence each has weight a constant times (n)_m
2^(−dn); their sum is O_K(n^(3K+2)2^(−2(K+1)n)). The support tail is at most
U_M=O_K(n^(M+1)2^(−2Mn/3)), which has the claimed order. Retained types all
have support<=3K, so none is mistakenly subtracted from the tail. In particular
each omitted remainder is exponentially smaller than the last retained layer.
The constants depend on K; this is not a uniform theorem for K growing with n.

In the all-3-cycle sector, relative to the naive independent-Poisson term
δ_3^k/k!, the exact factor is

    (n)_(3k)/( (n)_3 )^k · 64^binom(k,2).

The first factor is finite-population depletion. The second is a genuine
pairwise enhancement: each pair of 3-cycles saves six deficit bits. The
independent-Poisson picture gets the first correction but not the higher ones.
For fixed k the depletion tends to1, so the factor64 survives at all large n.
This is a precise interacting cycle gas, not an informal independence claim.
The original e_μ retains gcd interactions for arbitrary cycle lengths.

## 4. Four prediction methods, plus scale and observer controls

**Identity/log-gamma.** Compute

    L0=binom(n,2)log10(2)−logGamma(n+1)/log(10).

This already gives excellent relative accuracy beyond250, certified at the
model level by Section2. It is the preferred giant-index forecast.

**Stirling hierarchy.** For integer n>=1 and p>=0,

    log(n!)=(n+1/2)log n−n+log(2π)/2
       +sum_(k=1)^p B_(2k)/(2k(2k−1)n^(2k−1)) + E_p,

where E_p has the sign of the first omitted term and its absolute value is
at most that term's magnitude. For completeness, Binet’s integral gives the correction to the first three
terms as 2∫_0^∞ arctan(t/n)/(exp(2πt)−1) dt. Expand arctan(x) by integrating
the finite geometric expansion of 1/(1+x²). The remainder has the next sign
and absolute value at most x^(2p+1)/(2p+1) for every x>=0. Integrating against
the positive kernel gives the displayed Bernoulli coefficients and bounds
the remainder by the next coefficient. We use this classical identity and
bound, not a novel literature result. At p=4,

    |E_4| <= 1/(1188 n^9).

Divide by log10's denominator log(10) for the log10-count error. This is
algebraic error, vastly larger than the exponentially tiny symmetry correction
at huge n, but still tiny for ordinary significant-digit prediction. The
script compares Stirling and logGamma numerically as two evaluation routes.

**Finite-support improvement.** Retain every type through support12, and
record the omitted-support envelope. Evaluate the tiny correction separately
in logarithmic form, rather than adding it to1 at inadequate precision.

**Defect-polynomial improvement.** For precision on the correction itself,
use the grouped P_k. This reveals all-3 dominance inside each fixed layer and
which mixed gcd types enter next. A support envelope remains available by
noting every omitted defect after k=K has support at least2K+3, so U_(2K+3)
is a valid (possibly double-counting and loose) upper bound for that remainder.

**Adaptive stopping.** Scan s=2,…,24 and stop when log U_(s+1) is below the
requested relative tolerance. Return failure honestly if the cap is exhausted.
For the requested50 relative decimal places, the identity term alone passes
at every default index. At n=30 this conservative procedure does not pass
within the cap; `validation.json` preserves that negative control. This is a
bound-dependent stopping rule, not an estimate that silently assumes accuracy.

**Local ratio and entropy scaling.** Exactly,

    b(n+1)/b(n)=2^n/(n+1),
    a(n+1)/a(n)=[2^n/(n+1)](1+δ_(n+1))/(1+δ_n).

The relative ratio error against its identity predictor is at most
δ_n+δ_(n+1), so the same envelope controls local growth. The log entropy has

    log2 a(n)=n(n−1)/2−n log2 n+(log2 e)n
              −(1/2)log2(2πn)−1/(12n ln2)
              +O(n^(−3))+log2(1+δ_n).
The dominant normalized limit is log2 a(n)/n²→1/2. The logarithmic output size
is O(log n) digits in the exponent, while printing a(n) requires Θ(n²) digits.
In precision chosen to resolve its mantissa, we retain about2log10(n) guard-
independent decimal places in the log arithmetic. We do not claim constant
bit complexity for arbitrarily long input indices.

## 5. What the bound does and does not classify

For the uniform **labelled** random tournament,

    E(|Aut(T)|−1)=δ_n.

Proof: sum over nonidentity permutations the probability they fix T; this is
exactly the normalized Burnside correction. Thus the probability of any
nontrivial automorphism is at most δ_n. It is not generally equal to δ_n.

A very different structural event is a two-vertex module (the two vertices
have identical relations to every outside vertex). Each fixed pair is a
module with probability2^(−(n−2)), so its expected number is

    binom(n,2)2^(−(n−2)).

This is an exact expectation, not an independence assertion or an ordinary-
class count. It is exponentially larger than δ_n, of order n²2^(−n) versus
n³2^(−2n). Indeed the two vertices cannot be exchanged by a tournament
automorphism because their internal edge has a direction. Therefore
“almost all tournaments have trivial symmetry” does not replace a recursive
module inventory, and astronomically accurate census forecasts do not
classify rare recursive structures. The support/gcd coordinate is about
permutation symmetry; modular decomposition needs separate quotient and
substitution data. The parent synthesis supplies that separate lane.

## 6. Numerical interpretation and reproducibility

`predict.py` uses mpmath1.3.0 and decimal precision2*digits(n)+90. It repeats
the identity-log computation at double precision and records the discrepancy.
This is a numerical stability diagnostic, not a proof of rounding correctness.
Printed mantissas have50 significant digits; logs have separate formatting
budgets. Their final digits include ordinary output rounding. Neither the
model envelope nor precision doubling establishes exact leading integer digits
when an interval could cross a decimal boundary. No exact digit counts are
asserted for a(n); the exponent is that of the numerical identity predictor.

The symbolic Stirling bound is separate from symmetry/model error and both
are separate from floating-point evaluation. For n=10^1000, a60-significant-
digit rendering of a correction logarithm is a coarse *absolute* rendering,
although highly precise relatively; use the executable at a larger formatting
budget if absolute accuracy in that correction's giant exponent is wanted.

Commands from repository root:

    python3 05-knowledge/results/session271_predictive_tournaments_20261111/forecast/predict.py --output 05-knowledge/results/session271_predictive_tournaments_20261111/forecast/forecasts.jsonl
    python3 05-knowledge/results/session271_predictive_tournaments_20261111/forecast/validate.py

Validation counts pair orbits directly for159 cycle-type/fixed-point controls,
checks the deficit inequality and all-3 leading coefficients through k=6,
and compares the independent Stirling/logGamma routes. It never computes
an exact complete A000568 count or queries the inherited cache.
