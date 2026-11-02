# Session138 synthesis: sharp stability and actual-context barriers

## Result map and inheritance

Session137 proved the nonempty determinant bound, its weighted extension, and an all-order variance coefficient asymptotic to `24/n^2`. It also separated width-three local universality from exact realization and proved finite-base normal forms for bounded actual sources. Checkpoint one strengthened quantitative stability to the correct `1/n` order and proved exact constants through seven. The continuation proves sharp order eight and all-parameter sign chambers for three actual mixed-profile seeds. The independently audited all-order continuation sharpens the asymptotic constant to eight with error `O(n^(-2))`. Its exact finite-order bound and probability/matching consequences have separate source pins from checkpoint one.

Here **PROVED** denotes an argument for the entire stated domain. A finite certificate can support such a theorem when an explicit reduction covers the full continuous domain or every order; it is not a sampled inequality. **FINITE-EXACT** labels the complete finite classifications themselves. Audit and replay scopes are in [VERIFICATION.md](VERIFICATION.md).

## 1. Edge rounding retains both the binary gap and a variance correction

For a real weighted tournament A, let

```text
u=det(I+A),  v=(-1)^n det A,
delta=H-u,   Delta=H-u-|v|,
q_e=A_ij A_ji,   V=sum_e q_e.
```

For one uncertain edge e with probability p, its two endpoint orientations `A^0,A^1`, and pair-deleted remainder R, each signed margin obeys exact interpolation. Taking the smaller margin gives

```text
Delta(A) >= (1-p)Delta(A^0)+p Delta(A^1)
             +p(1-p)[u(R)-|v(R)|].
```

For `n>=3`, R is nonempty. Its strictly positive determinant margin permits successive rounding without losing the product-law expectation of the terminal binary tournament T. Order two is excluded: there `Delta=0` even with positive variance.

The uncertainty graph records exactly the nonbinary edges. For a graph J on r vertices, let `i_even(J),i_odd(J)` count independent vertex subsets by parity, including the empty set in the even count. Define

```text
mu(J)=2^(1-r) min(i_even(J),i_odd(J)),    r>=1.
```

A remainder whose uncertainty is contained in J has `u-|v|>=mu(J)`. The proof retains principal Pfaffian squares on independent sets: binary even-order tournament Pfaffians are odd integers, and a binary observer supplies the opposite parity count.

For any deterministic edge order, deleting the current edge's endpoints from the still-uncertain graph gives graphs `J_t` and the certificate

```text
Delta(A) >= E Delta(T)+sum_t q_(e_t) mu(J_t).
```

The order may depend on original probabilities, but not on future rounding outcomes. Given a vertex cover of size c, the parity counts have an exact `O(2^c poly(n))` evaluator by enumerating independent subsets inside that cover. This is a fixed-parameter certificate algorithm, not a polynomial-time general independent-set or Hamiltonian counter. The graph certificate is orientation-independent; its optimum is not claimed to equal the optimal tournament inequality.

[Uncertainty-graph theorem](rounding_graph/PROOF.md)

## 2. Checkpoint-one multistar rounding gives the correct all-order scale

If the uncertain edges of a nonempty r-vertex remainder have a vertex cover of size c, the Pfaffian argument gives `u-|v|>=2^(-c)`. The zero-size binary-core case uses a separate uniform determinant bound, not an invalid parity-subset count.

Choose k ordered centers, first round all edges outside them, then round their stars in order. Average the ordered centers uniformly. The correction coefficient of every edge is the same, yielding

```text
Delta(A) >= E Delta(T)+gamma_(n,k)V,        n>=3, 1<=k<=n,
gamma_(n,k)=4[n-k+2-(n+k+2)/2^k]/[n(n-1)].
```

At k=1 this is `2/n`; optimizing k near `log_2 n` gives `max_k gamma_(n,k)~4/n`. A convenient explicit lower coefficient for `n>=4` is

```text
4[n-ceil(log_2 n)]/[n(n-1)].
```

Fair-pair examples show that the optimal unrestricted coefficient is at most asymptotic to `8/n`. This first-checkpoint lower bound settled the dimensional order and left a factor-two asymptotic gap relative to the fair-pair upper family. Its scope is preserved; the separate continuation in Section 8 closes that asymptotic gap. A dimension-free positive coefficient is impossible, and exact finite-order optimality is a separate question.

Any separately proved nonnegative binary lower bound can also be transferred through `E Delta(T)`; this does not assume an unproved weighted recurrence.

[Multistar proof](stability/MULTISTAR_STABILITY.md)

### The full matching-uncertainty class

Suppose the uncertain edges form a matching of r edges. All remaining orientations may be arbitrary binary ones. The independence polynomial is `(1+2z)^r(1+z)^t`; random edge order gives

```text
Delta(A) >= E Delta(T)+a_r V,       n>2r,
a_r=(4/r)[1-(3/4)^r],

Delta(A) >= E Delta(T)+b_r V,       n=2r>=4,
b_r=(4/r)[1-(3/4)^r-(1/4)^r].
```

These sharp matching-support coefficients hold beyond the H=1 locus. For n>=4, positive-variance equality in the unreserved matching-support inequality forces H=1. With at least two uncertain pairs they must all be fair; with one pair any fractional probability is allowed under the H=1 condition. Order three remains exceptional because Delta=V holds for every weighted matrix. No equality classification is asserted for the stronger inequality retaining the binary baseline. None of these matching-support conclusions is exported to arbitrary uncertainty graphs.

## 3. H=1 has a complete structure and exact extremal ratio

Since `H(A)=E H(T)` and every binary tournament has `H>=1`, equality `H(A)=1` means every supported binary completion is transitive. Two fractional edges sharing a vertex would permit a directed triangle. A fractional pair seen nonuniformly by another vertex would do the same. Therefore

> H=1 exactly for ordinal sums of genuinely fractional pairs and singletons.

For pair variances `q_1,...,q_r`,

```text
H=1,  u=product_i(1-q_i),  V=sum_i q_i,
|v|=product_i q_i if there are no singletons, and 0 otherwise.
```

For at least two positive pair variances, `Delta/V` strictly decreases as each variance increases, so the minimum requires fair pairs. The sharp H=1 constants are b_r at even order `2r>=4` and a_r at odd order `2r+1>=3`. Positive-variance equality at order at least four requires the maximum possible number of fair pairs, with one singleton exactly at odd order. At order three the ratio is identically one on this locus.

These examples supply global upper bounds on the unknown optimal stability coefficient; the structure theorem alone is not a global lower-bound proof.

[H=1 theorem](stability/H1_STABILITY.md)

## 4. Sharp unrestricted constants through order eight

The exact coefficients in `Delta>=bV` are

```text
n:  3    4    5    6     7      8
b:  1   3/4  7/8  3/4  37/48  87/128.
```

The first two are inherited. Checkpoint one established orders five through seven; the separate continuation establishes order eight. For positive variance at orders four through eight, equality is exactly the ordinal fair-pair construction above. Exact equality with the fair-pair candidate at higher orders remains open.

Let tau be the sum of directed-triangle weights, C_j the sum of directed j-cycle weights, and M33 the sum over unordered disjoint triangle pairs. Cycle rotations and collection multiplicities follow the source definitions.

### Triangle reserves and exact joint regions

Retaining different path lengths in the root recurrences gives the inherited-coefficient refinements

```text
delta >= tau+c_n V,    n>=2;
Delta >= tau+beta_n V, n>=4.
```

The strong domain matters: binary C3 has `Delta=0`, `tau=1`. In particular a nearest binary rounding with a directed triangle forces `Delta>=1/8` for `n>=4`.

At order five the exact set of real coefficient pairs for `Delta>=a*tau+b*V` is

```text
a<=1,   b<=7/8,   a/4+b<=17/16.
```

One corner is analytic; the other has a complete independently reconstructed Bernstein certificate. Their convex hull and downward closure give sufficiency; explicit tournament families force all three facets. The slanted facet has a continuous equality family, not only two isolated witnesses.

At order six an analytic signed-cycle expansion proves `Delta>=tau+3V/4`. Together with the extremal witnesses, the exact coefficient region is the rectangle `a<=1,b<=3/4`.

[Low-order proof and equality cases](stability/SHARP_LOW_ORDER_STABILITY.md)

### Graph certificates with a full-cube reduction

The local triangle inequality gives

```text
tau >= W(q):=sum_(incident unordered edge pairs e,f) q_e q_f
    >= 2V^2/n-V/4.
```

Its disjoint-edge lift controls `P21(q)`, the weighted sum of a two-edge wedge and a vertex-disjoint edge. Separately affine residual polynomials attain their minimum over `q_e in [0,1/4]` at graph endpoints. Thus complete undirected graph certificates, combined with the analytic cycle reductions, prove

```text
n=6: delta >= tau+(37/48)V+(3/4)C4+C5+C6+3M33;
n=7: Delta >= (37/48)V+(7/24)tau+C4/4+C5+C6+3M33.
```

The complete graph universes have 32,768 and 2,097,152 labeled graphs respectively, with no filters. Independent edge-toggle counting checks the same universes using a different kernel. Positive reserves force tau=0 in the unreserved equality cases, returning to the H=1 classification.

The wedge-plus-disjoint-edge term is essential to this reduction: omitting it admits the graph `C4 disjoint-union K3` with negative residual -41. That graph is a counterexample to a proposed proof reduction, not to the weighted tournament inequality.

[Triangle variance](stability/TRIANGLE_VARIANCE.md) · [Weak six](stability/SHARP_WEAK_SIX.md) · [Strong seven](stability/SHARP_SEVEN_GRAPH_CERTIFICATE.md)

### The order-eight continuation

The fully expanded cycle-cover ledger groups potentially negative four-cycle terms using complementary four-vertex determinants; the factor two on pairs of disjoint four-cycles is essential. Other complementary-cycle estimates leave

```text
Delta >= (87/128)V+tau/16+C4/8+C5+C6/2+C7+3M33.
```

For the remaining variance polynomial, endpoints reduce to the graph inequality

```text
41m-32m2+8m3-4m4+8p21 >= 0,
```

where m_j counts j-edge matchings and p21 counts a wedge plus a disjoint edge. This is the stronger certificate **without** an extra wedge term. Complete coverage comes from 1,044 seven-vertex graph orbits and all 128 eighth-vertex neighborhoods, giving 133,632 augmentations. An independent orbit/vertex-subset counting route verifies the same full scope. Multiaffine interpolation supplies the continuous-cube conclusion.

Four fair ordinal pairs have `H=1`, `u=81/256`, `|v|=1/256`, `V=1`, hence `Delta=87/128`. The positive triangle reserve forces H=1 at unreserved equality, and the H=1 classification forces exactly those four fair pairs.

The alternating finite-order coefficients are not fitted numerology. At odd order a singleton kills the full determinant, giving `(4/r)[1-(3/4)^r]`; at even order a full pair cover subtracts `(4/r)(1/4)^r`. This explains the fair-pair candidate's parity dependence without proving its optimality at unreviewed orders.

[Order-eight proof](stability/SHARP_EIGHT_GRAPH_CERTIFICATE.md)

## 5. One actual exterior edge reflects the module's Hamilton count

For a nonempty tournament U, write `h=H(U)` and let E count ordered two-path covers allowing empty paths. Two explicit strong one-port contexts, differing in one exterior edge, satisfy

```text
H(R_minus(U))=6E-h,
H(R_plus(U)) =6E+h.
```

The module is unchanged. Its scalar gap is exactly `2H(U)`, and both resulting tournaments are strong even when U is reducible. The complete path-cover polynomial is retained: for
`F(z)=sum_c p_U(c)(z)_c`, the operator difference is

```text
F_Rplus(z)-F_Rminus(z)
 = z[F(z+1)-F(z-1)]
   -z^2[F(z+1)+F(z-1)-2F(z)].
```

These signed polynomial formulas evaluate actual graph operations; they do not authorize subtracting tournaments or accepting unrealizable states. For a transitive four-vertex U, H increases by two while the unordered two-path-cover count decreases by twenty-two. Scalar positivity therefore does not imply profile monotonicity.

Ordinal padding multiplies E by powers of two while preserving h. Every attained half-gap h consequently has unbounded strong one-edge-flip realizations. The transitive special case gives the neighboring odd pair

```text
6*2^k-1,  6*2^k+1,   k>=1.
```

This is an infinite dyadic family, not an interval theorem or an induced extension. Four exterior vertices are the minimum for a nontrivial unit isolator. At that minimum, all 48 oriented low-to-high one-edge isolators have the same two coefficient rows and form two reversal-related relabeling orbits.

A source-certified inverse enumerates actual nonsingleton SCC cores and their E values. It never accepts a formal `(h,E)` pair. The complete reflector image through 200 has sixteen values, each with an actual witness; omitted values are omitted from this constructor only.

[Reflector, full operators, inverse, and classification](exact_resources/PROOF.md)

## 6. Bounded actual resources fail even local universality

Session137's finite-base theorem puts every bounded height/interior source bank inside finitely many exceptional points and affine-dyadic rays. For J points and K rays,

```text
|outputs mod (2^r-1)| <= J+Kr.
```

This is eventually less than `2^r-1`, so the bank misses a residue at an odd modulus. It is not even profinitely universal among odd integers.

More explicitly, for height at most H and interior at most I, let K(H) be the inherited exact nonsingleton-core bound. A complete overcount of bases has parent order at most

```text
N(H,0)=K(H)+1,
N(H,I)=K(H)+3+floor(log_2 I),     I>=1,
```

and size at most `C_N=sum_(n=1)^N 2^(n(n+1)/2)`. Any r with `C_N(r+1)<2^r-1` certifies an omitted residue. Universality at these moduli requires `N>=(sqrt(2)-o(1))*sqrt(r)`.

This separates source resources from substitution width: the width-three all-modulus controller must leave every bounded `(h,L)` bank, although its quotient width remains three. No new unrestricted spectral omission or sufficiency theorem for growing resources follows.

[Explicit budget theorem](exact_resources/RESOURCE_BUDGET.md)

## 7. Actual kernel identities and restricted zero-rigidity

For a homogeneous pair adjoined to an arbitrary tournament D, compare cyclic and middle completions with outside cuts r,s. With the actual ear polynomial F_D,

```text
H(cyclic)-H(middle)
 = 2F_D(r)+(s-r) dot grad F_D(r)
 >= 2H(D)>0.
```

This includes empty D. It gives exact twin-capacity differences; a homogeneous triangle also has internal capacities `2H` after deleting one triangle vertex. An exact quotient run-word compiler reconstructs all remaining capacities from first derivatives of the actual transition polynomial. These are identities between actual kernels, not free positivity conditions.

A strengthened formal-network perturbation preserves the twin identities, internal-triangle identities, positivity, divisibility, and layer means, yet still creates additional zeros. It violates the cross-block derivative tensor of the actual parent. It is not a tournament counterexample; it identifies information a successful general proof must retain.

### Ordinal modules: all zeros and full-profile comparison

For `T=TT_k[C3,...,C3]`, every scalar H-restoring completion of a deleted-triangle card is canonical or has one specified adjacent-block defect and the middle-pair flip. The labeled restoring-cut count is one for k=1, four at an endpoint root for k>=2, and seven at an interior root.

More strongly, every completion U is entrywise comparable with T in its ordered path-cover profile:

- canonical U has the identical profile;
- a lower profile has strictly smaller H;
- an upper noncanonical profile has larger H, or equal H and strictly larger second coordinate.

The all-order reduction uses fifteen explicitly listed cores, whose whole profiles have independent literal-permutation checks. It is not extrapolation from small k.

Therefore `C3[T,A,B]`, for arbitrary nonempty A,B, is zero-rigid under star edits inside T **when all orientations to A and B remain fixed**. A two-run context coefficient detects the noncanonical scalar zeros. Changing those outside attachments destroys the hypothesis and remains open. An actual strong-quotient hostile has decreasing H but increasing higher profile coordinates, so the ordinal profile dichotomy cannot be exported to arbitrary strong quotients.

[Actual kernels](zero_rigidity/README.md) · [Ordinal zero classification](zero_rigidity/ORDINAL_ZERO_CLASSIFICATION.md) · [Profile dichotomy and restricted strong theorem](zero_rigidity/PROFILE_DICHOTOMY.md)

## 8. Asymptotic closure: constant eight and a quadratic error

**PROVED, independently audited.** This continuation strengthens, rather than silently rewriting, the first-checkpoint result.

Let `K_n=inf_(V>0) Delta/V`, and let kappa_n be the sharp H=1 fair-pair coefficients from Section 3. For every `n>=4`, the theorem is

```text
kappa_n/[1+n*kappa_n^2/(2*c0)] <= K_n <= kappa_n,
c0=(1-exp(-4))/4,
K_n=8/n+O(n^(-2)).
```

The upper family is already actual. The lower bound combines three mechanisms:

1. **Multiplicative Pfaffian margins.** Varying one skew edge in either signed principal-minor sum gives a quadratic whose completed-square residual dominates its leading coefficient. Iterated endpoint interpolation yields
   `u-|v|>=product_e(1+4q_e)^(-1)>=exp(-4V)` for every nonempty remainder.
2. **Nonlinear rounding.** Sum the one-edge corrections in a fixed order; a decreasing-exponential integral gives
   `Delta>=B+(1-exp(-4V))/4`, where `B=E Delta(T)`.
3. **Extract a weighted matching but retain B.** A greedy maximum-weight matching of total variance M satisfies `W>=(V-M)^2/(2n)`. The first-checkpoint matching-support theorem and binary proper-cycle bound `B>=tau>=W` then give
   `Delta>=kappa_n V+B-kappa_n sqrt(2nB)>=kappa_n V-n*kappa_n^2/2`.

The last additive error, combined with the nonlinear bound, gives the positive finite-order ratio bound above, uniformly over all variance scales. Discarding B too early loses the improved rate. The statement is not exact finite-n equality `K_n=kappa_n`, and the lower bound is not advertised as sharp at small orders.

### Near-candidate structure is quantitative, not exact rigidity

For every matrix with `Delta/V<=kappa_n`, the proved consequences are

```text
B <= 2n*kappa_n^2,       W<=tau<=B,
V-M <= 2n*kappa_n,
Pr(T nontransitive) <= 2n*kappa_n^2.
```

Thus the binary expected gap is `O(1/n)`, and only `O(1)` **variance mass**, not necessarily O(1) edges, lies outside the greedy matching. Random independent rounding is transitive with probability `1-O(1/n)`. These bounds do not force the original cross edges to be binary, every pair to be fair, or every supported completion to be transitive at a fixed finite n.

A second analytic mechanism uses resolvent norms rather than Pfaffian coefficient comparison. It proves `u-|v|>=exp(-12V)` and `Delta>=B+(1-exp(-12V))/12`. These weaker constants still give the same asymptotic rate through baseline-preserving matching extraction; the displayed sharper finite-order bound uses the Pfaffian constant four.

[Asymptotic theorem](asymptotic_stability/PROOF.md) · [Independent audit](independent/asymptotic_audit/AUDIT.md) · [Resolvent cross-check](stability/RESOLVENT_ROBUSTNESS_CROSSCHECK.md)

## 9. Exact sign chambers for actual mixed-profile seeds

Take the three fixed twelve-vertex star completions M_867, M_995, M_996 of the strong module `M=Q4[C3,C3,C3,C3]`. Let `A_a=TT_a[C3]` for positive integers a, and measure the change inside the actual cyclic context `C3[M,A_a,A_b]`, retaining the original outside attachments.

The all-parameter classification is:

| Seed | Thin boundary `min(a,b)=1` | Thick quadrant `a,b>=2` |
|---|---|---|
| 867, 995 | Strictly positive | Strictly negative |
| 996 | Strictly positive | Strictly positive |

None gives zero at any positive integer pair. A sign change from multiplicity one to two is not an intermediate-value argument on an integer parameter set.

The 996 result is stronger: for arbitrary nonempty side tournaments A,B,

```text
H(C3[M_996,A,B])-H(C3[M,A,B]) >= 135618 H(A)H(B)>0.
```

Singleton sides attain the coefficient. For 867/995, thin positivity also allows one side C3 and the other arbitrary nonempty. The negative chamber extends to actual sides with colored profile polynomials `P_A=p^2 R_A`, `P_B=p^2 R_B`, where `p(m)=m^3+2m` and both R have nonzero nonnegative monomial coefficients. This includes ordinal sums containing two triangle components and arbitrary additional tournament components; it is sufficient, not a necessary classification of all sides.

### Why this is all-parameter, not a scan

For ordered profile difference d, exact cyclic run-word contraction gives a tridiagonal kernel

```text
D(A,B)=f_A^T K f_B,
K_(r,r)=d_(r-1)+3d_r+d_(r+1),
K_(r,r+1)=d_r+d_(r+1).
```

For seed 996 this matrix is nonnegative and its first diagonal entry is 135618. For the other two seeds, its rows have a single positive-to-negative transition. The base two-triangle profile is positive and log-concave. Multiplication by a nonnegative monomial polynomial makes its binomial-profile likelihood ratios nondecreasing; a one-crossing comparison propagates the exact negative base response in both side variables. The local identities and row signs, not finite parameter samples, carry the proof.

[All-parameter sign theorem](zero_balance_later/NO_ZERO_THEOREM.md)

## What is not concluded

The unrestricted H-spectrum, global boundedness of exact spectral width, arbitrary strong-quotient scalar zero-rigidity, and LRC(14) remain open. Exact equality `K_n=kappa_n` beyond eight remains open. The new no-zero theorem closes one explicit unbounded construction family, not every mixed-profile seed or every actual side context. The asymptotic constant is now proved; this is not an exact finite-order equality theorem beyond eight.
