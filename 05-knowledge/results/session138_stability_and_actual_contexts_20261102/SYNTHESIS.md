# Session138 synthesis: exact stability and actual-context barriers

## Result map and inheritance

Session137 proved the nonempty determinant bound, its weighted extension, and an all-order variance coefficient asymptotic to `24/n^2`. It also separated width-three local universality from exact realization and proved finite-base normal forms for bounded actual sources. Checkpoint one strengthened quantitative stability to the correct `1/n` order and proved exact constants through seven. The continuation proves sharp order eight and all-parameter sign chambers for three actual mixed-profile seeds. The independently audited all-order continuation sharpens the asymptotic constant to eight with error `O(n^(-2))`. Its exact finite-order bound and probability/matching consequences have separate source pins from checkpoint one. A third independently audited continuation proves exact fair-pair optimality and the complete positive-variance equality classification for every n>=1500. The final finite-order continuation adds sharp order nine. The established exact domain is n=3,...,9 and n>=1500; orders 10 through 1499 remain unresolved.

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

Any separately proved nonnegative binary lower bound can also be transferred through `E Delta(T)`; this does not assume an unproved weighted recurrence. For cycle weights, this applies to directed simple cycles of length at least three (even cycles therefore have length at least four). Weighted two-cycle terms `q_e=p_e(1-p_e)` are variance corrections, not independent binary cycle events: the two opposite arcs of one unordered edge are correlated.

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

These examples supply global upper bounds on the optimal stability coefficient; the structure theorem alone is not a global lower-bound proof.

[H=1 theorem](stability/H1_STABILITY.md)

## 4. Earlier finite-order certificates: sharp constants through eight

The exact coefficients in `Delta>=bV` are

```text
n:  3    4    5    6     7      8
b:  1   3/4  7/8  3/4  37/48  87/128.
```

The first two are inherited. Checkpoint one established orders five through seven; the separate continuation establishes order eight. For positive variance at orders four through eight, equality is exactly the ordinal fair-pair construction above. These first two checkpoints cover orders through eight. Section 13 adds order nine; the separate argument in Section 10 establishes all orders n>=1500. Orders 10 through 1499 remain open.

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

The uniform formula is `kappa_n=(4/r)[1-(3/4)^r-1_(n even)*(1/4)^r]`, with `r=floor(n/2)`. The parity term is the full-matching determinant: an odd singleton kills det A. This is structural, not fitted numerology, and is not by itself a global optimality proof.

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

The last additive error, combined with the nonlinear bound, gives the positive finite-order ratio bound above, uniformly over all variance scales. Discarding B too early loses the improved rate. This asymptotic argument alone does not establish exact finite-n equality `K_n=kappa_n`, and its lower bound is not advertised as sharp at small orders. Section 10 supplies a distinct exact argument for n>=1500.

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

## 10. Exact fair-pair optimality for every n>=1500

**PROVED, independently audited, including positive-variance equality.** For all integer n>=1500, `K_n=kappa_n`. This proof depends directly on the already accepted uniform rounding, matching-support, triangle-variance, and H=1 results; it does not require a historical larger-threshold candidate.

Assume `V>0` and `Delta<=kappa_n V`, and keep the expected binary baseline `B=E Delta(T)`. Let D be the variance discarded by descending-weight greedy matching. The earlier bounds yield

```text
V>49(n-1)/1600,       V<=n/7,
B<=kappa_n D,         D<=2n*kappa_n.
```

Mark a selected edge light if its variance is at most `eta=1/44`, and let L be the number of light selected edges. If L>=1, delete those edges and apply the sharper nonperfect-matching coefficient to the retained heavy matching. Its retained variance exceeds n/120, so it has more than n/30 edges. The tail `(3/4)^(n/30)<=1/(1000n)` follows from an exact endpoint comparison and monotonicity, not testing many orders.

Using `floor(n/2)` rather than n/2 avoids a parity error. The retained-matching estimate then forces `L<(80/3)D`. Greedy edge groups have disjoint wedge contributions, giving

```text
D_heavy<=44B,
D_light<=2*sqrt(LB).
```

Together with `B<=kappa_n D`, these would imply

```text
1 < 44*kappa_n+2*sqrt((80/3)*kappa_n)
  <= 4/17+sqrt(320/561) < 1,
```

an exact contradiction: the final squared slack is `137/9537>0`. Therefore there are no light selected edges. Every selected variance is greater than eta and hence greater than kappa_n. Any discarded positive mass would then force `B>kappa_n D`, another contradiction. Thus `D=B=0`.

The uncertainty graph is consequently a matching, proving the universal lower bound; the fair ordinal family attains it. At equality, B=0 forces every supported binary completion to be transitive. The H=1 scalar minimization then gives **exactly a maximal ordinal collection of fair pairs, with one singleton at odd order**.

At these orders, V=0 equality is exactly the deterministic transitive case by the inherited binary theorem. This is distinct from the n=3 exception, where Delta=V for every weighted triangle.

The threshold 1500 is sufficient, not claimed optimal. The arithmetic controls verify endpoint and squared inequalities exactly; all-order scope comes from the monotonicity and matching argument. This is neither a tournament census at order 1500 nor an extension of an eight-vertex graph scan.

[Exact theorem](threshold_refinement/PROOF_N1500.md) · [Independent audit](independent/threshold1500_audit/AUDIT.md)

## 11. Reusable sufficient cyclic-response criteria and a minimal actual hostile

For any actual nonempty same-order modules M0,M1, let `d=f_M1-f_M0` be their ordered path-cover profile difference. The exact finite tridiagonal kernel of Section 9 applies without requiring the modules to share a deletion card:

```text
H(C3[M1,A,B])-H(C3[M0,A,B])=f_A^T K f_B
```

for all nonempty actual side tournaments A,B. Coordinates above the module order plus one vanish exactly; this is not a profile truncation approximation.

**Certificate P.** Entrywise `K>=0` implies nonnegative response for all actual sides. If `K11>0`, the sharp lower bound is `K11 H(A)H(B)>0`, with equality at singleton sides. Within this entrywise-nonnegative class, K11>0 is necessary for strict positivity over all nonempty sides. Entrywise positivity is not claimed necessary over realizable side profiles.

**Certificate N.** Supply an actual base V with positive log-concave profile v, a pivot c within its support, and `v^T K v<0`. Require every row before c to be entrywise nonnegative and every row after c to be entrywise nonpositive; the pivot row may be mixed. Then response is strictly negative for every actual pair whose colored profiles satisfy

```text
P_A=P_V R_A,      P_B=P_V R_B,
```

where the R polynomials are nonzero and monomial-nonnegative. The likelihood-ratio and one-crossing argument extends the finite base certificate to the whole stated chamber. Ordinal extra tournament factors supply actual examples. For an ordinal-power family, a supplied base `V=TT_k[F]` can certify all multiplicities at least k; the theorem does not assert that any successful k exists.

The exact certifier validates literal tournament masks before computing actual profiles. Its lower-level formal-profile interfaces do not imply realizability. A failed sufficient test returns **NOT_CERTIFIED**, not a negative theorem or a decision that no sign certificate exists.

### The row condition is indispensable even for actual strong modules

Two strong seven-vertex completions of the same six-vertex card have difference

```text
d=(-12,0,240,480,240,0,0).
```

For actual base K1 the profile `(1)` is log-concave and the response is −36. But with both sides TT2 the response is +828; with both sides C3 it is +67932. The multipliers are respectively m and `m^2+2`, both nonnegative. Thus removing the row condition makes the negative-propagation conclusion false for actual tournaments, not merely formal vectors. Here `K22=228>0`, explicitly violating the required high-row sign for the only possible base pivot c=1.

**Minimum module order seven is FINITE-EXACT.** A complete census through order six covers 33,867 labeled tournaments, 39 order-specific profiles, and all 308 distinct same-order profile-pair kernels. Every such kernel is entrywise one-signed. Therefore a negative base response cannot turn positive under any admissible multipliers at those smaller module orders. The producer reconstructs profiles from 23,717,429 literal vertex words; an independent implementation uses Hamilton-subset and ordered-cover dynamic programming instead. These are genuinely different finite kernels, not two independent literal-word traversals.

The short non-log-concave vector `(2,1,2)` remains a **formal** boundary control for a different hypothesis. No actual non-log-concave tournament profile is claimed. The accepted seven-vertex hostile retains log-concavity and isolates only the omitted row-crossing condition; it gives a sign reversal, not a scalar zero.

[General criteria and certifier](zero_framework_later/GENERAL_CRITERION.md) · [Independent audit](independent_zero_framework_later/AUDIT.md)

## 12. An explicit false-lemma boundary: odd does not mean unit

A proposed termwise Pfaffian shortcut is **REFUTED**, not used as a premise. A binary complementary Pfaffian is odd, but need not have magnitude one. The retained exact six-vertex example has full Pfaffian `3x-2` and complementary coefficient three. At x=2/3 the full square vanishes despite `x^2=4/9>0`, refuting the proposed individual-square lower bound.

The full weighted determinant margin is nevertheless `u-|v|=17/9`, with `u=32/9` and `v=-5/3`. Thus the example does **not** refute the proposed total matching-support margin, let alone an accepted stability theorem: other principal terms carry excess energy. The unit-coefficient version survives when the complementary Pfaffian is independently known to be ±1; a general compensating total-margin refinement remains open.

The accepted independent-set proof retains entirely binary principal subsets. The accepted multiplicative robustness proof uses the full Gram/Cauchy sum and its actual coefficient energy. Neither assumes that an arbitrary odd complementary Pfaffian is a unit. Independent full antisymmetrization and Leibniz sums verify the counterexample, separately from the producer's recursive Pfaffian and Gaussian determinant calculation.

[Exact hostile and scoped survivor](stability/PFAFFIAN_UNIT_COEFFICIENT_HOSTILE.md) · [Independent audit](independent_pfaffian_hostile/AUDIT.md)

## 13. Sharp order nine: the final finite continuation

**PROVED, independently audited.** The exact variance coefficient is 175/256, with the stronger reserve

```text
Delta >= (175/256)V+(105/256)tau+C4/16+(5/8)C5+C6/4+C7+C8.
```

The complete signed cycle-cover ledger includes the full covers 3+3+3 and 2+2+2+3 without choosing a sign for det A. Four-cycle grouping uses complementary five-vertex signed margins and the necessary factor one-half on complementary four-cycles. A six-vertex variance estimate controls the triangle coupling, yielding

```text
Delta >= P+(181/256)tau+M23+C4/16+(5/8)C5+C6/4+C7+C8,
P=V-M22+M222-M2222.
```

The separately affine residual is

```text
F=P+(19/64)W+P21-(175/256)V,
1024F=81m-64m2+16m3-4m4+19w+16p21
```

at graph endpoints. Complete augmentation coverage starts with all 1,044 seven-vertex orbits, covers 133,632 eight-vertex augmentations and exactly 12,346 eight-vertex isomorphism classes, then adds all 256 ninth-vertex neighborhoods. The 3,160,576 final augmentations cover every nine-vertex graph up to isomorphism; this is not a literal scan of all 2^36 masks.

The independent implementation uses adjacent-transposition BFS, equitable-refinement canonicalization, vertex-addition matching/wedge recurrences, and 344 separate full-DP checks. These differ from the producer's permutation orbit construction, degree-cell canonicalization, and fixed pattern-mask counts. Exact minimum is zero. Separate affinity extends the graph result to the entire real variance cube; the local triangle and disjoint-edge lifts complete the theorem.

The wedge coefficient 19/64 is optimal **for this graph reduction**, forced by P3 disjoint-union three edges. It does not assert that 105/256 is the optimal actual triangle reserve. The fourteen zero augmentations represent three graph types, and some are joined by genuine continuous zero segments of the graph polynomial. Thus graph-zero geometry alone cannot force the tournament equality locus.

Instead the strictly positive triangle reserve forces tau=0 at positive-variance equality. Every supported binary completion is then transitive; the accepted H=1 classification gives exactly **four fair ordinal pairs and one singleton**. That family has H=1, u=81/256, v=0, V=1 and Delta=175/256.

[Order-nine theorem](stability/SHARP_NINE_GRAPH_CERTIFICATE.md) · [Independent graph audit](independent/n9_audit/AUDIT.md)

## What is not concluded

The exact stability coefficient is proved for n=3,...,9 and every n>=1500. **Orders 10 through 1499 remain unresolved.** The unrestricted H-spectrum, global boundedness of exact spectral width, arbitrary strong-quotient scalar zero-rigidity, and LRC(14) also remain open. The no-zero theorem closes one explicit unbounded construction family, not every mixed-profile seed or every actual side context.
