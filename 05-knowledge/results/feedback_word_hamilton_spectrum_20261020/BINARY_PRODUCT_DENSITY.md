# Binary-product density and an all-cutoff obstruction to intercyclic coverage

**Scope: ordinary Hamilton-path counts. Proof-based asymptotic theorem with small independent controls. The global forbidden-set conjecture remains OPEN.**

This is a separate addendum. It does not modify the three previously frozen packets. The product approximation and exact dyadic-phase cancellation were developed jointly during independent proof review. No external novelty or priority claim is made.

Write `S_1` for the set of ordinary Hamilton-path counts of tournaments that become transitive after deleting at most one vertex, and `S_ic` for the corresponding spectrum of tournaments having no two vertex-disjoint directed cycles. All logarithms without a subscript are natural logarithms.

## 1. Main conclusions

For a binary word `b=b_1...b_m`, let `H(b)` count Hamilton paths in the tournament consisting of a transitive chain `v_1 -> ... -> v_m` and an apex `x`, where `b_i=1` means `x -> v_i`. Call a word **reduced** if `m>=2`, `b_1=1`, and `b_m=0`.

Let

`W(X) = #{reduced binary words b, of any length, with H(b)<=X}`.

Then, as real `X -> infinity`,

`W(X) = X log 2 + O(sqrt X)`.                                      (1)

This counts words, not distinct numerical values or isomorphism classes. Complement-reversal identifies two words with the same H, except at explicitly counted fixed words. Consequently,

`|S_1 intersect [1,X]| <= (log 2 / 2) X + O(sqrt X)`,               (2)

and the exact intercyclic-spectrum theorem from `intercyclic_rays.md` gives

`|S_ic intersect [1,X]| <= (log 2 / 2) X + O(sqrt X)`.              (3)

There are `X/2+O(1)` positive odd integers up to X. Thus the **upper relative asymptotic density** of `S_ic` among positive odd integers is at most

`log 2 = 0.6931471805...`,

and the lower relative asymptotic density of its complement among positive odd integers is at least

`1-log 2 = 0.3068528194...`.

These bounds hold at **all cutoffs**, not only dyadic cutoffs. They do not assert existence of a density for the distinct-count spectrum. Nor do they assert that the excluded odd values are globally unattainable: any such value that is attainable requires at least two vertex-disjoint triangles.

## 2. Exact binary-product identity

The weighted-inversion formula proved directly from Hamilton path words in `proofs.md`, Section 2, is

`H(b)=1+2 sum_{i<j} b_i(1-b_j) w(j-i)`,

where `w(1)=1` and `w(d)=2^(d-2)` for `d>=2`.

Define

`N(b)=sum_{i=1}^m b_i 2^(m-i)`,

`R(b)=sum_{i=1}^m (1-b_i) 2^(i-1)`,

`A(b)=#{i: b_i=1 and b_{i+1}=0}`,

`B(b)=sum_{i>j} b_i(1-b_j) 2^(j-i-1)`.

Expanding the product gives

`N(b)R(b)/2^m = sum_{i,j} b_i(1-b_j) 2^(j-i-1)`.

Terms with `i=j` vanish. Terms with `i<j-1` agree exactly with their coefficients in `H(b)-1`. An adjacent inversion has coefficient one in the product and coefficient two in `H(b)-1`. All backward terms are exactly B. Therefore

`H(b) - N(b)R(b)/2^m = 1 + A(b) - B(b)`.                         (4)

Adjacent `10` occurrences cannot overlap, so `0<=A(b)<=floor(m/2)`. Moreover,

`0<=B(b)<=sum_{i=1}^m sum_{d=1}^{i-1} 2^(-d-1)<=m/2`.

It follows that

`|H(b) - N(b)R(b)/2^m| <= 1+m/2`.                               (5)

In particular, for `U_m=N/2^m` and `V_m=R/2^m`,

`|H(b)/2^m-U_m V_m| <= (1+m/2)2^(-m)`.                          (6)

The estimate holds for every binary word, not only reduced words.

## 3. A uniform product limit, with an explicit error

Sample a reduced word of length m uniformly. Its `m-2` interior bits are independent unbiased bits. Let

`s=floor((m-2)/2)`.

The first s free binary digits of `U_m` are `b_2,...,b_{s+1}`. The first s free binary digits of `V_m` are `1-b_{m-1},...,1-b_{m-s}`. These index sets are disjoint because `2s<=m-2`. Both variables have first binary digit one.

Append fresh independent unbiased tails to these two prefixes. The resulting variables U and V are independent and uniform on `[1/2,1]`. They can be coupled to the original word so that, deterministically,

`|U_m-U|<=2^(-s-1)` and `|V_m-V|<=2^(-s-1)`.

All four variables lie in `[1/2,1]`, apart from irrelevant endpoint conventions. Thus

`|U_m V_m-UV|<=2^(-s)`,

and (6) gives

`|H(b)/2^m-UV|<=delta_m`,

`delta_m=2^(-s)+(1+m/2)2^(-m)`.                                 (7)

This argument does **not** claim that `U_m` and `V_m` are independent. They are deterministic functions of the same finite word. Only their selected disjoint prefixes, and the limiting coupled variables, are independent.

The distribution function `F(z)=P(UV<=z)` is obtained by integrating the constant density four over the square `[1/2,1]^2`:

```text
F(z) = 0,                              z <= 1/4;
       4z log(4z)-4z+1,                1/4 <= z <= 1/2;
       4z-3-4z log z,                  1/2 <= z <= 1;
       1,                              z >= 1.
```

For the first interior interval, integrate `4(z/u-1/2)` from `u=1/2` to `u=2z`. For the second, the full strip `1/2<=u<=z` is allowed, followed by the partial strip `z<=u<=1`. These integrations give the displayed formulas.

The derivative is `4 log(4z)` on `(1/4,1/2)` and `-4 log z` on `(1/2,1)`. Hence F is globally Lipschitz with constant

`K=4 log 2`.

Write `C_m(z)=P(H(b)/2^m<=z)` for the finite reduced-word distribution function. The deterministic coupling implies

`F(z-delta_m)<=C_m(z)<=F(z+delta_m)`,

and therefore the quantitative **uniform** CDF estimate

`sup_z |C_m(z)-F(z)| <= K delta_m`.                              (8)

Its right side is `O(2^(-m/2)+m2^(-m))`. Uniformity is essential below because the dyadic phase changes with X.

### Moment cross-check

For uniformly sampled reduced words of length `m>=3`, the exact mean is

`E H = 9*2^(m-4)+1/2`.

To verify this directly, put `l=m-2`. The elementary identity

`sum_{d=1}^{l-1}(l-d)w(d)=2^(l-1)-1` for `l>=1`

follows by induction, since the increment at l is `sum_{d=1}^l w(d)=2^(l-1)`. In the weighted-inversion formula, the outer pair contributes `2^(m-2)`, each group of endpoint/interior pairs has expected contribution `2^(m-3)`, and the interior pairs contribute `(2^(m-3)-1)/2`. Including the constant one gives the stated mean.

The coupling also gives convergence of the first two moments, because both normalized H and UV lie in `[0,1]` and their difference tends to zero uniformly. Since `E U=3/4` and `E U^2=7/12`,

`E(H/2^m) -> 9/16`,

`Var(H/2^m) -> (7/12)^2-(9/16)^2=55/2304`.

These distributional checks are consistent with, but are not substitutes for, the uniform CDF argument.

## 4. Only two word lengths are partially counted

Every reduced word satisfies

`2^(m-2)<H(b)<2^m`.                                             (9)

For `m>=3`, the outer inversion `(1,m)` alone gives the lower bound `H>=1+2^(m-2)`. At `m=2`, the only reduced word is `10`, with `H=3`, so the same strict bound holds.

For the upper bound, deleting x from a Hamilton path leaves an increasing chain segment before x and an increasing segment after x. Thus a path is uniquely specified by the subset of chain vertices before x, giving at most `2^m` possibilities. The subset `{v_1}` is invalid because `x -> v_1`, whereas the corresponding path would require `v_1 -> x`. Hence `H<=2^m-1`.

Let `X>=4`, `L=floor(log_2 X)`, and `t=X/2^L`, so `1<=t<2`. By (9), all reduced words of length at most L are counted, no words of length at least `L+3` are counted, and only lengths `L+1,L+2` are partial. There are `2^(m-2)` reduced words of length m. Therefore the following is an **exact** identity:

`W(X)=2^(L-1)-1`

`     +2^(L-1) C_(L+1)(t/2)+2^L C_(L+2)(t/4)`.                (10)

Applying (8) yields

`W(X)=2^(L-1)[1+F(t/2)+2F(t/4)]-1+E(X)`,

where

`|E(X)|<=K[2^floor(L/2)+2^floor((L+1)/2)+(2L+7)/8]`.             (11)

Indeed, multiplying `delta_m` by `2^(m-2)` gives

`2^floor((m-1)/2)+(m+2)/8`.

The bound in (11) is `O(sqrt X)` uniformly over the phase `1<=t<2`.

## 5. Exact cancellation of the dyadic phase

For `1<=t<2`, the arguments `t/4` and `t/2` lie in the first and second interior branches of F respectively. Substituting their formulas gives

`F(t/2)=2t-3-2t log t+2t log 2`,

`2F(t/4)=2t log t-2t+2`.

Consequently,

`1+F(t/2)+2F(t/4)=2t log 2`.                                   (12)

Combining (10)–(12) proves the sharper explicit form

`W(X)=X log 2-1+E(X)`,

with E bounded by (11), and hence proves (1). No limiting subsequence or restriction to powers of two is used.

## 6. Passing from words to the ordinary-H spectrum

Deleting leading zeros and trailing ones removes only source and sink vertices and preserves H. Thus every nontransitive value in `S_1` is represented by a reduced word. Conversely every such word supplies a one-apex tournament.

The involution

`b_1...b_m -> (1-b_m)...(1-b_1)`

preserves reduced words and H, by tournament duality and reversal of the chain order. The sublevel set counted by W(X) is invariant under this involution.

A fixed word requires `b_i=1-b_(m+1-i)`, which is impossible for odd m. At even m there are exactly `2^(m/2-1)` fixed reduced words. As (9) bounds relevant lengths by `L+2`, the number `J(X)` of fixed words in the sublevel set satisfies

`J(X)<=sum_{1<=r<=floor((L+2)/2)}2^(r-1)`

`     =2^floor((L+2)/2)-1=O(sqrt X)`.

There are `(W(X)+J(X))/2` involution orbits, and H is constant on each orbit. Different orbits may still give the same count, so this proves only an upper bound, which is exactly what is needed:

`|S_1 intersect [1,X]|<=1+(W(X)+J(X))/2`.

The extra one accounts for the transitive value `H=1`. Combining the bounds gives (2), or explicitly for `X>=4`,

`|S_1 intersect [1,X]| <= (log 2/2)X`

`  +2 log 2[2^floor(L/2)+2^floor((L+1)/2)+(2L+7)/8]`

`  +2^(floor((L+2)/2)-1)`.                                      (13)

The independently audited exact intercyclic theorem states

`S_ic = S_1 union {3*2^r+1:r>=2} union {9*2^r-3:r>=1}`.

Each added ray has only `O(log X)` terms up to X, proving (3). Since every directed cycle in a tournament contains a directed triangle, the complement of the intercyclic class is exactly the class containing two vertex-disjoint triangles.

The consequence is a stronger obstruction to the **intersecting-cycle construction strategy**, not a resolution or refutation of the conjecture that the only forbidden positive odd ordinary Hamilton-path counts are 7 and 21.

## 7. Verification and provenance

Reproduce the small controls from the parent workspace:

```text
python session16-hamilton-spectrum/verify_binary_product.py
python -O session16-hamilton-spectrum/verify_binary_product.py
```

The checks cover every reduced word through length twelve, with independent subset-DP checks through length eight; the exact integer-scaled identity (4), extremal cutoff (9), complement-reversal, disjoint-prefix independence, finite CDF bounds, non-dyadic phase cancellation, and the exact two-length formula (10). Floating-point evaluations are used only to test the analytic logarithmic formulas, not to establish the theorem or any exact combinatorial equality. The program uses explicit exceptions, not assertions, so optimized mode retains every check.

Dependencies: the weighted formula and elementary path encoding in `proofs.md`; the independently audited exact spectrum in `intercyclic_rays.md`. Its cited Bang-Jensen–Bessy feedback theorem and Irving–Omar/Grinberg–Stanley odd-cycle formula enter only through that prior intercyclic theorem. The present binary-product asymptotic is proved here using elementary bit independence, integration and orbit counting. No large tournament enumeration is used.
