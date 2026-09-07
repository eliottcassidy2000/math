# Independent referee: two-direction deletion cost and adaptive repair

**Status: PASS — analytic proof and independent exact controls.**
The repaired proof, including the stronger four-line-length constant below, is accepted. The only requested wording repair was to restrict the bound `4L(L-2)` to line lengths `L>=3`; lines of length below three contribute zero. The parent applied that clarification. No mathematical counterexample was found.

Primary proof: [primary proof](continuing10_20260907_no3_repair.md). This referee did not read, import, or execute its producer. Its independent verifier is entirely standard-library rational and combinatorial code. The inherited row-conditioned mean and permutation concentration are used in their already audited scopes; row order is arbitrary and only columns are uniform.

## 1. Deterministic object and adaptive consumer

The rotated difference/sum graph retains each original cell as a distinct unit-capacity edge. Capacity two at each endpoint encodes occupancy at most two in both selected diagonal directions. An integral maximum flow therefore retains exactly the largest subset meeting those two constraints. Edge capacity one is needed to prevent duplication of a cell. No condition on other integer slopes is silently added to this potential.

Let tau be the minimum number of deletions for these two directions. Each overfull difference line needs at least its occupancy minus two deletions, with total F. Points lying on difference lines of occupancy at most two are disjoint from that first set. At least `(Z_s-2)_+` such points must separately be deleted on sum line s. Thus the primary inequality `tau>=F+sum_s(Z_s-2)_+>=F+I` is valid even if other deletions are also necessary. Each line counted by I has three points on three different singleton difference lines.

For equal-cardinality boards B,B', intersection of an optimal retained subset with B' loses at most `|B\B'|` points. Applying this in both directions gives the stated replacement Lipschitz bound. A transposition of two row or two column labels on a saturated board changes at most four original cells.

For an arbitrary final no-three-in-line set T, the retained set `B intersect T` satisfies both selected diagonal constraints. Therefore `|B\T|>=tau(B)` deterministically. The algorithm may inspect the start and choose its later operations adaptively. Its success with an original-point loss budget k is contained in the event `tau<=k`; no independence of its decisions or endpoint sampling is required. Counts of local operations are merely ways to upper-bound that loss budget. The statement supplies no time bound for unrestricted search.

## 2. Uniform permutation avoidance

For a forbidden N-by-N board of maximum row and column degree Delta and m cells, the primary bound on ordered k-tuples is correct. Choose a pair of positions. Once the first forbidden entry has been selected, at most `2 Delta` partners share a row or column. The other k-2 positions are arbitrary. Union over the `binom(k,2)` pairs gives

```
0 <= m^k-Q_k <= k(k-1) Delta m^(k-1),    k>=2.
```

Repeated entries are included in this overcount. The cases m=0 and k=0,1 are harmless and treated separately. For every fixed k, `Q_k/(k!(N)_k)` differs uniformly from `(m/N)^k/k!` by `O_(k,Delta)(1/N)`. The range `0<=m/N<=Delta` is compact. Odd/even Bonferroni bounds and uniform exponential Taylor remainders give

```
P(no forbidden permutation entries)=exp(-m/N)+o(1),
```

uniformly over all such boards. The order of limits is correct: choose a finite Bonferroni order for the desired uniform error, then N sufficiently large. A finite truncation is not being claimed uniformly in an unbounded order.

## 3. Actual eight-choice event and uniform discarded count

The skeleton is simple and 2-regular, so every source column has exactly two distinct row neighbours. On a fixed sum line, three distinct selected row labels determine distinct target column labels. If two selected rows shared a source column, each would be the other's companion row, which eligibility excludes. Hence all eight choices of their source neighbours use three distinct source columns.

The companion of the assigned cell `(r,s-r)` is `(t,s-r)`, where t is the other row incident with that source column. It cannot lie on the selected sum line because t differs from r. It lies on the difference line of the selected row r_j exactly when `r+t=2r_j`. The eligibility condition excludes this equality for each other selected row; for r_j=r it would again imply t=r. Thus every forced companion lies outside the full four-line union U.

The eight intended assignment events are pairwise disjoint: a fixed target label cannot be assigned to both distinct source columns that could create the same target cell. Conditional on any one intended assignment, the remaining permutation is uniform on N=n-3 labels. Each free source column has two row neighbours, each meeting U in at most four cells, so it forbids at most eight labels. A target label meets U in at most four rows, each with two source neighbours, so its forbidden column degree is also at most eight. This remains true after the forced columns and labels are removed.

Avoiding all remaining forbidden entries makes the selected sum line have exactly the selected three cells, with each selected difference line occupied exactly once. It therefore produces an actual I event. Distinct selected triples on the same sum line produce disjoint such events because the occupancy is exactly three. Events on different sum lines need not be independent; expectation only sums their indicators.

Each row r has two companion rows t and at most two midpoint labels `(r+t)/2`, hence at most four forbidden other row addresses. For `L>=3`, the number of unordered ineligible row triples is bounded above by `4L(L-2)`, with any overcount harmless. Summing over the at most 2n-1 lines gives `O(n^3)`, uniformly in the skeleton and row labelling. The full selected triple count is exactly

```
sum_s binom(L_s,3)=binom(n,3)+2binom(n,4)=n^4/12+O(n^3).
```

After multiplication by `8/(n)_3`, all discarded triples cost only O(1) in the expectation. The uniform avoidance o(1), multiplied by the entire event count, costs o(n). Neither estimate depends on the fixed row ordering.

## 4. Stronger four-line mean constant

The count m of forbidden source-column/target-label pairs satisfies

```
m <= 2|U| <= 2(L_sum + L_diff,1 + L_diff,2 + L_diff,3).
```

Indeed each grid cell of U forbids at most the two source columns incident with its row. Duplicate forbidden entries and removed assigned labels only decrease m. Retaining the actual lengths therefore improves `exp(-8)` to the exponential of this geometric envelope divided by n-3.

Write the normalized sum-line length as theta in [0,1]. Its row parameter is uniform on an interval of length theta; the difference coordinate `d=s-2r` ranges over `[-theta,theta]`. A difference line at d has normalized length `1-|d|`. Hence for three selected rows the limiting exponential weight is

```
exp[-2(theta+3-|d1|-|d2|-|d3|)].
```

The unordered-triple Riemann measure on one line contributes `theta^3/6`. There are asymptotically two copies of each length, and eight disjoint intended assignments. Removing coincident and ineligible triples changes only the already bounded `O(n^3)` lower-order terms. The continuous limiting integrand is bounded, including at theta=0. Thus, uniformly over the discrete skeleton and fixed rows,

```
liminf E I/n >= epsilon,
epsilon=(8/3)exp(-6) integral_0^1 theta^3 exp(-2theta)
              [ (exp(2theta)-1)/(2theta) ]^3 dtheta
       = exp(-2)/12-exp(-4)/2+5exp(-6)/4+exp(-8)/6.
```

The product inside this deterministic triple integral is the integration of a product over a cube of three row coordinates. It is not an assertion that diagonal occupancy events or permutation assignments are independent.

Independent simplification gives `epsilon=exp(-6)/3` times the integral of `exp(4theta)-3exp(2theta)+3-exp(-2theta)`, confirming every coefficient. A separate Jensen check uses the theta^3-weighted mean `E theta=4/5`, giving mean exponent cost `6-4/5=26/5` and the valid lower bound `(2/3)exp(-26/5)`. In particular the new constant strictly improves the earlier conservative `(2/3)exp(-8)`.

With `alpha=1-5exp(-2)` and `gamma=alpha+epsilon`, the inherited uniform conditional bound `E F>=alpha n-21` and `tau>=F+I` give the stated uniform lower asymptotic mean gamma. No explicit finite threshold for this improved asymptotic constant follows from the argument.

## 5. Concentration, mixtures, and counting

Expose n-1 column labels. Two possible choices at a reveal are coupled by a transposition of the remaining labels. The conditional mean range is at most four by the deterministic Lipschitz bound. Hoeffding's bounded-range argument gives

```
P(tau<=k | G,rho) <= exp[-(E tau-k)_+^2/(8(n-1))].
```

The independent row-conditioned mean must be used here; an arbitrary mixture mean is not substituted into this conditional estimate. The inherited common finite lower bound b_n remains usable because tau>=F. Uniformity of the new asymptotic mean then gives, for `k_n<=kappa n` and `kappa<gamma`, the uniform rate `-(gamma-kappa)^2/8`. At k=0 this bounds full no-three success, since that event implies zero deletion potential.

Mixtures over G and row labellings preserve a uniform conditional bound when columns remain conditionally uniform. The inherited exact counting model has that property, so the ratio of successful uncolored saturated boards to all such boards inherits rate `-gamma^2/8`. No result is asserted for a changed starting distribution lacking this conditional uniformity, or for nonexistence of extremal boards.

## 6. Independent finite controls

The verifier uses exhaustive retained-edge subsets rather than the primary's integral-flow implementation. It checks the exact deletion inequality on 30 fresh saturated boards, every safe retained subset of those boards, and 210 transpositions. A separate finite forbidden-board engine checks all permutation outcomes for 24 bounded-degree matrices through N=6, exact rook factorial moments, the collision estimate and both Bonferroni directions.

For n=4 through 8 and two frozen row orders, it checks all 754 selected native row triples, finds 126 eligible triples, verifies all 1,008 intended source assignments, and checks both degree-eight bounds and the sharper actual four-line-length count. It independently completes 4,800 conditional permutations and verifies the literal occupancy event, not just a forbidden-matrix proxy. Exact rational bookkeeping recovers the new exponential integral and the Jensen cost.

The finite controls are positive and hostile checks of the proof mechanisms. They are not used to infer an asymptotic statement from a finite board census. No producer is imported or run.

Reproduction:

```
python3 -B 04-computation/continuing10_20260907_no3_repair_audit.py
python3 -B -O 04-computation/continuing10_20260907_no3_repair_audit.py
```

The parent owns final proof filing and promotion. The retained-point coordinate is the faithful link from this potential to arbitrary adaptive repair; it is stronger than a statement about independent restarts while remaining distinct from a universal search-complexity lower bound.

Exact alternating Taylor enclosures independently verify the additional numerical statement `epsilon>1/200`. The updated primary proof Sections 3-4, including (4a)-(4b), were reread after the stronger constant was added.

Both frozen runs pass **22,186 always-active exact gates**, with byte-identical actual LF stdout.

- Referee source SHA256: `4f9a9bf9e5a8ef99a047764e477871519f327c34e7453f3d399dc249e3061f8b`.
- Normal and optimized stdout SHA256: `617c410c8721a2517af9f3ca517a965dcfd15156358cc77348c3da6e37680971`.
