# Independent audit of the uniform diagonal-density theorem

**PASS: full analytic audit and separate exact controls. No proof repair
required.** The theorem in
[the producer report](overnight8_20260906_no3line_diagonal_density.md)
is accepted uniformly over all simple bipartite 2-regular source graphs.
Its joint-Poisson companion is accepted with the stated fixed-number-of-
targets condition and the repaired Holder domination.

The audited producer source SHA256 is
`ead04c2f82c7c77a4f94c5b2ca379c7cc56dcdb37e252f4276a3826b2972c200`;
the frozen output SHA256 is
`87e151596665aac1ec36d10f348a9555bc68e981d43ad317cb017bfdb431e553`.
No producer files were edited by this audit.

## Accepted result and actual target

For independent uniform injective shore labels of any simple 2-regular
bipartite graph with n vertices per shore, let S count triples on all
slope-one grid diagonals and let X count all nonaxis collinear triples.
The accepted formulas are

```text
E S = (2n-5)/3,                         n>=3,
Var S = (40/9)n+O(1),                   uniformly in the cycle type,
E[(S/n-2/3)^2] = 40/(9n)+O(n^-2),
limsup_(n->infinity) n sup_G P_G(X=0) <=10.
```

All coordinate determinants are ordinary integer determinants. Every
slope-one triple is one actual collinear triple, with no multiplicity,
so `0<=S<=X`. Degree two forbids axis triples. Consequently X=0 is the
actual no-three-in-line event in this saturated model, and
`P(X=0)<=P(S=0)` is the correct direction. S=0 is not sufficient for
success; both the producer and this audit retain a literal hostile.

The exact mean is independent of cycle type, and the variance remainder
is uniform. Thus mixtures inherit the same formulas, without a between-
type variance term. Uniformly sampling all distinct saturated boards is
one such mixture, with orbit-size weights. The conclusion concerns their
density. It proves neither eventual nonexistence nor a sharp decay rate,
and supplies no exponential bound or CLT.

## Uniform overlap and covariance audit

The probability normalization is exact. After expanding the factorial
moments, the selected grid cells form a graph H whose row and column
vertices are equality classes of the specified coordinates. Its event
probability is `inj_shore(H,G)/((n)_r(n)_s)`. The cells and their positions
are named; an additional automorphism divisor would be wrong. Non-induced
containment, rather than induced subgraph probability, is essential.

For a fixed colored pattern, each connected component has at most n
target placements: choose one root coordinate and follow the prescribed
partial-matching colors. Every further vertex is determined, and cycle
closure can only reject the placement. In the source graph, choosing
the component roots and a spanning tree gives at most
`n^c 2^(v-c)` maps. Requiring global injectivity can only reduce this.
For n>=2k, the injection denominator is bounded below by a fixed
multiple of n^v. Thus the contribution is at most a constant depending
only on k times `n^(2c-v)`. There are only finitely many equality patterns
for fixed k. None of these constants contains a cycle count or component
size of G.

Each component has at least two vertices. The only order-one pattern is
a set of isolated edges; the only order-1/n pattern is one two-edge path
plus isolated edges. Every other simple pattern has v-2c>=2. In particular
a four-cycle already costs two powers, even when G has linearly many
four-cycles. Degree-forbidden patterns have zero contribution. This is
the required **uniform**, rather than merely pointwise, O(n^-2) bound.

The first matching correction is also uniform. An edge has itself and
exactly two intersecting neighbors, so one bad position-pair has exactly
6n edge assignments. The intersection of two distinct position-pair
constraints has at most k-2 free roots, even if those constraints share
a position. Bounded degree controls its remaining choices. The first two
Bonferroni bounds therefore give the claimed second-order remainder in
the ordered matching count. Expansion of the actual falling denominators
then yields

```text
p_k=(2/n)^k(1+k(k-1)/(4n)+O_k(n^-2)).
```

For a prescribed two-edge path, its named center and ordered neighbors
have exactly2n embeddings. Adding k-2 independent edges gives leading
count `(2n)^(k-1)`; intercomponent collisions lose a free root. Dividing
by `(n)_(k-1)(n)_k` gives the stated leading path probability. This
normalization is the same for either shore as center.

For two edge-disjoint target matchings, the one-overlap tuple count is
`ab(R+C)L^(a-1)M^(b-1)+O(n^(a+b-2))`. A fixed row-overlap pair has
exactly `R(L-1)_(a-1)(M-1)_(b-1)` preliminary choices. Additional
collisions lose a second free choice; the identical argument holds on
columns. The same pair cannot share both coordinates because the target
cells are disjoint. These are polynomial bounds in L,M<=n, so the
remainder remains uniform when a target is short or empty.

Subtracting the product of the two exact single-matching moments leaves

```text
ab*2^(a+b-1)/n *
 [theta^a eta^b-theta^(a-1)eta^(b-1)(R+C)/n]+O(n^-2).
```

The first term is the difference between the matching correction at
a+b and those at a,b; the second is the fact that the path probability
is half the leading matching probability. Both terms are necessary.
There is no independent-occupancy assumption in this calculation.

## Variance constant and summation audit

The exact third matching count and the sum of diagonal triple counts
give the reported mean, including its -5/3 boundary correction. The
falling-factorial square identity produces the single-line variance
kernel `8theta^5+8theta^4+(4/3)theta^3` with a uniform O(1/n) error.
Its aggregate leading coefficient is98/15.

For two distinct diagonals the row-plus-column overlap is exactly twice
the minimum of the two lengths on the same sign side, and twice the
positive part of their length sum minus n on opposite sides. The formulas
agree at the main diagonal. The corresponding compact kernel is piecewise
polynomial and globally Lipschitz, so its two-dimensional Riemann-sum
error is O(1/n), with a geometric constant independent of G. Inserting or
removing the diagonal pairs changes the covariance sum by only O(1).

The referee independently obtains

```text
integral x^2 y^2 min(x,y) =1/14,
integral x^2 y^2(x+y-1)_+ =71/1260,
cross coefficient =2-32(1/14+71/1260)=-94/45,
98/15-94/45=40/9.
```

There are O(n^2) uniform O(n^-2) covariance remainders and O(n) uniform
O(n^-1) single-line remainders. Their sum is O(1); this is the needed
justification for using a growing family of diagonals. The squared bias
is25/(9n^2), so the stated L2 formula follows. Finally Chebyshev gives
`Var(S)/(E S)^2`, whose leading coefficient is
`(40/9)/(2/3)^2=10`. All inequalities apply to the actual zero event.

## Joint-Poisson companion and the corrected domination

Conditioning on the row labels, a single target matching gives each
source column vertex at most two allowable labels, each contributing at
most one hit. Counting k nonattacking assignments proves
`E[(Y_T)_k|rho]<=2^k`, including the zero case k>n. Its factorial
generating series gives the valid one-target exponential bound.

A union of targets can give a single assignment multiple hits. The
referee independently recomputed the n7 hostile through a weighted
permutation matrix W. Each row and column sums to4, and each row's
squared-entry sum is6. Therefore

```text
E(Y)_2 =
 [(sum W)^2-sum(row sums)^2-sum(column sums)^2+sum W_ij^2]
       /[n(n-1)] + [sum W_ij^2-sum W]/n
 =18-10/(n-1).
```

At n7 this is49/3>16. The discarded binary-union domination is false.
Holder correctly replaces it by
`E A^(sum Y_i)<=exp(2(A^K-1))` for fixed K and A>=1.

The factorial Vandermonde identity bounds the total degree-k generating
terms by `E(sum Y_i)_k R^k/k!`. Evaluating the valid bound at1+2R gives
a geometric2^-k tail on each bounded complex polydisc, uniformly in n,G.
Termwise application of the fixed-order overlap lemma is consequently
legal. Passing from the factorial to ordinary generating variables and
taking Cauchy coefficients yields the stated joint Poisson limit,
including zero means. This argument keeps K fixed. The main density
theorem uses the finite-order uniform covariance proof, not a growing-K
Poisson assertion.

## All-C4 companion and primary-source scope

The block count in equation(16) correctly distinguishes blocks used
twice from those used once, with factors2 and4 and the two factorial
divisors. At k3 it gives exactly the stated two-block and three-block
containment probabilities and their fractions of the common full mean.
The referee requested making `n>=4` explicit beside these auxiliary
three-matching probabilities: n2 has no grid3-matching and their displayed
denominator is zero. This is a scope clarification, with no change to
the main theorem, producer formulas, or computed universe.

The referee independently checked the primary
[Guy–Kelly paper](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B126DA7E4957722BAC70AC7B7F6E1FA2/S0008439500056770a.pdf/the-no-three-in-line-problem.pdf):
its theorem gives the grid-triple asymptotic, while its later independence
step is explicitly assumed. Removing axis triples changes only O(n^4).
Multiplying the proved count by the two exact containment probabilities
therefore gives the reported logarithmic and n-log-n means. The
[Voutier v2 primary record](https://arxiv.org/abs/2603.00215v2)
explicitly concerns a corrected heuristic. Neither source supplies or
is used as a proof of the new zero-event estimate. No external priority
is asserted by this audit.

## Independent exact controls and frozen files

The [referee source](../../04-computation/overnight8_20260906_no3line_independent_audit.py)
uses only the standard library and imports no producer. It enumerates
**distinct boards directly**, by selecting each row's two columns and
maintaining exact column degree two. It generates all6 boards at n3,
90 at n4, and2,040 at n5, then classifies the resulting bipartite cycle
components. This is a separate sample construction from shore-label
enumeration. The orbit count `n!^2/product_l((2l)^m_l m_l!)` is checked
to verify that the two sample spaces give the same conditional laws.

All five producer mean/variance/diagonal-zero rows match exactly, as do
the n4 full-zero probabilities1/36 and1/2. Additional finite controls
give n5 full-zero probabilities1/45 on C10 (32 of1,440 boards) and0
on C4+C6 (zero of600). These are bounded controls, not asymptotic claims.
A distinct S=0,X=1 board is retained.

The source also evaluates both overlap integrals by elementary rational
polynomial arithmetic. A second route reconstructs the exact discrete
diagonal-kernel sums, whose numerator degrees are at most8 and6 after
splitting into integer triangles. Independent extra integer arguments
verify the resulting polynomials. They have leading coefficients
`-47/180` and `98/5`, which recover the covariance and local constants
after the respective factors8 and1/3. These kernels are the leading
covariance approximation, not claimed to be the full finite variance.
The weighted-union hostile is checked by the separate matrix identity
above, without a permutation census.

Normal and optimized runs have byte-identical LF output and pass
**4,357 explicit gates**:

```text
python -B 04-computation/overnight8_20260906_no3line_independent_audit.py
python -B -O 04-computation/overnight8_20260906_no3line_independent_audit.py

source 563739706a1e5aa3dcdbd06112f713e167947eb6886c80ed3198d97d6d2807cf
output b51ca52a7e37226531eda9cb39746cc73cf02e8386f40fc1d6cc48093a0fe874
semantic bank 557c73eb2b0d8a0f4576418c76ba9301e06dba9c47a1bb00e7b77c25fc7e34f3
```

All files are outside the repository for parent-managed integration.
No incoming mathematical file, shared navigation, or Git state was
changed by this audit.

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
