# A uniform diagonal-density obstruction for randomly labeled saturated boards

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent referee](overnight8_20260906_no3line_independent_audit.md)
accepts the uniform remainders, variance constant, probability consequence,
and repaired joint-Poisson companion. Its separate board generator and
finite-diagonal polynomial path reproduce all constants and small controls.
This report proves a probability bound uniform over every degree-two cycle
type. It gives no nonexistence theorem and no external priority claim.
Files are staged outside the repository pending the parent agent's review.

## 1. Inheritance, source boundary, and the chosen coordinate

The closest proved mechanism is
`05-knowledge/results/overnight6_20260906_no3line_palettes.md`, together with
its independent audit: keep actual injective row and column labels, use
whole event unions, and only then average. The earlier non-induced copy
normalization is in
`05-knowledge/results/overnight_20260906_moments_pairprofile_theorem.md`.
The finite hostile is the n=4 pair C8 and 2C4: both have full triple mean 2,
but zero-event probabilities 1/36 and 1/2. The corrected near miss is the
nonpositive formal type-weight functional in
`05-knowledge/results/overnight4_20260906_no3line_cycle_defect.md`.

The least-used sidecar here is the **complete colored overlap graph of
points chosen from specified target matchings**, with its row and column
injections. This is smaller than the full line-event census, but retains
exactly the collision information needed for the first covariance term.

Concept board: injective labels; actual palettes; colored graph overlaps;
bounded degree; line occupancy factorial moments; a one-sided zero-event
gate. Source is the positive injective-label probability space; target is
the vector of line occupancies; map is intersection with specified grid
matchings. It preserves the implication that an occupied line with at
least three points makes the board fail. It discards other directions and
most cycle information. The needed sidecar is the whole colored union,
not a product of scalar cycle responses. A finite board with S=0 but X>0
below is the cheapest hostile to reversing the implication.

The Anchor remains LRC(14), open; the Niche is the independent Laurent/jet
portfolio; this Wildcard uses the same retained-coordinate principle but
does not import a theorem between those subjects.

The user-supplied [Guy--Kelly paper](https://doi.org/10.4153/CMB-1968-062-3)
and [Voutier's correction, v2](https://arxiv.org/abs/2603.00215v2) are the
historical context. The latter explicitly discusses a corrected heuristic,
not a proved asymptotic upper bound. Their current primary records were
checked again. The proof below is self-contained and does not assume their
independence heuristic or use a Poisson approximation as a hypothesis.

## 2. Model and theorem

Let G be any simple bipartite 2-regular graph with n vertices on each shore.
Assign independent uniform bijections rho and sigma from its two shores to
{0,...,n-1}, and put B={(rho(u),sigma(v)):uv in E(G)}. Thus B contains
exactly two points in each row and column. Collinearity throughout means
zero **integer** determinant, not congruence modulo a prime.

For 1-n <= d <= n-1 let

    D_d={(r,c):c-r=d},  L_d=n-|d|,
    Y_d=|B intersect D_d|,  Z_d=binom(Y_d,3),  S=sum_d Z_d.

Let X count all nonaxis collinear triples. Distinct slope-one diagonals are
disjoint, so S counts actual triples without multiplicity and 0<=S<=X.

**Theorem.** For every such G and n>=3,

    E_G S = (2n-5)/3.                                      (1)

Uniformly over all choices of G, as n tends to infinity,

    Var_G S = (40/9)n + O(1),                             (2)
    E_G[(S/n-2/3)^2] = 40/(9n) + O(n^-2),                (3)
    limsup_n n sup_G P_G(X=0) <= 10.                     (4)

Uniformity in (2) means that an absolute constant bounds the remainder
for every sufficiently large n and every cycle type, including a number
of C4 components proportional to n. Consequently (4) also holds for any
mixture of the cycle types, in particular the uniform distribution on all
simple saturated boards. The fraction of those boards having no three in
line is O(1/n). This does not rule out rare successful labelings at any n.
The number 10 is the coefficient supplied by this variance argument; no
optimality claim is made.

## 3. The uniform colored-overlap lemma

Let T_1,...,T_K be pairwise edge-disjoint partial matchings in the n by n
grid, and set Y_i=|B intersect T_i|. Fix nonnegative integers k_i, with
k=sum k_i. Expanding product_i (Y_i)_(k_i) chooses ordered, distinct edges
from each T_i; disjointness ensures all k grid cells are distinct.

Partition their row positions by equality and separately their column
positions. The resulting simple bipartite graph H has k labeled edges,
r row vertices, s column vertices, v=r+s vertices and c components; it has
no isolated vertices. For fixed k only finitely many such patterns exist.
For each pattern the probability that its specified grid copy lies in B is

    inj_shore(H,G) / [(n)_r (n)_s].                       (5)

This is non-induced containment. Inverse row and column labels are uniform
injections, proving (5); there is no additional automorphism factor.

There are at most n^c placements of this fixed colored pattern into the
T_i: in each component choose one root coordinate; following an edge of
its specified matching color determines the neighbor uniquely. Extra
equations, including a cycle closure, can only reject that choice.
There are at most n^c 2^(v-c) homomorphisms of H into G: choose each root
and follow a spanning tree, with at most two choices at each new vertex.
Injective maps form a subset. For n>=2k, the contribution of this pattern
to the factorial moment is therefore O_k(n^(2c-v)), uniformly in G and
the T_i. The same argument works for fixed maximum degree Delta, replacing
2 by Delta.

Every component has at least two vertices. The exponent is zero only
when every component is a single edge. It equals -1 only for one two-edge
path together with isolated edges. Every other pattern contributes
O_k(n^-2). In particular a cycle component has at least four vertices,
so is already in this remainder: cycles have not been omitted or treated
as trees. Patterns forbidden by the degree-two condition contribute zero.

The all-single-edge pattern also gives the useful uniform moment statement

    E product_i (Y_i)_(k_i)
       = product_i (2|T_i|/n)^(k_i) + O_k(n^-1).         (6)

Indeed, nonmatching target tuples number O_k(n^(k-1)), while the number of
ordered k-matchings of G is (2n)^k+O_k(n^(k-1)); denominator (n)_k^2
completes the normalization. The more precise first correction follows.

## 4. The first covariance correction, retaining the two-edge path

Let p_k be the probability that a specified grid k-matching is contained
in B. Then p_k=k!m_k(G)/(n)_k^2, where m_k counts edge matchings. Uniformly
over G, for every fixed k,

    p_k=(2/n)^k [1+k(k-1)/(4n)+O_k(n^-2)].              (7)

To prove this without a cycle-profile assumption, start with all (2n)^k
ordered edge tuples. For a selected pair of positions, exactly 6n ordered
edge choices intersect: each first edge has itself and its two neighboring
edges as the three possibilities. Thus the first exclusion term is
binom(k,2)*6n*(2n)^(k-2). Intersections of two distinct pair constraints are
O_k(n^(k-2)): their constraint graph on positions has at most k-2 connected
components, and the edge-intersection relation has bounded degree 3.
The first two Bonferroni inequalities give the stated uniform remainder.
Dividing by
(n)_k^2=n^(2k)[1-k(k-1)/n+O_k(n^-2)] gives (7).

Let q_k be the probability for a specified two-edge path plus k-2 disjoint
edges, all vertices otherwise distinct. It is the same to leading order
for either shore as the path's center. Choosing the center and its ordered
neighbors gives exactly 2n possibilities; each other edge gives 2n.
Coincidences between these components cost one free root and contribute
O_k(n^(k-2)). The denominator is (n)_(k-1)(n)_k. Hence

    q_k=2^(k-1)n^-k+O_k(n^(-k-1)).                     (8)

Now take two edge-disjoint grid matchings T,U of lengths L,M. Let R and C
be the sizes of their common row sets and common column sets. Fix a,b>=1,
k=a+b, theta=L/n, eta=M/n and omega=(R+C)/n. The number of target ordered
tuples whose union is exactly one two-edge path plus isolated edges is

    D=ab(R+C)L^(a-1)M^(b-1)+O_(a,b)(n^(k-2)).          (9)

For example, specifying the two positions that share a row gives exactly
R(L-1)_(a-1)(M-1)_(b-1) choices before excluding other collisions. The
column version contributes C. Two simultaneous overlap constraints cost
two free choices, giving the remainder. Sharing both row and column with
the same pair is impossible because T and U have disjoint cells.

The all-single-edge count is (L)_a(M)_b-D+O(n^(k-2)). The remaining patterns
contribute O(n^-2) by Section 3, even if their embedding probability is
larger than the matching probability. Combining (7)-(9) and subtracting
E(Y_T)_a E(Y_U)_b gives

    Cov((Y_T)_a,(Y_U)_b)
      = [ab*2^(a+b-1)/n]
        [theta^a eta^b - theta^(a-1)eta^(b-1)omega]
        + O_(a,b)(n^-2).                               (10)

All errors are uniform even when L or M is small, zero, or depends on n.
Only bounded degrees and a fixed number a+b of chosen positions were used.
This uniformity licenses summing (10) over a growing number of line pairs;
a fixed-K limit alone would not license that step.

## 5. Summing the slope-one diagonals

The inherited exact matching count is
m_3=2n(n-2)(2n-5)/3. Also

    sum_d binom(L_d,3)=binom(n,3)+2binom(n,4)
                      =n(n-1)^2(n-2)/12.

Multiplication by p_3=4(2n-5)/[n(n-1)^2(n-2)] proves (1).

For a single line of normalized length theta, (6) and

    (Y)_3^2=(Y)_6+9(Y)_5+18(Y)_4+6(Y)_3

give, uniformly in its length,

    Var binom(Y,3)=8theta^5+8theta^4+(4/3)theta^3+O(n^-1).

The sum of these individual variances, divided by n, tends with error
O(n^-1) to

    2 integral_0^1 [8x^5+8x^4+(4/3)x^3] dx = 98/15.   (11)

For distinct diagonals d,e, (10) at a=b=3 becomes

    Cov(Z_d,Z_e)
      =(8/n)[theta_d^3 theta_e^3
          -theta_d^2 theta_e^2 omega_(d,e)]+O(n^-2).   (12)

Writing t=d/n, s=e/n, their row and column interval overlaps are exact:

    omega=2 min(theta_d,theta_e)                 if d,e have the same sign;
    omega=2 max(theta_d+theta_e-1,0)             if they have opposite signs.

Zero endpoints affect only O(n) of the O(n^2) pairs and cause O(1) total
error in the covariance sum. The formulas also agree where a sign is zero.
The bounded piecewise polynomial kernel has a two-dimensional Riemann-sum
error O(n^-1). Including the d=e terms in this kernel changes the answer
by only O(1), so their exclusion is harmless for the leading coefficient.

The covariance contribution divided by n is

    8 [integral_-1^1 (1-|t|)^3 dt]^2
      -32(J_same+J_opposite)
    =2-32(1/14+71/1260)=-94/45,                        (13)

where

    J_same=integral_[0,1]^2 x^2 y^2 min(x,y) dx dy=1/14,
    J_opposite=integral_[0,1]^2 x^2 y^2 (x+y-1)_+ dx dy=71/1260.

There are O(n^2) uniform O(n^-2) remainders in (12) and O(n) uniform
O(n^-1) remainders in (11); both sum to O(1). Thus
98/15-94/45=40/9 proves (2). Equations (1)-(2) prove (3). Finally,

    P(X=0)<=P(S=0)<=Var(S)/(E S)^2
        <=[(40/9)n+C]/[(2n-5)/3]^2,

for one absolute C and all sufficiently large n, proving (4). The negative
aggregate covariance in (13) is essential: deleting it would replace the
coefficient 10 by 147/10. No sign is assigned separately to arbitrary pairs.

## 6. Joint Poisson companion and a repaired conditional shortcut

For any fixed K and edge-disjoint partial grid matchings T_i with
|T_i|/n tending to theta_i, their occupancies converge jointly, uniformly
over cycle types, to independent Poisson variables of means 2theta_i.
This is a conclusion of (6), with the following domination proof.

Condition on rho. For one target matching T, each source column vertex
has at most two allowed labels. Because T is a matching, a particular
column assignment can contribute at most one hit. The count of k
nonattacking allowed placements is at most binom(n,k)2^k, so

    E[(Y_T)_k | rho] <= 2^k,
    E[(1+u)^Y_T | rho] <= exp(2u),  u>=0.             (14)

One must NOT apply the same binary-hit argument to the union of K target
matchings: one column assignment can then contribute several edges.
Instead Holder's inequality gives, for A>=1,

    E[A^(sum_i Y_i)] <= product_i E[A^(K Y_i)]^(1/K)
                      <= exp(2(A^K-1)).              (15)

The expansion of E product_i (1+z_i)^Y_i is absolutely and uniformly
convergent on each bounded complex polydisc: bound its absolute degree-k
terms by E[(sum Y_i)_k] R^k/k!, and use (15) at 1+2R to dominate the tail
by a geometric factor. Termwise use of (6) yields
exp(2 sum_i theta_i z_i). Coefficients of the ordinary probability
generating function converge by integration on a fixed product of circles,
proving the joint Poisson statement, including zero means.

**Exact hostile to the discarded domination.** Let n=7, fix rho=id, and
take G and the target union both equal to the cyclic matchings
T0={(i,i)} and T1={(i,i+1 mod n)}. A uniform sigma gives a weighted
permutation board with weight 2 on its diagonal and weight 1 on the two
neighbor diagonals. Direct counting gives

    E[(Y0+Y1)_2 | rho=id]=18-10/(n-1)=49/3 > 16.

Thus the tempting bound (2K)^k fails already for K=k=2 in this control.
The single-target bound and Holder survive. The main density theorem uses
only the finite-order colored-union proof, not this domination argument.

As another consequence of (6), a fixed number of long parallel diagonals
would already force P(X=0) to tend to zero by taking that fixed number
arbitrarily large after n tends to infinity. That weaker argument alone
does not give (4), nor justify growing K in a Poisson limit.

## 7. All-C4 companion: which block occupancy carries the mean

For n=2m and G=mC4, a k-matching uses j blocks twice and k-2j blocks once.
Counting block choices and the two perfect matchings/four single edges in
each C4 gives exactly

    m_(k,j)=(m)_(k-j) 2^j 4^(k-2j) / [j! (k-2j)!].   (16)

For n>=4 and a prescribed grid triple, the two possible component-occupancy
probabilities are therefore

    P_(2,1)=12/[n(n-2)(n-1)^2],
    P_(1,1,1)=8(n-4)/[n(n-2)(n-1)^2].                 (17)

For n>=4, among the expected full collinear triples their exact fractions are
3/(2n-5) and (2n-8)/(2n-5). By the inherited Guy--Kelly grid-triple count
T3(n)=(3/pi^2)n^4 log n+O(n^4), the two-block contribution has mean
(36/pi^2)log n+O(1); the full mean is (24/pi^2)n log n+O(n).
Thus a two-block-only local approximation discards the leading full mean.
These expectation facts do not enter (1)-(4), and do not themselves give
a zero-event bound. They identify a separate lost coordinate in proposals
that truncate events by how many cycle blocks they meet.

## 8. Exact controls, reproduction, and scope

The standalone verifier uses Python and SymPy, no imports from a producer
inherited from an earlier round. Its exact universes include all 21 cycle
types for 2<=n<=8 (edge subsets versus cycle-rook polynomials through degree
six); full shore-label spaces for C6, C8, 2C4, C10 and C4+C6; exact symbolic
covariance coefficients for 1<=a,b<=3; both overlap integrals; the genuine
falling-factorial square; conditional Holder gates; and all 7! column
permutations in the weighted-union hostile. All checks use integers or
rationals and explicit exceptions that remain active under -O.
An independent selected-edge-subset/type-copy calculation reconstructs
the second moments for every displayed carrier, including shore-preserving
automorphism factors and ordered event-pair multiplicity. It additionally
gives Var S=1335953/88200 for C16 and 228169/14700 for 4C4 at n=8.

| G | labelings | E S | Var S | P(S=0) |
|---|---:|---:|---:|---:|
| C6 | 36 | 1/3 | 2/9 | 2/3 |
| C8 | 576 | 1 | 49/36 | 25/72 |
| 2C4 | 576 | 1 | 19/9 | 5/9 |
| C10 | 14400 | 5/3 | 461/120 | 101/360 |
| C4+C6 | 14400 | 5/3 | 1969/450 | 47/200 |

The different finite variances are compatible with a universal leading
term and a bounded cycle-dependent remainder. At n=4, the C8 board

    {(0,0),(0,2),(1,2),(1,3),(2,1),(2,3),(3,0),(3,1)}

has S=0 but X=1 by a literal determinant calculation. The full n=4 X
histograms independently recover the inherited zero probabilities 1/36
and 1/2, preserving the old hostile. For n=4 the latter means 9 good
distinct 2C4 boards out of 18, or 18 good named-palette pairs out of 36;
these two sample spaces must not be conflated.

Reproduction from this outside staging directory:

    python3 -B 04-computation/overnight8_20260906_no3line_diagonal_density.py
    python3 -B -O 04-computation/overnight8_20260906_no3line_diagonal_density.py

The exact computation verifies controls and finite identities; the uniform
remainders and asymptotic theorem are proved in Sections 3-5, not inferred
from the finite bank. There is no CLT claim, exponential tail claim,
asymptotic independence premise, modular-collinearity substitution, or
claim that O(1/n) is the true order of the successful-board probability.

The mechanism beyond the earlier conditional-cumulant theorem is a
uniform expansion graded by v-2c of an actual colored overlap graph.
An isolated edge gives the independent Poisson term, one two-edge path
gives the covariance, and every further connected identification loses
another power of n. Actual injection denominators make this grading a
probability calculation; the unlabelled formal type product did not.

**Frozen reproduction:** 23,124 active gates. Normal and optimized outputs
are byte-identical, with LF line endings. Producer source SHA256:
`ead04c2f82c7c77a4f94c5b2ca379c7cc56dcdb37e252f4276a3826b2972c200`.
Output SHA256:
`87e151596665aac1ec36d10f348a9555bc68e981d43ad317cb017bfdb431e553`.
Semantic payload SHA256:
`8ffd51855c5588de5979cfc6f8aac3e201e24372ec6dc4a3cdb47e4fe8c82264`.

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
