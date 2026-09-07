# Two diagonal directions give a repair barrier and a stronger probability rate

Status: **PROVED ANALYTICALLY + INDEPENDENTLY AUDITED**; the stated finite
universes are separate **FINITE-EXACT** controls. No extremal no-three-in-line
conjecture or external priority claim is made.

## Inheritance and the changed object

The closest proved mechanism is the uniform fixed-row mean and permutation
concentration in [overnight11](overnight11_20260906_no3line_rowfreeze.md).
[Overnight12](overnight12_20260906_no3line_count_restart.md) transfers it to
counts and conditionally uniform restarts, but expressly leaves purposeful
within-board changes open. Its marginally-uniform adaptive restart hostile
prevents silently substituting independence. Here the algorithm may inspect
the entire start and choose every later operation adaptively.

The least-used coordinate is the **number of original points retained**.
Rotate the representation, not the random measure: a point (r,c) becomes an
edge from its difference d=c-r to its sum s=c+r. The inverse is
r=(s-d)/2,c=(s+d)/2, so the edge, parity and geometry are retained. This new
bipartite graph is generally not regular. Its capacity-two matching number
has a geometric meaning unavailable from the original cycle type.

Concept board: frozen row labels; exact permutation avoidance; two diagonal
directions; capacitated matching; retained-point distance; adaptive repair.
The Guy--Kelly paper proves the triple-count asymptotic and then explicitly
assumes independence on printed p.530. We use no such assumption. The
[original paper](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B126DA7E4957722BAC70AC7B7F6E1FA2/S0008439500056770a.pdf/nothreeinline_problem.pdf)
was reread for this distinction; its heuristic is not a proved dependency.

## 1. Exact deletion potential

For any finite board B with no repeated cells, define tau(B) as the fewest
points to delete so that both slope +1 and slope -1 lines have occupancy at
most two. Write Y_d for difference-line occupancy and

    F(B)=sum_d (Y_d-2)_+.

The difference/sum graph has one unit-capacity edge per point and capacity
two at every vertex. A maximum integral flow through this bipartite network
retains exactly the largest subset satisfying both diagonal constraints.
Thus tau=|B|-maximum flow. Unit capacity on each original edge is necessary:
cloning both endpoints without an edge-capacity sidecar can retain the same
cell twice. This is the ordinary integral augmenting-path theorem, not a new
matching theorem.

Let Z_s count points on the sum line s whose difference line has Y_d<=2.
Then the following deterministic certificate holds:

    tau(B) >= F(B) + sum_s (Z_s-2)_+.                    (1)

Indeed every successful deletion removes at least F points from the disjoint
difference lines of occupancy >2. Independently, on each sum line at least
(Z_s-2)_+ of its points on the remaining difference lines must be removed.
These two sets of required deletions are disjoint. Equality is not asserted.

In particular let I(B) count sum lines containing exactly three points,
each of which lies on a difference line of occupancy exactly one. Such a
line contributes one to the second term, so tau>=F+I.

If |B|=|B'| and k=|B\B'|, then

    |tau(B)-tau(B')|<=k.                               (2)

An admissible retained subset of B loses at most k points on intersection
with B'; reverse the argument. A column or row transposition of a saturated
board changes at most four points, hence changes tau by at most four.

## 2. Uniform bounded-degree permutation avoidance

**Lemma.** Fix Delta. For every forbidden N by N zero-one matrix with at
most Delta forbidden entries in each row and column, put m equal to its
number of forbidden entries. A uniform permutation avoids them all with
probability

    exp(-m/N)+o(1),                                   (3)

uniformly over these matrices as N tends to infinity.

Proof. If W is the number hit, the number Q_k of ordered k-tuples of
forbidden entries having distinct rows and columns satisfies, for fixed k,

    0 <= m^k-Q_k <= k(k-1) Delta m^(k-1).

For each chosen pair of positions, an entry has at most 2 Delta forbidden
partners sharing a row or column; take a union bound. Interpret k=0,1
separately. Therefore E binom(W,k)=Q_k/[k!(N)_k]
=(m/N)^k/k!+O_(k,Delta)(1/N), uniformly, since 0<=m/N<=Delta.
The odd/even Bonferroni sums bound P(W=0). First fix their order and let
N grow; then let the order grow. The exponential Taylor remainders converge
uniformly on [0,Delta], proving (3). In particular the lower bound is
exp(-Delta)-o(1). No independence between forbidden cells is used.

## 3. An additional linear deletion cost, uniformly after freezing rows

Fix any simple bipartite 2-regular skeleton G on n+n vertices and **any**
row labelling rho. Only the column labelling sigma is uniformly random.
Let alpha=1-5 exp(-2) and

    epsilon=exp(-2)/12-exp(-4)/2+5 exp(-6)/4+exp(-8)/6,
    gamma=alpha+epsilon.

**Theorem.** Uniformly over G and rho,

    liminf_(n->infinity) inf_(G,rho) E_sigma tau(B)/n >= gamma.  (4)

Proof of the new term. On one sum line of length L select three distinct
target rows r1,r2,r3, with target columns ci=s-ri. Each row has two source
column neighbours, giving eight choices of three intended source edges.
Only lines with L>=3 enter this count. Call the row triple eligible if,
for every selected row r and either of
its source-column neighbours with other row t, no other selected row is t
or (r+t)/2. Each r excludes at most four other row addresses. Hence at most
4L(L-2) triples are ineligible (an overcount is harmless); lines with L<3
contribute zero triples.

For every eligible triple all eight choices have distinct source columns.
Their other three forced cells (t,s-r) lie on neither the selected sum line
nor any of the three difference lines: equality with the difference line
through (rj,s-rj) would say r+t=2rj. Let U be the union of these four full
grid lines. Conditional on an intended assignment, the remaining n-3
source columns map uniformly to the remaining n-3 labels. For each source
column, its two row neighbours permit at most eight forbidden labels in U;
each target label likewise forbids at most eight source columns. Lemma (3)
therefore gives conditional probability at least exp(-8)-o(1) that no
additional occupied cell lies in U, uniformly in all choices.

On this event the selected sum line is counted by I. Each intended
assignment has probability 1/(n)_3. The eight assignments are disjoint.
Different triples on the same sum line also give disjoint events because
its occupancy is exactly three. The sum-line lengths are n once and every
1,...,n-1 twice, so

    sum_s binom(L_s,3)=binom(n,3)+2 binom(n,4)=n^4/12+O(n^3).

The total discarded ineligible triples is O(n^3), uniformly in G,rho.
Consequently

    E I >= [8/(n)_3]*(n^4/12+O(n^3))*(exp(-8)-o(1)),
    liminf inf_(G,rho) E I/n >= (2/3)exp(-8).

This already gives a positive improvement. Retaining the actual lengths
in the same forbidden board strengthens it. If the four lengths are
L_s,L_d1,L_d2,L_d3, its forbidden-entry count is at most twice their sum:
each grid cell forbids at most the two source columns adjacent to its row.
With N=n-3, (3) therefore gives the uniform lower bound

    exp[-2(L_s+L_d1+L_d2+L_d3)/(n-3)]-o(1).             (4a)

Overlaps of grid lines or forbidden entries only improve this bound. The
omitted O(n^3) triples have bounded weight, so they do not affect the linear
coefficient after multiplication by 8/(n)_3. The deterministic Riemann sum
of (4a), over every target triple and sum line, can be evaluated explicitly.
For either half of the square put theta=L_s/n. The three rescaled selected
row addresses range over [0,theta], and the corresponding difference-line
length fractions are 1-|theta-2r_j|. Equal row addresses and the unique
middle sum line affect only lower-order terms. Thus

    liminf inf_(G,rho) E I/n
      >= (8/3) exp(-6) integral_0^1 exp(-2theta)
           [integral_0^theta exp(2|theta-2r|) dr]^3 dtheta
       = (1/3) exp(-6) integral_0^1 exp(-2theta)
           (exp(2theta)-1)^3 dtheta
       = epsilon.                                     (4b)

The factor 8/3 is eight edge choices times two halves of the square divided
by 3! for unordered triples. This product integral sums deterministic row
addresses; it does not assume probabilistic independence of line events.
The positive integral proves epsilon>0; exact Taylor enclosures in the
certificate additionally give epsilon>1/200. Combine tau>=F+I with the
inherited **uniform conditional** bound E F>=alpha n-21 to prove (4).
The constant is conservative; neither optimality nor an explicit finite
threshold for (4) is claimed.

## 4. Concentration, adaptive repair, and counts

Put mu_tau(G,rho)=E_sigma tau. Exposing n-1 column labels and coupling two
completions by a transposition gives conditional increment interval length
at most four by (2). The inherited elementary exponential lemma therefore
gives, for every k>=0,

    P(tau<=k | G,rho)
      <= exp(-(mu_tau-k)_+^2/[8(n-1)]).                 (5)

For n>=4 a completely explicit finite version replaces mu_tau by

    b_n=max(0, alpha n-21, 2(n-9)/15+2/[n(n-1)]).

For every fixed 0<=kappa<gamma and every sequence of integer budgets
k_n<=kappa n, the new asymptotic consequence is

    limsup (1/n) log sup_(G,rho) P(tau<=k_n | G,rho)
      <= -(gamma-kappa)^2/8.                           (6)

Let an arbitrary algorithm inspect the full initial board and use arbitrary
randomness, intermediate states, new points, cycle changes, and adaptively
chosen moves. If its final set T has no three collinear points, then

    |B\T| >= tau(B).                                  (7)

Indeed B intersect T satisfies both selected diagonal conditions. Thus (5)
and (6) apply to every procedure allowed to lose at most k original points.
For example q row/column transpositions permit at most 4q original-point
losses; q degree-preserving two-switches permit at most 2q; and q single
point relocations permit at most q. Additional deletions add their count.
These are operation-budget consequences, not running-time bounds for
unrestricted search or impossibility of reaching a solution eventually.

Full saturated no-three-in-line success implies tau=0. Conditional mixing
over any G,rho keeps these bounds when the columns remain uniform. For
uniform uncolored saturated boards, with the exact inherited total B_n,

    limsup (1/n) log(N_n/B_n) <= -gamma^2/8.             (8)

This strictly improves the inherited alpha^2/8 rate by using two directions
without an independence assumption. The scale remains exp(-constant*n),
whereas B_n has log B_n=2n log n-2n+O(log n); no extremal nonexistence follows.

## 5. Exact controls and reproduction

The standalone source reconstructs every **70,087** saturated board for
n=2,...,6, with no geometric or cycle-type prefilter. Its recursion prunes
only impossible residual column capacities. Every board checks the full
integral-flow potential and (1); a separate literal subset-deletion search
checks its exact minimum through n=5. Column transpositions are exhaustive
through n=4 and explicitly sampled thereafter. All 1,152 eligible triples
on cyclic skeletons n=6,...,12 check all eight source assignments, every
forced companion, both forbidden degrees and the retained four-line count.
Thirty explicit forbidden matrices check the collision estimate and every
Bonferroni truncation against complete permutation distributions.

The mean flow deletion costs for n=2,...,6 are exactly

    0, 2/3, 62/45, 28/15, 26438/11325.

An isolated-line positive control at n=5 has rows
(01,34,02,23,14): F=2,G=1,I=1,tau=3. The board
(01,34,24,13,02) has F=1,G=1,tau=3, so (1) need not be equality.
The certificate also retains a literal no-three-in-line failure with tau=0;
two-direction clearance remains a necessary condition. At n=5 there are
92 two-direction-zero boards, versus the inherited 32 full successes.
The one-cell flow control checks the unit edge capacity that endpoint
cloning alone loses. Exact rational Taylor bounds validate the displayed
new constant; the asymptotic theorem relies on the analytic proof above,
not extrapolation from these finite means.

Reproduce with `python -B 04-computation/continuing10_20260907_no3_repair.py`
and the same command with `-O`. The JSON certificate is written beside the
source. All gates raise explicitly under both modes. The independently
written audit reconstructs its own native objects without importing this
producer.

## 6. Decisive next object

The exact flow potential can be much larger than the elementary F+I
certificate. A sharper uniform mean should retain the actual sizes and
intersections of the four-line forbidden board in (3), or a richer flow
dual, rather than multiply marginal diagonal probabilities. The retained
initial points connect any such improvement immediately to purposeful
repair and to zero-event counting. Changing the initial sampling law still
requires its own argument; adaptivity after sampling does not.

The [independent referee](continuing10_20260907_no3_repair_audit.md)
accepts the full analytic proof, including the retained-length integral,
and passes 22,186 exact gates by a separate edge-subset engine. The primary
passes 122,799 gates. Both outputs are raw LF and unchanged under `-O`;
the regenerated JSON is byte-identical.
