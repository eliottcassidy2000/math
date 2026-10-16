# Native survivor packets and a smaller residual invoice

2026-10-15. **PROVED:** the conformal packet construction, its minimum
rowwise residual, and the conditional iteration/energy implications below.
**FINITE-EXACT:** the declared first-reset controllers and hostile schedule.
**UPDATE:** the [extended census](collatz_reset_debt_20261015f.md) refutes the
unrestricted all-scale one-reset `7X` bound at `k=23`; discounted and paired
repairs pass the declared `k=8..23` census. **OPEN:** the required `k>=38`
tail, unconditional amplitude nonvanishing, and universal ROOT entry.
Self-audited; no independent referee
or historical novelty is claimed for signed matching or norm inequalities.

[Program](../../04-computation/experiments/collatz_survivor_packets_20261015e.py)
and [exact output](collatz_survivor_packets_20261015e.json).

## 1. Inheritance and the changed question

The closest proved mechanism is the original-source signed-fibre invoice in
[Gram/slack](collatz_gram_slack_20261015b.md), followed by the first odd-ray
half-payment in [O2, odd-ray overlap](collatz_odd_ray_overlap_20261015d.md).
The canonical hostile is that note's second-block sign reversal: the same
terminal column changes sign between physical rows. Its orthogonal mod-3
channel also has zero first-block payment. The corrected near miss is to
iterate the first cone contraction without transporting its hypotheses.
The least-used coordinate is **which literal source mass remains unmatched
inside each signed terminal fibre**, rather than just the terminal's sign.

Anchor: a native packet with a small residual. Niche: first-reset fibres and
their exact 18-step history-preserving lift. Wildcard: the critical boundary
of the incoming infinite-depth Fourier operator, developed in
[the Wiener continuation, W2](collatz_conductor_wiener_20261015d.md#5-the-critical-boundary-pays-one-conductor-moment).
The live board is original row weights, ordered native clocks, unmatched
mass, reset length, conductor tail, and cross-stage debt. The method cards
used are **Compute the repair quotient before testing the residual defect**
and **Find the hidden second coordinate in a nearly true theorem**.

Incoming commits b56bef2502 and e62673a9f7 supply the completed 2^40 census,
the failed quiet-phase prediction, a larger *conditional* amplitude budget,
and explicit infinite-depth eigen-distributions. Their source/status and
exact output were read, not independently recomputed at that large scale.
No stationarity or random-sign assumption enters the constructions here.
The subsequent audit commit 0015ac67b8 was also read: MISTAKE-621 removes
certificate nesting as evidence for the energy premise and clarifies the
finite-depth spectral scope. Neither interpretation is used below; the
energy seeds and safe conditional budget remain the stated inputs.

The important change is that a contracting packet need not be another copy
of the original annular cone. It can be the mass already authenticated to
cancel, together with a precisely limited amount of unmatched mass. This
gives an explicit packet for **every finite native controller**, with an
explicit residual price. What remains open is an all-scale bound on that
price, not the existence or interpretation of the packet.

## 2. A minimum residual on the actual source coordinates

Let a finite native refinement A send each current terminal e to one checked
terminal f with sign s_e in {+1,-1}. Store its literal integer path length;
ROOT is marked and fixed. This section excludes zero-deletion columns for
simplicity. The tested graphs have no nonROOT cycles. A separately verified
zero-deletion extension can only improve the resulting norm bound.

For one real row b write u=||b||_1 and v=||bA||_1, so 0<=v<=u. A split
b=c+r is **conformal** if c_e and r_e have the sign of b_e (or vanish).
The original coefficient at every coordinate remains exactly c_e+r_e.
This permits fractional allocation of a source weight, not a change in that
weight, its source, its endpoint, or its clock.

**Proposition P1 (PROVED).** For 0<rho<1 there is a conformal split with

    ||r||_1=||rA||_1=d_rho(u,v):=(v-rho*u)_+/(1-rho),
    ||c||_1=u-d_rho(u,v),
    ||cA||_1=v-d_rho(u,v)<=rho*||c||_1.                 (1)

Moreover d_rho(u,v) is the smallest possible residual norm among conformal
splits satisfying this rowwise contraction condition.

Proof. In each target fibre, multiply coefficients by their checked signs,
then match positive and negative mass. Retain every cancelled pair in c.
The unmatched mass totals v and has a single sign in each target fibre.
Select any d=d_rho(u,v) units of that unmatched mass as r. This is possible
because d<=v, a consequence of v<=u. Pull the selected pieces back to their
original coordinates with the original signs. No new edge or guessed label
is used. The three norm identities in (1) follow. If v<=rho*u, choose r=0;
otherwise v-d=rho*(u-d). For minimality, any conformal split with residual
norm r0 has ||c||_1=u-r0 and, by the reverse triangle inequality,
||cA||_1>=v-r0. Hence v-r0<=rho*(u-r0), giving r0>=d. QED.

The matching is deterministic after ordering sources, but its order does
not affect d. Its coordinate list is a receipt; the two scalars u,v suffice
to compute its price. The minimum assertion is **rowwise and conformal**:
it is not a minimum over arbitrary signed decompositions, all common
controllers, or a joint quadratic optimization over multiple rows.

For the rational choice rho=7/10,

    d(u,v)=(10v-7u)_+/3.                               (2)

This is slightly stronger than the earlier 1/sqrt(2) contraction. It makes
the packet, its residual weights, and acceptance checks exactly rational.
An injective signed map has v=u and forces r=b, c=0: adding a name or a
phase register to a terminal is not itself cancellation.

## 3. The residual is no longer amplified by the Pell constant

Keep the original bins and positive denominators d_i. Write

    Q(B)^2=sum_i ||B_i||_1^2/d_i.

Apply P1 to each row, with u_i=||B_i||_1, v_i=||(BA)_i||_1 and residual
r_i=d(u_i,v_i). Then, coordinatewise,

    v_i <= (7/10)u_i+(3/10)r_i,
    Q(BA) <= (7/10)Q(B)+(3/10)Q(R).                    (3)

The second statement is the ordinary weighted Euclidean triangle
inequality on the vector of row norms. Also Q(CA)<=7Q(C)/10. All rows use
the **same native A**; their source allocations may differ.

This improves the earlier generic cone-with-defect estimate because the
residual is a conformal subset of the authenticated survivors. The earlier
estimate allowed cancellations between C and R before applying A and paid
for them twice. No Q--Pell stochastic covariance law is imported here.

For changing controllers A_j, iteration gives the more flexible bound

    Q(B_m) <= (7/10)^m Q(B_0)
       +(3/10) sum_(j=0..m-1) (7/10)^(m-1-j) Q(R_j).   (4)

Thus a large early residual can be paid later. Uniform stage bounds are
sufficient; a bound on this **discounted sum** is also sufficient. This is
an actual debt-transport identity, not an assumption of independence.

**Corollary P2 (PROVED implication; premise OPEN).** For X=2^k, the original
square-root bins, and the halved initial rows, suppose k+8 finite native
refinements have the P1 residual Q(R_j)<=7X at every stage. Then

    Q(B_(k+8)) < (15/2)X,
    ||G_(k+8)||_entry1 <= Q(B_(k+8))^2 < (225/4)X^2.   (5)

Indeed Q(B_0)<=8X^(3/2), and

    (7/10)^2<1/2,           8(7/10)^8<1/2.

The initial term in (4), divided by X, is at most
8(7/10)^8((7/10)sqrt(2))^k<1/2. The residual geometric sum is at most 7X.
The earlier generic estimate was 83X/2; the new one is 15X/2.

If these controllers exist at every k>=38, the inherited fine-bin cost 4/3
gives D_X<(691/12)X^2. Using statement 3 of
[THM-4658, two-octave-memory-energy-descent-critical-lines](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md),
this implies E_X<=60X^2 for every k>=38. The incoming exact seeds are

    E_(2^38)=70317261195513340128655,
    E_(2^39)=23439502132608158412453.

Both are below 60X^2. The induction seam is already paid when X>=512:
sqrt(3/X)<=(sqrt(60)-sqrt(691/12))/2, using
3*48^2*60/29^2<512. This energy statement remains conditional.

The new incoming certificate has certified budget margin
14348.206382850238 X^2 beyond 2^40 (the reported bisection value, not an
exact optimal constant). A less demanding sufficient target therefore also
exists: **replace 7X by 119X** in P2. It gives coarse budget (239/2)^2 X^2,
fine budget (171379/12)X^2, and E_X<=14300X^2 by the same induction from
k=38. Here 14300-171379/12=221/12, and the seam holds for X>=8192 since
3*48^2*14300/221^2<8192. The safe constant 14300 is below the incoming
certified margin. This enlarges the allowable residual without asserting
that either 7X or 119X has been proved at all scales. Amplitude nonvanishing
would still not be universal ROOT coverage.

## 4. A universally defined reset controller, and its exact finite test

Start with original sources 1..8X and every halving edge in that interval.
Close the observed graph, retaining source, endpoint, sign, and integer
clock. At each stage, from every current nonROOT terminal follow actual odd
blocks until the first valuation a=v_2(3n+1)>=2, or ROOT. Retain **all**
edges, then close the complete union again. This reuses an old suffix when
a new path enters it; it does not discard its clock.

Every individual reset extension is finite without assuming Collatz. For
an odd x, each valuation-one step satisfies U(x)+1=3(x+1)/2 and decreases
v_2(x+1) by exactly one. After v_2(x+1)-1 such steps the next valuation is
at least two. This is a local termination counter, not a rank decreasing
across successive resets. Every finite number of these stages is therefore
defined for every positive source.

For each row record u_i before the stage and v_i after it. The remaining
all-scale question is the explicit, original-weight inequality

    sum_i (10v_i-7u_i)_+^2 <= 441*b*X^2,                 (6)

where b=2^floor(k/2). This is exactly Q(R)^2<=49X^2, not a proxy counting
merges or unresolved sources. If Q(B)<7X, any remaining stages can be
identities: choose C=0 and R=B. Their clocks are zero and their residuals
are not omitted from the bound.

**FINITE-EXACT universe:** every k=8,9,...,18, every source 1..8*2^k, all
original bins, and reset stages until Q<7X or k+8 stages. Each of these
eleven cases reaches Q<7X, in 1..8 nonidentity stages. In every measured
nonidentity stage the residual is below 7X. Identity padding supplies the
remaining stages with the same source weights and terminal labels. This
certifies the finite cases of (5); it does not extrapolate (6).

At k=18 the reset controller needs eight nonidentity stages, and its final
Q/X is about 6.70507. Its largest nonidentity residual is about 2.35911X.
The saved integer numerators, not these decimal displays, decide acceptance.
There are 2,397 terminal columns, 15,890 unresolved original sources, and
3,295,799 observed edges; the largest compiled native clock is 300. Thus
this small measurement certificate deliberately leaves ROOT membership
unresolved for some of its sources. Eight reset stages are not eight
integer Collatz steps.

The hostile uses the **same original sources and closure rule**, but just
one odd block from each boundary terminal per stage. At k=18 its stage-three
residual already exceeds 18X. Thus one-block extension does not satisfy the
7X target, even after full observed-graph glue. This is a failed schedule,
not a counterexample to universal Collatz or to the existence of better
native packets. Its whole k+8-stage history also tests the discounted rule
(4), separately from a uniform per-stage rule.
In fact its certified discounted bound is below 5.00459X, although some
individual residuals exceed 18.66X. Its actual final Q is below 4.79743X.
So the discounted alternative passes this finite case that the uniform
7X test rejects. The code rounds each residual norm upward to an integer
before forming that exact rational bound; no floating square root decides it.

Validation includes 1,875 exhaustive four-coordinate rational matching
checks; an injective hostile; independent literal c_I(n) versus event-sweep
rows and full packet reconstruction at k=8,9; authentication of every added
native edge; exact local clock recurrences at every compiled source; and
independent literal T replays of selected clocks. Sums of squares use Python
integers. An explicit 32bX bound protects the int64 event sweep. No eventual
ROOT label oracle, fitted sign, random independence, or presumed convergence
is used in constructing or accepting a controller.

## 5. Recovering 7/8 and 18/19 without losing their predicates

The [natural-number index](natural_number_roles_20261004.md) already separates
the many appearances of these numbers. The strongest relevant inheritance
is [atom-route memory, section 4](collatz_atom_route_memory_20261003.md),
itself crediting the earlier finite-prefix realization and trunk isometry.
It already proves the general internal lift, including

    4n+1,       64n+35,       262144n+184471.             (7)

We use these as reset-fibre coordinates, not as new discoveries. For p
initial valuation-one steps followed by a reset of valuation a>=2, the
endpoint is

    z=(3^(p+1)n+3^(p+1)-2^(p+1))/2^(p+a).

Keeping that head and endpoint, increase a by h=2*3^p. The source changes by

    n' = 2^h n +(2^h-1)(1-(2/3)^(p+1)).                (8)

The multiplier order ord_(3^(p+1))(2)=2*3^p makes the intercept integral.
The old source is recovered exactly from the affine lift; the endpoint and
the first p valuations are unchanged. The native clock grows by h, which
is even, so its sign is unchanged. At p=2 the anchor is -19/27 and h=18.
This is a concrete, history-preserving connection between 18 and 19.

The older [cubic field tower](cyclotomic_depth_towers_20261004.md) has a
different map using the **same exponent-order law**: Frobenius on primitive
27th roots gives degree 18, while

    2^18-1=27*7*19*73,         2^9+1=27*19.

The cubic-root-plus-one generator misses precisely the order-19 factor at
this level; its additional phase register restores the missing coordinate.
That finite-field norm is not a native Collatz cancellation map. Its useful
instruction here is to retain the missing fibre coordinate instead of
declaring the old generator complete.

For 7/8, the genuine inherited odd-ray payment uses a factor-four dilation
and a three-phase mod-3 weight: 1-4^(-3)=63/64. Separately F8 has seven
nonzero elements and lies in F64. These explain the common arithmetic
factor 7 in two specified constructions. The allowance **7X** and the
**eight spare stages** in P2 are convenient sufficient constants, not
invariants forced by those field orders. The choice 7/10 is rational and
satisfies its explicit two inequalities; it does not assert a new
finite-field interpretation of the norm.

The lift (8) carries a limitation directly relevant to the target. Sources
on one of its rays have the **same** native sign at their common endpoint.
Their collision does not pay unless their original weighted row
coefficients have opposite signs. The missing sidecar is precisely the
unmatched mass recorded by P1. This is why the larger reset controller is
tested against (6) rather than against a count of ray coincidences.

## 6. Remaining proof targets

The original concrete target is (6) for the explicit first-reset
controller at all k>=38 and all required stages. It remains OPEN; the
[later counterexample at k=23](collatz_reset_debt_20261015f.md) is below
that tail and refutes only its unrestricted all-scale strengthening.
Discounted debt and paired-stage estimates are the preferred repaired
targets. The original controller has no unspecified
packet dictionary or unknown terminal-label choice. A counterexample would
identify an actual row-weight debt, allowing a longer reset batch or a
context-selected native join to be tested on the same source rows.

Two weaker alternatives remain useful: bound the discounted sum in (4),
or use the larger 119X allowance when the goal is the incoming amplitude
certificate. The Fourier continuation supplies strong stability at a
critical regularity boundary and prices its inverse by one conductor
moment; no intertwining map has yet converted that moment into the native
survivor bound (6). That conversion must preserve the marked source and
integer clock, not merely the formal eigen-equation.

Reproduce with `python3 04-computation/experiments/collatz_survivor_packets_20261015e.py`
and the same command with `python3 -O`. Normal and optimized runs have identical SHA256
`041166f15a7a39a3498aebe054630be4b3e69a3beae38fce38f69c1ae2d787d3`.
