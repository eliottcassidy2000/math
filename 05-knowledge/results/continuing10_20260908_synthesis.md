# Continuing synthesis: actual paths, repair distance, joint moments, and component poles

**CURRENT CHECKPOINT, 2026-09-08. Status: PROVED scoped results with
FINITE-EXACT certificates and independent audits.** General LRC(14),
planar JC, the original-response sign problem, and the untouched wall-core
no-return problem remain **OPEN**. No external priority claim is made.

The [previous synthesis](continuing9_20260907_synthesis.md) records the
fixed-moment, ratio-tree and genus inputs. The first two proof checkpoints
of this session are `4441760e1` and `e291fcc94`; `0c099746a` publishes the
complete clock and all-m chart theorem after incoming research through
`1d64840fe`. The [manifest](continuing10_20260908_manifest.json)
pins the exact artifacts and normal/optimized verification. Older frozen
arrays, sources and certificates are retained as provenance.

## The portfolio and inherited objects

The anchor was the fifteen connected-complement clock-7200 survivors.
The niche was the joint Newton-circuit/interlacer image. The wildcard
revisited the no-three-in-line paper and global polynomial carriers.
Incoming work on labelled component responses then supplied a further
geometric connection. Five concepts stayed live: actual path composition,
distance to a feasible configuration, joint moment positivity, complete
special-fibre poles, and the exact category of a proposed operation.

The closest proved mechanisms were the complete 126-profile minimum-tree
consumer, fixed-first-two-moment circuit surjectivity, row-uniform diagonal
defect concentration, and rational integration on a fixed invariant fibre.
The hostile tree was jointly realizable and had a common phase attaining
its edge minima; C-only circuit feasibility did not pay the D moments;
single-direction diagonal defect missed independent repair costs; and a
globally smooth function with generic affine-line fibres still lacked a
global conjugate. The underused sidecars were actual endpoint pairs,
capacity-one cell edges, the next moment row, and separately labelled
special-fibre components.

## 1. LRC: a complete clock closure and twenty further exclusions

**PROVED.** If thirteen distinct positive integer speeds form a primitive
row, some selected six-label set A has gcd(A)=7200, and the **actual**
strict-atlas graph on its complementary seven labels is connected, then
the full row has a time with weak clearance at least 1/14.
The [complete proof composition](continuing10_20260907_lrc_clock7200.md)
and [final audit](continuing10_20260907_lrc_last_b_audit.md) state the
exact atlas and physical lift. There is no height box, unit speed,
connectivity on A, decoder-span equality, or unit-gcd assumption on B.

The inherited minimum-tree computation covers all 76,814 admissible words
and leaves fifteen. Complete native wedges with margins (4,24,18),
(16,24,18), and (24,18,30), combined with connectedness, remove thirteen.
The remaining words are (1,9,16,18,24,32,60), with excess 116, and
(5,8,9,30,32,36,48), with excess 103. Forced edges of credit 114 and 102
leave only a two- or one-point budget. Sheet divisibility forces the next
arms to have zero credit. All 70 products for the first word give endpoint
credit 142; all 180 and 918 products for the two second-word paths give
at least 112 and 140. The partition 13+1+1 exhausts the actual residual.

The substantive new operation was to alternate **native path closure and
forced graph structure**. The initial relaxation erased actual rational
edge labels. Recovering the labels gave actual endpoint pairs; forbidding
their bad wedges forced more of the graph; the forced graph then made
the remaining native alphabets much smaller. The proof does not require
the relaxed pair minima to be jointly attained.

The [bounded upper-clock continuation](continuing10_20260908_lrc_upper_twenty.md)
uses the existing minimum-tree theorem to exclude exactly the twenty
largest remaining clocks. Its [independent C++ audit](continuing10_20260908_lrc_upper_twenty_audit.md)
checks 183,164 unpruned multisets, all 14,890 word/clock evaluations,
614 complete pair tables and 661,538 compatible atlas edges. The least
margin is 133, at clock 11952 and word (1,6,8,8,9,36,48).

The necessary array in this connected-complement reduction now has
**7,625 clocks, maximum 11,935**, with canonical JSON SHA256
`51c00db52c18732ec6128e4d9e1ed197813ab107bfd3dec7dc7ace9c5c3c67c7`.
Only 7200 and those twenty clocks were removed in this session.
The inherited bound gcd(B)<=90 retains its reduction scope.
Arbitrary entry and disconnected complements remain **OPEN**.

**Failure boundary.** The same ratio 578:801 gives credit 114 on the
effective 3600-grid with sheet 2, but 112 on the effective 1800-grid with
sheet 4. A scale transfer must recompute its quotient clock and saturated
arm alphabet. Minimum trees use possible edges; a maximum tree is allowed
only after its edges are known to be actual. A partial minimum tree is
not a monotone pruning statistic when another vertex is added.

## 2. No-three-in-line: a repair-distance obstruction

The [two-direction repair theorem](continuing10_20260907_no3_repair.md)
returns to the independence assumption in the original
[Guy--Kelly paper](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B126DA7E4957722BAC70AC7B7F6E1FA2/S0008439500056770a.pdf/nothreeinline_problem.pdf).
Its new argument uses a uniform column permutation of **any fixed** simple
2-regular row/column skeleton and any fixed row order. It does not assume
that the diagonal occupancies are independent.

Rotate a cell (r,c) to an edge joining c-r to c+r. Give every diagonal
vertex capacity two and every cell edge capacity one. If B is the board,
the exact minimum number of deletions clearing these two directions is

    tau(B) = |B| - maximum integral flow.

This is a distance to a necessary feasible set. Every final no-three-in-line
configuration T satisfies |B minus T| >= tau(B), including arbitrary
adaptive repair procedures with insertions and temporary rearrangements.
The claim bounds original-point loss, not algorithmic running time.

The strengthened mean bound is

    liminf E[tau(B)]/n >= gamma,
    gamma = 1 - 5 exp(-2) + exp(-2)/12 - exp(-4)/2
            + 5 exp(-6)/4 + exp(-8)/6
          = 0.3285980553011038... .

It is uniform over the fixed skeleton and row order. The additional term
counts antidiagonal triples whose points have individually safe opposite
diagonals, using bounded-degree permutation avoidance with the actual
forbidden cells retained. The transposition change bound then gives
the probability rate gamma^2/8 = 0.01349708524345841... and exponentially
small probabilities for repairs retaining all but kappa*n original
points when kappa<gamma. Explicit finite-n mean bounds and tail formulas
are in the proof. The audit uses independent flow/deletion and exact
permutation controls; the complete board universe for n=2 through 6 contains
70,087 boards.

**Failure boundary.** Clearing these two directions is necessary, not
sufficient, for full no-three-in-line. The rate is exponential in n,
whereas the number of all row/column boards grows on the n log n scale;
this is no extremal nonexistence theorem. Allowing cell capacity above
one would compute the wrong repair problem.

## 3. Newton circuits: the full closed nonpositive octant is excluded

The earlier fixed-moment result realizes every circuit sign word; all
27 anchored C-only words survive. Restoring the D moment row changes the
answer. The [joint-moment theorem](continuing10_20260907_nonpositive_circuits.md)
proves that there is no positive-product point in the native anchored
B/C/D model with all three circuit ratios at most one and both J3>=0
and J4>=0 for D. This includes ties and the whole unbounded z tail.

J3 confines the feasible x,y region to a compact polynomial chart.
At its minimal allowed z, exact tensor Bernstein coefficients prove
J4<0 and its z derivative negative; its second derivative is -6.
The compact certificate therefore controls every larger z. There are
133 and 40 strictly positive coefficients in the two required tensors,
independently recovered by the auditor.

This is stronger than a finite scan or the earlier open-neighbourhood
exclusion. Every full anchored positive-product model point must have
at least one circuit ratio above one. Two-negative configurations still
survive the degree-six D conditions: (95,68,1) and (86,50,9) are exact
controls, and both fail the next degree-eight D matrix. The full
two-negative classification, original-response transport and wall core
remain **OPEN**. A prose coefficient in the draft sixth moment was
repaired to -915y/7 before promotion; the executable and determinants
already used that value.

## 4. Global carriers: smooth fibres still retain component obstructions

The [complete source-linear classification](continuing10_20260907_dg_linear_carrier.md)
describes every global function on the graph-complement surface W of
source t-degree at most one. For F=A(x)t+B(x), A!=0, every polynomial g with
J(F,g) in C[x] is h(F)+q(x), and the image intersection is A*C[x].
This all-degree descent excludes a global conjugate for every nonconstant
member of the source-linear layer.

The stronger hostile is the globally critical-free fourfold-root family

    F=a(x-h)^4*t + a(x^2-4hx) + B0,   a*h != 0.

Its generic fibres on W are affine lines and it has the rational primitive
1/[3a(x-h)^3]. Yet the exceptional fibre has two disjoint reduced components.
Cancelling the primitive's pole on one component by a function of F
necessarily creates a pole on the other. For the distinct-double-root
family, a nonzero residue already prevents rational integration.

Incoming research now excludes every nonconstant global output-pencil
member of source degree at most three, using the full quadratic and
cubic canon. Thus “try degree two” is superseded as a search direction.
The source-linear theorem remains a useful explicit classification and
the smooth fourfold example retains its component obstruction.

### A new all-m realization: repair order, chart change, and ramification

The incoming component-response mechanism turns that hostile into an
[unbounded family on the actual surfaces W_m](continuing10_20260908_dg_unbounded_torsion.md),
accepted by an [independent analytic and exact audit](continuing10_20260908_all_m_torsion_audit.md).
For each m>=1 and a*h!=0, let Q_m be the complete polynomial part of
(x-h)^(2m)/x^m and set

    F=a(x-h)^(2m)*t+a*Q_m(x)+B0,   f0=F(h,t),   e=2m-1.

F is globally smooth on W_m. Its generic fibres there are affine lines;
the special fibre consists of two disjoint reduced affine lines. In the
original source module C[x,t]/D_F C[x,t], the unit class has **complete
annihilator ideal ((F-f0)^e)**. Its j-th canonical connection derivative
has annihilator ((F-f0)^(e+j)). There is no bound on a prospective mate's
degree in this conclusion. Unboundedness varies m, not just F on W2.

The geometric mechanism is explicit. If E1 is the component x=h and
u=1/(x-h), then W_m minus E1 is another actual affine plane with
coordinates (F,u). The rational primitive is G0=u^e/(a*e), so (F,G0)
becomes a finite flat power map of degree e on this new chart. The old
boundary D is u=0, with ramification index e when m>=2. On the original
source, however, G0 has its pole on E1. The two affine charts exchange
exactly the divisors responsible for the obstruction. The actual form
is omega=a^-1*u^(e-1)*dF wedge du; dropping it would confuse this
ramified power map with a polynomial Jacobian-one pair.

The equality between repair exponent and power degree has a sharp
boundary. At h=0 and m>=2 the power degree is still 2m-1, but the omitted
fibre component has multiplicity m, F is critical, and the exact unit
annihilator drops to ((F-B0)^2). Component multiplicity is indispensable.
At m=1 the admitted family has exponent one and an unramified alternative
chart map. These controls explain the mechanism, rather than merely
checking a finite list of examples. The [incoming connection review](continuing10_20260908_incoming_connections.md)
records the source, target, maps and stopping boundaries of the other bridges.

### The full source-linear classification and an order-two gap

The [complete classification for every m](continuing10_20260908_dg_linear_all_m.md)
now closes the surrounding rationally integrable, globally smooth linear
layer. Its [independent audit](continuing10_20260908_linear_classification_audit.md)
proves necessity at every source degree and checks all parameter boundaries.
A reciprocal polynomial has a rational primitive only when it is constant
or a single pure power: the orders of the primitive's finite poles cannot
pay its zero at infinity if there is more than one distinct root.

For m>=2, every admitted nonconstant first function has A=a(x-h)^n,
a*h!=0, and B equal to the polynomial part of A/x^m up to a free
constant. The exact allowed exponents are m+1<=n<=2m-2 or n=2m.
The missing n=2m-1 has one unavoidable critical point on the boundary.
The positive unit orders on each such surface are

    {m,...,2m-3} union {2m-1},

with an empty interval omitted. Across all m the positive **unit** orders
are exactly {1,3,4,5,...}; order two is absent in this classified layer.
This does not exclude order two for other response classes, connection
derivatives, or higher source-degree functions. For n<2m the generic
global fibres are punctured lines and a separate special fibre contains D;
the affine-line fibre theorem above concerns n=2m.

There is also an actual map exhibiting loss under restriction. On W1,
F=a*t+B0 has zero unit class in the original C[x,t] response module,
but its unit has exact order one in O(W1)/D_F O(W1). The global witness
is -x*t=1+r*b. Restriction preserves the derivation while killing this
nonzero global class because it removes the pole component D. This is
a precise reason to retain the coefficient ring when importing a
component-response obstruction.

## 5. Connections that survive the incoming work

| Source and target | Map and preserved predicate | Lost information, required sidecar, next test |
|---|---|---|
| Actual LRC graph -> endpoint pair | Compose native rational edge labels; preserve an actual pair and uniform overlap | Margins alone lose prime valuations and quotient clock. Recompute both before testing another surviving clock. |
| Diagonal occupancy -> repair distance | Capacity-two matching on unit-capacity cell edges; preserve the necessary feasible set | A count loses shared-cell competition. Keep exact edge capacities; test additional slope families against the existing strict flow duality gap. |
| Circuit signs -> joint interlacer feasibility | Retain the actual factorial D moments and their principal minors | C-only feasibility loses the next moment row. Attack the two-negative cells at degree eight; do not infer their emptiness from the present octant theorem. |
| Rational primitive -> polynomial response | Compare principal parts on separately labelled fibre components modulo a common function of F | Generic fibre type loses exceptional components. Supply an actual derivation-preserving map before transferring the obstruction to moving-source responses. |
| Genus obstruction -> another operation | Retain both the invariant and the commuting vector field | A rational primitive is not a commuting scalar-time selfmap. Incoming high-genus rational mates are compulsory hostiles to that transfer. |

The incoming [full cusp-ideal theorem](planar_jc48_sep08_cusp_ideal.md)
closes nonzero rational scalar times for every invariant in
C+Delta*C[p,y]. Earlier coefficient-by-coefficient genus extensions are
therefore saturated. The genus-one case needs the actual pole divisor
of its vector field. A composition using different invariants remains
outside that theorem; inverse cancellation is its first hostile test.

The incoming [leading exactness theorem](planar_jc48_sep08_leading_exactness.md)
places dx/A^(1/n) in the actual radical coefficient field for a rational
mate of A(x)t^n+... . It generalizes the source-linear residue test but
discards lower coefficients. The [complete quadratic pencil classification](planar_jc48_sep08_exact_pencils.md)
and current [planar board](planar_jc48_sep06_board.md) retain those lower
rows. The subsequently audited 6+2 closure leaves only boundary partitions
7+1 and 5+2+1 in the proved DG quartic table. Pending further audits are
not dependencies. A leading exact differential alone does not close these
lower-section problems.

The [new incoming review](continuing10_20260908_incoming_tail.md) also
records the stronger all-finite 7+1 hostile: every formal inverse row is
exact and a formal Laurent mate exists, but no rational mate exists.
Normalized trace rules out an algebraic mate as well. Extending the same
formal coefficient hierarchy is saturated on that locus; rational
realization or the actual compact-fibre differential must be restored.

The incoming [THM-4458, one-sided adverse-leak budget](../../01-canon/theorems/THM-4458-lrc-one-sided-adverse-leak-budget.md)
offers a concrete cross-thread test. For the same real observer R and
physical target C=R+L, charge only adverse signed products of L against
r=R-mean(R), minimizing the irrelevant constant gauge by a weighted
median. An upper bound B_R(L)<=b gives Var(C)>=(E-b)_+^2/E, E=mean(r^2).
Our actual forest profiles may supply information about L only after an
explicit map to the present-Q owner/word row is written. A nonpositive
forest lower profile still allows C=0 and exact cancellation B_R=E, so
the variance theorem cannot create entry from that lower bound alone.
Favorable error is free; changing the observer or squaring a loose upper
bound above E is not allowed. This is a new decisive test, not an
excluded LRC profile or an imported PDE theorem.

The September 8 affine-pole correction reinforces the same mechanism:
repairing a pole at a boundary point can create poles at ordinary affine
points over the same denominator zero. The next candidate must retain
every such point before using a global pole-space dimension count.

## Next decisive tests and stopping rules

1. For another LRC clock, first run the existing complete minimum-tree
   consumer. If it survives, recover a native low-credit graph and alternate
   endpoint closure with forced topology. Stop scale analogies at the first
   changed quotient or saturated alphabet; disconnected entry needs its
   own theorem.
2. For no-three-in-line, use the exact repair dual rather than simply
   summing incompatible directional defects. A third slope is useful only
   if its shared-point cost survives a cheap exhaustive hostile.
3. For the anchored two-negative circuit cells, retain the degree-eight
   D moment condition. Try a finite chart plus concave-tail certificate
   only after proving the chart covers the entire proposed region.
4. For quartic global carriers, use the two surviving weighted boundary
   partitions and the actual nonconstant lower section. Any rational
   primitive must pass every affine and boundary pole, not only the
   leading radical differential.

These are distinct research obligations. None supplies the missing entry
for a general conjecture merely by analogy with a successful local result.
