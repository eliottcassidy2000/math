---
id: THM-4460
title: "Adverse-budget continuation states and free-shear composition collapse"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.
  The exact fixed-observer continuation quotient, optional partial transport,
  finite scanner, and shear subadditive-collapse obstruction are proved.
  All 4096 auxiliary density activation states separate by depth seven,
  sharply; their current adverse budgets have only 19 values.
  Canonical LRC owner realization and Euler trajectory applications remain OPEN.
source: cross-concepts-20260908
depends_on:
  - THM-4458-lrc-one-sided-adverse-leak-budget
  - THM-4457-euler-sharp-transverse-shear-distance-budget
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-840-hamming-five-continuation-congruence-boundary
  - THM-853-closed-cf-return-semigroup
  - THM-4459-pure-shear-mixtures-have-exact-nuclear-moment-cost
script: 04-computation/cross_concepts_catalysis_20260908.py
output: 05-knowledge/results/cross_concepts_catalysis_20260908.out
script_sha256: b597ca9410c12f0d18d4c43147321bb4b04caad3eb5b4b3e1ecb94591ea4c087
output_sha256: 42215e9a0f4512a4b33a60a28a3012edf8d4a3bb27be2073fb9c509253240ef1
hash_basis: raw LF bytes
audit: >
  Independent root review of optional-transport duality, fixed-observer
  continuation equivalence, scanner edge isolation, p13 normalization,
  complete-mask cancellation, full-future separation, and sharp depth-seven
  boundary. Hinge and transport agree on 7020 profile pairs; Moore and direct
  activation-subset response partitions agree on all 4096 masks. Eight-shear
  reconstruction is checked on 6561 trace-free matrices; the two-shear
  stretching hostile is independently reproduced by THM-4459's script.
  Normal and optimized outputs agree with all explicit checks enabled.
---

# THM-4460 -- continuation states for adverse budgets, and the zero-cost shear trap

**Status: PROVED ELEMENTARY + FINITE-EXACT in the specified censuses;
INDEPENDENTLY AUDITED and promoted on 2026-09-08.** LRC(14),
canonical owner realization, and Euler continuation/singularity applications
remain OPEN. No external literature priority claim or Lean claim is made.

## Inheritance and the two genuine transfers

The recovered mechanism is the universal continuation profile in
[THM-2176, Gordian continuation profile](../../01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md).
It refines a present observable to the coarsest quotient supporting every
named future operation. Its metric companion is
[THM-2191, catalytic localization](../../01-canon/theorems/THM-2191-catalytic-localization-of-the-gordian-metric.md):
the common-context envelope is the largest invariant pseudometric below an
actual nonexpansive metric. These are existing mechanisms, not new theorems
being rediscovered here.

The controls are
[THM-2183, order-join metric product](../../01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md),
where common contexts have no metric contraction, and
[THM-840, Hamming-five continuation boundary](../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md),
where a present statistic fails the operation-kernel test. The least-used
sidecar is the actual operation word, as quantified by
[THM-853, closed CF return monoid](../../01-canon/theorems/THM-853-closed-cf-return-semigroup.md):
144 present classes there become 20,419 future classes. The other THM-853
slug is not used. MISTAKE-249 repairs right-versus-two-sided congruence;
our additions and mask unions are commutative. MISTAKE-547 forbids changing
the observer while keeping its cancellation invoice.

The Anchor is the live adverse-leak gate of
[THM-4458](../../01-canon/theorems/THM-4458-lrc-one-sided-adverse-leak-budget.md).
The Niche is its minimal continuation state and a finite density activation
monoid. The Wildcard is whether treating pure Euler shears as free operations
admits the common-context construction. It fails before that construction's
metric hypothesis is reached.

| Board object | Predicate / operation | Sidecar or obstruction |
|---|---|---|
| One-sided budget | Charge adverse leakage after a fixed observer | Present zero cost is not a removable direction |
| Continuation congruence | Add the same future leak or activate the same label | Visible leak coordinates or full activation mask |
| Partial transport | Match positive-role capacity to negative-role capacity | Retain labels, signs and capacities |
| Catalytic metric envelope | Remove common nonexpansive contexts | A genuine pseudometric is required first |
| Pure shear cone | Add locally shear-like gradients | Cone is not addition closed; subadditive closure collapses |

The first map sends a fixed-role leak vector to its adverse-budget function
on future additions. It preserves every future budget and hence every budget
variance certificate. It discards only a constant on visible coordinates and
coordinates where the fixed centered role vanishes. The needed sidecar is a
visible coordinate class, and, for physical transport, the canonical target/
word realization remains an independent obligation. The cheapest hostile is
`B_r(0)=B_r(r)=0`, split by the common future `-r`.

The second map attempts to send an Euler gradient to the additive cost of
its distance from the pure shear cone. It destroys the obstruction entirely:
every trace-free matrix is a sum of zero-cost pure shears. The corrected
carrier must retain how gradients combine, or fluctuation/mixing data. This
is a matrix statement, not a construction of an Euler trajectory.

## 1. Fixed-observer budget as optional partial transport

Fix `r in R^p`, nonzero and centered. Write `P={i:r_i>0}`,
`N={j:r_j<0}`, `V=P union N`, and `a_i=r_i`, `b_j=-r_j`.
The totals of the two sets of capacities agree. THM-4458 defines

```
B_r(L) = min_b (1/p) sum_s (-r_s(L_s-b))_+.
```

**PROVED (re-expression of THM-4458's dual).**

```
B_r(L) = (1/p) max_f sum_(i in P,j in N) f_ij (L_j-L_i),
f_ij>=0, sum_j f_ij<=a_i, sum_i f_ij<=b_j.               (1)
```

Transport is optional: no prescribed total mass is imposed. To prove (1),
take any feasible dual vector `z` of THM-4458. Its selected positive and
negative masses agree. Coupling those masses gives a feasible `f`; conversely
the marginals of `f`, divided by their capacities, give such a `z`. The
objective agrees because the cost is the difference of two node values.

A direct exact compiler sorts positive-role sources by increasing `L_i`
and negative-role sinks by decreasing `L_j`, fills the best available pair,
and stops once its profit is nonpositive. For a fixed amount of transported
mass the optimum uses the cheapest source masses and most valuable sink
masses; pairing them does not affect the separable objective. The sorted
marginal profit is nonincreasing, proving the stopping rule.

In particular,

```
B_r(L)=0 iff max_(j in N) L_j <= min_(i in P) L_i.       (2)
```

This also follows directly by choosing a constant between these extrema.
The zero set is a favorable cone, not an equivalence class.

## 2. Exact continuation state, and a finite scanner

**PROVED.** For fixed `r`, define `L ~ L'` when

```
B_r(L+K)=B_r(L'+K) for every real vector K.
```

Then

```
L ~ L' iff L-L' is constant on V.                       (3)
```

There is no restriction on their differences outside `V`. Indeed, put
`K=-L` and `K=-L'`. Then both `B_r(L-L')` and `B_r(L'-L)` vanish.
Applying (2) in both directions forces every visible coordinate to agree.
Conversely such differences are a constant gauge on all nonzero summands,
so every translated budget agrees. This explicitly computes the universal
continuation congruence of THM-2176 for the new budget; it does not invent a
new generic congruence theorem.

The symmetric part is a norm on this quotient:

```
N_r(L) = B_r(L)+B_r(-L)
       = min_b (1/p) sum_s |r_s| |L_s-b|,
B_r(L)-B_r(-L) = -(1/p) sum_s r_s L_s.                  (4)
```

The two formulas follow from THM-4458's median identity. Positive
definiteness on the quotient follows from (2); the triangle inequality is
subadditivity of both terms. Thus `d([L],[L'])=N_r(L-L')` is already
translation invariant. Its THM-2191 common-context localization is itself:
common additive contexts create no catalytic contraction in this metric.
Savings in one present scalar budget have a different meaning.

The minimal hostile has `p=2`, `r=(-1,1)`: both `0` and `r` cost zero, but
the common continuation `-r` gives costs `1` and `0`. The displayed exact
audit also uses `p=3`, `r=(-1,0,1)`, giving `2/3` and `0`, to keep a zero-role
coordinate visible in the audit. That coordinate is harmless only for the
fixed observer: for `r'=(-1,1,0)`, the old invisible vector `-e_2` costs
`1/3`. This is the precise boundary of observer changes.

There is an explicit finite scanner on any bounded box. Suppose
`|L_s|<=M` for `s in V` and choose `T>2M`. For an edge `(i,j) in P x N`,
let the future context be

```
K_i=-T, K_j=T,
K_k=2T for k in P\{i}, K_k=-2T for k in N\{j}.
```

Every transport edge except `(i,j)` has negative profit, while that edge
has positive profit. Formula (1) gives

```
B_r(L+K) = min(r_i,-r_j) (2T+L_j-L_i)/p.                (5)
```

Consequently `|V|-1` such queries along a bipartite spanning tree reconstruct
all visible differences. This proves a finite explicit recovery mechanism;
no minimal-query theorem for arbitrary nonlinear encodings is claimed.
The contexts in (5) are unrestricted real additions, not physical density
updates. The next section supplies an independently lawful finite operation
class on auxiliary density indicators.

## 3. A finite physical-density monoid: 19 present states, 4096 future states

Use the exact auxiliary circle densities in THM-4458. Set
`D={y:||y||<1/14}`, `A_s={y:||y-s/13||<1/14}`, and
`R_s=int 1_D 1_(A_s)`. For each labelled subset `F of {1,...,12}`, put
`w_s=1_(A_s)` when `s in F`, and `w_s=1_D` otherwise, including `s=0`.
Then `C_s=int w_s 1_(A_s)`, `L=C-R`, and every density has mass `1/7`.
The lawful operation `T_j` switches density `j` on, i.e. `F -> F union {j}`.
These maps are commuting idempotents; they are not unrestricted additions.

The inherited exact intervals give

```
R_0=1/7, R_1=R_12=6/91, R_s=0 otherwise,
r_0=144/1183, r_1=r_12=53/1183, r_s=-25/1183 otherwise,
E(R)=2508/1399489.
```

Let `a=|F intersection {1,12}|` and `m=|F intersection {2,...,11}|`.
The leak is `1/13` at activated adjacent labels and `1/7` at activated
negative-role labels. **PROVED by the weighted-median formula:**

```
B_r(L_F) = min(325m, 1750-371a+150m)/1399489.            (6)
```

For completeness, the weights at leak values `0,1/13,1/7` are respectively
`500-53a-25m`, `53a`, `25m`, divided by `1183`. Half of the total is
`250/1183`. A minimizing median is therefore `0` or `1/13`; the third
candidate only ties at `m=10`. Evaluating the two hinges yields (6).

Moreover `B_r(L_F)=E(R)` iff `F` is the full mask. The first numerator
`325m` cannot equal `2508`; the second can equal it only when `a=2,m=10`
(check the three possibilities for `a`). That full mask has constant
`C_s=1/7`; every other mask has a nonconstant target and is certified by
THM-4458's exact adverse gate. This is a complete iff only for this finite
auxiliary family, not an iff for general profiles.

The full future congruence is exactly equality of labelled masks. If `F,G`
differ at label `j`, activate every other label. One future state is the
full mask and the other is the mask missing just `j`. Their budgets differ
by the preceding paragraph. Thus the scalar's current permutation symmetry
does not survive individually labelled operations.

**FINITE-EXACT:** in all `2^12=4096` masks there are 19 distinct present
budgets. Equality under all continuations of length at most `d` has counts

```
d       0     1     2     3     4     5     6     7
classes 19  2434  3148  3688  3973  4073  4094  4096.
```

These counts were obtained independently by Moore successor refinement and
by direct response partitioning over every activation subset of size at
most `d`. Depth seven is sharp. For `F=empty`, `G={1}`, every context with
at most six activations either merges the states or leaves both on the
`325m` branch of (6). The context `{2,3,4,5,6,7,12}` of size seven gives
`1950/1399489` and `1908/1399489`. The analytic branch check independently
explains why the first six levels miss this hidden adjacent label.

The old norm gate `E(L)<E(R)` certifies 34 of the 4096 masks, whereas the
exact adverse gate certifies all 4095 nonconstant ones. If only an upper
bound on `B` is available, THM-4458's one-sided/enclosure rule must still be
used; these exact-value counts do not turn a loose upper bound above `E`
into a certificate.

This finite experiment is a genuine transfer of the old CF-return/continuation
mechanism into the new budget, with actual indicator-density operations.
It remains outside the canonical scalar-cover owner action. The next live
question is which smaller observation is closed under the *actual* owner/
word operation; applying arbitrary scanner contexts there would be a type
error. It supplies no exclusion of the 165 valuation profiles or LRC(14).

## 4. Why additive completion destroys the pure-shear distance

Let `N={q p^T:p dot q=0}` include zero, inside real trace-free `3x3`
matrices. THM-4457's observer is `d_F(M)=dist_F(M,N)`.

**PROVED.** If `h` is any nonnegative subadditive function on trace-free
matrices and `h(M)<=d_F(M)` everywhere, then `h` is identically zero.

Every pure shear has zero `h`. Every trace-free matrix is a sum of at most
eight such shears: put

```
S12=[[1,1,0],[-1,-1,0],[0,0,0]],
S23=[[0,0,0],[0,1,1],[0,-1,-1]].
```

If `M` has diagonal `(x,y,-x-y)`, then `x S12+(x+y) S23` has the same
diagonal. The remaining matrix is the sum of its six off-diagonal multiples
of `E_ij`. Each of these eight pieces has rank at most one and trace zero.
Subadditivity now gives `0<=h(M)<=0`. QED.

Thus `d_F(X-Y)` fails even the triangle inequality; it cannot be fed into
THM-2191 as if it were a metric. The greatest nonnegative subadditive
minorant of the root cost is zero. This is a precondition failure, not an
instance contradicting the valid catalytic-localization theorem.

Two shear summands already give the smallest possible witness:

```
S=E12,
T=(1,1,1)^T(1,0,-1),
M=S+T, w(M)=(1,-2,0), alpha(M)=-3/5.
```

Both have zero shear distance, while `M` has nonzero vorticity and rank two.
THM-4457 gives `d_F(M)>=3 sqrt(6/7)/5>0`. The first failed implication is
that two zero-cost operations have a zero-cost combined gradient. Retaining
the full gradient is the strongest immediate repair. A mixture or convex
relaxation needs a fluctuation/covariance invoice rather than only distance
of the mean from a set of elementary gradients. The separate
[THM-4459, pure-shear nuclear moment cost](../../01-canon/theorems/THM-4459-pure-shear-mixtures-have-exact-nuclear-moment-cost.md)
sharpens that repair; it is not assumed in this proof. No PDE trajectory is realized by this matrix example.

## Reproduction, exact universe, and stopping boundary

Run, from the repository root:

```
python 04-computation/cross_concepts_catalysis_20260908.py
python -O 04-computation/cross_concepts_catalysis_20260908.py
```

[Source](../../04-computation/cross_concepts_catalysis_20260908.py) and
[matching output](../../05-knowledge/results/cross_concepts_catalysis_20260908.out).
Normal and optimized output bytes agree. Raw LF SHA-256 pins:

```
script b597ca9410c12f0d18d4c43147321bb4b04caad3eb5b4b3e1ecb94591ea4c087
output 42215e9a0f4512a4b33a60a28a3012edf8d4a3bb27be2073fb9c509253240ef1
```

All checks use explicit exceptions, including under `-O`. The controls are
7020 pairs `R,L in {-1,0,1}^p`, `p=2,3,4`, excluding constant `R`; 18198
scanner queries; all 4096 labelled auxiliary density masks; and all `3^8`
trace-free integer matrices whose first eight entries belong to `{-1,0,1}`.
Hinge and optional transport agree; Moore and direct-context refinements
agree. Positive controls include gauge invariance, exact cancellation and
zero shear pieces; hostiles include favorable-budget continuation, observer
change, and two-shear stretching. No finite computation is extrapolated to
canonical LRC rows or physical Euler evolution.

Searches used exact budget, continuation, cone, partial transport, and
localization statements across canon and result notes, together with their
named correction mechanisms. Existing THM-4458 duality and THM-2176 generic
continuation are credited. The source/replay additions are the exact
budget quotient, scanner, density activation classification, and shear
subadditive-collapse obstruction. META-PATTERNS used: localize a genuine
metric before hunting catalysts; search the statement before the method;
separate observer and continuation. No new method card is promoted here.
