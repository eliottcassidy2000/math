# A quantitative feedback-arc budget gate for interval coverage

**Status: PROVED scoped counting theorem, with FINITE-EXACT portable controls.**
The target that every positive odd integer except 7 and 21 is a tournament
Hamiltonian-path count remains **OPEN**. This packet supplies a necessary
complexity budget for a potential constructor, not that constructor. No
literature-priority claim is made.

## 1. Inheritance and scope

`H(T)` is the raw number of directed Hamiltonian vertex words of a tournament;
set `H(empty)=1`. `fa(T)` is its minimum number of backward arcs over all
linear vertex orders. Deleting every singleton strongly connected component
(SCC) gives an induced core with the same H and no larger fa.

Closest proved mechanisms:

* Session133 `core_reduction/THEORY.md` gives the exact unrestricted core
  order budget `K(X)=max{3a+4d:3^a 5^d<=X}`, equivalently restricting `d<=2`.
* The October20 `feedback_word_hamilton_spectrum_20261020/run_density.md`
  already proves density zero for fixed feedback size and fixed run number.
  Fixed backward-arc budget is a consequence of that theorem; we do **not**
  claim the qualitative fixed-budget obstruction as new.
* THM-4748, `two-reversed-arcs-exact-path-counts-...`, already classifies
  two-backward-arc families. Those formulas are not rediscovered here.
* Overnight packet g gives the six-backward-arc family
  `H=2^(n+1)-9`, `n>=9`, and every odd count in
  `[2^k-235,2^k-1]` for every `k>=14`, on at most `k+1` vertices (118
  deficits). Packet h proves that deficit nine has minimum dyadic source
  slope 512, so that source enters only at exponents `k>=10`. Neither
  licenses increasing the target exponent to evade a source-entry cost.

The addition here is an **exact SCC-count-capped core budget**, its finite
binomial gate, and a uniform growing-budget consequence with leading
constant one. The prime lane is retained explicitly. Literal directed odd
cycles, not just their supports, are the objects in a cycle conflict graph;
an arbitrary graph independence polynomial is not a tournament constructor.

### Live board

| Lane | Object / operation | Preserved predicate | Missing information or boundary |
|---|---|---|---|
| Anchor | Odd interval / sparse-arc presentations | Actual H and actual backarc budget | A fixed budget supplies too few scalar values |
| Niche | SCC floors / exact multiplicative knapsack | Maximum core order at bounded H and SCC count | This is not the exact maximum at bounded fa |
| Wildcard | Literal odd cycles / conflict graph | H through the corrected cycle formula | Collapsing equal supports changes multiplicity |

## 2. Exact core order with a cap on nonsingleton SCCs

The only external mathematical input is the exact strong-tournament minimum:

```
f(3k)=3*5^(k-1),  f(3k+1)=5^k,  f(3k+2)=9*5^(k-1),  k>=1.
```

It is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong
Tournaments*, EJC 13 (2006), N3, Theorem 1 and Corollary 1:
[primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).
The paper was reopened for this packet. It also records SCC multiplicativity.
Our use of that input is inherited through the pinned Session133 theorem.

For integers `X>=1`, `b>=0`, let `J_b(X)` be the maximum number of vertices
of a core `C` with `H(C)<=X` and at most b nonsingleton SCCs. The empty core
is permitted.

**Theorem 1 (exact capped-SCC budget).** Include zero in the maximum. Then

```
J_b(X) = max { r+2u+3v :
              1<=r<=b, 0<=u<=2r, v>=max(0,r-u),
              3^u*5^v<=X }.
```

Equivalently, loop over `1<=r<=min(b,floor(log_3 X))` and `0<=u<=2r`,
set `v=floor(log_5(floor(X/3^u)))`, and retain the term if its lower bound
on v holds. Integer division/powers evaluate every floor exactly.

**Proof.** A strong component of each congruence class can be represented
numerically by a base order 3, 4 or 5, followed by `t>=0` increments of
three. Its minimum count is respectively `3*5^t`, `5^(t+1)` or `9*5^t`.
Suppose there are a bases of order 3, d of order 4 and c of order 5, with
r=a+d+c components and t total increments. Put

```
u=a+2c,   v=d+t.
```

The total order is

```
3a+4d+5c+3t = r+2u+3v,
```

and the product of strong minimum counts is `3^u*5^v`. Necessarily
`0<=u<=2r` and `v>=max(0,r-u)`. For an arbitrary core, replacing each
component by a strong minimum of the same order can only decrease H.
Thus its order is included in the displayed upper bound.

Conversely, for every feasible `(r,u,v)`, choose

```
c=max(0,u-r), a=u-2c, d=max(0,r-u), t=v-d.
```

All four are nonnegative integers and `a+c+d=r`. Distribute t increments
among the r components. Take actual strong tournaments attaining the
cited minima, and form their ordinal join. It has the specified order,
exact count `3^u*5^v<=X`, and r nonsingleton SCCs. Thus the maximum is
attained. QED.

**One-dimensional form.** For `b>=1`, put
`v(u)=floor(log_5(floor(X/3^u)))`. Then

```
J_b(X)=max_{0<=u<=min(2b,floor(log_3 X))}
        [2u+3v(u)+min(b,u+v(u))].
```

For fixed u,v, the largest admissible number of components is
`r=min(b,u+v)`: its lower bound `r>=ceil(u/2)` follows from `u<=2b`.
The objective increases with v, so choose its largest permitted value.
The case u=v=0 accounts for the empty core. For `b=0` set `J_0(X)=0`
separately; forgetting that branch would incorrectly allow increments
without any component. This form needs only O(b) integer-log trials.

This is an exact extremal theorem for an **SCC cap**. Its extremizer may
have much more than b feedback arcs. We use it only as an upper bound on
the order of a bounded-fa witness, never as an equivalence of constraints.

Useful consequences are

```
J_b(X)<=K(X)<=3 log_3 X,
J_b(X)<=b+floor(log_5(X^3)).
```

The second follows from `H(component)>=5^((order-1)/3)` and at most b
components. Neither inequality replaces the exact formula in finite tests.

## 3. An exact finite counting gate

Let `A_b(X)` be the set of positive integers `h<=X` for which some actual
finite tournament T has `H(T)=h` and `fa(T)<=b`. Put

```
N=J_b(X),  M=N(N-1)/2,
B_b(X)=sum_{j=0}^{min(b,M)} binom(M,j).
```

**Theorem 2 (budget gate).** `|A_b(X)|<=B_b(X)`.

**Proof.** In every linear presentation each nonsingleton SCC contains a
backward arc, since it is not transitive. Distinct SCCs require distinct
such arcs. A tournament with fa at most b therefore has at most b
nonsingleton SCCs. Delete singleton SCCs, preserving H and not increasing
fa. The remaining order is at most N by Theorem 1.

Pad with source singleton SCCs to exactly N vertices. This preserves H
and the old arc budget: place the new sources before an order of the old
core that has at most b backwards arcs. Relabel this total order as
`0,...,N-1`. Its presentation is specified by a subset of at most b of
the M unordered pairs. There are exactly B_b(X) such subsets. Distinct
H-values cannot have the same presentation, so their number is at most
this count. When N=0 the sole empty presentation accounts for H=1. QED.

The padding step avoids a gratuitous sum over all smaller orders. The
formula counts presentations, not isomorphism classes, and so is only
an upper bound on the number of values.

**Interval corollary.** If every odd integer in an interval `I subset[1,X]`
has a realization using at most b backward arcs, then

```
#{odd integers in I} <= B_b(X).
```

For `X=2^k-1`, the open dyadic interval `(2^(k-1),2^k)` contains
`2^(k-2)` odd integers. The first positive pigeonhole certificates from
this particular bound are:

| Arc budget b | First k for this counting bound | J_b(2^k-1) |
|---:|---:|---:|
|1|8|11|
|2|17|23|
|3|29|40|
|4|40|55|
|5|52|72|
|6|64|88|
|10|117|161|

These are **not** the first actual failures of each construction class;
the bound is deliberately coarse. For b=6 at k=64, precisely computed,

```
B_6(2^64-1) = 4,359,919,199,804,281,268
          < 4,611,686,018,427,387,904 = 2^62.
```

Therefore at least `251,766,818,623,106,636` odd integers in that dyadic
interval have no six-backward-arc realization. This says nothing about
whether they have realizations with more backward arcs.

At k=128 the same upper bound is `13,088,468,648,540,145,585,103`, using
N=171, out of `2^126` odd targets. The fixed six-arc deficit-nine family
remains valid: one value at each scale and a fixed-width window are
consistent with this strong sparsity.

## 4. Uniform growing-budget obstruction

All logarithms in this section are natural. Fix `0<epsilon<1` and let

```
b(X)=floor((1-epsilon) log X / log log X).
```

**Theorem 3 (near-logarithmic arc budget is necessary for dense coverage).**

```
|A_{b(X)}(X)| <= X^(1-epsilon+o(1)).
```

Thus the proportion of odd targets below X having any realization with
this low arc budget tends to zero. In particular a constructor covering
all sufficiently large allowed odd integers cannot have a uniform arc
budget bounded by this expression.

**Proof.** Since `N<=K(X)=O(log X)`, we have `M=O((log X)^2)`. Also
`N>=J_1(X)>=1+3 floor(log_5 X)` for X>=5, by the actual strong minima
of order `3v+1` and count `5^v`. Thus `M=Omega((log X)^2)`. For large X,
`b<=M/2`, and the standard product bound on binomial coefficients gives

```
B_b(X) <= (b+1) (eM/b)^b.
```

Write `L=log X`. Then `log M<=2 log L+O(1)` and
`log b=log L-log log L+O(1)`. Consequently

```
log B_b(X)
 <= log(b+1)+b(1+log M-log b)
 <= (1-epsilon)L * (1 + (log log L+O(1))/log L) + o(L)
  = (1-epsilon+o(1))L.
```

Exponentiate. QED.

Equivalently, for every fixed epsilon, density-one odd integers n have
**no** realization with
`fa(T)<=(1-epsilon) log n/log log n`. This is a necessary condition on
any realization; it is not an assertion that a realization exists.
The n-dependent formulation follows because the threshold is eventually
increasing and the finite initial range has density zero.

### Prime and 7-times-prime targets cannot be rescued by products

The same bad-budget set has relative density zero among primes and among
numbers `7p` with p an odd prime. No prime number theorem is required:

```
4^m/(2m+1) <= binom(2m,m) <= (2m)^pi(2m).
```

For the upper bound, each prime contributes at most 2m, because each
floor-difference in its binomial valuation is zero or one. Taking
logarithms proves `pi(X)>=c X/log X` for some fixed c>0 at all sufficiently
large X. This lower bound dominates `X^(1-epsilon+o(1))`; replace X by X/7
for the second target set.

Moreover a tournament with prime count has just one nonsingleton SCC.
The same holds for count 7p, using the inherited permanent omission of 7:
a nontrivial factorization into SCC counts would require a factor 7.
This retains the source-specific strong-component obligation rather than
pretending ordinal products fill the prime gaps. One may tighten N to
`J_1(X)` when counting these targets; the asymptotic leading constant is
unchanged.

## 5. Checks, hostiles, and what this does not prove

`verify_arc_budget.py` uses only the Python standard library and explicit
exceptions (never `assert`). Its deterministic output agrees normally and
under `python -O`.

* The capped-core formula is compared with a direct multiplicative knapsack
  over all strong-order minimum counts, for every `1<=X<=5000`, `0<=b<=6`:
  35,000 exact comparisons, also checking the one-dimensional form. This
  tests the arithmetic formula, not the
  externally cited strong minimum theorem.
* All 1,100 labelled tournaments of orders 0 through 5 are counted both by
  subset DP and literal permutations. Minimum fa is checked by subset DP
  and all vertex orders. SCC deletion, SCC cap, induced fa monotonicity,
  order bounds, and selected source-padding operations are checked on
  actual adjacency matrices.
* Actual interval-chain strong-minimum witnesses of every order 3 through
  15 are counted independently. The all-order attainment input remains
  the cited theorem, not this finite range.
* The actual inherited six-arc tournaments at orders 9 through 15 have
  H=`2^(n+1)-9` and minimum fa exactly 6 by two separate subset DPs.
* A transitive module in a directed triangle has H=`2^t+1`, not 3, when
  its size is t. Arbitrary twin contraction would invalidate the kernel;
  only singleton SCC deletion is used.
* A five-vertex hostile has seven literal odd directed cycles but only
  five supports: raw H=15 whereas collapsing supports gives I(2)=11.
  This reproduces the corrected THM-209 distinction and blocks an
  unauthenticated conflict-graph construction.

The numerical gate does not widen the near-dyadic window, synthesize
arbitrary deficit sources, or settle the two-apex/unbounded-word lane.
A fixed feedback-vertex set may have an unbounded number of backward arcs;
Theorem 3 is not a prohibition on that strategy. Multiplicative closure
of a bounded-arc family can increase its arc budget, but cannot add a new
prime count. A useful next positive step must grow the actual arc budget
or escape this presentation count by a different complexity carrier.
