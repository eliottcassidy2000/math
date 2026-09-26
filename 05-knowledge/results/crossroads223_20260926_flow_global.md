# A charge that counts actual flipped pairs: a strict two-step density bound

Status: **DRAFT PROOF + FINITE-EXACT independent checks**. This concerns
global two-step descent in the pairing family. It does not prove Collatz,
an asymptotic optimum for two-step pairings, or a bound for longer horizons.

## Inheritance and the missing coordinate

The [preceding flow note](crossroads223_20260926_flow.md) proved an exact
`9/8` penalty for a source-weighted executed-cost functional. That functional
counts a shared flip again for each source, so it is not global flip density.
The source `7,11` control showed the first failure: one actual flip can rescue
both while one source abandons its private optimum.

The current [pairing-transitions note, section 2.4](procgen_brackets_20260924_pairings_transitions.md)
already gives the exact global two-step clauses and the lower bound `1/4`
from their multiplicative matching. It reports exact finite window optima
near `0.2908`, not an infinite optimum. The
[pairpeak note](procgen_pairpeak_20260926_pairing_peak_price.md) separates
private prices from a single globally consistent pairing.

[THM-2300](../../01-canon/theorems/THM-2300-small-owner-multipliers-force-same-character-relation-multiples.md)
offers a useful methodological comparison: separately forced source and
image atoms need not share a multiplier. Here separately optimal source
certificates need not share consistent pair bits. This is a comparison of
lost information; no LRC theorem is a dependency of the results below.

The active objects are: actual pair ownership; signed landing clauses;
the canonical multiplicative matching; its unmatched bank; and the new
`UUDD` option at horizon four. The first decisive probe sought a prefix
where the signed clauses cost more than the positive matching. It found
pair cutoff 23. The second probe tried extending that new clause beyond
two steps and produced an explicit global four-step counterexample.

## 1. Exact global objective and its tree

For pair `P_i={2i-1,2i}`, the bit `epsilon_i` changes the map as follows:

```
epsilon_i=0:  2i-1 -> 3i-1,   2i -> i;
epsilon_i=1:  2i-1 -> i-1,    2i -> 3i.
```

Let `F={i:epsilon_i=1}`. Property `P_L` means every integer `n>=3` has an
iterate strictly below `n` within `L` steps. The global objective is the
number of distinct indices in `F`, with each index charged once.

Write `U(n)=ceil(3n/2)` and `D(n)=floor(n/2)`. A first down step descends;
`UD(n)<n` for `n>=2`; and two up steps do not. Thus `P_2` is exactly the
following inherited clauses, for every source-pair index `i>=2`:

```
i=2k:    epsilon_(2k)+epsilon_(3k)=1;
i=2k+1:  epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).
```

**Proposition 1 (exact prefix optimization).** Retain the clauses whose
two endpoints are at most `X`. Their constraint graph on `2,...,X` is a
tree. Every feasible assignment extends to an infinite `P_2` member.
Consequently minimizing `|F intersect [1,X]|` over all global `P_2`
members is an exact linear-time binary tree dynamic program.

Indeed every `j>=3` has exactly one parent:

```
j=0 mod3: parent(j)=2j/3;
j=1 mod3: parent(j)=(2j+1)/3;
j=2 mod3: parent(j)=(2j-1)/3.
```

The parent is smaller than `j` and at least two. Its edge is respectively
an opposite-bit constraint, a child-at-most-parent constraint, or a
child-at-least-parent constraint. Each new child can always be assigned,
so no consistency condition arrives later at an old vertex.

For completeness the elementary LP relaxation is also integral. Give the
root its proposed Bernoulli marginal. On a complement edge the child is
forced; on an edge with child marginal `q<=p` choose the child to be one
with probability `q/p` when its parent is one, and zero otherwise. On an
edge with `q>=p`, set the child to one when its parent is one, and use
probability `(q-p)/(1-p)` otherwise. Zero-probability cases are harmless.
This constructs a distribution on valid binary assignments with any
feasible fractional marginals, proving the convex-hull assertion.

The DP and its reconstructed bit assignments reproduce the earlier exact
window values without a solver: 874 flips at `X=3000` and 29,080 at
`X=100000`. At `X=223` the exact value is 66. Different prefix optima need
not be compatible with each other; no infinite optimal-density limit is
asserted here.

## 2. A new cut in the globally unmatched bank

The canonical matching is

```
M = {{2r,3r}: r>=1 and v3(r) is even}.
```

Its edges are disjoint, and a global `P_2` member chooses exactly one
flipped endpoint on each. The number with both endpoints at most `X` is

```
B(X)=sum_(j>=0)[floor(X/(3*9^j))-floor(X/(9*9^j))]
    =X/4+O(log(X+2)).
```

Every odd index with even 3-adic valuation is globally unmatched. We now
force additional flips in that bank, so these charges cannot be recycled
onto the matching endpoints.

**Proposition 2 (extra disjoint clauses).** Every global `P_2` member obeys

```
epsilon_(36s+17)+epsilon_(54s+23)>=1       for every integer s>=0.
```

Set

```
i=16s+7, j=24s+10, h=24s+11, l=36s+15,
a=36s+17, b=54s+23.
```

Here `j` is even and `i,h,l` are odd. The exact landing clauses give

```
epsilon_j+epsilon_l=1,
epsilon_j<=epsilon_i<=epsilon_h<=epsilon_a,
epsilon_l<=epsilon_b.
```

Adding proves the cut. Both `a,b` are odd and equal to two modulo three,
so neither belongs to any edge of `M`. Moreover the progressions
`36s+17` and `54t+23` never meet: an intersection would require
`gcd(36,54)=18` to divide six. Thus the extra edges are mutually disjoint
and also disjoint from the canonical matching, for all parameters.

**Corollary (genuine global density lower bound).** For every global `P_2`
member and every integer `X>=1`,

```
|F intersect [1,X]| >= B(X)+max(0,1+floor((X-23)/54)).
```

In particular, without assuming the density exists,

```
liminf_(X->infinity) |F intersect [1,X]|/X
    >= 1/4+1/54 = 29/108 = 0.2685185185... .
```

This is strictly stronger than the inherited matching bound. It is weaker
than the observed finite window optima. The improvement counts actual
distinct flipped pairs, not private costs, source multiplicities, or a
formal-cell measure.

The smallest prefix where the exact signed-clause optimum exceeds `B(X)`
is `X=23`: the values are six and five. This minimality statement refers
only to the declared prefix universe; it is checked by the exact DP for
all smaller cutoffs. The six-index witness is `7,10,11,15,17,23`.

## 3. Why the new cut does not transfer to horizon four

There is no new first descent at step three. If the first two steps have
not descended, they were `UU`, and
`UUD(n)>=floor(9n/8)>=n`. Thus `P_3=P_2`. At step four the new option is
`UUDD`, which can descend without satisfying the two-step landing clauses.

Define the inherited free-zero global `P_2` assignment recursively:

```
b(1)=b(2)=0;
b(3k)=1-b(2k);
b(3k+1)=0;
b(3k+2)=b(2k+1)              (indices >=3).
```

Each argument on the right is smaller. The three clauses in section 1
verify `P_2` immediately. In this member `b(17)=0` and `b(23)=1`.
Now change only `b(23)` to zero.

**Proposition 3 (global hostile control).** The edited assignment is a
global `P_4` member but violates the extra clause at `s=0`, since both
bits 17 and 23 are zero. Its only source requiring more than two steps is

```
30 -> 45 -> 68 -> 34 -> 17.
```

Only map values at 45 and 46 change. A source `n>46` that encounters either
point has already descended, so its original two-step certificate remains
valid. Direct exact verification of the 44 sources `3<=n<=46` finds the
displayed four-step path and no other delay. The mechanism is also local:
the change at 46 makes its step downward, while the change at 45 delays
the former landing of source 30. Source 45 itself still descends in two
steps, `45 -> 68 -> 34`.

This counterexample invalidates an extension of Proposition 2 to all
`P_4` members. It does not show that the numerical bound `29/108` is false
for `P_4`: the one-bit change preserves this example's original density.
The precise failed implication is the requirement that an upward landing
must immediately move down.

## 4. Reproduction and stopping boundary

Run `python3 04-computation/experiments/crossroads223_20260926_flow_global.py`
or the same command with `-O`. Output is retained in
[crossroads223_20260926_flow_global.out](crossroads223_20260926_flow_global.out).
Both modes give identical output; all checks are explicit.

Independent controls include 65,534 exhaustive Boolean assignments through
cutoff 16 against the tree DP; exact optimal reconstruction at the reported
prefixes; all 64 assignments of the six-bit local implication; 1,000 affine
instances; direct matching/extra-edge disjointness at cutoff 100,000; and
the finite set sufficient for the global `P_4` surgery proof. The reported
100,000 optimum is an exact dynamic program, not a time-limited solver.

The useful transfer from the first lane is to charge a forced choice to a
bank disjoint from already charged pairs. The sidecar is the actual signed
incidence and an injective assignment to unmatched pair indices. At longer
horizons the newly legal `UUDD` alternative destroys the two-step implication
tree. A corresponding longer-horizon inequality must retain these path
alternatives and their shared bits. This note stops at that exact boundary;
it does not extrapolate the trace-cost penalty or the matching tree to
global private-price gluing.
