# An exact convergent series for the two-step prefix price

Status: **PROVED + independently audited by root and geometry; FINITE-EXACT coefficients; THM-4491.**
This note concerns the pairing family with global two-step descent. It does
not prove Collatz or a longer-horizon price result. The limit below concerns
finite-prefix optima and the minimum lower natural density; natural-density
attainment is a separate question.

## 1. Inheritance and the new invariant

The [global flow note](crossroads223_20260926_flow_global.md) proved that the
exact two-step clauses form a rooted tree on pair indices `i>=2`. Let

```
C(X)=min{|F intersect [1,X]|: F defines a global P_2 pairing}.
```

Every feasible finite prefix extends globally, so `C(X)` is the exact tree
DP optimum. Its numerical values fluctuate: `C(10000)/10000=0.2907`, while
`C(100000)/100000=0.2908`. Those finite ratios alone do not prove a limit.

The useful new object is the *local cost increment* at each tree vertex.
These nonnegative increments telescope the full optimization cost and
separate its different depth scales. The carry offsets remain bounded
after rescaling a descendant back to its source. This is what turns the
finite DP into a convergent series, without assuming independent bit choices.

## 2. Full-depth differences and local cost increments

For the full tree of descendants through depth `h` below an index `i`, let
`F_h^b(i)` be the minimum number of flipped vertices when the root bit is
fixed to `b`. Put

```
Delta_h(i)=F_h^1(i)-F_h^0(i),
tau_h(i)=min_b F_h^b(i)-sum_children min_b F_(h-1)^b(child).
```

At depth zero, `Delta_0=1` and `tau_0=0`. For `h>=1`, the exact recurrence
uses the following child differences from depth `h-1`:

```
i even, child j=3i/2, x=Delta_(h-1)(j):
    Delta_h(i)=1-x,
    tau_h(i)=1{x>=1}.

i odd, children j=(3i-1)/2, k=(3i+1)/2,
x=Delta_(h-1)(j), y=Delta_(h-1)(k):
    Delta_h(i)=1+min(x,0)+max(y,0),
    tau_h(i)=min(max(-x,0),1+max(y,0)).
```

These follow directly by fixing the root to zero or one in the inherited
complement and monotonicity clauses. Induction gives

```
|Delta_h(i)|<=h+1,        0<=tau_h(i)<=h.
```

The same bound in terms of maximum subtree depth holds for a tree cut at
an arbitrary index `X`, including nodes with only one surviving child.
The residue `i mod 2^h` determines the full-depth recurrence. Define the
nonnegative integers

```
A_h=sum_(i=0)^(2^h-1) tau_h(i).
```

Residue zero here is a formal parity type for the finite-depth recurrence;
it does not introduce the integer zero into the positive pairing tree.

As a further exact control, `sum_(i mod 2^h) Delta_h(i)=2^h`. Indeed the
even children run through all residues once, and both odd child maps do
the same. The average even-root difference is `(1-E Delta)/2`; the odd
contribution is `(1+E Delta)/2`.

## 3. The limit and its exact series

**Proposition.** The prefix price has a limit

```
lim_(X->infinity) C(X)/X = alpha
    =sum_(h>=1) A_h/3^(h+1).
```

Writing `alpha_H` for the sum through depth `H`, an explicit error bound is

```
0<=alpha-alpha_H<=(H+3)(2/3)^(H+1).
```

**Proof.** For a tree cut at `X`, define the local increments using the
actual child subtrees. Summation telescopes: `C(X)` is the sum of these
increments over the vertices `2,...,X`. They are nonnegative.

Every descendant `v` at depth `h` below `i` satisfies

```
(3/2)^h(i-1)+1 <= v <= (3/2)^h(i+1)-1,
```

because each child is `(3/2)i` plus an error of absolute value at most
one half. Therefore, for fixed `h`, throughout the interval

```
1+(2/3)^(h+1)(X-1) < i <= (2/3)^h(X+1)-1,
```

all descendants through depth `h` are present and all depth `h+1`
descendants are absent. The difference from the macroscopic interval
`((2/3)^(h+1)X,(2/3)^h X]` involves only a bounded number of indices for
each fixed `h`. On this full-depth band, the local increment is the
periodic function `tau_h(i)`.

For fixed `h`, division by `X` and averaging its finitely many residues
gives the contribution

```
[(1/3)(2/3)^h] [A_h/2^h] = A_h/3^(h+1).
```

To control the remaining small-index region uniformly, a vertex `i>=2`
has depth at most `log_(3/2)((X-1)/(i-1))`. Its increment is bounded by
that depth plus one. Thus the normalized total from `i<=epsilon X` is
`O(epsilon(1+log(1/epsilon)))+o(1)`, by the integral bound for
`sum log(X/(i-1))`. This tends to zero with epsilon. Fixed finite depth
bands therefore exhaust the normalized DP cost and prove the limit.

Finally `A_h<=h 2^h`, so the series is absolutely convergent and

```
sum_(h>H) A_h/3^(h+1)
 <=(1/3)sum_(h>H)h(2/3)^h
 =(H+3)(2/3)^(H+1).
```

## 4. What this says about infinite pairings

Every global `P_2` pairing satisfies `|F intersect [1,X]|>=C(X)`. Hence

```
liminf_(X->infinity)|F intersect [1,X]|/X >= alpha.
```

This lower density is attainable by an infinite member, although the
construction does not establish the existence of its natural density.
To see this, fix any feasible prefix through `Y`. At a later cutoff `X`,
forcing that prefix changes the unrestricted optimum by at most
`O(Y(1+log X))`. There are at most `2Y` child subtrees on the boundary of
the old prefix; forcing a child bit costs at most its absolute root
difference, bounded by its depth plus one. Changing the old prefix itself
costs at most `Y`. These subtrees are disjoint, so the costs add.

Choose a sequence of cutoffs with
`X_(j-1) log(X_j)/X_j -> 0`, for example successive squares starting above
two. Extend the previously fixed prefix to a conditionally optimal prefix
at each cutoff. Their union is a global `P_2` assignment. Its density
along these cutoffs tends to `alpha`; the universal lower bound gives
lower natural density exactly `alpha`. Thus `alpha` is the minimum lower
natural density over global two-step pairings.

This statement must not be silently replaced by an assertion that one
member has natural density `alpha`, or that a separately defined
upper-density objective has the same value.

## 5. Exact coefficients and reproduction

The first twenty coefficients are

```
1, 1, 4, 7, 16, 31, 62, 125, 242, 504,
985, 2014, 4004, 8106, 16183, 32437, 65039,
129976, 260488, 521302.
```

They give the exact interval

```
3041388727/10460353203 <= alpha <= 3089623223/10460353203,
0.2907539227382626... <= alpha <= 0.2953650955221956... .
```

In particular, the lower endpoint is a new asymptotic lower bound for
every global two-step member. It is not a lower bound on every finite
prefix ratio; the finite value at 10,000 illustrates that distinction.

Run `python3 04-computation/experiments/crossroads223_20260926_flow_series.py`
or the same command with `-O`. The retained
[output](crossroads223_20260926_flow_series.out) agrees in both modes.
The integer recurrence is independently compared with direct full-tree
optimization in all 254 residue cases through depth seven. The universal
depth bounds and mean-difference identity are checked through depth twenty.
The printed decimal interval is a rendering of the exact rational bounds.

The averaged local increments in this range are at most one half, but no
uniform theorem about that sharper bound is claimed. A sharper tail or a
proof of natural-density attainment would be additional work. The guarded
actual-incidence tree is essential: the global horizon-four counterexample
in the preceding note invalidates direct reuse of these two-step clauses.
