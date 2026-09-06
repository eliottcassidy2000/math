# Exact second weights in every cycle layer, with a sharp perturbation-type transition

Status: **PROVED ANALYTICALLY + INDEPENDENTLY AUDITED**. All finite controls
below are **FINITE-EXACT**. Root and the independent `nc2_seed` referee
audited the full proof; the latter also independently checked the refined
n>=16 benchmark transport. This extends the audited Hamilton second-weight
result to every cycle length; it does not infer equality classes from an
asymptotic bound.

## 1. Statement, inheritance and the retained root

Let n>=16 and `3<=k<=n`. Put

```text
e=e_(n,k)=(n-2)!/(n-k)!,
a_(n,k)=[2-2/(n-2)]e,
b_(n,k)=[2-4(k-3)/((n-2)(n-3))]e.
```

As in the [audited cross-scale theorem](overnight_hexagon_sep05_crossscale.md),
remove the full zero set Z_k and first-weight equality set E_k. The second
distinct positive negative-k-cycle weight is exactly

```text
d2(n,k)=min(a_(n,k),b_(n,k)).                         (A1)
```

Its equality classes are:

| Scale condition | Negative-edge representative |
|---|---|
| n>2k-3 | two adjacent edges |
| n<2k-3 | two disjoint edges |
| n=2k-3 | both types |

Global edge negation supplies additional equality classes **exactly when
k is even**. Thus the labelled multiplicity is `3 binom(n,3)` in the
adjacent regime, `3 binom(n,4)` in the disjoint regime, and their sum at
the transition, doubled when k is even. The number of relabelling orbits
is respectively one, one, or two, likewise doubled when k is even.

The closest proved mechanisms are the normalized Hamilton dominance in
the cross-scale theorem and the complete R2 classification in
[the audited exact Hamilton second-weight proof](overnight_hexagon_sep05_hamilton_second.md).
The new sidecar is the **unmatched marked vertex** after contracting an
additional negative edge. Dropping that vertex would substitute the wrong
signed matching count, especially in the short layers.

The hostiles are retained explicitly: the signed Hamilton values
`W_3(1)=-1`, `W_5(1)=0`, `W_6(3)=-4`, and the negative-C5/apex minimum
failure on K6. The live board is: normalized transposition rigidity;
induced matching distribution; unmatched marked vertex; short sign tables;
parity Bonferroni; adjacent/disjoint scale. No new large switching-class
census is used. The n>=16 cutoff is sufficient, not claimed minimal.

## 2. The same finite structural carrier suffices uniformly

Let H satisfy `c_k(H)<=min(a_(n,k),b_(n,k))`. The cross-scale theorem
gives, for every pair's interior disagreement count r,

```text
c_k(H)/e >= T_(n,n)(r)/(2(n-2)!).
```

To exclude `3<=r<=n-5` at the same cutoff n>=16, retain the benchmark's
dependence on k as well as the transposition count. Put

```text
R_k=T_(n,k)(r)/(2e_(n,k)),
U_d(k)=b_(n,k)/e_(n,k),
q=r(n-2-r),  D4=(n-2)(n-3)(n-4)(n-5).
```

For k>=6 the exact cross-scale formula gives

```text
(R_k-U_d(k))-(R_n-U_d(n))
 =2(n-k)[q(q-n+3)-2(n-4)(n-5)]/D4.                 (T)
```

In the present range `q>=3(n-5)`, both factors in `q(q-n+3)` increase
with q, and

```text
q(q-n+3)-2(n-4)(n-5) >=4(n-5)(n-7)>0.
```

The independently audited Hamilton n>=16 proof gives
`R_n>U_d(n)` for every such disagreement. Therefore (T) gives
`R_k>U_d(k)>=min(a_(n,k),b_(n,k))/e` at every k>=6.

For k=3,4,5, the normalized transposition value is q/(n-2), and

```text
q/(n-2)-a_(n,k)/e >=(n-9)/(n-2)>0.
```

These two arguments exclude all R3 disagreements at n>=16. The earlier
proof using a single k-independent adjacent-edge benchmark required n>=18;
keeping the benchmark coordinate removes that unnecessary loss. Thus every
row pair has

```text
r in {0,1,2,n-4,n-3,n-2}.                            (A2)
```

The proved R2 classification, valid already at n>=13, now says H is
switching-equivalent, possibly after global negation, to an arbitrary
matching or one of

```text
P3, P4, C3, C4, C5,
```

with all unused vertices isolated. Pj denotes the path on j vertices.
That classification retains all sign/complement cases and is independent
of the layer k. We now evaluate this finite list of types, keeping the
unbounded matching-size parameter.

## 3. All-layer signed matching sums are positive at n>=16

Let W_l(j) be the signed Hamilton sum on K_l with j negative disjoint
edges. The exact identity and sufficient sign theorem from the Hamilton
proof are

```text
W_l(j)=(1/2) sum_(i=0)^j binom(j,i)(-4)^i(l-i-1)!,
l-j-1>=4 => W_l(j)>0.                               (A3)
```

The second statement follows from the actual integral
`(1/2) integral exp(-t)t^(l-j-1)(t-4)^j dt` by pairing t=4-x and
t=4+x, not from a scalar analogy. The complete short tables are

| l | W_l(j), in order j=0,...,floor(l/2) |
|---|---|
| 3 | 1, -1 |
| 4 | 3, -1, 3 |
| 5 | 12, 0, 4 |
| 6 | 60, 12, 12, -4 |
| 7 | 360, 120, 72, 24 |
| 8 | 2520, 1080, 600, 312, 216 |

Let `W_(n,k)(m)` denote the signed k-cycle sum in K_n with m negative
matching edges. For a uniformly chosen k-subset S, let J be the number
of matching edges wholly inside S. Every k-cycle lies in exactly one
such subset, so

```text
W_(n,k)(m)=binom(n,k) E[W_k(J)],
E J=m k(k-1)/[n(n-1)] <= k(k-1)/[2(n-1)].            (A4)
```

For n>=16, this signed sum is strictly positive at **every** length:

- For k>=9, `k-J-1>=ceil(k/2)-1>=4`, so (A3) applies termwise.
- For k=7,8, every entry in the corresponding short table is positive.
- For k=6, the table gives `E W_6(J)>=12-16 Pr(J=3)`.
  Since `Pr(J=3)<=E J/3<=5/(n-1)`, this is at least
  `12-80/(n-1)>0`.
- For k=5, all entries are nonnegative, and some five-subset has J=0:
  a maximum subset choosing at most one endpoint from each matching edge
  has size n-m>=ceil(n/2)>=8. Thus a strictly positive term occurs.
- For k=4, `E W_4(J)=3-4 Pr(J=1)>=3-4 E J>=3-24/(n-1)>0`.
- For k=3, `E W_3(J)=1-2 E J>=1-6/(n-1)>0`.

Consequently

```text
W_(n,k)(m)>0 for n>=16, 3<=k<=n, 0<=m<=floor(n/2).    (A5)
```

This averaging argument is why the negative small Hamilton entries are
not counterexamples to the large-ambient-order statement.

## 4. Strict matching monotonicity requires the unmatched root

Write w_(n,k)(m) for the negative-cycle weight of the matching. Add a new
negative edge disjoint from the existing m edges. The weight change is
the old signed sum over cycles containing that edge. For k>=4, contract
the edge. Each contracted cycle has two expansions, and the contracted
vertex is necessarily present and is **unmatched**. Hence

```text
w_(n,k)(m+1)-w_(n,k)(m)=2 W^*_(n-1,k-1)(m),          (A6)
```

where the star means cycles containing that fixed unmatched vertex.
For k=3 the increment is directly n-2, since no old matching edge can
occur in a triangle that contains the new disjoint edge.

For n>=16 the right side of (A6) is strictly positive. Choose k-2 vertices
among the N=n-2 non-root vertices of the contracted graph. If J is their
number of full matching edges, the signed Hamilton count on the resulting
k-1 vertices is W_(k-1)(J). Here `J<=floor((k-2)/2)`, not the larger
unmarked maximum `floor((k-1)/2)`.

- For k>=9, the integral exponent is
  `(k-1)-J-1>=ceil(k/2)-1>=4`, giving positivity termwise.
- For k=8, the W_7 table is positive. For k=7, only entries
  `W_6(0),W_6(1),W_6(2)` can occur; the negative `W_6(3)` is impossible
  because the marked vertex is unmatched.
- For k=6, the W_5 table is nonnegative. A no-hit four-subset exists,
  because N-m>=ceil(N/2)>=7, and supplies a positive term.
- For k=5, the three chosen non-root vertices contain at most one matching
  edge. Its probability is `6m/[N(N-1)]<=3/(N-1)`, so the average W_4
  is at least `3-12/(N-1)>0`.
- For k=4, the two chosen vertices form a matching edge with probability
  `m/binom(N,2)<=1/(N-1)`. The average W_3 is therefore at least
  `1-2/(N-1)>0`.

Thus

```text
w_(n,k)(m+1)>w_(n,k)(m)                               (A7)
```

whenever n>=16, `3<=k<=n`, and `m+1<=floor(n/2)`. In particular, among
matching representatives beyond one edge, precisely the two-edge matching
minimizes the weight, with value b_(n,k).

## 5. The other finite motifs cost more than twice the edge

For a signing with f negative edges, a cycle containing j of them obeys

```text
1_(j odd) >= j-2 binom(j,2).
```

The numbers of k-cycles containing two specified edges are

```text
I_adj=(n-3)!/(n-k)!,
I_disj=2(k-3)(n-4)!/(n-k)!.
```

Both are at most `2e/(n-2)` for k<=n. Summing the elementary parity
inequality therefore gives

```text
c_k(H)/e >= f-2f(f-1)/(n-2).                         (A8)
```

For f=3,4,5 this is strictly greater than two at n>=16: at n=16 the
three lower values are `15/7,16/7,15/7`, and they increase with n.
Thus P4,C3,C4,C5 all lie strictly above both two-edge competitors.
Their full-cycle subsets and short-layer exceptions cannot invalidate
this bound, because (A8) is a literal per-cycle inequality.

The remaining P3 is exactly the adjacent two-edge signing and has weight
a_(n,k). This proves the comparison between all **unnegated** classes in
the structural list.

## 6. Global negation, the transition and complete equality

For even k, global negation preserves every k-cycle sign, so it supplies
the stated extra equality families. For odd k it replaces w by N_(n,k)-w.
For a matching, (A5) gives w<N_(n,k)/2, so its global negative has weight

```text
> N_(n,k)/2 = [n(n-1)/(4k)]e >=[(n-1)/4]e>2e.
```

Every finite motif has at most five negative edges, so w<=5e by the union
bound. Its global negative has weight at least

```text
[n(n-1)/(2k)-5]e >=[(n-1)/2-5]e>2e                 for n>=16.
```

This accounts for every odd global-negative case, including the balanced
and single-edge classes hidden inside the matching list. Therefore only
the adjacent and disjoint two-edge types can attain the second weight.
Their exact difference is

```text
a_(n,k)-b_(n,k)=2(2k-n-3)e/[(n-2)(n-3)],            (A9)
```

which gives the three regimes in (A1), with **both** families retained at
the equality boundary. Notice that the parity here is the cycle length k,
not necessarily the ambient order n.

For multiplicity, choose the centre and two leaves for an adjacent pair,
giving `n binom(n-1,2)=3 binom(n,3)` possibilities. The disjoint count is
`3 binom(n,4)`. Distinct two-edge signings, even of different shapes, cannot
be switching-equivalent: their symmetric difference has at most four edges,
whereas a nonempty cut has at least n-1>4. Their global negatives are
separate, because the original negative-triangle counts are `2(n-3)` or
`2(n-2)`, and the complementary counts are larger than both for n>=16.
These counts also distinguish the shape and global-sign relabelling orbits.

## 7. Exact controls and remaining boundary

```bash
python3 04-computation/overnight_hexagon_sep05_crossscale_second.py
python3 -O 04-computation/overnight_hexagon_sep05_crossscale_second.py
```

One checker path enumerates every literal simple cycle on K5 through K9,
at all lengths and all admissible matching sizes. A separate path chooses
induced subsets by their numbers of full matching edges, single endpoints,
and unmatched vertices, then uses exact signed Hamilton sums. It checks
the complete subset counts before comparing weights. The unmatched-root
contraction identity is tested separately at every length and matching
size for n=16,17,18,19,20,21,30,40, including the n=19 and n=21 transition
layers and all short negative/zero sign controls. No full switching-class
enumeration at any of these larger orders is performed.

The refined exact all-layer theorem now starts at the same order 16 as
the Hamilton theorem. No conclusion for all shorter ambient orders is
inferred. The formerly open cumulative second weights are now proved for
n>=16 in [the cumulative second-weight theorem](overnight_hexagon_sep05_cumulative_second.md):
mixing parity retains an antibalanced third competitor, and the sharp
worst-truncation isolation constant is six. Higher weights and general
edit-distance estimates remain open. No LRC(14), tournament, or Boolean
quotient result is asserted.

Normal and optimized runs are byte-identical: 17,478 explicit failure
gates and 72,238 literal cycles. Raw-LF SHA-256:

```text
5cf7d1f6178385117c31b2d3d2940cc851970468af8acb7312e9c1b0127545c0  04-computation/overnight_hexagon_sep05_crossscale_second.py
2cdf775ff42ac830add86b97215566a4db6be61dcc5e68f537937767c7c47840  05-knowledge/results/overnight_hexagon_sep05_crossscale_second.out
```
