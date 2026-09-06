# Exact second Hamilton weight from bounded disagreement and matching positivity

Status: **PROVED ANALYTICALLY + INDEPENDENTLY AUDITED**, with separately
labelled **FINITE-EXACT** controls. Root and the independent `nc2_seed`
referee checked the threshold, full structural classification, scoped
matching integral sign, recurrence, all motif/complement cases and equality
multiplicity. The parent session has reserved THM-4433 for this result and
the preceding stability bracket. Reservation is not a proved dependency.

## Inheritance, hostile and statement

The closest mechanism is the exact transposition difference and row-type
classification in [the all-cumulative gap proof](overnight_hexagon_sep05_d7_d8_gap.md),
promoted separately as THM-4427. The first quantitative continuation was
[the asymptotic second-weight bracket](overnight_hexagon_sep05_hamilton_stability.md).
The new operation is to keep one more allowed row-disagreement value and
classify the resulting bounded-degree graph. Its previously unused sidecar
is a positive integral for the exact signed matching count.

Two hostiles are load-bearing. The negative C5 plus a positive apex on K6
has Hamilton weight 20, below the single-edge weight 24. Also, a perfect
three-edge negative matching on K6 has signed Hamilton sum **-4**; the
matching integral is not universally positive. We prove and retain the
precise sufficient positivity condition rather than assuming an analogy.

Let n>=16 and let H be an edge signing of K_n, modulo cut switching. Let
`w(H)=c_n(H)` count unoriented negative Hamilton cycles, and put

```text
e=(n-2)!,       b=2(n-4)(n-3)!.
```

Remove the zero classes (balanced only for odd n; balanced and antibalanced
for even n) and all first-weight classes (single edges for odd n; single
edges and their global negatives for even n). Then

```text
the second distinct positive weight is exactly b.                 (E1)
```

Its equality classes are exactly two disjoint negative edges, allowing
global edge negation when n is even. There are `3 binom(n,4)` labelled
classes for odd n and `6 binom(n,4)` for even n. They form one relabelling
orbit for odd n and two for even n. The second distinct positive eigenvalue
of the unnormalized multiplicity-weighted Hamilton-cycle Laplacian is `2b`.

This is an all-order theorem from n=16, not a finite census. No sharpness
of the cutoff is claimed. The full spectrum, higher weights and cumulative
second eigenvalues remain outside its scope.

## 1. A low-weight signing has only two exceptional disagreements per row pair

For two vertices u,v write

```text
r(u,v)=#{x outside {u,v}: sign(ux)!=sign(vx)},
q=r(n-2-r).
```

Transposing u,v and multiplying edge signs gives a negative K_(2,r).
Its exact Hamilton weight and the xor triangle inequality give

```text
2w(H) >= T_n(r)=2q[(n-2)(n-3)-2q](n-5)!.
```

Suppose `3<=r<=n-5`. Then `3(n-5)<=q<=(n-2)^2/4`. The quadratic
`F(q)=q[(n-2)(n-3)-2q]` is concave. To show `T_n(r)/2>b`, compare its
two real endpoint values with `b/(n-5)!=2(n-4)^2(n-3)`. The excesses are

```text
n^3-26n^2+193n-444,
(n-4)(n^3-22n^2+124n-200)/8.
```

On writing n=16+t these become

```text
t^3+22t^2+129t+84,
(t^4+38t^3+500t^2+2504t+2976)/8,
```

strictly positive for t>=0. Thus

```text
w(H)<=b => r(u,v) in {0,1,2,n-4,n-3,n-2} for every u!=v.           (E2)
```

Switching at any collection of vertices either preserves a particular
pair's r or replaces it by n-2-r, according to whether the switches at
u,v agree. Global edge negation preserves r. Hence the predicate (E2)
survives both operations used below.

## 2. Complete structural classification of (E2) for n>=13

We prove the classification with the slightly stronger range n>=13.
Root-gauge the signing so all edges at vertex zero are positive. Let G
be the negative-edge graph on the other m=n-1 vertices. Applying (E2)
to (0,v) says every degree is at most two or at least m-3. Call these
vertices low and high, and let their numbers be a,b respectively.

There are at most 2a low/high edges. A high vertex misses at most two
vertices in total, so there are at least b(a-2) low/high edges. Consequently

```text
b(a-2)<=2a,          (a-2)(b-2)<=4.
```

If a,b>=3, the positive integers a-2,b-2 with product at most four have
sum at most five, forcing a+b<=9. This contradicts m>=12. Thus one type
has size at most two. Complement G if necessary, which means global edge
negation followed by restoration of the root gauge. We may assume b<=2.

### Switch the at-most-two high vertices

Now switch exactly those b high vertices in the **full** signed graph,
without subsequently restoring the root gauge. Call the new negative-edge
graph J, on all n vertices. The old root has degree b<=2. A low vertex of
old degree d<=2 with q high neighbours now has degree `d+b-2q<=b+2<=4`.

For a high vertex, let h be its number of neighbours among the other high
vertices. Its number of missed low vertices is at most
`2-[(b-1)-h]=3-b+h`. After the switch its degree is at most

```text
(3-b+h)+h+1 <= b+2 <=4,
```

where the last +1 is the now-negative edge to the old root. Hence J has
maximum degree at most four. Every pair then has r<=8. Since n>=13 implies
n-4>=9, the large alternatives in (E2) are impossible. Thus

```text
r_J(u,v)<=2 for every pair.                                      (E3)
```

### Degrees four and three are impossible

Suppose a vertex u has degree four. A non-neighbour v satisfies

```text
4+deg(v)-2|N(u) intersect N(v)| <=2.
```

It must therefore have at least two common neighbours with u. The number
of length-two incidences through N(u), excluding returns to u, is at most
`4*(4-1)=12`. There are n-5 non-neighbours, each consuming at least two
of these incidences, so `2(n-5)<=12`, or n<=11, a contradiction.

The maximum degree is now at most three. If u has degree three, each
non-neighbour v must have at least one common neighbour with u, again by
(E3). At most `3*(3-1)=6` length-two incidences are available, giving
`n-4<=6`, or n<=10. This too is impossible. Hence J has maximum degree two.

### The remaining graph shapes

If every degree is at most one, J is an arbitrary matching. Otherwise
it has a degree-two vertex. There cannot be another nontrivial connected
component: a degree-two vertex in the first component and a positive-degree
vertex in the other would have r>=3.

Thus there is a single nontrivial component, a path or cycle. A path on
at least five vertices has an endpoint and a vertex at distance three
with disjoint neighbour sets of respective sizes one and two, contradicting
(E3). A cycle on at least six vertices has two vertices at distance three
with disjoint two-element neighbour sets, giving r=4. Therefore the only
nonmatching possibilities are

```text
P3, P4, C3, C4, C5, with all other vertices isolated.              (E4)
```

Here Pj denotes the path on j vertices. We have proved that every signing
satisfying (E2), for n>=13, is switching-equivalent to one of the graphs
(E4) or a matching, possibly after global negation. Conversely all these
graphs have (E3), so the classification is exact, not merely necessary.

## 3. Exact matching counts, the sign boundary and monotonicity

Let `W_n(m)` be the signed Hamilton sum for a negative matching of m edges,
where `0<=m<=floor(n/2)`. The number of Hamilton cycles containing a chosen
set of j disjoint edges is

```text
2^(j-1)(n-j-1)!.
```

For j=0 this means `(n-1)!/2`; otherwise contract the j edges and retain
their independent orientations. Expanding the product of the cycle signs
gives the exact integer identity

```text
W_n(m) = (1/2) sum_(j=0)^m binom(m,j)(-4)^j(n-j-1)!
       = (1/2) integral_0^infinity exp(-t)t^(n-m-1)(t-4)^m dt.       (E5)
```

The second equality follows by expanding the polynomial and applying the
elementary factorial integral. No limiting probabilistic approximation or
numerical quadrature enters.

### Precisely scoped positivity

Put a=n-m-1. If m is even, the integrand is nonnegative and not identically
zero. If m is odd and a>=4, pair the negative contribution at t=4-x,
0<x<4, with the positive contribution at t=4+x. Their positive/negative
absolute-value ratio is

```text
exp(-2x) [(4+x)/(4-x)]^a.
```

For 0<x<4,

```text
log((4+x)/(4-x)) >= x/2,
```

as follows, for example, by differentiating the difference, which vanishes
at x=0 and has positive derivative thereafter. Thus the ratio is at least
one when a>=4. The remaining tail t>8 is strictly positive. Consequently

```text
n-m-1>=4 => W_n(m)>0.                                             (E6)
```

This is a sufficient condition, not an iff. In particular the exact small
values `W_3(1)=-1`, `W_5(1)=0`, and `W_6(3)=-4` prevent a universal-sign
overclaim. The last is the perfect matching on K6.

Let w_n(m) be the negative-cycle weight of the matching. For n>=4 and
m+1<=floor(n/2), from (E5), or
by conditioning on one new disjoint edge and contracting it,

```text
W_n(m)-W_n(m+1)=4W_(n-1)(m),
w_n(m+1)-w_n(m)=2W_(n-1)(m).                                     (E7)
```

For n>=16 and m+1<=floor(n/2), the exponent for this last signed sum is
`n-m-2>=ceil(n/2)-1>=7`, so (E6) makes every increment strictly positive.
Thus among matching classes other than zero and a single edge, the unique
smallest matching size is two. Its weight is exactly

```text
w_n(2)=2e-4(n-3)!=b.                                              (E8)
```

## 4. All finite motifs and all global-negative cases

Write `a=(n-3)!`, `c=(n-4)!`, and `d=(n-5)!`. Literal edge-subset
inclusion-exclusion gives the following complete list. A linear forest
with j edges and h nontrivial path components is contained in
`2^(h-1)(n-j-1)!` Hamilton cycles; a proper cycle component is contained
in none. This accounts for every subset, including the full C3/C4/C5
subsets which cannot be extended to a Hamilton cycle at n>=16.

| Negative-edge motif | Exact weight | Excess above b |
|---|---|---|
| P3 | 2e-2a | 2a |
| P4 | 3e-8a+4c | c(n^2-9n+22) |
| C3 | 3e-6a | (n-4)a |
| C4 | 4e-16a+16c | 2c(n^2-11n+32) |
| C5 | 5e-30a+60c-40d | d(3n^3-53n^2+320n-664) |

Every excess is strictly positive for n>=16. The two displayed quadratics
are positive by completing squares, or direct translation n=16+t. The
last cubic becomes `3t^3+91t^2+928t+3176`, also positive. At n=6 the C5
formula instead gives 20, so the inherited hostile is explicitly retained.

It remains to undo the possible global negation in the classification.
For even n it preserves all Hamilton signs, so both globally signed
two-edge matching families give equality. For odd n it replaces w by
`N-w`, where `N=(n-1)!/2` is the total cycle count.

For a matching, (E6) gives `w<N/2`, so its global negative has weight
`>N/2=(n-1)e/4>2e>b`. Every finite motif has at most five edges, so the
union bound gives `w<=5e`, and its global negative has weight at least
`N-5e=[(n-1)/2-5]e>=5e/2>b` for n>=16. This includes all sign cases;
no odd antibalanced or globally negative edge class is silently removed.

Combining (E2)--(E8), the only classes with nonzero weight at most b are
the single-edge minimizers and the stated two-edge equality classes.
This proves (E1) and its full equality classification.

## 5. Multiplicity and remaining boundary

There are three perfect matchings on each four-set, hence `3 binom(n,4)`
labelled two-edge signings. Two such signings cannot be switching-equivalent:
their symmetric difference has at most four edges, whereas a nonempty cut
of K_n has at least n-1>4 edges. Their global negatives are distinct from
them for n>=16, as the negative-triangle counts are respectively
`2(n-2)` and `binom(n,3)-2(n-2)`. These two counts are different. Thus
global negation doubles the equality multiplicity exactly when n is even,
and the different triangle counts also distinguish the two relabelling
orbits in that case.

The proof deliberately uses a conservative cutoff. It does not settle
the exact second weight at orders 9 through 15. The previous asymptotic
bracket remains valid there and is not silently replaced. The new result
also does not settle second weights for fixed short cycle layers: adjacent
and disjoint two-edge perturbations reverse order with the ambient scale,
as recorded in the preceding stability note.

## Reproduction

```bash
python3 04-computation/overnight_hexagon_sep05_hamilton_second.py
python3 -O 04-computation/overnight_hexagon_sep05_hamilton_second.py
```

The independent controls enumerate all Hamilton cycles on K3 through K9,
all admissible matching sizes there, and every listed finite motif. They
verify the inclusion-exclusion and recurrence by exact integers, retain
negative/zero matching signed-sum hostiles, and audit the all-height
endpoint and motif polynomials. No large switching-class census, numerical
integration or imported Walsh source is used. The structural classification
is proved above; finite examples and local capacity checks do not replace it.

Normal and optimized runs have identical output: 15,347 explicit failure
gates, 23,116 literal Hamilton cycles, and exact matching recurrences and
scoped sign checks through n=100. The claimed equality family on K16 is
enumerated separately (5,460 representatives; 10,920 after even global
negation), without enumerating the other `2^105` switching classes.
Raw-LF SHA-256:

```text
6327ac1bb9097b80d45ab9d9be1a0171e70512c8e32de872ab4fffe5d1c9b386  04-computation/overnight_hexagon_sep05_hamilton_second.py
93d5db2082fdf705b3a2f3a5f60e3e508d722e48b25e1fddaa333f44bc9aa001  05-knowledge/results/overnight_hexagon_sep05_hamilton_second.out
```
