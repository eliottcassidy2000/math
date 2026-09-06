# Exact cumulative second cycle weights: the antibalanced competitor survives mixing parity

Status: **PROVED ANALYTICALLY + INDEPENDENTLY AUDITED**. Root and the
independent `nc2_seed` referee audited the complete proof, including every
global-sign case. The latter separately audited the sharp constant-six
and central-winner corollaries below. The finite structural controls are
**FINITE-EXACT**, not a substitute for the all-order proof. This extends
THM-4433 without a new theorem identifier.

## 1. Statement and inheritance

Let H be a signing of the labelled complete graph K_n, modulo cut
switching, and let c_k(H) count its negative unoriented simple k-cycles.
For n>=16 and 3<=D<=n-1, set L=D+1 and

```text
S_(n,D)(H)=sum_(k=3)^L c_k(H),
e_k=(n-2)!/(n-k)!, N_k=n!/[2k(n-k)!],
E=sum_(k=3)^L e_k, D2=(n-2)(n-3),
A=[2-2/(n-2)]E,
B=sum_(k=3)^L [2-4(k-3)/D2]e_k,
C=sum_(3<=k<=L, k odd) N_k.                         (1)
```

The zero class is balanced; the first positive weight is E, attained
exactly by the single-negative-edge classes. After removing these zero
and first classes, the exact minimum is

```text
d2_cumulative(n,D)=min(A,B,C).                       (2)
```

The equality classes are precisely the following competitors whose value
attains the minimum in (2): two adjacent negative edges (value A), two
disjoint negative edges (value B), and the antibalanced class (value C).
Thus their labelled multiplicities are respectively `3 binom(n,3)`,
`3 binom(n,4)`, and `1`, added if there are ties. Each listed type is one
relabelling orbit. No claim that ties never occur is needed or made.

Unlike a single even layer, the cumulative problem does **not** give the
global negatives of both two-edge types as extra equality classes. Odd
layers retain the distinction. The antibalanced competitor can genuinely
win: at n=16, D=9 (L=10),

```text
E=140807044,
C=234801776 < B=260444884 < A=261498796.              (3)
```

For odd L, C>2E and cannot win. The minimum formula retains the exact
even-L boundary without claiming a simpler all-order phase diagram.

The closest proved mechanisms are
[THM-4427, all cumulative signed-cycle gaps by transposition rigidity](../../01-canon/theorems/THM-4427-all-cumulative-signed-cycle-gaps-by-transposition-rigidity.md)
and [THM-4433, signed Hamilton second minimum and cross-scale stability](../../01-canon/theorems/THM-4433-signed-hamilton-second-minimum-and-cross-scale-stability.md),
with the fully audited extension
[exact second weights in every cycle layer](overnight_hexagon_sep05_crossscale_second.md).
The detailed R2 carrier proof is in
[the Hamilton second-weight proof](overnight_hexagon_sep05_hamilton_second.md).
The new sidecar is the **odd/even weighted split under global negation**.
Passing from each individual second minimum to their sum loses that
coordinate and does not prove (2).

The live board is: normalized transposition rigidity; complete R2 carrier;
strict matching monotonicity; adjacent/disjoint scale reversal;
odd/even pairing; the antibalanced class. The corrected near miss is the
two-competitor formula refuted by (3). The inherited hostile controls
`W_3(1)=-1`, `W_5(1)=0`, `W_6(3)=-4`, and the negative-C5/apex example
remain relevant: matching positivity uses the proved large ambient order,
not a false assertion for every small Hamilton graph. Targeted inheritance
recovers these mechanisms instead of rerunning a larger switching census.

## 2. One transposition again forces the complete finite carrier

For a pair u,v let r be the number of vertices x outside the pair on
which the ux and vx signs disagree. Write q=r(n-2-r), and let
T_(n,k)(r) be the negative-cycle weight of H plus its u,v transposition
(addition means symmetric difference of negative edge sets). The XOR
triangle inequality gives

```text
c_k(H)>=T_(n,k)(r)/2.
```

For 3<=r<=n-5 the exact benchmark-aware transport from the all-layer
second-weight proof gives, at n>=16 and k>=6,

```text
T_(n,k)(r)/(2e_k)>2-4(k-3)/D2.
```

For k=3,4,5 its left side is q/(n-2), which is at least
`3(n-5)/(n-2)>2` and hence also strictly above the displayed benchmark.
Summing therefore proves `S(H)>B` for every such pair.

Any H whose cumulative weight is at most min(A,B,C) must consequently
have every pair disagreement in

```text
{0,1,2,n-4,n-3,n-2}.                                (4)
```

The complete, independently audited R2 classification applies already
for n>=13: after switching and possibly global edge negation, H has an
arbitrary negative matching, or one of

```text
P3, P4, C3, C4, C5,
```

with unused vertices isolated. This is a structural theorem for all
signings, not an assumption that extremizers have sparse representatives.
We evaluate all of its global-sign cases below.

## 3. The representatives before global negation

The all-layer theorem proves that the negative-k-cycle weight of a
matching is strictly increasing with its number m of edges, for every
n>=16 and 3<=k<=n. Its proof retains an unmatched marked vertex after
contraction and checks all short signed Hamilton tables. Hence its
cumulative weight is strictly increasing in m. The matching m=0 is
balanced, m=1 is the first class, and m=2 has weight B; all m>=3 exceed B.

The path P3 has weight A. For each remaining motif P4,C3,C4,C5, let f be
its number of negative edges, so f=3,4,5. The parity Bonferroni estimate
from the all-layer proof is

```text
c_k(H)/e_k >= f-2f(f-1)/(n-2)>2                    (5)
```

throughout n>=16. It applies also when k is shorter than the motif or
when the motif is itself a complete cycle. Thus these motifs all have
cumulative weight >2E, whereas A<2E. The only unnegated second candidates
are the adjacent and disjoint pairs.

## 4. Global negatives: an explicit odd/even audit

Global edge negation sends c_k to N_k-c_k for odd k and leaves it
unchanged for even k. Pair each odd length x with the next even length
x+1, beginning with (3,4). If L is odd there is one final unpaired odd
length. For a pair put a=e_x. Then

```text
e_(x+1)=(n-x)a,
N_x/a=n(n-1)/(2x).                                  (6)
```

We show every globally negated carrier type except the empty matching
has weight >2E. The globally negated empty matching is precisely the
antibalanced class and has weight C; it must not be discarded.

### A single negative edge before negation

Its pair contribution after negation is `N_x-a+(n-x)a`. This exceeds
`2(a+(n-x)a)` if

```text
n(n-1)>2x(n+3-x).
```

The right side is at most `(n+3)^2/2`, and
`n(n-1)-(n+3)^2/2=(n-9)(n+1)/2>0`. The unpaired odd contribution, if
present, is N_L-e_L>2e_L since N_L/e_L>=(n-1)/2>=15/2.

### Two adjacent negative edges before negation

Its odd weight is at most 2a and its even weight is exactly
`alpha e_(x+1)`, where `alpha=2-2/(n-2)`. The pair contribution is
therefore at least `N_x-2a+alpha(n-x)a`. To exceed twice the pair's
first weight it suffices that

```text
N_x/a > 4+2(n-x)/(n-2).
```

The left side is at least n/2>=8 because x<=n-1; the right side is
strictly below 6 because x>=3. The unpaired odd contribution is at least
N_L-2e_L>2e_L.

### A matching with at least two negative edges before negation

The proved signed matching positivity gives odd negated weight strictly
greater than N_x/2. Strict matching monotonicity gives even weight at
least the disjoint-pair weight, which is at least `beta e_(x+1)`, with
`beta=2-4/(n-2)`. Thus a sufficient paired inequality is

```text
N_x/2+beta(n-x)a>2(a+(n-x)a),
```

equivalently

```text
n(n-1)(n-2)>8x(3n-2-2x).                            (7)
```

The right side is at most `(3n-2)^2`, since their difference is
`(4x-3n+2)^2`. The remaining gap is

```text
n^3-12n^2+14n-4;
at n=16+t this equals t^3+36t^2+398t+1244>0.
```

An unpaired odd contribution is >N_L/2, while
`N_L/(2e_L)>=(n-1)/4>2`. This excludes every globally negated matching
with m>=2, without an upper bound on the matching size independent of n.

### The motifs with three to five negative edges before negation

Their even weights exceed 2e_k by (5). Each odd negated weight is at
least `N_k-5e_k`, since the union bound gives c_k<=5e_k before negation.
But `N_k/e_k>=(n-1)/2>=15/2`, so this lower bound exceeds 2e_k as well.
This completes the global-sign classification.

## 5. Equality, the surviving antibalanced scale, and uniform isolation

Sections 2--4 show that every signing with weight at most min(A,B,C)
belongs to exactly one of the three stated candidate types, and each
candidate has its stated value. The labelled multiplicities follow by
choosing the three vertices and middle vertex of P3, or the four vertices
and one of their three perfect matchings. Distinct representatives with
at most two edges cannot be switching-equivalent: their symmetric
difference has at most four edges, whereas a nonzero K_n cut has at
least n-1>4. Nor can an antibalanced signing have such a representative:
it would require every triangle negative, while two edges meet at most
2(n-2)<binom(n,3) triangles. Thus no counted classes alias.

For an odd/even pair, (6) gives

```text
N_x/(e_x+e_(x+1))
 =n(n-1)/[2x(n+1-x)] >= gamma_n,
gamma_n=2n(n-1)/(n+1)^2.                            (8)
```

The unpaired odd ratio, if any, is at least `(n-1)/2>gamma_n`. Hence
C>=gamma_n E. The disjoint weight satisfies B>=beta E, and

```text
beta-gamma_n=2(n^2-9n-4)/[(n-2)(n+1)^2]>0          (9)
```

for n>=16; A>beta E. Thus all three candidate values exceed E, as
required for the genuine second distinct positive weight, and

```text
gamma_n <= d2_cumulative(n,D)/E < 2,
sup_(3<=D<=n-1) |d2_cumulative(n,D)/E-2|
 <= (6n+2)/(n+1)^2 -> 0.                           (10)
```

This is uniform even when the truncation D changes with n.

Finally, if L is odd then C>=N_L, while

```text
E/e_L <= sum_(j=0)^(L-3) 1/j! < 3,
N_L/e_L >= (n-1)/2 >= 15/2.
```

Consequently C>2E. The antibalanced win in (3) is therefore genuinely a
mixed-parity, even-top phenomenon, not an omitted global symmetry of an
individual layer. An exact simpler description of its even-top winning
window remains a possible refinement of the explicit formula (2).

### The worst-scale asymptotic constant is exactly six

Set `delta_n=min_(3<=D<=n-1) d2_cumulative(n,D)/E`. Choose an even top L
such that `|L-1-(n+1)/2|<=1`; such a legal top exists for every n>=16.
Write t=n-L, so t>=(n-5)/2. For odd j>=5 with j<=L-1,

```text
N_(j-2)/N_j=j/[(j-2)(n-j+1)(n-j+2)]
 <=3/[(t+2)(t+3)] <=12/(n^2-1).
```

Summing the finite geometric tail of odd total-cycle counts gives
`C<=N_(L-1)/(1-12/(n^2-1))`. Meanwhile E is at least the first-weight sum
of the last odd/even pair, and its ratio to N_(L-1) obeys (8). The distance
of L-1 from (n+1)/2 is at most one, giving

```text
gamma_n <= delta_n <= C/E <= U_n,
U_n=gamma_n/[(1-4/(n+1)^2)(1-12/(n^2-1))]
   =2n(n^2-1)/[(n+3)(n^2-13)].                     (11)
```

All denominators are positive. The lower and upper bounds both have
expansion `2-6/n+O(n^(-2))`, so

```text
lim_(n->infinity) n(2-delta_n)=6.                    (12)
```

Thus the worst cumulative truncation has a different sharp first-order
isolation loss from the Hamilton layer's disjoint-edge value
`2-4/(n-2)`. Mixing parity changes which configuration controls stability.

In fact the chosen central even top has the antibalanced class as its
**unique** second minimizer for every n>=16. Because L<=(n+5)/2, both
pair weights have normalized value at least

```text
eta_n=2-2(n-1)/[(n-2)(n-3)].
```

For B this follows by bounding its positive weighted average of decreasing
layer ratios by its top ratio; A/E>eta_n directly. On the other hand,

```text
eta_n-U_n
 =2(2n^4-29n^3+55n^2+149n-273)
  /[(n-3)(n-2)(n+3)(n^2-13)]>0.                    (13)
```

The numerator at n=16+v is
`4v^4+198v^3+3470v^2+24810v+56958`, which is positive for v>=0.
Hence `C/E<=U_n<eta_n<=min(A/E,B/E)`. This proves an all-order central
winning region, without asserting the exact endpoints of the whole
antibalanced winning window.

## 6. Reproducible exact controls and scope

Run from the repository root:

```bash
python3 04-computation/overnight_hexagon_sep05_cumulative_second.py
python3 -O 04-computation/overnight_hexagon_sep05_cumulative_second.py
```

The explicit finite universe is every n=16,...,60, every top L=4,...,n,
every matching size 0,...,floor(n/2), and all five finite motifs, with
both global signs. Only the balanced and first single-edge classes are
removed. The checker independently evaluates literal edge-subset cycle
inclusion-exclusion, including the proper/full-cycle boundary, rather
than reusing the proof's induced-subset matching averages or inequalities.
It compares the complete carrier's minimum and full equality-type set
with (2), not just the three proposed witnesses. It also checks first
separation, (10), every odd-top exclusion, the hostile (3), and exact
symbolic identities used in the paired bounds. It performs no switching
class census and does not assume that a finite carrier check proves the
all-order R2 reduction.

Normal and optimized runs agree: 1,575 (n,L) cases, 186 strict
antibalanced winners, 56 central even-top witnesses for (11)--(13), no
multiple-type ties in this finite box, and 7,199 explicit failure gates.
The sharp constant-six limit and all-order central strict margin are also
checked symbolically. The stored output is
[the exact control transcript](overnight_hexagon_sep05_cumulative_second.out).
The absence of ties in that box is not promoted to an unbounded claim.

Connection map: the source is the proved individual-cycle R2 mechanism;
the target is a mixed-scale cumulative second minimum; the map is a
positive sum of cycle weights and their transposition lower bounds. It
preserves the common carrier and matching monotonicity but destroys a
single parity's global-negation symmetry. The retained odd/even sidecar
repairs that loss. The cheapest decisive test was the explicit n=16,
L=10 antibalanced hostile, followed by full carrier controls. The new
uniform lower bound then comes from the same pairing geometry.
