# The next Hamilton spectral scale is asymptotically twice the edge minimum

Status: **PROOF CANDIDATE**, independent all-order audit pending. This is a
quantitative consequence of the transposition mechanism, not a larger
switching-class census and not an exact finite-order second-minimum theorem.

## Inheritance and precise carrier

The closest mechanism is Section 2A of
[the all-cumulative transposition-rigidity proof](overnight_hexagon_sep05_d7_d8_gap.md),
prepared for the separately reserved THM-4427. The relevant object is the
binary code obtained by evaluating a switching class on every unoriented
Hamilton cycle of K_n. Its Hamming weight is `w(H)=c_n(H)`. Retain the
zero-class sidecar: for even n the balanced and antibalanced classes both
have weight zero; for odd n only the balanced class does.

The canonical hostile remains the negative C5 with a positive apex on K6,
whose weight 20 is below the edge weight 24. The corrected near miss here
would be claiming the two-disjoint-edge construction is the second minimum:
we prove an asymptotically sharp bracket, not that exact equality statement.
Targeted searches of the current theorem/results corpus found no earlier
signed-Hamilton second-minimum or quantitative transposition-isolation law;
no external priority claim is made.

The concept board is: code Hamming weight; row disagreement; transposition
differences; exact integer endpoints; two-edge intersection; ambient versus
cycle-length scale. The least-used sidecar is that the maximum possible
product `r(n-2-r)` is an **integer** endpoint, with a parity correction.

## 1. An explicit all-order isolation bracket

Let n>=9, let `e_n=(n-2)!`, and remove the zero classes and every class with
weight e_n. Denote the minimum remaining weight by `d2(n)`. Thus this is the
second distinct positive weight, or equivalently half the second distinct
positive eigenvalue of the unnormalized multiplicity-weighted Hamilton
cycle Laplacian. For even n the first-weight classes include both a single
negative edge and its global negative. For odd n they include only the
single-negative-edge classes.

Then

```text
36/35 <= d2(9)/e_9 <= 10/7,

2(n^2-9n+22)/[(n-2)(n-3)] <= d2(n)/e_n
                             <= 2(n-4)/(n-2)       for n>=10.   (S1)
```

In particular

```text
lim_(n->infinity) d2(n)/e_n = 2.                                (S2)
```

Both parity subsequences have this limit. The bracket width in the second
line of (S1) is exactly `4(n-5)/[(n-2)(n-3)]`.

### Lower bound

The row-disagreement classification from the cited research proof says:
unless H is balanced, antibalanced, a single-edge class, or a global negative
of such a class, some pair of vertices has

```text
2 <= r <= n-4.
```

Its transposition difference has the exact weight

```text
T_n(r)=2q[(n-2)(n-3)-2q](n-5)!,    q=r(n-2-r),
```

and the xor triangle inequality gives `w(H)>=T_n(r)/2`. The possible q
range lies between `2(n-4)` and `floor((n-2)^2/4)`. The function
`F(q)=q[(n-2)(n-3)-2q]` is concave, so its minimum is at an endpoint.

For even n, upper-endpoint minus lower-endpoint values factor as

```text
F((n-2)^2/4)-F(2(n-4)) = (n-10)(n-6)^2(n-4)/8.
```

For odd n, the exact integer upper endpoint gives

```text
F(((n-2)^2-1)/4)-F(2(n-4))
  = (n-7)(n-5)(n^2-14n+41)/8.
```

The first expression is nonnegative for even n>=10. In the second, writing
n=11+t turns the last factor into `t^2+8t+8`, strictly positive for t>=0.
Thus for every n>=10 the minimum is at r=2 (or n-4), giving the lower
ratio in (S1). At n=9 the two q endpoints are 10 and 12, with F values
220 and 216, so the lower ratio is `216*4!/7!=36/35`.

For odd n the excluded structural families also contain the antibalanced
class and globally negated single edges, but those have weights
`(n-1)e_n/2` and `(n-3)e_n/2`, respectively. They exceed `2e_n` for n>=9
and therefore do not compromise either lower bound for d2(n).

### Upper bound and sharp asymptotic scale

Use two disjoint negative edges. Each edge lies in e_n Hamilton cycles.
Exactly `2(n-3)!` cycles contain both: contract both edges and account for
their two independent orientations. A Hamilton cycle is negative precisely
when it contains one of the two edges, so this class has weight

```text
2e_n-4(n-3)! = [2(n-4)/(n-2)]e_n.                    (S3)
```

It is not a zero or first-weight class. For example, its number of negative
triangles is `2(n-2)`, whereas those classes have zero, `n-2`,
`binom(n,3)` or `binom(n,3)-(n-2)` negative triangles, all different for
n>=9. Thus it is a legitimate upper witness for d2(n). Taking limits of
the two rational bounds proves (S2).

## 2. A threshold form of stability

For every fixed `0<epsilon<1`, if

```text
n > 3+8/epsilon
```

then every nonzero Hamilton codeword satisfying

```text
w(H) <= (2-epsilon)e_n
```

is a first-weight single-edge class (allowing global negation when n is
even). Indeed the lower bound equals

```text
2 - 8(n-4)/[(n-2)(n-3)] > 2 - 8/(n-3) > 2-epsilon.
```

This is exact switching-class rigidity below a quantitative threshold, not
a claim about the number of negative edges in an arbitrarily chosen gauge.
No general edit-distance estimate above the threshold is asserted.

## 3. The near-minimizer operation changes with the scale

For a fixed k-cycle layer on K_n, k>=3 and n>=max(k,4), let
`e_(n,k)=(n-2)!/(n-k)!`. Two negative edges sharing a vertex have weight

```text
2e_(n,k) - 2(n-3)!/(n-k)!.
```

For two disjoint edges, the weight is

```text
2e_(n,k) - 4(k-3)(n-4)!/(n-k)!,
```

where the second term is zero at k=3. These are exact inclusion-exclusion
counts, since there are only two negative edges. Their difference is

```text
w(adjacent pair)-w(disjoint pair)
 = [2(2k-n-3)/((n-2)(n-3))] e_(n,k).                 (S4)
```

Hence the disjoint pair is the better of these two competitors in the
Hamilton regime k=n>=9, but the adjacent pair is better when `n>2k-3`.
At `n=2k-3` they tie. This identifies an actual change in the perturbation
type under the operation of enlarging the ambient graph while holding
cycle length fixed. It does not identify the global second minimum in
either regime. Forgetting cycle length would lose the coordinate that
determines this reversal.

## Reproduction and boundary

```bash
python3 04-computation/overnight_hexagon_sep05_hamilton_stability.py
python3 -O 04-computation/overnight_hexagon_sep05_hamilton_stability.py
```

The source checks exact polynomial factorizations, all-height positivity
after translation, the squeeze limits and the scale transition. A separate
literal enumeration of all simple cycles at every length on K5 through K9
checks both two-edge intersection formulas and their parity weights. It
uses no large switching-class census and does not import the earlier
Walsh implementations. The all-order lower bound rests on the analytic
transposition and classification proof; finite checks are only controls.

Still open here: the exact finite-order second minimum and its equality
classes, an edit-distance estimate beyond the one-edge threshold, and
stability for cumulative or fixed-short-layer operators. No LRC(14),
tournament Hamilton-path, or Booleanized quotient claim follows.

Normal and optimized runs produced identical 85-gate output, including
72,238 directly enumerated simple cycles. Raw-LF SHA-256:

```text
8dc0ed295a591ecfb03c425738c243f738582bbe1b6bce5fff5ca2a6e54cd3b2  04-computation/overnight_hexagon_sep05_hamilton_stability.py
a0f586ff797272e083e5860def9d427a279bec6f6609ca6dd97c0a1984b26afe  05-knowledge/results/overnight_hexagon_sep05_hamilton_stability.out
```
