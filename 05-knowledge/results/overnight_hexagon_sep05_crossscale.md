# Hamilton transposition differences dominate every normalized cycle layer

Status: **PROVED ANALYTICALLY + INDEPENDENTLY AUDITED**; all stated finite
controls are **FINITE-EXACT**. Root and the independent `lrc_components`
referee separately reconstructed the incidence counts, normalized
comparison, parity/equality cases and uniform second-weight quantifiers.
The new result is uniform in cycle length, not an inference from a few
additional profile vectors.

## Inheritance and the missing stratum

The source is the transposition difference from
[THM-4427, all-cumulative signed-cycle gaps by transposition rigidity](../../01-canon/theorems/THM-4427-all-cumulative-signed-cycle-gaps-by-transposition-rigidity.md).
The closest continuation is the quantitative isolation bracket in
[the Hamilton stability note](overnight_hexagon_sep05_hamilton_stability.md),
followed by the n>=16 exact second Hamilton weight in
[the matching-rigidity note](overnight_hexagon_sep05_hamilton_second.md).

The canonical hostile is still the negative C5/apex signing on K6, so an
all-order claim below the stated cutoff is forbidden. The new cheapest
hostile is simpler: in a triangle layer, a negative K_(2,r) has **all** its
negative cycles in the stratum containing exactly one of the two distinguished
vertices. A Hamilton-only enumeration would count zero there.

The live board is: cycle length versus ambient order; exact one-vertex
stratum; two-star difference; normalized code weight; row-type rigidity;
adjacent/disjoint two-edge perturbations. The map keeps the cycle and its
incidence with the two transposed vertices. Forgetting that incidence loses
the stratum that makes the all-length formula correct. The preserved
predicate is exact cycle-sign parity, not a surrogate connectivity statistic.

## 1. Exact transposition count at arbitrary cycle length

Let n>=5, `3<=k<=n`, and consider the signing whose negative edges are
exactly those of K_(2,r) between vertices u,v and a set R of r of the other
n-2 vertices. Put

```text
s=n-2-r,       q=rs,
e_(n,k)=(n-2)!/(n-k)!,
T_(n,k)(r)=number of negative unoriented simple k-cycles.
```

For k=3,4,5,

```text
T_(n,k)(r)/(2e_(n,k)) = q/(n-2).                     (C1)
```

For 6<=k<=n,

```text
T_(n,k)(r)/(2e_(n,k))
 = q[(n-k+2)(n-4)(n-5)
        +(k-5)(n^2-7n+14-2q)]
     /[(n-2)(n-3)(n-4)(n-5)].                       (C2)
```

### A complete incidence partition

Write `(a)_j=a(a-1)...(a-j+1)`, with `(a)_0=1` and `(a)_j=0` when
the nonnegative integer a is smaller than j.

Cycles containing neither u nor v are positive. Cycles containing exactly
one contribute

```text
2q(n-4)_(k-3).                                      (C3)
```

Indeed, choose the distinguished vertex in two ways. Root an oriented cycle
there, choose its first and last neighbours in opposite R/S classes in
2rs ways, order the remaining k-3 vertices, and divide by reversal two.
This term vanishes precisely at the Hamilton boundary k=n (apart from
the trivial q=0 cases).

For cycles containing both u and v, root an oriented cycle at u and retain
the position of v. When they are adjacent, division by reversal gives
`2q(n-4)_(k-4)`. When their cyclic distance is two and k>=5, the common
neighbour cancels in the sign product. Summing that common neighbour's
choices gives a second identical contribution.

When both cyclic distances are at least three, k>=6, there are k-5 possible
positions for v. The four distinct neighbour slots have an odd number of
R entries in exactly

```text
4q[(r-1)(r-2)+(s-1)(s-2)]
```

ways. After ordering the remaining vertices and dividing by reversal,
this contribution is

```text
2q[(r-1)(r-2)+(s-1)(s-2)](k-5)(n-6)_(k-6).          (C4)
```

At k=4 the opposite position has repeated neighbours and is always positive,
so only the adjacent contribution `2q` remains among the cycles containing
both vertices. At k=3 that entire stratum is positive. These two small
boundary cases and the k=5 adjacent/distance-two count give (C1).

For k>=6, adding (C3), both short-distance terms and (C4), then using

```text
(r-1)(r-2)+(s-1)(s-2)=n^2-7n+14-2q,
```

gives (C2). Thus neither the one-vertex stratum nor a repeated-neighbour
configuration has been suppressed.

## 2. The normalized Hamilton layer is the worst one

For this comparison assume n>=6. Let
`R_(n,k)(q)=T_(n,k)(r)/(2e_(n,k))`. At k=n, (C2) gives

```text
R_(n,n)(q)=q[(n-2)(n-3)-2q]/[(n-2)(n-3)(n-4)].
```

For k>=6, direct subtraction yields

```text
R_(n,k)-R_(n,n)
 = 2q(n-k)(q-n+3)/[(n-2)(n-3)(n-4)(n-5)].          (C5)
```

For k=3,4,5, the difference is

```text
2q(q-n+3)/[(n-2)(n-3)(n-4)].                        (C6)
```

Therefore, whenever `q>=n-3`,

```text
R_(n,k)(q)>=R_(n,n)(q) for every 3<=k<=n.            (C7)
```

In particular the interior row-disagreement range `2<=r<=n-4` has
`q>=2(n-4)>=n-3` throughout this range. For these witnesses, increasing cycle length
can only weaken the normalized transposition bound. This is the exact
cross-scale comparison; a raw unnormalized count would not have this form.

## 3. Every individual cycle layer has the edge minimum once n>=9

For every n>=9 and every `3<=k<=n`, let

```text
Z_k={balanced}                         for odd k,
Z_k={balanced,antibalanced}            for even k,

E_k={single-edge switching classes}    for odd k,
E_k={single-edge classes and their global negatives} for even k.
```

Then the full minimum and equality statement is

```text
c_k(H)>=e_(n,k) for H outside Z_k,
equality exactly on E_k.                              (C8)
```

To prove it, transposing u,v in any H again gives the K_(2,r) difference,
and xor parity implies `c_k(H)>=T_(n,k)(r)/2`. If an interior disagreement
exists, (C7) and the strict Hamilton estimate in THM-4427 imply
`c_k(H)>e_(n,k)`. If no interior disagreement exists, the proved row-type
classification gives precisely balanced, antibalanced, single-edge and
globally negative single-edge classes.

The parity audit completes the proof. For even k the two zero families
have weight zero and both signed-edge families have weight e_(n,k).
For odd k, let `N_(n,k)=n!/[2k(n-k)!]`. The globally negative families have
weights N_(n,k) and N_(n,k)-e_(n,k). Their normalized total is

```text
N_(n,k)/e_(n,k)=n(n-1)/(2k)>=(n-1)/2>=4,
```

so neither can be a zero or first-weight class. This proves (C8), including
the complete zero set, with no finite base or deletion induction in n>=9.
In particular the short individual layers k=3 through seven are covered
uniformly as well as every longer layer.

Together with the existing antibalanced odd/even pairing, (C8) also gives
a direct analytic proof of every cumulative gap for n>=9. The older
finite bases and deletion proof remain needed for the smaller ambient
orders and retain their historical scope.

## 4. Uniform second-weight isolation across all cycle lengths

Let d2(n,k) be the minimum c_k after excluding Z_k and E_k, and set

```text
lambda_9=36/35,
lambda_n=2(n^2-9n+22)/[(n-2)(n-3)] for n>=10.
```

For every n>=9 and every `3<=k<=n`,

```text
lambda_n <= d2(n,k)/e_(n,k)
 <= min(2-2/(n-2), 2-4(k-3)/[(n-2)(n-3)]).          (C9)
```

For the lower bound, any class outside the four row-rigid families has
an interior disagreement, and (C7) transfers exactly the Hamilton endpoint
lower bound from the stability note. For odd k, the two additional globally
negative rigid families have weight at least 3e_(n,k), hence exceed this
bound. For the upper bound, use respectively adjacent and disjoint pairs
of negative edges, whose exact weights were counted in that same note.
Both have weight greater than e_(n,k) when n>=9, so neither is an excluded
zero or first-weight class.

Consequently the convergence to two is **uniform** over all cycle lengths:

```text
sup_(3<=k<=n) |d2(n,k)/e_(n,k)-2|
 <= 8(n-4)/[(n-2)(n-3)] ->0                 for n>=10.             (C10)
```

This includes fixed k, Hamilton k=n, and every joint scaling k=k(n).
For each fixed `0<epsilon<1`, the condition `n>3+8/epsilon` forces every
nonzero k-cycle codeword with weight at most `(2-epsilon)e_(n,k)` to lie
in E_k, simultaneously for all k.

The exact finite-order second weight outside the Hamilton regime remains
open here. In particular this argument does not erase the adjacent/disjoint
perturbation reversal at `n=2k-3`, and does not transport the exact Hamilton
equality classification to all k.

## Reproduction and explicit universe

```bash
python3 04-computation/overnight_hexagon_sep05_crossscale.py
python3 -O 04-computation/overnight_hexagon_sep05_crossscale.py
```

Every simple unoriented cycle on K5 through K9, at every length and every
split r, is enumerated literally. Counts are compared separately in the
one- and two-distinguished-vertex strata before comparing the total formula.
The nonzero triangle one-vertex stratum is an explicit hostile to carrying
the Hamilton-only count into the shorter layers. Rational symbolic identities
verify (C5), (C6), the k=5 and Hamilton boundaries, and the endpoint formula;
additional exact comparisons through n=150 corroborate uniform dominance.
The all-order proof is the displayed algebraic identity, not the finite grid.

Normal and optimized runs have identical output: 33,510 explicit gates and
72,238 directly enumerated cycles. Raw-LF SHA-256:

```text
ccff310e8a7b5cadd309cd77b9ca9f13a8ac95c1b73b4990d77326c6fe309b3c  04-computation/overnight_hexagon_sep05_crossscale.py
20f8ab0614eba4992719676a73f869c9656caf07fc60cbd89126a5f811da158d  05-knowledge/results/overnight_hexagon_sep05_crossscale.out
```
