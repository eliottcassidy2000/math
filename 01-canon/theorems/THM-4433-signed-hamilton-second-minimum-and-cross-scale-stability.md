---
id: THM-4433
title: "Signed Hamilton second minimum and cross-scale stability"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 second research wave
---

# THM-4433 -- Signed Hamilton second minimum and cross-scale stability

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**

For an edge signing of K_n let c_k be the number of negative unoriented
simple k-cycles, and e_(n,k)=(n-2)!/(n-k)! the single-edge value. Work
modulo vertex switching. For odd k the zero class is balanced; for even k
the balanced and antibalanced classes are zero.

## Cross-scale first minimum and uniform stability

For every n>=9 and every3<=k<=n, the first positive value of c_k is exactly
e_(n,k), attained only by a single negative edge, together with its global
negative for even k. The
[all-k proof and exact stratum counts](../../05-knowledge/results/overnight_hexagon_sep05_crossscale.md)
extend the transposition mechanism of THM-4427 without a K8 census premise.

Writing d2(n,k) for the next distinct positive value, for n>=10,

```text
lambda_n <= d2(n,k)/e_(n,k)
 <= min(2-2/(n-2), 2-4(k-3)/((n-2)(n-3))),
lambda_n=2(n^2-9n+22)/((n-2)(n-3)).
```

For n=9 the uniform lower constant is36/35. Thus the normalized second
value tends uniformly to two across all cycle lengths. For0<epsilon<1
and n>3+8/epsilon, a positive value at most(2-epsilon)e_(n,k) is forced
to be a first-minimum signing.

A vertex transposition produces a negative K_(2,r). Its normalized negative
k-cycle count is at least its normalized Hamilton count; the exact excess
retains q=r(n-2-r). This makes the Hamilton scale the worst obstruction.
It does not make the two candidate equality shapes scale-independent:
their difference changes sign at n=2k-3.

## Exact Hamilton second minimum

For every n>=16,

```text
d2(n,n)=2(n-4)(n-3)!.
```

Equality consists exactly of two disjoint negative edges, and their global
negative when n is even. There are3*binom(n,4) labelled switching classes
for odd n and6*binom(n,4) for even n, respectively one and two relabelling
orbits. The second signed-cycle operator gap is twice this value.

The [full classification, matching-integral positivity and equality proof](../../05-knowledge/results/overnight_hexagon_sep05_hamilton_second.md)
are part of this theorem. Transposition rigidity reduces signings to
matchings or one of P3,P4,C3,C4,C5. The exact signed matching sum is

```text
W_n(m)=1/2 sum_j binom(m,j)(-4)^j(n-j-1)!
      =1/2 integral_0^infinity exp(-t)t^(n-m-1)(t-4)^m dt.
```

Its positivity and contraction recurrence order the matching weights;
the five remaining motifs have explicit strict excess. Small hostile
values W_6(3)=-4 and W_5(1)=0 are retained, not extrapolated away.

## Exact second minimum at every cycle scale

The [full all-layer proof and independent marked-root audit](../../05-knowledge/results/overnight_hexagon_sep05_crossscale_second.md)
strengthen the Hamilton conclusion to every n>=16 and3<=k<=n:

```text
d2(n,k)=e_(n,k)*min(2-2/(n-2),2-4(k-3)/((n-2)(n-3))).
```

Equality is exactly two adjacent negative edges if n>2k-3, exactly two
disjoint negative edges if n<2k-3, and both types at equality. Global
negation adds classes precisely when **k**, not necessarily n, is even.
The labelled multiplicities are3*binom(n,3),3*binom(n,4), or their sum,
respectively, doubled for even k. Relabelling gives one, one or two orbits,
likewise doubled for even k.

The improved threshold compares each transposition defect with the actual
disjoint-edge benchmark before transporting from Hamilton scale. The
same R2 classification follows at n>=16. Matching monotonicity contracts
the added edge to an **unmatched marked vertex**; induced-subset averaging
and complete short sign tables prove positivity. Forgetting the marked
vertex would admit a false small-layer sign. A literal parity Bonferroni
inequality puts the other three-to-five-edge motifs above2e.

Normal/optimized replay passes17,478 gates and72,238 literal small cycles;
root and independent agent reviews passed the proof and cutoff improvement.
## Exact mixed-parity cumulative seconds

The [full cumulative proof and independent global-sign audit](../../05-knowledge/results/overnight_hexagon_sep05_cumulative_second.md)
close the separate mixed-parity problem. For n>=16 and3<=D<=n-1 put
L=D+1, E=sum_(k=3)^L e_(n,k), and N_k=n!/[2k(n-k)!]. Then

```text
d2_cum=min(A,B,C),
A=(2-2/(n-2))*E,
B=sum_(k=3)^L(2-4(k-3)/((n-2)(n-3)))*e_(n,k),
C=sum_(3<=k<=L,k odd)N_k.
```

Equality consists exactly of the minimizing adjacent-edge pair, disjoint-edge
pair, or antibalanced types, respectively, retaining every type if tied.
Their labelled multiplicities are3*binom(n,3),3*binom(n,4),1, added on
ties. No extra global-negative edge-pair equality survives the odd layers.
The second operator gap is twice this weight.

The third competitor is necessary: n=16,L=10 gives
C=234801776 < B=260444884 < A=261498796, with E=140807044.
For every n>=16 an even L with |L-1-(n+1)/2|<=1 has the antibalanced
class as its unique second minimizer. With delta_n=min_D d2_cum/E,

```text
2n(n-1)/(n+1)^2 <= delta_n <= 2n(n^2-1)/((n+3)(n^2-13)),
lim_(n->infinity)n*(2-delta_n)=6.
```

The proof retains the odd/even split under global negation and explicitly
excludes every other R2 carrier. Normal/optimized finite controls compare
the full classified carrier in1,575(n,L) cases with7,199 gates; absence
of ties in that box is not claimed for all orders. This cumulative result
is a separate proof, not a sum of the individual second-minimum formulas.

The [earlier stability proof](../../05-knowledge/results/overnight_hexagon_sep05_hamilton_stability.md),
exact-second15,347-gate audit, and cross-scale33,510-gate audit are mutually
consistent and independently reviewed. Exact second values at arbitrary
ambient orders below16 and higher-weight classification remain OPEN here.
