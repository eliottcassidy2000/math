# Silent square births are exactly two Pell clocks

Status: PROVED, independently checked in graph-audit/PELL_RESET_ADDENDUM.md. This refines the proved square-sum arrival formula; it is not a Hamiltonian-existence criterion.

Let G_n have vertices1..n and edges between distinct vertices whose sum is a square. The new vertex n has neighbors s²−n for n<s²<2n, so its birth degree is

d(n)=floor(sqrt(2n−1))−floor(sqrt(n)).

For n>=2 both floor increments are0 or1. The lower floor jumps exactly at n=k². The upper floor jumps when n=ceil((m²+1)/2), with m the newly admitted square root. A square birth is 'silent' in degree when these two increments cancel.

## Exact classification

For k>=2, d(k²)=d(k²−1) if and only if some integer m satisfies

m²−2k²=−1 or −2.

Indeed the upper threshold condition is 2k²−2<m²+1<=2k², equivalent to 2k²−3<m²<=2k²−1. An integer square in this interval must be 2k²−1 or2k²−2. Conversely either equality produces the upper jump. There cannot be two such m because consecutive positive squares differ by at least3.

The −1 family is the negative Pell equation, with positive solutions

m+k sqrt2=(1+sqrt2)(3+2sqrt2)^r, r>=0.

The −2 family has m=2a and k²−2a²=1, with

k+a sqrt2=(3+2sqrt2)^r, r>=0.

The r=0 cases give k=1, excluded from the transition claim because d(0) is not defined by the graph convention. All genuine silent square births are therefore the union of

k=5,29,169,... and k=3,17,99,577,...,

squared. The first values are9,25,289,841,9801,28561,... . These are all-order classifications from the standard Pell solution formula, not extrapolations from a search.

## What happens at25

At n=24 the birth neighbors are1 and12, from squares25 and36. At n=25 they are11 and24, from squares36 and49. Thus degree2 persists while one square channel is removed and a new one appears. Even a surviving channel changes its actual neighbor label by−1 when n increases; unchanged channel count is not unchanged incidence. The cancellation is the exact negative-Pell identity7²−2·5²=−1.

At n=9 the analogous positive-Pell identity4²−2·3²=−2 already produces a silent square birth, but G_9 has no Hamiltonian path. Therefore the Pell property alone cannot explain why25 is the permanent Hamiltonian threshold. The required additional structure at25 is the marked ear11–25–24 and the endpoint/port compatibility proved in the companion graph notes.

## Density and scope

Let S(X) count genuine silent square births at most X. With λ=3+2sqrt2, each family grows on the n-scale as a positive constant times λ^(2r). Consequently

S(X)=log(X)/log(λ)+O(1).

They form a zero-density subset even among the square birth times. Away from these Pell clocks, a perfect-square birth lowers the birth degree by one. The complete vertex-arrival graph still changes at every step; 'silent' refers only to this degree difference, not to an unchanged graph.

The governing norm m²−2k² is an indefinite quadratic form. This is an explicit arithmetic mechanism, not a claimed connection to the MRB CM point or Monster representation theory. In particular Q(sqrt2) here is a real quadratic field, distinct from the imaginary quadratic CM field Q(sqrt(−2)).

## Every fixed birth degree has an explicit finite generating polynomial

For any integer j>=0, split by r=floor(sqrt n)>=1 and define

L_(r,j)=max(r², ceil(((r+j)²+1)/2)),
U_(r,j)=min((r+1)²−1, floor((r+j+1)²/2)).

Then d(n)=j exactly when n belongs to one of the nonempty integer intervals [L_(r,j),U_(r,j)]. The intervals for different r are disjoint. Only

1<=r<=floor((1+sqrt2)(j+1))

can contribute: n>=r² and n<=U imply2r²<=(r+j+1)². Conversely the interval formulas themselves impose both floor equalities exactly.

Therefore the complete arrival-degree generating polynomial is

D_j(z)=Σ_(r:L<=U) z^L (1−z^(U−L+1))/(1−z).

At z=1 each summand has removable value U−L+1. This is an explicit finite polynomial for every j, not a fitted recurrence. The first counts D_j(1), j=0..9, are

3,12,24,34,47,59,68,84,92,105.

Thus even the tempting rule12j suggested by degrees1 and2 fails at degree3. The degree-one births are exactly3,5..12,16..18. Degree-two births are exactly13..15,19..32,36..40,49..50, agreeing with the companion graph theorem. This makes the early rises and reversals instances of one exact all-scale interval law.
