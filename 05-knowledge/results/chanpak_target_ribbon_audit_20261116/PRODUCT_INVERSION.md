# Target inversion: a growing-seam obstruction that survives changing the cut

Use the notation and exact formula (★) of [PRODUCT.md in the previous seam packet](https://github.com/eliottcassidy2000/math/blob/a872455acffc668b3f577c1fb3941bc70e1cfd50/05-knowledge/results/chanpak_signed_seam_bank_20261116/PRODUCT.md). This note only concerns that transitive two-block architecture.

For width w and 0≤z<2^(w−1), the reversal satisfies 0≤rev_w(z)≤2^w−2, and 0≤J(z)≤z≤2^(w−1)−1. (For each occupied upper bit i≥1, the geometric contributions from all lower positions sum to at most 2^i−1, and adjacency adds at most one.) Therefore
```
5 − 3·2^(r−1) − 2^(k−1) ≤ C ≤ 2^k − 1.
```
These one-sided bounds give an exact target-first finite algorithm in the family: enumerate the positive integers P with
```
(N − (2^k−1))/2^t ≤ P ≤ (N − (5−3·2^(r−1)−2^(k−1)))/2^t.
```
For each such P, use the divisor-and-digital-remainder test in that note, requiring C=N−2^t P as an integer equality.

**Uniform dyadic counterexample.** For every m≥1, there is no realization of 2^m+3 in this family with t≥max(k,r)+2, regardless of the cut, total tape order, or widths. In fact |C−3|<2^(t−1) by the displayed bounds. If m<t then 2^m is at most 2^(t−1), whereas 2^t P+C−3>2^(t−1) since P≥1. If m≥t, then 2^m−2^t P=C−3 is divisible by 2^t and so must vanish. The left factor in (2^(k−1),2^k] of a dyadic product must be 2^k, and the right factor in [2^r,3·2^(r−1)) must be 2^r. This forces b=a=0 and C=1, contradicting C=3. In contrast 2^m+1 has the empty-mark realization at every such split with k+r+t=m. More generally an odd shift s with |C−s|<2^(t−1) obeys the same cancellation, and among them only s=1 can be present. The leading divisor obstruction is exact even though all widths and the cut may change.

**Density obstruction for moderately growing spacers.** Write n=k+r+t. The left factor is at least 2^(k−1)+1, and the right factor at least 2^r, so H≥2^(n−1)+(2^t−3/2)2^r−2^(k−1)+5≥(3/8)2^n, using t≥1 and k≤n−2. Also H<(3/2)2^n+2^k≤(7/4)2^n. Consequently any H≤X has n≤log₂(8X/3). At order n, for any L≥1, all t≥L together offer at most
```
Σ_{t=L}^{n−2} (n−t−1) 2^(n−t−2) ≤ n·2^(n−L−1)
```
mask assignments, even if counted with multiplicity. Suppose t≥(1+ε)log₂ n for fixed ε>0. Orders n≤½log₂(8X/3) contribute O(X^(1/2)log X) using L=1, and the larger orders contribute O(X/(log X)^ε), summing geometrically. Hence this entire changing-width growing-spacer slice has density zero, despite exact all-height repair banks. This is a necessary budget restriction for this family, not a statement about the other seam architectures or the full tournament spectrum. In particular a density-one route within (★) would have to exploit thin spacers, not just enlarge these exact banks.

## Removing the widths from the inverse problem

There is a finite catalogue per odd part, not just eventual stabilization for a fixed aspect ratio. For odd positive d,e define k₀=ceil(log₂d), r₀=floor(log₂e). Permit a width-zero all-zero mask, with reversal and energy both zero. Set b₀=2^k₀−d and a₀=e−2^r₀, and admit the pair iff e<3·2^(r₀−1). Define C₀(d,e) by the remainder in (★) evaluated at these minimal widths and masks.

Then the **whole** spectrum of (★), allowing every positive k,r,t, is exactly the union of the dyadic rays
```
de·2^s + C₀(d,e),
s ≥ 1 + [k₀=0] + [r₀=0].
```
Indeed the two high factors must be d·2^u and e·2^v for u,v nonnegative. Their normalized interval conditions uniquely force k=k₀+u, r=r₀+v, with u≥[k₀=0] and v≥[r₀=0]. The actual masks are b₀·2^u and a₀·2^v, which preserve both reversal and energy. The exponent is s=t+u+v. Conversely put u=[k₀=0], v=[r₀=0], t=s−u−v; every ray point is physical. The tape length for the ray point is k₀+r₀+s, regardless of how its exponent is distributed. Thus distributing the dyadic budget is pure redundancy for target inversion, not a new repair coordinate. At fixed odd part M only divisors d of M contribute, and taking e=1 proves that every odd positive M has a ray, but only at its prescribed remainder. This does not anchor arbitrary targets.

## The inverse covers every guarded transitive tape with at least three old vertices

In the guarded head model, before appending the sentinel and the head keep the two extreme vertices of the transitive tape unmarked. If any interior vertex is unmarked, declare it the spacer of width one: this gives (★) directly. If all n−2 interior positions are marked, their palindrome makes E−L=0 and J=2^(n−2)−2 for n≥3, yielding H=3·2^(n−2)+3. The singleton mask at old index 1 has E=2, L=2^(n−2), J=0 and gives the *same* H. At n≥4 this leaves another unmarked interior position. At n=3 the empty mask gives the same count 9 instead. Thus the width-free decoder solves membership in the whole guarded-head count family at n≥3, not just those input masks that visibly leave an interior hole. A two-vertex tape is the small exceptional guarded case H=5. This scope does not include all orientations of a transitive tape and an apex with no unmarked penultimate guard.

The replay `invert_product.py` tests the exact target-first algorithm against complete enumeration of 64 small seam dimensions, recovers all three masks in the earlier growing window from their counts, and checks width-free decoding on those small feasible inputs.
