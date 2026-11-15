# When two fixed constructions agree after an arbitrary prefix

**Scope.** These statements concern comparing fixed pairs of cores and fixed suffixes, allowing a common prefix or a fixed dyadic rescaling of that prefix on either arm. Changing a core with the absolute target is not ruled out.

Let a>0 be the value of the prefix, written without leading zeros, and use the exact binary apex recurrence
```
p(0)=0,    p(2z+1)=p(z),    p(2z)=p(z)+z+(z mod 2).
```
This p is the augmentation in the [earlier suffix decoder](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md). By reverse-complement it is also the augmentation of the apex over the ordinary binary word of value a with leading bit 1. For an arm with fixed suffix v, fixed profile pair and constants d,t,c,D from [ABSOLUTE_CARRY](ABSOLUTE_CARRY.md), the exact readout is
```
N(a) = (D−t)a + t p(a)+c + t[a odd and v starts 0],
```
with the last term zero for an empty suffix. This is the seam correction for one adjacent inversion: all non-adjacent cross inversions give the product a(2^k−1−rev_k(b)). Thus each prescribed difference between two fixed arms, after fixing prefix parity, reduces to an affine line in the (a,p(a)) plane. A fixed dyadic rescaling on either arm obeys the same statement, because the recurrence peels any appended zeros.

## 1. Exact line search and a sparse obstruction

Let ρ(a) be the bit-reversal fraction 0.bits_lsb_first(a), and L the bit length. From the recurrence, summing over zero digits gives
```
|p(a)−a(1−ρ(a))| ≤ L+2.
```
Each replaced floor and adjacency carry costs at most one, and the discarded geometric tail costs at most one. In particular for any constants η,C, a solution of p(a)−η a=C in the width-L shell must have
```
|ρ(a)−(1−η)| ≤ (L+2+|C|)/2^(L−1).
```
There are O_C(L) grid points of spacing 2^(−L) in that interval. Reverse each of those numerators and test the exact equality. This yields O_C(log² X) candidates for a≤X, with an explicit decoder. If the p-coefficient vanishes in the original difference, one has instead a purely linear equation in a. Only cancellation of both prefix coefficients can give a persistent difference on a positive-density set of prefix addresses for fixed arms. This deduction does not require any stochastic assumption about carries.

## 2. Dyadic slopes: every long solution is on a finite set of exact dyadic rays

Assume 1−η=j/q with q=2^r, r≥1, 0≤j≤q. If η lies outside [0,1] there are only finitely many solutions by 0≤p(a)≤a. For large widths where q(L+2+|C|)/2^(L−1)<1, the integer floor(qρ(a)) must be j or j−1. Write a=qz+s. Then floor(qρ(a))=rev_r(s), and the exact split of the recurrence says
```
p(a)−η a = p(s)−η s + p(z) + Δ_s(z mod 2) + (j−1−rev_r(s))z,
Δ_s=1−highest_bit(s).
```
Thus the two allowed width-r low words force one of the two equations
```
z−p(z)−Δ_s(z mod 2) = constant,     or
p(z)+Δ_s(z mod 2) = constant.
```
Endpoints j=0,q just omit an out-of-range branch. Nonintegral right-hand sides cannot have a solution.

The first equation has only finitely many odd solutions, since p(2y+1)=p(y)≤y. For even z>0 write z=2^ell u with u odd and ell≥1. Then z−p(z)=u−p(u)−1≥(u−1)/2, independent of ell, so finitely many u suffice and each accepted u fills its entire dyadic ray. Handle z=0 separately. In the second equation even z has p(z)≥z/2. Every odd z has z=2^ell(u+1)−1 for a unique even u≥0 and ell≥1, found by stripping trailing ones. Then p(z)=p(u), again leaving only finitely many u and entire rays. This classifies all long solutions exactly as finitely many rays a=A·2^ell+B, plus finite exceptional points, with explicit finite searches. Fixed parity can only discard rays or initial terms. It covers dyadic comparison slopes, not arbitrary rational slopes more sharply than §1.

## 3. A shifted physical +2 pair obtained backward from a nonzero defect

For even prefixes define the defect E=a−p(a). This is unchanged when appending further zeros. In particular E=2 holds on both a=3·2^ell and a=5·2^ell, ell≥1: their odd roots have p(3)=0,p(5)=2. Use:
```
arm 1: U=TT3, W=TT1, empty suffix: (D,t,c)=(15,1,9)
arm 2: U=TT1, W=C3,  suffix 0:    (D,t,c)=(30,3,15)
```
Put prefix 2a on arm 1 and a on arm 2. Both are even, so the adjacency correction vanishes. Both have the same nominal 30a, and their *exact* low corrections are 9−E=7 and 15−3E=9. Thus these two actual tournament constructions give
```
90·2^ell+{7,9},     150·2^ell+{7,9},     ell≥1.
```
A target in those families selects a physical prefix from its exact value by dividing N−7 or N−9 by 90 or 150 and checking a power of two. This uses a nonzero deficit rather than the old prefix 1 followed by zeros. At ell=1,2,3 the correction 9 already generates a nonzero carry at the indicated scale. The all-order result is an exact recentering on two specific shifted rays, not arbitrary carry coverage.
