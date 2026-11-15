# Two interacting tails: divisor-volume growth and an absolute inverse

The [two-sided recontact](TWO_SIDED_RECONTACT.md) has a different size budget from the one-sided windows. It retains two independent physical binary addresses x,y; no approximation to the count is needed. For any fixed branch, write
```
N(x,y)=c+Qxy+t x+s y+h[p(x)+p(y)]+f(x mod 2)+e(y mod 2).
```
In the displayed pair and triple all coefficients are positive. This note distinguishes *address volume*, *finite congruence coverage*, and *distinct absolute coverage*.

## Address volume exceeds linear

For x,y≥1 and a sufficiently large fixed integer B, the exact inequalities 0≤p(z)≤z imply
```
Qxy ≤ N(x,y) ≤ Q(x+B)(y+B).
```
The elementary divisor summation gives the number of pairs of positive integers with product at most Y as Y log Y+O(Y): sum floor(Y/x) for x≤Y and use the harmonic estimate. Discarding the first B rows or columns costs O_B(Y). It follows that the number of distinct **address pairs** with N(x,y)≤X on one branch is
```
(X/Q) log X + O(X)
```
as X→∞, with natural log. The left-only and right-only axes add O(X). These are presentations, not counts; the product term makes collisions a serious separate question. This evades a linear *presentation* bottleneck but not by itself an absolute density bound.

## Exact absolute-target decoder

Given N, at least one of x,y is at most sqrt(N/Q) when both are positive. Handle either zero axis and try the possible short address on each side. For example, at fixed x the unknown y solves
```
Y=A y+h p(y)+e(y mod 2),   A=Qx+s,
Y=N−c−t x−h p(x)−f(x mod 2).
```
If Y<0 the branch fails. Otherwise split y=2^r u+b, 0≤b<2^r, with r≥1. Use the exact bit reversal split in the [binary suffix decoder](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md), absorbing the parity term into the known low b. At a fixed b, the possible interval of u has real length O(N/(A 2^r)^2+1), with an absolute constant depending only on the branch. Choose 2^r comparable to max(2,sqrt(N)/A); this leaves O(1) candidate u per b and O(1+sqrt(N)/(1+x)) candidate words in total. Check the original integer equality to return literal addresses, not merely residues. Apply the same procedure with x and y swapped when y is short. Summation gives O(sqrt(N) log N) candidate tests overall, with O(log N) elementary integer operations per test. This is a constructive search bound, not an existence theorem.

## Why the actual order-nine window is congruence-complete as a bank

Take x=y=1 in the matched triple, so the left singleton is marked 1 and right singleton 0. The three literal order-nine witnesses have N=999,1001,1003 and the same continuation state h=45,S=2(432+129)=1122,d=0. With the end guards the directed Hamiltonian-cycle count is h: deleting the apex from a directed cycle leaves a Hamilton path of the ordinal spine, and the left/right guards allow exactly the required endpoints. Appending any further common word ending 0 keeps this count at 45 and keeps the three Hamilton-path values consecutive.

The [fixed-cycle congruence theorem](../chanpak_cross_problem_theorems_20261114/FIXED_CYCLE_CONGRUENCES.md) applies to these exact guarded cores. Since gcd(45,1122)=3, each arm covers the odd-compatible classes congruent to its own initial N modulo gcd(3,M), for every finite modulus M. The three initial counts cover all residues modulo 3, so their union covers **every odd-compatible class** modulo every M, all at cycle count 45. This is a literal relation with order nine and a matched common tail, not a nine-count window and not cofinite absolute coverage.
