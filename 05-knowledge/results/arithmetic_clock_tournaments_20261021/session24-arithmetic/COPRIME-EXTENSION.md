# Axis reversal is prime for coprime odd cyclic factors

Status: elementary proof, extending the maximal-module argument. Independent finite controls are supplementary, not the proof. No literature-priority claim.

Let A=Z/m and B=Z/n, with m,n>1 odd and coprime. Let alpha:A\{0}->+/-1 and beta:B\{0}->+/-1 be any skew Cayley tournament signs, so alpha(-a)=-alpha(a), beta(-b)=-beta(b). Define the axis-reversed product on A x B by

- s(0,b)=beta(b);
- s(a,0)=-alpha(a);
- s(a,b)=alpha(a) when a,b are both nonzero.

The claim is that the resulting tournament is prime, even if either factor tournament is decomposable.

## Proof

The carrier A x B is cyclic by coprimality. Every cyclic Cayley tournament is strong, since one of the two generator differences1,-1 is present and gives a directed Hamiltonian cycle.

In a strong tournament, maximal proper modules are disjoint. Intersecting modules have a module union. Two overlapping maximal proper modules would therefore cover all vertices, and their nonempty regions X=M\N, I=M intersect N, Y=N\M would be oriented X->I->Y and X->Y, or all reversed, contradicting strength. The maximal proper modules consequently form a canonical partition. Vertex transitivity makes their sizes equal. If the tournament is decomposable, their size is nontrivial, and translation invariance forces the block H through0 to be a nontrivial proper additive subgroup: for every h in H, the block H+h meets H at h and hence equals H.

Every subgroup H of A x B splits as H_A x H_B. To see this directly, choose an integer k with k=1 modulo m and k=0 modulo n. Multiplying (a,b) in H by k gives (a,0) in H, and subtraction gives (0,b) in H. Therefore all combinations of the coordinate projections occur in H.

No nontrivial proper rectangular subgroup is a module:

1. If H_A is nonzero and H_B is proper, choose nonzero a in H_A and b outside H_B. The outsider z=(0,b) sees (a,0) and (-a,0) with the opposite mixed-coordinate signs alpha(a), alpha(-a).
2. If H_B=B and H_A is proper, choose a outside H_A and nonzero b in B. The outsider z=(a,0) sees (0,0) and (0,b) with opposite signs -alpha(-a), alpha(-a).
3. The remaining nontrivial possibility is H_A={0} and H_B nonzero proper. Choose any nonzero a in A and nonzero b in H_B; the same pair as in case2 distinguishes H.

This contradicts the subgroup-module forced by decomposability. Hence the product is prime.

## Sharp scope warning and exact controls

Coprimality is sufficient; it is not asserted necessary for every factor pair. It cannot be omitted uniformly: for m=n=3 with alpha(1)=beta(1)=+1, the diagonal {(0,0),(1,1),(2,2)} is a proper module. This carrier is C3 x C3, not cyclic C9, so no false CRT identification is made.

`coprime_extension.py` independently computes every anchored pair's minimum module closure:

- orders9 and5: all16*4=64 skew-factor pairs give prime tournaments on45 vertices;
- orders3 and5: all2*4=8 pairs give prime tournaments on15 vertices;
- orders3 and3: all2*2=4 pairs are decomposable, with size3 proper closures; the stated diagonal is checked directly in the positive-sign case.

These finite cases are controls on the proof, not a substitute for it. Normal and optimized Python runs give identical JSON receipts.
