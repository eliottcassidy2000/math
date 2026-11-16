# A growing alphabet of whole-tape determinants with a literal +2

There are Ω(X^((log_3 2)/4)) distinct positive whole-tape determinants from physical +2 repairs with both counts ≤X. This exceeds the necessary determinant-menu exponent 1−(3/2)log_3 2 from the earlier bounded boundary obstruction. We do not assert that the determinant choices yield distinct Hamilton counts or select arbitrary absolute counts.

## Exact gap

Put letters 2a+b on a transitive spine, with a,b recording that apex A,B beats the site. Let Δ=H_AB−H_BA. Let L=(L_a,L_b) be the forward binary addresses, R=(R_a,R_b) the reverse-complement addresses, and D(L,R)=L_b R_a−L_a R_b. The [arc contraction identity](../unit_shields_wheel_followup_20261116/apex/SHIELD.md) is
```
Δ(x y) = Δ(x)+Δ(y)+D(L(x),R(y))+λ(last(x),first(y)),
λ(x,y)=a_x−b_x+b_x a_y−a_x b_y.
```
In particular λ=0 when both letters are diagonal (0 or 3).

Take central patch p=0001120. Direct contraction gives L(p)=(2,12), R(p)=(95,103), Δ(p)=8 (weighted term 6 and adjacent term 2). Let u,v be diagonal words of equal width m, with common forward address T for u and common reverse-complement address S for v. Their direct cross term is zero. The internal patch has diagonal boundaries, so
```
Δ(u p v) = 8−8T+10S.
```
Choose T=5t+2 and S=4t+1 to make Δ=+2 at every allowable t. This is one change of the mutual apex arc with neither exterior edited.

## An exact missing-digit packing of whole determinants

For j≥1 put m=4j+7. Set
```
t = 25·16^j + Σ_{i=1}^{j−1} ε_i 16^i,      ε_i ∈ {0,1,2,3}.
```
The leading seven-bit chunks of T and S are 125 and 100, both starting in 11. The last nibbles are 2 and 1; each internal nibble is 5ε_i for T and 4ε_i for S, without carries. The left tape u therefore starts 33 and the right tape v ends 00.

Write r_m for reverse complement of an address in width m. For w=u p v,
```
L_a(w)=2^(m+7)T+2^m·2+r_m(S);    L_b(w)−L_a(w)=10·2^m;
R_a(w)=2^(m+7)S+2^m·95+r_m(T);   R_b(w)−R_a(w)=8·2^m.
```
Thus 10S−8T=−6 gives
```
D(L(w),R(w)) = 2^m [166·2^m+10 r_m(T)−8 r_m(S)] >158·4^m.
```
In particular none vanishes. The internal reversed nibbles contribute −10 rev_4(5ε)+8 rev_4(4ε), respectively 0,−84,−42,−126. The fixed high and low chunks then give
```
D(L(w),R(w)) = (339/2)·4^m −800·2^m −42·2^(m+7)·Z,
```
where Z runs through base-16 words with j−1 digits in {0,1,2,3}. Each digit is independently selected by the permutation 0,2,1,3 of ε and the digit order is reversed. This proves separated, not contiguous, determinant control.

Both tournaments are strong: both apices beat the first spine vertex, and the last beats both. For either apex order a Hamilton path assigns the n=2m+7 spine sites to three increasing runs, so H≤2·3^n. Fix the first two vertices 33 as the starts of runs 1/2, and the last two vertices 00 as the ends of runs 0/1. Each interior vertex independently chooses any run, giving H≥3^(n−4) in both orientations. For X large take the maximum j with 2·3^(8j+21)≤X; the 4^(j−1) determinants at this j have both counts below X.

## Arbitrary-gap variant: growing *cut* determinants

For any even d, the [programmable-gap construction](../unit_shields_wheel_followup_20261116/apex/PROGRAM.md) supplies P=s c rev(s) with diagonal s of width k, c=12 or 0123, and Δ(P)=d. Let e=L_b(P)−L_a(P), equal to 2^k or 2^(k+1), respectively. Wrap by a diagonal u starting 33 of arbitrary width m and its reverse complement rc(u). Its common forward address T runs through [3·2^(m−2), 2^m−1]. At the cut u|P rc(u), the cross determinant is exactly −eT, because reverse-and-swap symmetry gives R_a(P)−R_b(P)=−e and both coordinates of R(rc(u)) agree. The other cross term cancels it exactly, keeping Δ=d. This proves a complete arithmetic interval of nonzero cut determinants. A cut determinant is not the combined tape's own determinant; the explicit +2 family above addresses that intrinsic quantity.

## A Vitali analogue and its limitation

For a finite family of intervals, select a longest remaining one, discard those intersecting it, and repeat. The chosen intervals are disjoint and their concentric triple dilations cover the union: a discarded interval is no longer than the chosen one it intersects. The sum of disjoint lengths is at least one third of the union length. This controls a measure, not individual integer membership.

A bracket with two physical endpoints need not realize any values between them. E.g. values 1 mod 4 supply endpoints for all consecutive length-4 brackets, covering the odd line as a union of continuous intervals, yet miss every value 3 mod 4. Our +2 pairs certify both odd grid points of each chosen interval, but determinant packing alone gives no count placement at a given admissible N.
