# Repeatable common selective returns after full-cover dominance

Start with the two equal-H return pairs from
[Session 320](../session320_second_moments_and_native_contexts_20261113/README.md).
The [all-visit test](FOUR_VISIT.md) rules out common strong homogeneous
substitution of the whole returned tournament. Here we instead attach
a new module using a named cut of the old tournament. This yields a
repeatable, exact strong return with arbitrarily many added vertices.

## A bivariate itinerary calculation

Let U be the original nine- or ten-vertex source, before the Session
320 homogeneous return. Compress it to vertex 0 of its old quotient.
Add one new marked quotient vertex V. In the order-13 return use outcut
7167 on the old full labels (all but label 10); in the prime order-15
return use outcut 6144 (old labels 11,12). Both cuts are uniform on U.
After compression the refined outgoing-mask quotients are
```
R_13 = (24,9,35,4,14,27)
R_15 = (108,97,2,6,79,92,12).
```
Here U fills 0 and V fills the last vertex. Let F_i(U) count ordered
nonempty i-path covers and let H(U)=F_1(U). Every Hamilton path compresses
to a legal quotient word with each exterior singleton exactly once and
i,j runs in the two modules. Fill runs by independent ordered covers.
The exact quotient tables (rows i=1,2,3, columns j=1,2) are
```
R_13: [29,15] [24,12] [4,2]
R_15: [53,27] [48,24] [12,6].
```
There are no other entries. The replay counts these marked words by
a finite subset/counter recursion; the count is independent of the
sizes and internal structures of U,V. In particular the prime cut has
two out-neighbors yet no feasible three-run layout.

Write L(U) for the first-column linear form in F_1,F_2,F_3. Twice the
second column equals the first plus (1,0,0). Hence for either context,
and for **any nonempty tournaments U,V**,
```
H(R[U,V]) = ((L(U)+H(U))/2) Q_V(2) − H(U) H(V).
```
Thus the two named responses H(U),L(U) are sufficient for the entire
new reservoir, not just a fixed prefix size. The cross-interface term
is retained in the formula. Both quotients are strong, so all these
inflations are strong, even if V is not.

## Consecutive exact returns at unbounded sizes

For the old order-13 pair, L(U)=1068387 and H(U)=1047 on both sides.
For the old prime order-15 pair, L(U)=9189287 and H(U)=1915 on both
sides. Taking V=TT_m for m>=1 gives respectively
```
534717 · 2^m − 1047
4595601 · 2^m − 1915.
```
Each additional transitive vertex beats the same named old set and
loses to every earlier new vertex. Thus equality holds after every
addition along one nested chain, not merely at isolated sizes.
The second ray is 2 modulo 3 at every m; its intercept and slope
are coprime. No claim of arbitrary prime or odd coverage follows.

There is a separate route starting from the prime pair: first add
the unique out-neighbor cut {11} (mask2048), which gives 7655505
on both sides but **does not** persist if cloned immediately.
After that common step, the cut 38912 gives matched one- and two-run
counts (15316755,7659335), with no three- or four-run contributions.
Cloning this second cut gives the exact ray 7659335·2^m−1915.
The actual old-mask DP verifies each of these layout claims.

The new ray pairs remain distinct. With two optional path slots,
concatenate two independent layout words on complementary old sets.
Their run counts add, so the coefficient of a j-run filling is a
subset convolution. For the first two rays the differences of these
two-slot layout coefficients (old X minus old Y, starting at j=0) are
```
(-134400, -397920, -394080, -130560)
(-2053440, -6143040, -6134400, -2044800).
```
All higher entries vanish. For m>=1 at least the j=1 contribution
is strictly negative, so Q(2) differs at every size. This also proves
nonisomorphism up to converse. What repeats is a **named selective
interface**, not equality of the full slot polynomial.

## Exhaustive one-apex boundary on the original pairs

For the order-13 sources there is one common nontrivial cut equalizing
H (7167). For the prime sources there are three (2048,5120,6144),
giving 7655505,13775313,9189287. Among separately chosen nontrivial
cuts, there are respectively 51 and 25 ordered cut pairs with equal H,
but none with both equal H and equal Q(2). No such one-apex strong
addition keeps either original absolute H unchanged. These statements
are finite, for the specific labelled source pairs. The selective rays
above are proved for all lengths by the bivariate layout theorem.
