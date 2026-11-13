# Sparse target spectrum of the two-apex bivariate language

Here "one-apex" means a transitive tournament TT_n on its natural linear order plus one exceptional vertex x, allowed an arbitrary binary out-cut to that order. This includes all orientations, all orders, and arbitrary ordinal padding. This is an obstruction to that particular language, not to tournaments generally.

The inherited Session321 selective contexts have path count
```
T(U,V) = A(U) Q_V(2) − H(U) H(V),
A_13(U) = −3H(U)+6Q_U(2)+2Q_U(3),
A_15(U) = −3H(U)+6Q_U(2)+6Q_U(3).
```
These coefficients follow by substituting F2=Q(2)−2H and F3=Q(3)−3Q(2)+3H into the two exact itinerary matrices. Here Q(t) counts covers by t optional ordered slots.

## Exact correlation formula for a one-apex word

Write a_i=1 if x→i. Let I_d count pairs i<j with a_i=1,a_j=0 and j−i−1=d. Every directed cycle uses x. Such a cycle with extremes i,j uses an arbitrary increasing subset of the d intervening vertices; the cycle is odd exactly when the chosen subset has even size. No two such cycles can appear in one packing. Thus the odd-cycle packing formula yields
```
Q(z) = z^(n+1) + sum_d I_d z^(n−d−1) [(z+1)^d+(z−1)^d],
H    = 1 + 2 I_0 + sum_(d≥1) I_d 2^d.
```
The second formula is evaluation at z=1; the d=0 convention is the polynomial value 1+1.

In particular, Q(2) modulo 2^k (n+1≥k) depends only on farthest inversion distances d≥n−k+1: the bracket 3^d+1 is even, so all earlier terms vanish. If k≥2 this uses only the first and last k−2 bits, with the empty-width case interpreted literally. In contrast H modulo 2^k sees only short distances d<k, but throughout the entire interior. This is an exact boundary/interior split, not an equality of their response functions.

For an explicit separation, start with a long word having a single initial 1 and all subsequent bits zero. Change two consecutive interior zeros to 1's, far from the first and last k−2 positions and at least two positions away from the initial 1 and terminal bit. This adds one adjacent and two distance-two inversions, hence flips H modulo 4, while leaving Q(2) equal modulo 2^k. The output T(U,V) with V these two words differs by 2 modulo 4 for every common U and k≥2, since H(U) is odd. No truncation of the outer Q(2)-boundary, even at arbitrary prescribed precision, subsumes the interior endpoint response.

## Zero density in both variables, not just fixed rays

If the word has no inversion, the whole tournament is transitive. Otherwise its active span is from the first 1 to the last 0, of length r≥2. The leading zeros form a transitive ordinal prefix to the active apex core, and trailing ones form a transitive suffix. There are at most 2^(r−2) active words. For H and Q, the allocation of the outside padding between the two ends is immaterial: if its total length is ell then H is unchanged and Q(t) is multiplied by t^ell. This follows either by slot concatenation or absence of directed cycles across ordinal joins.

The one extreme inversion of every active word forces
```
Q_core(t) ≥ t [(t+1)^(r−2)+(t−1)^(r−2)]
```
for real t≥1; it is exactly that inversion's contribution to the packing expansion. Therefore Q_V(2)≥(2/9)·3^r·2^ell and Q_U(3)≥(3/16)·4^r·3^ell. For each positive Y, the number of *candidate response profiles* of V with Q_V(2)≤Y is O(Y^alpha), alpha=log_3(2): summing first over ell gives an upper bound by a constant times
```
sum_(ell≥0) (Y/2^ell)^alpha.
```
Indeed, for each ell, summing 2^(r−2) over r≤log_3(CY/2^ell) is bounded by that power. The wholly transitive case adds O(log Y).

For every nonempty U,V, Q(2)≥2H from the two choices of a single nonempty slot. Thus
```
T(U,V) ≥ [A(U)−H(U)/2] Q_V(2) ≥ 2 Q_U(3) Q_V(2).
```
Counting V for each U gives O(X^alpha) times the sum of Q_U(3)^(-alpha) over candidate U profiles. The transitive series converges, and the one-apex part is dominated by
```
const · sum_(r≥2) (2/4^alpha)^r · sum_(ell≥0) 3^(−alpha ell),
```
which converges since 4>3. This proves a uniform O(X^alpha) bound on **distinct absolute targets at most X** realized by either context when both source and reservoir vary freely over all one-apex words. It also holds after taking a finite union of such contexts. Since alpha<1, these physical families cannot cover a positive-density odd tail; they cannot cover all primes or all 7p targets either (even Chebyshev's lower bound gives more than O(X^alpha)). This does not require enumerating numerical target coverage and does not constrain more general multi-apex or heterogeneously changing sources.

Reproduce the word formula against literal subset-convolution ordered-slot counts with replay_one_apex.py; the asymptotic proof is independent of the finite replay.

## Stronger: arbitrarily many serially composed apex cores still miss almost every target

Allow both U and V to be arbitrary nonempty ordinal sums of arbitrary length, with letters either a singleton or an active one-apex core. (Every active core is strong: the apex reaches the first vertex, the first reaches the last along the linear order, and the last returns to the apex. Every intermediate vertex is on that directed cycle.) This includes a clock with unrestricted repetitions and all internal binary word choices at every letter.

Let L(s) sum Q_C(2)^(-s) over the singleton and all candidate active words C, counting duplicate isomorphism types if necessary. For an active word of r bits, set w(d)=2^(r−d−1)(3^d+1). Each interior position j=2,...,r−1 forces one additional inversion: from the first vertex if its bit is 0, or to the last if its bit is 1. These are all distinct pairs, and none is the extreme pair. Therefore
```
Q_C(2) ≥ B_r = 2^(r+1) + w(r−2) + sum_j min(w(j−2),w(r−j−1)).
```
For r=2,...,6 the five B_r are 12,32,84,216,548. The coarse extreme-pair lower bound Q_C(2)≥2·3^(r−2) handles the tail. Consequently
```
L(1) ≤ 1/2 + sum_(r=2..6) 2^(r−2)/B_r + (3/2)(2/3)^5
     = 1189691/1242864 < 24/25.
```
An explicit exponent below 1 follows from Hölder: L(3/4) < 7, because 2^(-3/4)<3/5 and 2/3^(3/4)<9/10, so the geometric estimate is less than 3/5+6. At s=199/200, interpolation gives
```
L(s) ≤ L(1)^(49/50) L(3/4)^(1/50)
     < (24/25)^(49/50)·7^(1/50) < 1,
```
the final inequality being the integer test 7·24^49 < 25^49.

Now Q(2) multiplies along ordinal sums. Their complete Dirichlet weight sum is bounded by the geometric series 1/(1−L(s)), including the empty sequence. Since Q(3)≥Q(2), the same is true in the U role. Counting pairs by T(U,V)≥2Q_U(3)Q_V(2) shows that the number of distinct targets ≤X is O(X^(199/200)) for the full two-module serial language. This is weaker numerically than the single-core exponent, but applies to arbitrary numbers of changing components. The proof also allows adding any further block alphabet whose sum of Q(2)^(-s) remains less than the available margin, with the same fixed bivariate context. In particular, the missing ingredient cannot be supplied merely by unrestricted serial reuse of the one-apex words in these two inputs.
