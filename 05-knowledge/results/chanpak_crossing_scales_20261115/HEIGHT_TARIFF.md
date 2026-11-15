# Sharp square-root capacity of separated high/low placement

The previous equal-Q reservoir and this session's new +2 crossings are physically exact at unbounded height. Nevertheless this note proves that the entire **carry-free high/low regime** in which they live is thin, even with a changing binary host, both core ratios independent, and infinitely many core types in the one-apex ordinal alphabet. The upper bound is O(sqrt X) through X and is sharp in exponent and scale. This does not cover binary glues allowed to carry their low part through the high place or cores outside the stated alphabet.

Use the exact [gap-separated formula](CHROMATIC_INTERFACE.md)
```
N=2^g D+c,
D=d·2^k+t(2^k−rev_k(b)),
c=G+d b+t(A(v)−1),
d=(Q_U−h_U)Q_W, t=h_U h_W, G=(Q_U−h_U)(Q_W−h_W)+2t.
```
Here v is the k-bit suffix of value b (0≤b<2^k), g≥1, and the host word is 1 0^g v. The no-carry condition is exactly 0<c<2^g, not an asymptotic approximation. The case k=0 and empty v uses b=rev=0 and A(v)=1.

The basic floor is G≥d/2, by Q_W≥2h_W. Hence D≥d·2^k and c≥(d/2)(b+1), for **all** nonempty cores. If N≤X and c<2^g, this forces
```
2^k ≤ X/(d·2^g),    b+1 < 2^(g+1)/d.
```
It also forces d<sqrt(2X), since D≥d and c≥d/2. For a fixed profile pair, put a=X/(d2^g), b_*=2^(g+1)/d. At a fixed g, the number of possible (k,v) is bounded by the sum of min(2^k,b_*) over 2^k≤a. When b_*≤1 it is zero. Otherwise it is at most b_*[3+log_2^+(a/b_*)] by splitting the geometric sum at b_*, and also at most 2a. For 2^g≤sqrt X use the first bound, summing O((2^g/d)(1+log_2 X−2g)) to obtain O(sqrt X/d). For 2^g>sqrt X use 2a and sum the geometric tail, again O(sqrt X/d). All constants are absolute. Distinct N are no more numerous than the presentations counted.

It follows that for **any** target-independent choice of core profile alphabets, the number of exact outputs through X in this grammar is at most
```
C sqrt X · sum over distinct profile pairs with d<sqrt(2X) of 1/d.
```
So any positive-density theorem based solely on these separated words needs that truncated reciprocal weight to grow as at least a constant times sqrt X. Merely making the high-quotient union infinite is not enough.

In particular, when both cores are arbitrary ordinal sums of one-apex blocks, the reciprocal sum is bounded. The [independent-core capacity estimate](../chanpak_heisenberg_target_20261115/HOST_CAPACITY.md) shows exp convergence of the multiset Euler product at exponent 1: there are at most 2^(n−3) strong cyclic one-apex block types at order n, and each contributes Q≥2·3^(n−3). Singleton copies contribute a convergent geometric factor. Since d≥Q_UQ_W/2, the required sum is at most twice the square of the convergent sum of 1/Q over distinct profiles. This is another demonstration of the 2-versus-3 tension: binary choice creates profiles, but the two-cover cost of a cyclic core suppresses its contribution.

For a matching lower bound, take U=W=TT1. Then d=2,t=1,G=3, and each of the 2^k suffixes gives a distinct D with 2·2^k<D≤3·2^k. Also c=2(b+1)+A(v)≤3·2^k. Choose g=k+2, so there is no carry and the N are distinct. They all lie below 15·4^k. For any large X choose the largest k with 15·4^k≤X; these 2^k distinct outputs give the matching lower bound c_0 sqrt X.

## Consequence for the next proof attempt

The chromatic list idea identifies correct interfaces, and the positive boundary calculations quantify which profiles are possible. But a prefix–suffix proof for all targets in this family cannot insist on separating the scales before choosing the source. It must either use actual carry with a matched checksum, or buy substantially higher weighted capacity from richer cores. Exact +2 examples and even a growing number of choices at one height do not pay that cost. This is an all-order counting obstruction, not inferred from any observed density.

Replay: run `python replay_height_tariff.py` for exact strong-word signatures, the lower construction, and validation of the per-profile counting inequalities at finite bounds.
