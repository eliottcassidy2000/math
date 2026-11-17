# Both signs of the sevenfold dyadic ray, and a five-arc route

**Proved:** every `7(2^k+1)` with integer k≥2 is a tournament Hamiltonian-path count (at k=1 the value is 21). The [minus ray](MERSENNE.md) occurs for every k≥3 (initial holes 7 and 21). More strongly, one fixed physical operation on an arbitrary nonempty module turns the count `128P₂+H` into exactly seven times that count, by five monotone single-arc reversals. This does not give universal sevenfold prime selection.

## An arbitrary-module plus reader

Let D have outgoing masks (270,140,264,272,167,463,31,333,82), marked q=8. Add B with all D−q beating B and B beating q. Write P=P₂(B) for ordered two-slot covers allowing empty slots and h=H(B). For the pointed D, let G=H(D), J=H(D−q), e be the full count ending at q, L=G−e, and W the number of ordered covers (X,Y) of D where Y starts at q and has at least two vertices, with X optional. The exact formula is

    H(D with B) = W P + (L+J−W)h.

Indeed B can be visited in at most two runs because every nonterminal B run must exit to q. One B-run at the end contributes L skeletons; one earlier run gives a split X|Y, with Y starting at q, contributing W+J (the cases Y={q} account for J). Two B-runs have the second terminal and force Y to extend beyond q, contributing W skeletons, each with P−2h choices. The replay evaluates the rooted profiles by independent subset path counting:

    (G,J,e,L,W) = (877,127,101,776,896).

Hence the output is `896P+7h=7(128P+h)`. For TT_m, P=2^m and h=1, giving all k=m+7≥8. For k=2,…,7 the respective counts 35,63,119,231,455,903 have literal witnesses in the replay. The main family is strong, as follows from its shared Hamiltonian cycle below. The order-ten case m=1 gives 7·257, improving the previous order-eleven construction without a minimality claim.

## The same block gives the unscaled count, five arc reversals away

Use the nine exterior vertices in order 4,5,7,6,0,1,2,3,8, with B immediately before 8. Start with the transitive order of the exterior, except reverse the first-to-last edge to 8→4. Excluding the exceptional edge, a Hamiltonian path has the unique exterior order and contributes h. Including it, the interior vertices (seven transitive exterior vertices followed by B) split between the two sides of 8→4, giving `2^7 P`. So this is a strong witness for the **unscaled** value `128P+h` with exactly the same vertices and B contacts as the final witness.

Now perform these five literal reversals, in order; none touches B:

| exterior pair reversed | resulting count |
|---|---|
| none | 128P+h |
| {1,7} | 256P+5h |
| {1,8} | 276P−9h |
| {3,4} | 475P+109h |
| {6,8} | 823P+13h |
| {4,6} | 896P+7h |

These affine reads use the same two-run formula and are checked from the exterior at each step. Each increase is strictly positive for every B, since P≥2h. All six tournaments retain the spanning directed cycle formed by any Hamiltonian path in B followed by 8→4→5→7→6→0→1→2→3 and the arc back to its start. Thus the sevenfold scaling is a five-arc path through strong tournaments, not an ordinal product with the forbidden count 7.

## Symmetric gap, and why it cannot be a prime twin construction

Feed the plus reader B, and feed the [old negative seam](MERSENNE_MODULES.md) B followed by one ordinal singleton. This leaves h unchanged and doubles P. The two strong counts are exactly

    F₋=7(128P−h),     F₊=7(128P+h),     F₊−F₋=14h.

Given the two counts, the physical block would have to realize `P=(F₊+F₋)/1792` and `h=(F₊−F₋)/14`. A transitive B gives exact gap 14 at every dyadic quotient center 2^k for k≥8, and the two rays together give the corresponding numerical gap also for 3≤k<8.

But P≡(−1)^|B|h (mod 3), by the [path-slot identity](OPPOSITE_SIGNS.md#why-the-paired-prime-examples-need-an-odd-block-at-odd-reader-levels). Since 128≡−1, even |B| makes the plus quotient divisible by 3 and odd |B| makes the minus quotient divisible by 3. This excludes simultaneous prime quotients greater than 3 in the matched construction. A blocker in a repair family is not a blocker in the whole spectrum.

There is an additional rigidity in the old negative scaffold. Its full pointed triple is `(64P+h,96P−2h,10P)`: the starting count also has at most two B-runs, and its readings at TT₁ and TT₂ are 20 and 40. With arbitrary right triple (d,e,f), the seam gives `(64d+96e−10f)P+(d−2e)h`. To retain the exact normalized coefficient −1 on h and an integral quotient slope, one needs 2e=d+7 and f≡0 (mod 7). Since e≥d, d is odd and 7 is forbidden, d≤5; then f≤d forces f=0. If the right port has a terminal Hamiltonian frame, f=0 places it in a last strong component with a nonempty ordinal prefix (a strong component admits every start by rotating a Hamiltonian cycle). The prefix multiplies the terminal two-cover count by its two-slot count, at least twice its Hamiltonian count. Thus e≥2d, so e=(d+7)/2 forces d≤7/3, and oddness gives d=1. This is exactly e=4,f=0 and slope 64. Hence changing only that framed right reader cannot keep coefficient −1 while introducing a slope divisible by 3; transitive padding only multiplies 64 by powers of two. A genuinely different contact geometry is needed for such a symmetric prime pair.

Replay: `node replay_plus_ray.js`. The computations certify the small witnesses and the affine coefficients, not the all-order counting argument.
