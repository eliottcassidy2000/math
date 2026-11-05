# Independent review: prime-cube rigidity

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197-orbit-wildcard/`:

- PRIME_CUBES.md: fdc49f17f3f568f18171e03aa0e1ebdf2197666e849c7b569986a11d33746cb9
- verify_prime_cubes.py: 7f2468fd19179e0e3678fc65ede2698987285cef306db307633580a53e02116c
- Normal/optimized output: d3b90153a122dc6fd642112bdf98628aa4658b7f879fc6bc3ef7566d05f26b7a

Both source replays reproduce the stated output hash. The import dependency is the already-reviewed wildcard verify.py, hash 9f354c23dfe8616f73f58cc15f56b6054a5f197bdb51fc944c6550b04d6ceba9.

## All-order proof

The free action of Aut(T) on Hamilton vertex words implies its order divides p^3. Cauchy's theorem then supplies a prime-order element whenever the group is nontrivial. A nontrivial vertex orbit of that element is a p-circulant, and is strong: any nonzero allowed step generates the additive group of prime order.

The Busch estimate is applied to this strong induced orbit, not the arbitrary ambient graph. The previously independently opened official Busch Theorem1 gives the stated weaker H>=5^((n-1)/3) bound. The integer base inequality and the monotonic comparison 5p^9>(p+1)^9 propagate the cutoff to every p>=17. This excludes any moving orbit there.

For p=11,13 the exhaustive connection-set minima are far above p^3. For p=5,7 the disjoint-path merging injection excludes two moving orbits. With exactly one moving orbit, every remaining vertex is fixed and sees it uniformly, so it is a module. If it is its entire SCC, the SCC Hamilton-product divisibility is impossible. Otherwise contraction inside its strong component produces a strong quotient. The contracted vertex lies in a directed triangle, which lifts to the claimed cyclic extension. Its count 3H+F2 gives 135>125 at p=5 and at least2163>343 at p=7. This covers unbounded singleton padding and arbitrary weak ambient tournaments without a global order cutoff.

The p=3 boundary is correctly excluded from the theorem: three ordered triangles have H=27 and a nontrivial triangle rotation. The labelled and path-pointed mass consequences follow from rigidity and make no existence or class-count claim.

## Independent finite evidence

The new audit program `prime-cubes-audit/walk_inclusion.cpp` does not use either source Hamilton subset recurrence. It counts arbitrary directed walks of length p-1 inside every vertex subset S, then applies inclusion-exclusion:

    H = sum_S (-1)^(p-|S|) * walks_(p-1)(T[S]).

A length p-1 walk that visits all p vertices visits each exactly once, so the surviving walks are precisely Hamilton vertex words. The computation uses signed 64-bit integers; at p=13 the maximum intermediate walk magnitudes and the complete subset sum remain well below the integer limit.

All108 connection-set presentations reproduce exactly the four source histograms, including minima92411 at p=11 and3669497 at p=13. For p=5,7 a separate literal path DFS on every induced subset independently gives F2=90 and1638. The result is recorded in `prime-cubes-audit/independent.json`.

Thus the large-prime argument is all-order, the exceptional orbit banks are complete finite certificates, and the module/ear reduction is the bridge between them. No conjectural circulant minimum or finite-order ambient census is used.
