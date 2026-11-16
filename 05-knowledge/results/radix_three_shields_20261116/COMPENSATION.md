# A width-four window with all its +2 moves beside a literal rooted triangle

Take a rooted tournament A+x, a transitive B after A (every arc A→B),
and choose the contact pattern between x and B. Write h=H(A),
f=H(A+x). Let R and S count even and odd rooted **start** paths
from x, weighted by the Hamilton count of the complement. Then
R counts odd-size A runs; S−h counts nonempty even-size A runs.

Four literal choices of A, with out-neighbor mask of x (root given
last), give the following profiles:

| outgoing rows of A | mask | f | h | R | S |
|--|--:|--:|--:|--:|--:|
| (42,12,1,4,15,30) |21|75|15|51|61|
| (34,16,3,7,13,30) |6|77|33|47|65|
| (34,12,1,5,15,30) |20|79|19|51|61|
| (48,13,1,5,14,30) |10|81|23|49|63|

All four obey R+S=112. Let B be a transitive tape of k≥2
vertices, with x beating every site but the last. The standalone
B+x shield has count 1+2^(k−1), hence a literal triangle when k=2.
The full count is

    H = f + (R+S)·2^(k−1).

Indeed for a rooted odd cycle wholly in B+x, choose an odd-size
subset of the first k−1 sites, then the last; this gives 2^(k−2)
cycles, weighted 2h. For a rooted cycle meeting both A and B,
its B run ends at the last site. Of the preceding sites, exactly
2^(k−2) subsets have the appropriate parity for any given A run.
Sum with the odd-cycle pack weight 2 and the weighted A response,
then include f for packs not using a root cycle in B. No new
cycle can use different A/B blocks without passing through x.
The four windows are unchanged at stride 2, e.g.
299,301,303,305 for k=2, 971,973,975,977 for k=4, and
1867,1869,1871,1873 for k=5. Both the shield arcs and every
incidence with it are identical across settings. The joined
tournaments are strong, since A and x communicate, A beats the
shield, and every shield site can reach x.

## All fixed shields preserving exactly these four low settings

Allow arbitrary nonempty B, still with every arc A→B. Set t=H(B).
Let C count weighted even paths in B starting at sites beaten by x
and ending at sites beating x. Let E and O count weighted even
and odd paths with just the latter endpoint condition. In all
cases the weight is H of the complement in B. The same cyclic
partition gives the exact output

    H = t f + 2(C−E) h + 2(O−E) R + 224E.

If this is to equal f plus a common translation on all four
seeds, the three difference vectors in (f,h,R) are
(2,18,−4), (4,4,0), (6,8,−2); they have nonzero determinant 96.
Thus necessarily t=1 and C=E=O. In particular B must be transitive.

In transitive order a vertex beating x at rank j contributes
equally to E,O if j≥1 (zero-based ranks), but rank zero gives
one odd path and no even one. Thus E=O means no incoming arc
to x at rank zero. For C=E, every even path ending in an
in-neighbor of x must start in an out-neighbor. Two distinct
in-neighbors already give a forbidden two-site path; conversely
zero or one in-neighbor at positive rank works. This classifies
every fixed shield in this ordinal-interface architecture for
these four specific seeds. It offers exact low-window
translations 0 and 112·2^j (j≥1), not a free adjacent +6 dial.
The independent two-dial bank uses a genuinely different
boundary response.

Source: the four seeds and constant-parity shield calculation
were identified in a companion research pass; the rank and
uniqueness arguments above are all-order. Odd-cycle packing:
[Grinberg and Stanley](https://arxiv.org/abs/2307.05569).
