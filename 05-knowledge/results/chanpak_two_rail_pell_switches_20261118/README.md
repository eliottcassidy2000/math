# Two rails, the 3–4–5 rotation, and the spoke-to-coupling transition

In the [Chan–Pak conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf), 7 and 21 are known holes; that these are the *only* holes is conjectural. The smaller problem of filling the odd-prime rail and the 7-times-prime rail (allowing the known exceptions, with the finite seeds 49, 63, 343) is equivalent by ordinal multiplication. Include 1 as the origin-side label, not as a prime. The first four odd columns have count pairs (1,7), (3,21), (5,35), (7,49): only the left endpoint, only the left endpoint, both, only the right endpoint, respectively. For 35 and 49, independent replay exhibits literal witnesses. Zero in the diagram is a coordinate origin, not a tournament count.

## The two arc directions are two readings, not two drawn edges

Fix every contact except one, e. Let C be the weight of odd-cycle packings not using either choice of e, and let W+ and W− sum H of the residual tournaments over odd directed cycles that use the chosen e in its two possible directions. By the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569),
`h+=C+2W+` and `h−=C+2W−`. All other arcs coincide. In rotated coordinates,
`M=(h++h−)/2=C+W++W−`, `D=(h+−h−)/2=W+−W−`.
Reversing the polarity reflects D. The actual local x7 test is the line `4D=3M` with D>0, or `W+−7W−=3C`. This gives a planar diagram of **observables**, not an embedding of the underlying complete graph.

At points (p,7p), the rotated legs are (4p,3p) and the hypotenuse is 5p. The two legs differ by p; after dividing out p, this is the primitive adjacent-leg triple (3,4,5). At p=3 the missing 21 is a rail endpoint; at p=7 the missing 7 is a rail endpoint and the earlier hole 21 now occurs as the shorter rotated leg. Geometry alone cannot tell if either endpoint is realizable.

The first three actual induced counts 1,3,5 are the beginning of `1+2^b`: reverse the extreme arc of a transitive tournament with b+2 vertices, and form the chain by inserting a new middle vertex, preserving every old arc. Starting at b=1 gives 3,5,9,17,…; a two-vertex directed edge supplies the start 1 on these same extremes. The excess doubles. Precisely those first two jumps are +2; 3→9 is a literal x3 as well as the start of the next multiple-of-3 column. See the [fixed induced-gap cutoff](../chanpak_log_induced_gap_cutoff_20261118/README.md).

## All-order classification: one spoke can multiply only by a power of two plus one

Keep a transitive spine fixed and change one incident arc of the apex. With positions numbered from zero, write w_i=1 for an outgoing arc and w_i=0 for an incoming arc. Every directed odd cycle passes through the apex; those with endpoints i<j and w_i=1,w_j=0 use an even subset of the j−i−1 intermediate vertices. Thus
`H(w)=1+sum_(i<j, w_i=1,w_j=0) (2^(j-i-1)+[j=i+1])`.
This is the inherited one-apex word law.

**Theorem.** In this fixed-spine spoke-flip family, every strictly increasing integral ratio is `1+2^b` with b≥1. The lower count is 1 or `1+2^a` with a≥1, and all these cases occur for every b. In particular, 7 cannot be the multiplier, at any height. This is about a spoke flip in a common transitive spine, not all arc reversals.

**Proof.** Since converse reverses the spine and complements all bits, assume the rising flip changes bit i from 0 to 1. If an earlier bit is 1, each gain from a zero to the right is at most the corresponding weight from that earlier 1 to that same zero (the inversion weight is nondecreasing in distance). Hence the increase is at most H−1, and there is no increasing integer ratio.

Otherwise write the suffix after i as `0^b s`, where s is empty or starts with 1. Let `Z=sum_(j:s_j=0)2^j` in zero-based positions of s, and h=H(s). Per zero, the contribution to h−1 is at most its weight 2^j, with equality precisely when all earlier bits are 1. If s starts with 1, the first bit alone supplies at least Z/2. In all cases `(Z+1)/2<h<=Z+1` (for Z=0 the inequality is literal). If b=0, the increase is Z (or zero for empty s), less than 2h, too small for an odd integer ratio ≥3. If b≥1, summing the new inversions, including the adjacency correction, gives the increase `2^b(1+Z)`. It has valuation b and lies in `[2^b h, 2^(b+1)h)`. If division by odd h is exact, the only integer with that valuation in the resulting interval is 2^b. Thus h=Z+1, which means s has at most one zero. Take s empty or s=`1^a0` to realize all cases.

**Adjacent-leg check.** An odd ratio k gives primitive rotated legs ((k−1)/2,(k+1)/2). The physical spoke ratios put a power of two on the *smaller* leg; the ratio 7=2³−1 instead puts 4 on the *larger* leg. In fact (3,4,5) is the only positive adjacent-leg integer Pythagorean triple with either leg a power of two. For a=2^j≥8 and b=a±1, hypothetical hypotenuse c satisfies `v2((c−1)(c+1))=j+1`, so `c ≡ ±1 mod a`. But `a+1<c<2a−1`. The cases a=1,2,4 give only 3,4,5. This proves a disjointness of **mechanism classes**, not a nonexistence of counts or ties at all Pell slopes.

## Two switches cross the obstruction, but still leave a placement defect

All arcs in the ordered blocks `0, A, B, 4, 5, 6` go forward except that B beats 0. A and B are transitive blocks of sizes a,b≥1. Reverse 0→4 (switch i) and 4→6 (switch j), independently. Put u=2^a,v=2^b. The exact count square is

    H00 = (u−1)(v−1)+2,
    H10 = (2u−2)v−u+3,
    H01 = 3 H00,
    H11 = (6u−6)v−3u+11.

To see why these are all-order formulas, both blocks can be visited at most twice. A can only be entered from 0 after leaving it. After leaving B, returning to B requires reaching A again, which can only be done through 0 (a first block in a path needs no predecessor). There is no other way back: the later triple has no arcs into A or B, even under either switch. Every quotient path word therefore contributes 1 or `2^a−2` for the A slots and 1 or `2^b−2` for the B slots, independently; two ordered nonempty subsets of a transitive block admit exactly one path order each. The displayed bilinear formulas now follow from the four (a,b) corners in {1,2}², computed by directed subset DP (see replay). The factor 3 on the 01 side also follows directly by ordinal join with the cyclic last triple.

The mixed correction `H11−H10−H01+H00=2(u−1)v+2` is positive: it is exactly what independent additive repairs would lose. At (a,b)=(2,1) the four counts are **5,11,15,35**. No one-arc ×7 pair exists through order 7 (labeled replay); the pair 5→35 at order 7 uses two coordinated reversals. The formulas give the exact general defect:

    H11−7 H00=(4−v)(u−1)−6.

Thus within this family equality occurs only at (2,1); at b=2 the defect is **−6 at every a**. These are actual counts `p_a=3·2^a−1` and `7p_a−6`, irrespective of whether p_a is prime. This supplies a bounded residual on an unbounded ray without freezing the small endpoint. It does **not** supply the +6 correction: at large height the induced-gap theorem rules out doing that by merely adding vertices while preserving all old arcs. The open construction is a third rewiring (or a different core) that adds six without changing the leading coefficient, followed by a way to choose the ray from the absolute prime target.

## Replay

```sh
python replay.py | diff -u replay.out -
g++ -O3 -std=c++17 -Wall -Wextra -Werror probe.cpp -o /tmp/rail-probe
/tmp/rail-probe | diff -u probe.out -
```

Integer tournament codes use lexicographically ordered unordered pairs (i,j), i<j; a set bit reverses i→j.
