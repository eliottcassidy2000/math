# Moving projections: three repairs on a fixed core, six on a changing core

**Proved constructions at every indicated guard height.** No extension of the Chan–Pak coverage range is claimed. A *route* here is a sequence of exact strong tournament witnesses, with a bounded batch of literal arc reversals carried out at the core/contact interface. Intermediate states within an arc batch need not belong to the route.

## The projection to hold fixed

Write a guarded old tape as TT_a → D → TT_b with a,b≥1. Add a sentinel t beaten by every tape vertex, and a head x with t→x, x→tape except for contacts V⊆D with V→x. The tournament is strong: every old vertex reaches t then x, and x reaches every old vertex via the left guard or directly. In a core D let h=H(D), q=Q₂(D), E count ordered optional two-slot covers whose first slot ends in V, L those whose second slot starts in V, and J those with both. The exact endpoint rule from the [computational tape session](https://github.com/eliottcassidy2000/math/blob/6d38e8ceb0ccc8a061270e52507820bd4ecf53e6/05-knowledge/results/chanpak_computational_tapes_20261115/HEAD.md) specializes to

`H=(h−J)+2^(a+b)q+2^a E−2^b L.` (1)

One can recover (1) by classifying the predecessor of x: no predecessor gives h; t immediately before x gives 2^(a+b)q−2^b L; a predecessor in V gives 2^a E−J, because a forbidden restart in V forces the complete left guard to have been used before x. In each latter case the ordinal guards supply an unconstrained subset split with 2^a or 2^b choices, as appropriate.

There is a complete bounded-difference criterion for any fixed finite bank of cores/masks and an affine guard schedule a=a₀+uk, b=b₀+vk, nonnegative integer slopes. Group like powers of 2^k in (1). Two profiles give bounded difference at all heights if and only if their coefficient of every strictly positive exponential agrees: a nonzero largest differing coefficient dominates the rest. The remaining difference is the signed difference of the constant terms.

Consequences:
- If u,v>0 are unequal, this schedule really does force q,E,L separately equal.
- If u=v>0, the protected high terms are only q and 2^a₀ E−2^b₀ L. The two large contact counts can move together.
- If u=0, v>0, the only protected term is 2^a₀ q−L (the common factor 2^b₀ being immaterial), and the exact difference comes from h+2^a₀ E−J. All of h,q,E,L may move. The symmetric rule holds for the left-growing schedule.

The criterion concerns a fixed bank and specified schedules. It does not prohibit additional correlations when the bank itself grows.

## One strong core, symmetric guards, three consecutive odds

Lexicographic forward-bit code 1265214 defines a strong seven-vertex D. Bit ij is 1 when vertex i beats vertex j for i<j. Mask bits index the selected vertices. Here h=95 and q=1276.

| Mask | E | L | paths ending in V | paths starting in V | J |
|---:|---:|---:|---:|---:|---:|
| 14 | 425 | 467 | 29 | 35 | 124 |
| 38 | 445 | 487 | 33 | 39 | 126 |
| 42 | 445 | 487 | 33 | 39 | 128 |

All three have E−L=−42. With a=b≥1 they produce the consecutive odds

`1276·4^a−42·2^a−33, −31, −29`

where the second and third expressions use the same high terms. At a=1 these are 4987, 4989, 4991. Run masks 42→38→14 to increase the count by exactly 2 each time; both transitions toggle exactly two arcs to x. No internal arc, guard, or secondary nontransitive block is changed. In the final step E and L each decrease by 20, so their height-dependent changes cancel. This is a sequence of two repairs, not a commuting Boolean two-control cube.

## One moving core, one growing guard, six consecutive odds

Fix a=1. The following six order-eight cores are individually strong. They all satisfy 2q−L=7617, and their odd residual h+2E−J runs 3617, 3619, ..., 3627.

| Core code | Mask | h | q | E | L | J |
|---:|---:|---:|---:|---:|---:|---:|
| 75823804 | 45 | 311 | 5000 | 2081 | 2383 | 856 |
| 47535549 | 98 | 245 | 4304 | 1877 | 991 | 380 |
| 207924374 | 106 | 309 | 5160 | 2145 | 2703 | 978 |
| 257190115 | 116 | 305 | 5168 | 2165 | 2719 | 1012 |
| 253178809 | 51 | 333 | 5376 | 2217 | 3135 | 1142 |
| 114827917 | 58 | 279 | 4776 | 2040 | 1935 | 732 |

Thus one core position, with no added nontransitive blocks, supplies `7617·2^b+3617+2j` for 0≤j≤5 at every b≥1. The first step alone changes q by −696, L by −1392, h by −66, E by −204, J by −476, leaving exactly +2. The first six-odd interval starts at 18851 (b=1). Coherent local relabelings are included in the replay, realizing the consecutive steps with 6, 7, 4, 6, 6 flips respectively, including both internal and head contact arcs. Nothing outside the eight core vertices and their head reads is changed.

An order-seven width-four version has an even cleaner physical route: `2360·2^b+1271+2j` for 0≤j≤3. Its core/mask sequence is (1564245,15), (1888791,82), (101621,83), (758958,28), and every +2 step uses exactly two internal reversals and three contact reversals with coherent names. Positive results here are uniform in height, not just a calculation at the first few heights.

## Search direction and scope

The right projection tuple is small: for equal clocks hash by (q,E−L) and order by h−J, or for a=1 hash by 2q−L and order by h+2E−J. Hashing the unprojected counts would discard the largest moves above. For a genuinely reusable cube, the two (or more) moves must also preserve the exact mixed differences, not just hit neighbouring scalar targets. That issue and prescribing an arbitrary absolute baseline remain open here. These sparse dyadic lanes have constant window widths. They do not by themselves yield arbitrary windows or a proof of coverage.

Replay with `python3 check_multishift_cube.py`. It calculates both endpoint-coordinate tables and independently runs subset Hamilton-path DP on the full strong tournaments at concrete heights, as well as verifying every carried relabeling and the batch-flip counts.
