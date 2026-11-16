# Exact repair menu on the transitive mirrored seam

**Proved, all orders.** This is the complete ordered +2 graph for the mirrored transitive core, not a finite search for examples.

Take two adjacent mirror copies of a transitive tournament of order n, with the same nonempty guards of length a on both sides. In the left core choose a mask V, and in the other core use its converse-reflected copy. The [signed-motion mirror identity](https://github.com/eliottcassidy2000/math/blob/c238b05d4b69f332d389ff8be0c40f6d4bfce189/05-knowledge/results/chanpak_signed_motion_20261115/MIRROR.md) says the resulting strong tournament has count

    H = 2^(2n+2a) + 1 − M(x),

where `x=Σ_{v∈V}2^(n−1−v)`, `M(x)=x²+(x mod 2)+2J(x)`, and

    J(x) = Σ_{i>j, bits i and j of x both set} (2^(i−j−1)+[i=j+1]).

Indeed the h of the transitive core is 1, its optional two-slot count is 2^n, L(V)=x, and T(V) is the last-vertex bit, now the low bit of x. For the J formula, take two marked vertices u<v in transitive order. If the first slot ends at v and second starts at u, the middle vertices may be assigned independently and all outer vertices are forced. For the opposite endpoints, all vertices are forced and the cover exists exactly when v=u+1. Their combined contribution is the displayed pair weight.

## The exact carry law

For 0≤x<2^n−1 write x=`2^(r+1)g+2^r−1`, where r=v2(x+1), so bit r is zero and bits below it are one. Then

    (M(x+1)−M(x))/2 = 1 + 3g + (g mod 2).

For r=0 adding the low bit increases J by g+[g odd], and the square and last-vertex contribution increases by 4g+2; this gives the formula. For r≥1 put B=2^(r+1)g for the unchanged high word. The old all-one tail has interaction 2^r−2: summing 2^(i−j−1) for 0≤j<i<r gives 2^r−r−1 and adds r−1 adjacent pairs. The new single set bit at r has no internal pair. For each high set bit i>r, its old pairs with the tail contribute 2^i−2^(i−r), whereas its new pair with bit r contributes 2^(i−r−1)+[i=r+1]. Thus the tail/high interaction changes by −B+3g+[g odd]. Adding the loss of the old tail and the square-plus-parity change 2x gives exactly 2+6g+2[g odd].

This law is in terms of the actual count, not a residue. The energy is strictly increasing, so two distinct masks have outputs differing by two iff they are neighbors in x order and g=0. In other words the only edges are at x=2^r−1→2^r for 0≤r<n. The first two edges (when n≥2) form the three-count window with energies 0,2,4; every other primitive edge is isolated. For n=1 there is one edge and no triple. No arbitrary-width bank can come from a mirror of two increasing transitive cores.

The half-gap samples precisely the positive integers coprime to 6 as g ranges over all nonnegative integers: it is 3g+1 for even g and 3g+2 for odd g. In particular the same arithmetic pattern that has appeared in the Collatz and 2/3 archive here has an exact physical meaning: after discarding the trailing-one carry, the untouched high word prices the next seam move. This is an identity for this repair family, not a claim of a dynamical conjugacy to Collatz.

## What this suggests

The old [mirror-complement address ledger THM-2079](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-2079-mirror-complement-safe-child-addresses.md) separates reversing a binary address from reversing the endpoint ownership. Here the two needed records are the carry word g and its bottom bit. Keeping just the longest trailing run is not enough to predict the signed cost. For a nontransitive core, one wants a bank whose pair interaction moves independently of this 1+3g tariff, rather than a longer copy of this special mirror.

## Odd core widths: the same carry tariff, but no triple

An independent background archive pass obtained the odd-width analogue. For a single transitive core on 2n+1 vertices, at equal guards the zero signed projection selects exactly reflection-symmetric masks: the head counts E,L are unsigned forward and backward binary addresses, so E=L iff the mask is palindromic. Write P for the binary address of occupied mirror pairs, nearest pair first, and t for the central contact. Let K(P) be the same pair expression as J above, but only among the chosen pairs on one half. The normalized residual J/2 is

    F(P,t) = P² + K(P) + t(P+[P odd]).

To verify: same-half interactions contribute K after division by two; opposite-half mirror pairs at distance 2j+2 supply the diagonal P² term, and pairs of different radii supply its cross terms. The central contact at distance j+1 contributes 2^j+[j=0] per selected radius in the normalized residual. The two P=0 states have identical output, and at fixed nonzero P the gap between t=0 and t=1 is P+[P odd]≥2.

Now let P=2^(r+1)g+2^r−1 as above. The gap from F(P,1) to F(P+1,0) is P+1−[P odd]+ΔK. The same carry expansion proved above gives *exactly* 1+3g+[g odd], including r=0. Hence the distinct energy states are strictly ordered lexicographically; the only actual +2 edges again occur at P=2^r−1, but every positive P has an intervening vertical gap >1, so none concatenate. There are exactly n different +2 edges at width 2n+1 and at width 2n, with at most two consecutive odds in the odd case and three in the even case (when n≥2). The large edge baselines separate geometrically: a new isolated radius at P=2^r has normalized energy 4^r in odd width, or 2^(2r−1) in even width for r≥1. This settles the symmetric zero-projection fiber, not arbitrary transitive masks.
