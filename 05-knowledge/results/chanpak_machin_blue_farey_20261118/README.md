# Native blue repair and the Machin/Farey tournament bracket

Follow-up to [the Machin and signed Collatz lift](../chanpak_machin_collatz_lift_20261118/README.md). This packet uses ratios of two independently realized Hamiltonian counts; it does not identify them with the single-tournament slot ratio P₂/(2H).

## 1. The exact repair is a blue metagraph line

With the marked path 0→1→⋯→6, these outgoing bit rows are mates:

    [46,68,8,82,103,78,5]     H=69
    [82,60,121,49,32,65,26]  H=71

The two tournaments are invariant under arc reversal followed by label reflection i→6−i, and both have the marked path. Flipping exactly the nonpath arcs takes one to the other. This is the native blue operation, not the one-arc 71↔73 pair of the previous packet. Both endpoints first occur at order seven (the order-six maximum is 45; transitive singleton padding preserves H). Among the 512 reflection-compatible order-seven frames there are seven frames with H=71; their mate counts, each with its frame counted, are 13,45,69,93,133,137,145. This is a frame count, not an isomorphism-class degree census.

## 2. Exact rational bounds and the first possible improvements

Two ratios made from actual tournament counts bracket π:

    311/99 < π < 355/113.

The two fractions have cross-determinant 2, as small as possible for unequal odd/odd fractions. Their primitive midpoint vector is (333,106), and that vector splits the interval into two cones, each of determinant 1.

Any reduced fraction strictly inside the lower cone has a unique positive integer expression
((p,q)=a(311,99)+b(333,106)); similarly use the endpoint (355,113) in the upper cone. Odd/odd parity forces a odd and b even. This proves there is no odd/odd improvement to either endpoint with denominator at most 310. The first lower improvement is 977/311, with (a,b)=(1,2). It is below π by exact rational interval arithmetic in Machin's formula.

The upper bound survives much farther. Write ((p,q)=a(333,106)+b(355,113)) with a positive even and b positive odd, and put δ=355/113−π>0. The distance down from 355/113 is exactly a/(113q), so a new upper bound needs 113qδ>a. For a=2 the smallest working odd b is 587, giving

    209051/66543.

Indeed b=585 gives 208341/66317<π, and b=587 gives 209051/66543>π. Machin intervals certify both signs and 113·66543·δ<4, excluding every a≥4 at a denominator this small. The conclusion holds for all odd/odd ratios, reduced or not (reduction only decreases the denominator), hence also for ratios of tournament counts. It is a least-denominator assertion, not a best-error estimate at larger denominators.

The improvement is itself physical: 209051=349·599 and 66543=3·41·541. Ordinal sums multiply the path counts, so small verified factors suffice on both sides. The replay contains literal witnesses for those factors, for the initial bracket, for 977, and for 223/71. The Archimedean upper 22/7 is impossible as a ratio of tournament counts even after common scaling: an odd numerator over an odd denominator stays odd/odd in lowest terms. No claim about the whole Chan–Pak spectrum follows.

## 3. Odd-root rays and Pell

For distinct positive odd roots r<s, consider 2^i r and 2^j s. If i=j the gap is a power of two exactly when s−r is a power of two. If i≠j, divide out the smaller power. The difference inside is odd, so the gap is a power of two exactly when the inside is ±1. The only possible direction of unequal height has the smaller root at higher height, with |2^d r−s|=1, d>0. For s=r+2 this reduces to (2^d−1)r in {1,3}. The complete list of (r,d) is (1,1),(3,1),(1,2). For roots 7,21,71, direct substitution gives no unequal-height solution and only the aligned 7 and 71 gap, 2^(i+6). This is a limitation on connecting odd roots by pure-power moves, not an obstruction to more general changes of tournament core.

Machin's 239 also lies in the exact identity (1+√2)^7=239+169√2, on the Pell sequence of real parts 1,7,41,239 of odd powers. The previous packet's Gaussian pair (119,120) has hypotenuse 169, so the Pell and Gaussian identities share a genuine norm equation, not just a numeral.

Run `node verify.js` for independent exact arithmetic, the order-six maximum, a complete blue order-seven frame check, and the literal ratio witnesses.
