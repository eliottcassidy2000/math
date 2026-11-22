# Perfect centers: exact primitive seams at four early exponents

**FINITE-EXACT.** The common-cut grammar yields a different minimum
problem from the two scalar birth problems. The following four
positive certificates land at the exact two integers \(E_p\pm1\),
\(E_p=2^{p-1}(2^p-1)\), without testing their residues or estimating
their distance. Both endpoints are strong and differ in one arc.

| Exponent | Lower, upper counts | Joint order | Minimum joint order proved? |
|---|---:|---:|---|
| 3 | 27, 29 | 6 | yes |
| 5 | 495, 497 | 8 | yes |
| 7 | 8127, 8129 | 10 | yes |
| 13 | 33550335, 33550337 | 24 | no |

For \(p=2\), a straddle is impossible because \(E_2=6\) and its
upper neighbor 7 is a [known global hole](../chanpak_gap_star_ports_20261122/README.md).
This is one necessary qualification to an every-exponent theorem.
The first three order minima follow from the preceding-order maxima
15, 189, 3357, respectively, as supplied by the complete labelled-to-
unlabelled census through order nine in this packet. Ordinal singleton
padding makes the maximum nondecreasing with order.

## Outgoing masks and birth incidences

Bits use vertex 0 at the low end. In each row, use the lower masks
and reverse just the marked unordered pair for the upper tournament.

| \(p\) | Lower outgoing masks | Reversed pair |
|---|---|---|
| 3 | 30,52,56,50,32,1 | 1,5 |
| 5 | 160,145,67,7,13,30,59,124 | 0,2 |
| 7 | 988,433,74,34,236,645,938,780,60,286 | 8,9 |

Subset-endpoint DP checks the two literal counts on each pair.
At \(p=3\), the first 27 and 29 fibers have 2 and 4 classes
respectively. There are **2 primitive pair classes** and both
lower classes contribute. Each has three eligible arcs, forming one automorphism orbit. At \(p=5\), the 495 and 497 first fibers
have 22 and 33 classes; there are **18 primitive pair classes**,
with only **14** of the 22 lower classes incident to any such pair.
This makes the marginal/joint distinction concrete: knowing the
two minimum fiber sizes does not specify whether the smallest
classes meet by one arc.

## A strong 24-vertex certificate at \(p=13\)

Use the [signed two-arm reader](../chanpak_gap_star_ports_20261122/README.md)
with no left arm and right address \(y=23955\). The lower deletion
core has masks \((0,33,27,211,163,141,55,71)\), tag mask 55,
and the flipped arc joins vertices 3 and 7. The respective profiles
\((h,Q,(e,o),(e',o'),k)\) are
\[
(111,2856,(633,744),(184,245),309),\qquad
(111,2856,(633,744),(186,246),311).
\]
For odd \(y\), the adjusted intercepts are \(420\) and \(422\).
With \(p(y)=5080\) and unsigned start coefficient 1377, these
give
\[
  420+1377\cdot23955+111\cdot5080=33550335
\]
and exactly two more. Appending the reverse-complement contact
word and the apex gives the two physical hosts on 24 vertices.
The replay constructs both, checks their strongness, and verifies
their labelled difference is exactly the chosen arc. This is an
absolute exact suffix inversion, not a proof that 24 is minimal.

Together with the earlier [exponents 17 and 19](../chanpak_gap_star_ports_20261122/PERFECT-TARGETS.md)
and [23 and 31](../chanpak_e31_physical_20261122/README.md), there
are exact physical examples at \(3,5,7,13,17,19,23,31\) on record
here. The Mersenne factor at 23 is composite; no primality is
needed for the exact count target. The general straddling
construction for unbounded exponents remains **OPEN**.

Run `python3 check_perfect_birth_seams.py`. To replay the class
multiplicities and minimality checks, pass the full order-nine TSV
emitted by the census program documented in [the birth note](BIRTH-MULTIPLICITIES.md).
