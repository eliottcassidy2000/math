# Thue–Morse gives the exact lonely chord on every Fermat grid

This extends [the dihedral chord-grid work](../fermat_dihedral_tournaments_and_lonely_chord_grids_20261127/README.md), [the \(9/4\) half-cycles](../abc_pythagorean_fermat_half_cycles_20261127/README.md) and the [period construction of Fermat polygons](../fermat_polygon_chord_energies_and_period_splits_20261127/README.md). Here a different speed pattern, consecutive powers of two, has one exactly determined maximizing orbit.

Put \(k=2^m,\ F_m=2^k+1\), \(m\ge1\). For speeds \(1,2,\ldots,2^{k-1}\), all consecutive powers of two, the best least circular distance over nonzero times \(a/F_m\) is
\[
h_m=M_m/F_m,\qquad M_m=\prod_{j=1}^{m-1}(2^{2^j}-1).
\]
There is exactly one maximizing orbit, the \(2k\) residues \(2^j M_m\bmod F_m\). The best shortest **squared chord** is \(4\sin^2(\pi h_m)\). Neither the formula nor uniqueness requires \(F_m\) to be prime.

Indeed \(2^k\equiv-1\bmod F_m\); the reduced denominator of any nonzero point stays odd and greater than 1. Its doubling period divides \(2k\) but not \(k\), hence is \(2k\), and the second half is the reflection of the first. The minimax problem is therefore the symmetric period-\(2k\) survivor problem for doubling. The known Thue–Morse periodic extremal-word theorem ([Allouche–Clarke–Sidorov, Propositions 2.15–2.16](https://arxiv.org/pdf/0803.4439); explicit symmetric-window formula in [Ermann, §3.2, equation (4)](https://arxiv.org/html/2611.28242v1#S3.SS2)) gives the unrestricted exact-period upper bound. For the transfer, any candidate word here is complement-closed under rotation, so its largest shift also controls both edges of the symmetric window. The extremal Thue–Morse word itself is a block followed by its complement, so it lies on this grid. In the explicit period formula with \(\nu=m+1\), cancel \(1-2^{-k}\) and multiply by \(F_m\) to get \(M_m\). For uniqueness, an orbit attaining \(M_m/F_m\) must contain \(M_m\) or \(-M_m\) in its numerator; these are on the same orbit.

## Nonsquares and the precise polygon resolution

If \(F_m\) is prime with \(m\ge2\), all the winners are nonsquares. Indeed \(M_m=\prod_{i=0}^{m-2}F_i^{m-1-i}\). By reciprocity and the Fermat product identity, the symbol \((F_i/F_m)\) equals the Jacobi symbol \((2/F_i)\) (even if an earlier \(F_i\) is composite). This is negative only at \(i=0,1\); the total contributing exponent is odd (one if \(m=2\), \(2m-3\) otherwise). And 2 is a square mod \(F_m\), so the entire orbit has this sign.

More precisely the polygon labels of all winners form the single coset \(M_m\langle2\rangle\). Its index is \(2^{k-m-1}\). In the dyadic period tree it is therefore **one full group at depth \(k-m-1\)**, not merely a subset selected by the quadratic character. A multiplicative character indexed by \(j\) relative to any primitive root survives averaging on the winning group precisely if \(2k\mid j\); the first informative such frequency is \(2k\), not necessarily the quadratic frequency. At 17 the required depth is 1, so the winner is exactly the nonsquare \(9/4\) half-cycle. At 257 and 65537 the depths are 4 and 11: the square/nonsquare test is a first cut, not the whole minimax solution.

| Grid | Best distance | Best on squares (finite exact) | Winning labels |
|---:|---:|---:|---:|
| 5 | \(1/5\) | \(1/5\) | 4 |
| 17 | \(3/17\) | \(1/17\) | 8 |
| 257 | \(45/257\) | \(25/257\) | 16 |
| 65537 | \(11475/65537\) | \(11061/65537\) | 32 |

The square-only maxima in this table are scans, not extrapolations. In particular at 257 and 65537 even the best strictly suboptimal orbit is still nonsquare in the scan.

## The Pollock carry and the grid-to-continuous gap

The same signed polynomial
\[
C_m(z)=\prod_{j=0}^{m-1}(1-z^{2^j})
\]
gives \(|C_m(2)|=M_m\) and vanishes to order \(m\) at 1. Its coefficients are the Thue–Morse signs, so the zero at 1 gives precisely the lower-moment cancellations of the [dyadic simplex carry](../fourth_simplex_parity_cube_and_matching_response_20261126/README.md). The value at 2 is the lonely-grid numerator, while the order at 1 is the moment cancellation. This does not assert an equivalence of the Pollock and runner spectra.

The continuous optimum is \(1/3\) for every \(k\ge2\): \(t=1/3\) works since the speeds alternate signs mod 3, and speeds 1 and 2 alone give the upper bound. The grid optimum decreases with \(m\), by the exact ratio \(h_{m+1}/h_m=1-2/F_{m+1}\), yet the grid **fails the conjectural threshold exactly for \(m=1,2\)**. Indeed these two values are \(1/5<1/3\) and \(3/17<1/5\). For \(m\ge3\), the period formula gives \(h_m=\prod_{j=1}^{m-1}(1-2^{-2^j})/[4(1+2^{-k})]\ge2/15>1/9\), using \(2^j\ge2j\) and \(\prod(1-a_j)\ge1-\sum a_j\). Now \(1/(k+1)\le1/9\). Thus a failure at 17 is a genuine grid obstruction despite an elementary continuous witness, but it does not propagate to the later grids of this particular family.

Run `node replay.js` for exhaustive integer scans through 65537 and independent BigInt checks of the predicted orbit through \(m=10\), including composite Fermat numbers.
