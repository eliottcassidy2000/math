# Thirteen fourth-power values, thirty-seven gaps, and the aliquot rays

The proposed pattern is exact once its domain is made explicit: **among the positive integers 1 through 50, thirteen are sums of four nonnegative fourth powers and thirty-seven are not.** Including zero adds one representable value, so the corresponding count on 0 through 50 is 14 versus 37. Zeros are allowed as summands throughout; “four powers” therefore includes representations padded with zeros.

This packet develops two precise mechanisms rather than identifying unrelated appearances of small integers:

- [Fourth powers and their gaps](FOURTH_POWERS.md): exact representations, the four-term budget, modular exclusions, and the 11/12/13 gap lengths.
- [Aliquot fibres and doubling rays](ALIQUOT.md): complete source sets, retained prime weights, and the exact affine appearance of 37 in the 25-ray.

## 1. The positive-domain 13/37 split

Below 81 the available fourth powers are 0, 1, 16. Up through 50 every representation therefore has the form

\[
16a+b,\qquad a,b\ge0,\qquad a+b\le4.
\]

The bands are 0–4,16–19,32–34,48–49, of sizes 5,4,3,2. Their total is 14; deleting zero gives 13. In the positive interval 1–50 the complementary gaps are 5–15,20–31,35–47, and 50, with lengths 11,12,13,1.

The distinction between congruence and resource budget is important. Every fourth power is 0 or 1 modulo 16, so four of them cannot sum to 37, whose residue is 5. The number 13 is also excluded by this modular test, although it counts the positive representable values in this particular interval. The endpoint 50 has an allowed residue 2, but would require three 16s and two 1s: five summands rather than four. The final singleton is an exact budget obstruction, not another modular obstruction.

## 2. The old aliquot association really is an inverse fibre

For positive source integers define \(s(n)=\sigma(n)-n\), the sum of proper divisors, and \(F(m)=\{n\ge1:s(n)=m\}\). The [prior Session69 report](https://github.com/eliottcassidy2000/math/blob/616f9d0a5f068ecaa6ed692621882a60833908e8/05-knowledge/results/session69_prime_ports_power_aliquot_20261026/session69-aliquot/REPORT.md) gives the complete fibres

\[
F(11)=\{21\},\quad F(12)=\{121\},\quad F(13)=\{27,35\}.
\]

Their source shapes differ: semiprime, prime square, and a collision between prime cube and semiprime. Completeness follows from the inherited bound \(n\le(m-1)^2\) for any source of a target \(m>1\), not merely from a convenient search cutoff. The [canonical source-shape theorem](https://github.com/eliottcassidy2000/math/blob/616f9d0a5f068ecaa6ed692621882a60833908e8/01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md) explains why these three shapes suffice throughout targets 2–14.

The [prior complete table through 45](https://github.com/eliottcassidy2000/math/blob/616f9d0a5f068ecaa6ed692621882a60833908e8/04-computation/experiments/aliquot_fibres_to_45_20261027.out) also gives

\[
F(37)=\{155,203,299,323\}.
\]

Thus 37 is not an aliquot gap. Neither is 13. Their fourth-power exclusions cannot be transported to aliquot untouchability.

## 3. A fresh exact 37 connection in the requested ray

The [earlier chain work](https://github.com/eliottcassidy2000/math/blob/616f9d0a5f068ecaa6ed692621882a60833908e8/05-knowledge/results/aliquot_power_124_20261026/chains/RESULT.md) preserves doubling/power addresses. Following 5→10 by squaring reaches 100; following 5→25 by two doublings, 25→50→100, reaches the same value: \(10^2=4\cdot25=100\).

These are operation edges, not aliquot edges. For an odd core m the fresh specialization of the divisor formula is

\[
s(2^a m)=2^a\bigl(m+2s(m)\bigr)-\sigma(m).
\]

Because \(s(25)=6\) and \(\sigma(25)=31\), it yields

\[
\boxed{s(2^a25)=37\cdot2^a-31}.
\]

The values at 25,50,100 are 6,43,117. This makes 37 a genuine affine coefficient of the requested family. It does **not** make that coefficient the cause of the fourth-power complement size.

## 4. Two exact observations restore the source

The carry is not only useful for prediction. For every positive integer n, put x=s(n), y=s(2n), c=y−2x and d=y−x. Then d>0 and

\[
n=2^{v_2(d)+1}c-d.
\]

The proof writes n=2^a m with m odd. The difference is d=2^a(2σ(m)−m), whose second factor is positive and odd, so its valuation recovers the original doubling exponent exactly. For 25, 50 and 100, the differences are 37, 74 and 148, recovering exponents 0, 1 and 2.

This works for a specified prime dilation as well, with the exact formula and necessary image guards in [ALIQUOT.md](ALIQUOT.md). It does not invert a single aliquot observation or prove termination of aliquot iteration. Nor can residues replace exact values: the sources 27 and 35 have identical entire doubling traces modulo 8, although their exact next outputs are 66 and 74.

## Scope and inheritance

Bounded source recovery found prior cyclotomic 13, exact 37 fibres, and operation-address laws. It did **not** identify an earlier paired 13/37 theorem explaining this fourth-power partition. The ray specialization is derived here as an elementary consequence, without a historical-priority claim. Search co-occurrences and numerical tables are not enough to establish such provenance.

Prime-base guards remain essential: \(s(p^3)=1+p+p^2\) holds for prime p, but substituting 4 falsely predicts \(s(64)=21\); the actual value is 63. Finally, 2 and 5 are proved untouchable, but these computations do not prove they are the only prime untouchables or settle any asymptotic density.
