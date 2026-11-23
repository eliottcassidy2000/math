# Kaprekar registers, decimal tournament counts, and exact digit shelves

**PROVED scoped / FINITE-EXACT.** This extends [the Heegner/absolute-target work](../chanpak_heegner_absolute_targets_20261123/README.md), [guarded slot synthesis](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md), and the repo's [earlier Kaprekar census](../kaprekar_cuboid_20260928.out). The ordinary Hamilton-path spectrum still has its open Chan–Pak coverage question.

- [Sorted subtraction and actual tournament images](KAPREKAR-IMAGES.md): exact memory from the last half of the digits *after* one step. Two-digit even-base dynamics gives precisely the doubling action from the five Heegner primes on the quotient by sign. Every periodic height at those five bases is realized, as is every odd decimal Kaprekar image at widths 2–5. The full nonzero binary image is realized except for \(7,21\).
- [Decimal congruences and height windows](DECIMAL-WINDOWS.md): guarded words now set an ordinary Hamilton count's last \(k\) digits to any odd ending while independently steering its real size to \(O_{10^k}(\sqrt X)\) of \(X\). The factor-seven parent and exits have their corresponding controlled windows. Independent joint \((H,Q_2)\) decimal phases have genuine 2-adic restrictions, and the full pair need not be height-balanced.
- [Physical decimal append and opposing arcs](DECIMAL-SHELF.md): five literal strong constructions realize \(10h+d\) for each odd digit \(d\) and every \(h=2^k+1,\ k\ge4\). A single arc crosses the even fixed point \(6174\), connecting actual counts \(6173\) and \(6175\) in order ten. This is a marked pair, not a tournament with even Hamilton count or a repeatable operation on decimal expansions.

## Why the decimal quotient helps, and where it stops

The factor-seven fork has center \(C=840Q_2+21H\), with offsets \(\pm2H\). Therefore the three readings reduce to
\[
(-H,\ H,\ 3H)\pmod{10}.
\]
Last digit 5 collapses the three phases to 5; each other odd digit gives three distinct phases. The missing absolute values \(7,21\) cannot be separated from realized counts by a last-digit obstruction (indeed their ordinary residues are realized modulo **any** fixed modulus). One must retain more than just the terminal phase.

The *arithmetic* maps \(x\mapsto3x-2,3x,3x+2\), on the half-index \(z=(x-1)/2\bmod5\), are \(z\mapsto3z,3z+1,3z+2\). They generate the group of all affine bijections of \(\mathbb F_5\): multiplication by 3 has order four; composing its inverse with \(3z+1\) gives a translation generator. This is not an iteration theorem for physical fork arcs. Nor is it the two-digit Kaprekar clock: that clock is multiplication by 2 modulo 11 **up to sign**, on the five odd representatives. Both are exact descriptions, but the state coordinates and realizability conditions are different.

A useful negative comparison is accelerated Collatz. For \(a=\nu_2(3x+1)\), computing \((3x+1)/2^a\bmod10^k\) from residues needs \(x\bmod 2^{k+a}5^k\) as an adequate input, not merely its last \(k\) decimal digits. Even knowing \(a\) need not repair the information lost: if \(1\le a<k\), \(x\) and \(x+10^k\) have the same \(a\), but their quotients differ by \(3\cdot10^k/2^a\), not zero modulo \(10^k\). Kaprekar's exact compression works only because of its special image, not because truncation generally commutes with digit dynamics.

The exact-lift criterion remains simple: to force \(H=N\) from equality mod \(10^k\), prove \(|H-N|<10^k\) (or produce a literal inverse). Fixed-precision window bounds are not estimates uniform in \(k\), and so do not close the general exact-target gap.

**Sharp arc itinerary.** [Actual Kaprekar walks](KAPREKAR-METAGRAPH.md) need seven edits up to relabelling and eight with fixed labels; both are attained. The direct five-cycle fails at all orders, and so does the [first jump of the base-16 orbit containing 15](BASE16-OBSTRUCTION.md). [General finite order caps](ITINERARY-CAPS.md) apply to open or labelled closed itineraries with prescribed checkpoint heights.
