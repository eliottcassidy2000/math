# Physical signed-quarter edges and reduced prime exception rays

The [signed arithmetic selector](CARRY_PARENT.md) decreases x by choosing the sign s such that 4 divides 3x+s, then taking the odd part. These notes identify a family with actual tournament witnesses at **both endpoints on both axes**, and sharpen which defective parents matter if the proof follows only prime quotients.

## An all-order family with literal endpoints

For k≥8 take either sign and put x=2^k±1. Its signed quarter parent is exactly y=3·2^(k−2)±1. The ordinary x and upper 7x both have the physical dyadic-plus/minus constructions in [the signed-height packet](../chanpak_six_signed_height_pairs_20261118/README.md). The positive y and 7y also have the mantissa-three positive readers from that packet. The negative y and 7y use the [paired negative reader](README.md), with slopes 96,672 and intercepts −1,−7, on transitive width k−7:
```
96·2^(k−7)−1 = 3·2^(k−2)−1.
```
Thus these signed arithmetic arrows connect four literal outputs, with explicit choices encoded by the absolute exponent. The construction of x need not mutate the earlier y witness; the four graphs are obtained from different contexts. This proves a real family of arithmetic edges with physical witness certificates, not an inverse for arbitrary sources.

## Which defective queries survive restriction to primes?

When using ordinal products to handle composites, it suffices to place prime atoms on the ordinary and sevenfold axes, with the already realized additional count seeds 49,63,343 (and the count-five and triangle bases). This does **not** assert that a prime target's signed parent is itself prime. A legal composite parent may be constructed separately by joining already attained smaller factors. For instance, for an upper request x whose parent y is composite and not divisible by 3, write y with its smaller prime factors. To get 7y, use either a sevenfold witness on a prime factor other than 7, or write the resulting pure power of 7 as a product of the count-49 and count-343 seeds. If y is only 7×a unit, use 49 directly. This accounts for the special multiplicative axis, rather than dividing an arbitrary path count.

The only bad parents for the signed selector are 7 on the ordinary axis and 1 on the upper quotient axis. For upper bad parent 1, x=(2^a−(−1)^a)/3, a≥2.

- For even a=2k≥6, the factorization (2^k−1)(2^k+1)/3 is composite: one of the factors contains 3, and even after its removal both factors exceed 1. At a=4 the prime is 5, with known partner 35.
- For odd composite a, a proper odd divisor d≥3 gives the divisor (2^d+1)/3 via the usual odd-exponent factorization of a sum of powers. Both factors exceed 1.
- Therefore the only additional prime-side upper defects are actual primes of the form (2^ell+1)/3 with ell itself an odd prime ≥5. This is a restriction on the exceptional exponents, not a claim that those expressions are generally prime.

For ordinary bad parent 7, the ray x=(7·2^a−(−1)^a)/3 cannot give prime targets when a≡2 modulo 3: the numerator is divisible by 9 (use −2≡7 modulo 9), and x≥9.

One route would therefore focus physical chart changes on these thin rays rather than try to repair every problematic integer with the same site. The all-target unknown is still the physical inverse for general signed parent edges; primality screening alone does not solve it. Nor does this reduction assume a relation between prime status and Collatz convergence.

## Actual chart changes at defects

The replay checks literal single-arc changes 17→19, 73→75, 297→299 on the ordinary side and 133→147, 287→301, 581→595 on the upper side. These have the exact repair gaps 2 and 14, respectively. We do not infer a universal repair from the examples.

There is a broader infinite bypass for some composite upper defects. Put R_k=(4^k−1)/3. If 7R_u is realized, then so is 7R_(u·2^j) for each j≥0, by the difference-of-squares identity R_(2k)=R_k(1+4^k); each new factor is a longest-arc reversal. Besides u=2, verified base values at u=3,5,7,9 are 147; 77·31; 301·127; and 147·3·19·73. Their factors are ordinary path counts: 31,127 are long-arc reversals, 3 is a triangle, 77 and 19 are in [the main replay](CARRY_PARENT.md), and the single-arc controls above give 147 and 301. The lower endpoint in the second single-arc ordinary control supplies 73. This enlarges the constructible defective family without treating the forbidden quotient as a valid parent.

Replay: `python3 signed_lane_replay.py`.
