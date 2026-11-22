# Minimum literal seams at exponents 4 and 6

**FINITE-EXACT.** With row masks indexed from zero, the following strong tournaments reach both endpoints at minimum order by the indicated single reversal:

| \(p\) | Order | Lower outgoing masks | Reversed arc pair | Counts |
|---:|---:|---|---|---|
| 4 | 8 | 192,13,1,5,47,207,30,94 | 0,1 | 119,121 |
| 6 | 9 | 266,436,121,194,201,345,3,101,220 | 4,8 | 2015,2017 |

The witnesses were located by fixed-seed random orientation search, but the counts and strong connectivity are independent of that search. [Reproduce](check_small_exponents.py). For \(p=4\), generate the 56 order-six classes by taking one per permutation orbit, then extend each by all 64 one-vertex incidences. Every order-seven class has such a presentation. None is even a lower endpoint of count 119, so no tournament of order at most seven can realize that count (pad by a dominated vertex to order seven). For \(p=6\), the preceding order-eight maximum is 661 from the [first-birth census](../chanpak_support_birth_grammar_20261122/BIRTH-MULTIPLICITIES.md), below 2015.

Together with the [previous exact p=3,5,7 minima](../chanpak_support_birth_grammar_20261122/PERFECT-BIRTH-SEAMS.md), every integer exponent from 3 through 7 has a literal minimum-order seam. This does not supply an induction across exponents.
