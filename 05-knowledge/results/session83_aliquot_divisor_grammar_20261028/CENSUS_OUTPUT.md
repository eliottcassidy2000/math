# Reproduced aliquot shape census: targets 2..2000

Independent divisor-addition sieve through n=3996001: PASS.

Full sorted fiber SHA-256: `01c598ee1831b6c1b255950803fd0271e66934fc408a7cec6c752663d83adfe9`

## Source incidences and singleton fibers

| Class | Source incidences | Singleton fibers | Square-sigma incidences |
|---|---:|---:|---:|
| I | 332 | 148 | 2 |
| II | 26647 | 87 | 283 |
| III | 1791 | 170 | 74 |

Total source incidences: 28770.

## Per-class multiplicity

| Class | Absent | Exactly one | At least two | Maximum | Maximizing targets |
|---|---:|---:|---:|---:|---|
| I | 1668 | 330 | 1 | 2 | [31] |
| II | 702 | 306 | 991 | 91 | [1891] |
| III | 871 | 696 | 432 | 6 | [1483, 1585] |
## Class presence masks

| Present classes | Target count |
|---|---:|
| none | 196 |
| I | 148 |
| II | 476 |
| III | 260 |
| I,II | 51 |
| I,III | 98 |
| II,III | 736 |
| I,II,III | 34 |

## First witnesses

- all three classes: t=40, classes={'III': [44], 'II': [74], 'I': [81]}.
- only III, nonempty: t=28, classes={'III': [28]}.
- I and III present: t=40, classes={'III': [44], 'II': [74], 'I': [81]}.
- II and III present: t=16, classes={'III': [12], 'II': [26]}.
- multiple I sources: t=31, classes={'I': [32, 125], 'II': [161, 209, 221]}.
- multiple II sources: t=17, classes={'II': [39, 55]}.
- multiple III sources: t=64, classes={'III': [56, 76], 'II': [122]}.
- Maximum full fiber size: 92; targets=[1891].
- Odd targets with only III: [].

## Exact source shape stratification

Exponents are an unordered decreasing tuple; the primes are distinct.

| Class | Exponent shape | Proper-divisor terms | omega | Incidences |
|---|---|---:|---:|---:|
| I | (2,) | 2 | 1 | 303 |
| I | (3,) | 3 | 1 | 14 |
| I | (4,) | 4 | 1 | 5 |
| I | (5,) | 5 | 1 | 3 |
| I | (6,) | 6 | 1 | 2 |
| I | (7,) | 7 | 1 | 2 |
| I | (8,) | 8 | 1 | 1 |
| I | (9,) | 9 | 1 | 1 |
| I | (10,) | 10 | 1 | 1 |
| II | (1, 1) | 3 | 2 | 26647 |
| III | (2, 1) | 5 | 2 | 486 |
| III | (3, 1) | 7 | 2 | 120 |
| III | (1, 1, 1) | 7 | 3 | 731 |
| III | (2, 2) | 8 | 2 | 18 |
| III | (4, 1) | 9 | 2 | 47 |
| III | (3, 2) | 11 | 2 | 11 |
| III | (5, 1) | 11 | 2 | 22 |
| III | (2, 1, 1) | 11 | 3 | 193 |
| III | (6, 1) | 13 | 2 | 10 |
| III | (4, 2) | 14 | 2 | 5 |
| III | (3, 3) | 15 | 2 | 2 |
| III | (7, 1) | 15 | 2 | 5 |
| III | (3, 1, 1) | 15 | 3 | 45 |
| III | (1, 1, 1, 1) | 15 | 4 | 35 |
| III | (5, 2) | 17 | 2 | 3 |
| III | (8, 1) | 17 | 2 | 2 |
| III | (2, 2, 1) | 17 | 3 | 18 |
| III | (4, 3) | 19 | 2 | 2 |
| III | (4, 1, 1) | 19 | 3 | 14 |
| III | (6, 2) | 20 | 2 | 1 |
| III | (5, 3) | 23 | 2 | 1 |
| III | (3, 2, 1) | 23 | 3 | 7 |
| III | (5, 1, 1) | 23 | 3 | 4 |
| III | (2, 1, 1, 1) | 23 | 4 | 7 |
| III | (2, 2, 2) | 26 | 3 | 1 |
| III | (4, 2, 1) | 29 | 3 | 1 |

## Parity stratification

| Class | Target parity | Source parity | Incidences |
|---|---|---|---:|
| I | even | odd | 307 |
| I | odd | even | 9 |
| I | odd | odd | 16 |
| II | even | even | 301 |
| II | odd | odd | 26346 |
| III | even | even | 779 |
| III | even | odd | 11 |
| III | odd | even | 40 |
| III | odd | odd | 961 |

## Earliest square-sigma source-target pairs

- I: [(40, 81), (57, 343)] (pairs are (target, source)).
- II: [(14, 22), (25, 119), (29, 115), (39, 217), (49, 527), (50, 94), (59, 265), (59, 517), (77, 1219), (79, 497), (89, 1207), (99, 1501), (103, 2201), (105, 679), (110, 214)] (pairs are (target, source)).
- III: [(74, 70), (78, 66), (154, 170), (191, 385), (219, 357), (231, 345), (254, 322), (266, 310), (294, 282), (361, 935), (366, 210), (373, 651), (420, 364), (428, 472), (501, 795)] (pairs are (target, source)).

PASS: complete fibers, independent oracle, factor re-encoding, shape, parity, gap, square-edge, and hostile controls.
