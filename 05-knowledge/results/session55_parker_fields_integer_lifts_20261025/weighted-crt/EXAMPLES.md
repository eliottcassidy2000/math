# Explicit complete example histograms

**FINITE-EXACT**, reproduced by both enumerations and CRT. Partitions are restricted-growth strings in row-major cell order: equal digits mean equal entries. Every partition not listed has count zero. No root multiplicities are counted.

## Modulus 35

| Exact partition | Collision mask | Arrays |
|---|---|---:|
| 012345678 | 0x0 | 24 |
| 012201120 | 0x61089888 | 38 |
| 012120201 | 0x8a054250 | 38 |
| 000000000 | 0xfffffffff | 12 |

Total: 112; nine-distinct labeled: 24; D4 position-orbits among distinct arrays: 3.

Local modulus 5: 012201120 has weight 2, 012120201 has weight 2, 000000000 has weight 3.

Local modulus 7: 012201120 has weight 6, 012120201 has weight 6, 000000000 has weight 4.

## Modulus 55

| Exact partition | Collision mask | Arrays |
|---|---|---:|
| 012345678 | 0x0 | 40 |
| 012201120 | 0x61089888 | 62 |
| 012120201 | 0x8a054250 | 62 |
| 000000000 | 0xfffffffff | 18 |

Total: 182; nine-distinct labeled: 40; D4 position-orbits among distinct arrays: 5.

Local modulus 11: 012201120 has weight 10, 012120201 has weight 10, 000000000 has weight 6.

Local modulus 5: 012201120 has weight 2, 012120201 has weight 2, 000000000 has weight 3.

## Modulus 201

| Exact partition | Collision mask | Arrays |
|---|---|---:|
| 012314012 | 0x10102420 | 264 |
| 012201120 | 0x61089888 | 1056 |
| 012120201 | 0x8a054250 | 1056 |
| 010222343 | 0x404600002 | 264 |
| 000000000 | 0xfffffffff | 68 |

Total: 2708; nine-distinct labeled: 0; D4 position-orbits among distinct arrays: 0.

Local modulus 3: 000000000 has weight 2.

Local modulus 67: 012314012 has weight 132, 012201120 has weight 528, 012120201 has weight 528, 010222343 has weight 132, 000000000 has weight 34.

## Why an existence antichain cannot count

For moduli 5,7,11 the only states are the two complementary three-block partitions A,B and the constant partition T. The local A,B,T weights are respectively (2,2,3), (6,6,4), (10,10,6). The product distinct count is 2·2·6=24 for 35 and 2·2·10=40 for 55. But the A weight at 35 is 2·6 + 2·4 + 3·6 = 38; dropping constant states before multiplication would incorrectly yield 12. At 55 the corresponding weight is 2·10 + 2·6 + 3·10 = 62. Constants are dominated for existence, yet carry indispensable counting mass.
