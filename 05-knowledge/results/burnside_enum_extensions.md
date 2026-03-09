# OEIS B-file Extensions via Burnside/Polya Partition Enumeration

## Summary of extensions computed by opus-2026-03-08-S48

### Core sequences (burnside_enum_v2.c)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000568 | Tournaments | 77 | 201 | +124 | Complete |
| A002785 | Self-comp oriented graphs | 100 | 300 | +200 | Complete |
| A000171 | Self-comp graphs | 100 | 439+ | +339+ | Running (to 500) |
| A000273 | Digraphs | 65 | 101 | +36 | Complete |
| A000595 | Binary relations | 51 | 100 | +49 | Complete |
| A001174 | Oriented graphs | 50 | 100 | +50 | Complete |
| A000088 | Simple graphs | 88 | 101 | +13 | Complete |
| A000666 | Symmetric relations | 81 | 113 | +32 | Extending (to 130) |
| A003086 | Self-comp digraphs | 80 | 100 | +20 | Complete |
| A005639 | Self-converse oriented graphs | 50 | 80 | +30 | Complete |
| A002499 | Self-converse digraphs | 50 | 80 | +30 | Complete |
| A002854 | Two-graphs / Euler graphs | 88 | 113 | +25 | Extending (to 130) |

### k-uniform hypergraph sequences

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000665 | 3-uniform hypergraphs | 29 | 81 | +52 | Complete |
| A051240 | 4-uniform hypergraphs | 19 | 77+ | +58+ | Gaps at 66-69 |
| A051249 | 5-uniform hypergraphs | 16 | 64+ | +48+ | Gaps at 41-47 |
| A309860 | 6-uniform hypergraphs | 15 | 60 | +45 | Complete |
| A309861 | 7-uniform hypergraphs | ~15 | 43 | +28 | Complete |
| A309862 | 8-uniform hypergraphs | ~15 | 32+ | +17+ | Running |
| A309863 | 9-uniform hypergraphs | ~14 | 30+ | +16+ | Running |
| A309864 | 10-uniform hypergraphs | ~14 | 30+ | +16+ | Running |
| A309865 | k-uniform triangle | 15 rows | 25 rows | +10 rows | Running |

### Triangle sequences (new enumerators)

| Sequence | Name | OEIS had | We computed | New entries | Status |
|----------|------|----------|-------------|-------------|--------|
| A008406 | Graphs by edges | rows 1-20 (1350) | rows 1-35+ | +many | Running (to 50) |
| A052283 | Digraphs by arcs | rows 0-20 (2681) | rows 1-30 (9020) | +6340 | Complete |

### Derived sequences via inverse Euler transform (euler_transform.py)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A001349 | Connected graphs | 76 | 100 | +24 | Complete |
| A003085 | Weakly connected digraphs | 64 | 100 | +36 | Complete |
| A051337 | Strongly connected tournaments | 50 | 200 | +150 | Complete |
| A086345 | Connected oriented graphs | 51 | 100 | +49 | Complete |
| A054919 | Connected binary relations | 51 | 98 | +47 | Complete |
| A054921 | Connected symmetric relations | 87 | 100 | +13 | Complete |
| A003049 | Connected Eulerian graphs | 88 | 103 | +15 | Complete |

### k-ary relation sequences (new enumerator, opus-S50)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000662 | 3-ary relations | 15 | 48 | +33 | Complete |
| A001377 | 4-ary relations | 7 | 26 | +19 | Complete |
| A051241 | 5-ary relations | 5 | 16 | +11 | Complete |
| (new) | 6-ary relations | 0 | 26 | +26 | Complete |
| (new) | 7-ary relations | 0 | 16 | +16 | Complete |
| (new) | 8-ary relations | 0 | 13 | +13 | Complete |
| (new) | 9-ary relations | 0 | 11 | +11 | Complete |
| (new) | 10-ary relations | 0 | 9 | +9 | Complete |

### Additional Burnside sequences (opus-S50)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A083670 | Antisymmetric relations | 51 | 81 | +30 | Complete |
| A101460 | Connected antisymmetric | 51 | 80 | +29 | Complete |

### Trivially derived sequences

| Sequence | Name | OEIS had | We computed | New terms | Formula |
|----------|------|----------|-------------|-----------|---------|
| A059735 | Complementary pairs tournaments | 50 | 200 | +150 | (A000568 + A002785)/2 |
| A334335 | InvEuler(A000568) | 76 | 200 | +124 | InvEuler(A000568) |
| A007869 | Graphs w/ even # edges | 50 | 100 | +50 | (A000088 + A000171)/2 |
| A054928 | Digraphs w/ even # arcs | 50 | 100 | +50 | (A000273 + A003086)/2 |
| A054934 | Oriented graphs up to arc reversal | 50 | 80 | +30 | (A001174 + A005639)/2 |
| A054960 | Graphs w/ odd # edges | 50 | 100 | +50 | (A000088 - A000171)/2 |
| A000250 | Symmetric reflexive relations | 40 | 113 | +73 | A000666/2 |
| A001173 | Half binary relations | 59 | 100 | +41 | A000595/2 |
| A047832 | Self-comp binary relations | 40 | 99 | +59 | A000171(4n+1) |

### k-multigraph sequences (new enumerators, opus-S50)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A004102 | 2-multigraphs (signed graphs) | 51 | 81 | +30 | Complete |
| A053400 | 3-multigraphs | 50 | 81 | +31 | Complete |
| A053420 | 4-multigraphs | 50 | 81 | +31 | Complete |
| A053421 | 5-multigraphs | 50 | 81 | +31 | Complete |
| A053465 | Connected 2-multigraphs | 51 | 80 | +29 | Complete |
| A053467 | Directed 2-multigraphs | 40 | 81 | +41 | Complete |
| A053468 | Directed 3-multigraphs | 40 | 81 | +41 | Complete |

### Hypergraph sequences (a000612_gmp.c)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000612 | Hypergraphs on n nodes | 13 | 31 (n=0..30) | +18 | Complete |
| A003180 | 2*A000612 | 12 | 21 | +9 | Complete |
| A055621 | A000612 differences | 13 | 20 | +7 | Complete |
| A323819 | Connected covering set-systems | 13 | 16 | +3 | Complete |

### New sequences (not yet in OEIS)

| Description | Terms | Derived from |
|-------------|-------|--------------|
| Connected 3-uniform hypergraphs | 80 | InvEuler(A000665) |
| Connected 4-uniform hypergraphs | 47 | InvEuler(A051240) |
| Connected 5-uniform hypergraphs | 40 | InvEuler(A051249) |
| Connected 6-uniform hypergraphs | 35 | InvEuler(A309860) |
| Connected 7-uniform hypergraphs | 28 | InvEuler(A309861) |
| Connected self-comp digraphs | 100 | InvEuler(A003086) |
| Connected self-converse oriented | 80 | InvEuler(A005639) |
| Connected self-converse digraphs | 80 | InvEuler(A002499) |
| Connected self-comp oriented | 299 | InvEuler(A002785) |
| Connected symmetric reflexive | 111 | InvEuler(A000250) |
| Connected 3-ary relations | 48 | InvEuler(A000662) |
| Connected 4-ary relations | 26 | InvEuler(A001377) |
| Connected 5-ary relations | 16 | InvEuler(A051241) |
| Connected graphs even edges | 100 | InvEuler(A007869) |
| Connected graphs odd edges | 100 | InvEuler(A054960) |
| Connected digraphs even arcs | 100 | InvEuler(A054928) |
| Connected half-binary-relations | 100 | InvEuler(A001173) |
| Connected comp-pairs tournaments | 200 | InvEuler(A059735) |
| Connected oriented, conn complement | 100 | 2*A086345 - A001174 |
| Connected antisym, conn complement | 80 | 2*A101460 - A083670 |
| Connected oriented up to arc reversal | 80 | InvEuler(A054934) |
| Connected digraphs up to arc reversal | 80 | InvEuler(A054933) |
| 6-ary relations on n nodes | 26 | C+GMP on [n]^6 |
| 7-ary relations on n nodes | 16 | C+GMP on [n]^7 |
| 8-ary relations on n nodes | 13 | C+GMP on [n]^8 |
| 9-ary relations on n nodes | 11 | C+GMP on [n]^9 |
| 10-ary relations on n nodes | 9 | C+GMP on [n]^10 |
| Connected hypergraphs | 17 | InvEuler(A000612) |
| Connected 3-multigraphs | 80 | InvEuler(A053400) |
| Connected 4-multigraphs | 80 | InvEuler(A053420) |
| Connected 5-multigraphs | 80 | InvEuler(A053421) |
| Connected directed 2-multigraphs | 80 | InvEuler(A053467) |
| Connected directed 3-multigraphs | 80 | InvEuler(A053468) |
| Directed 4-multigraphs | 81 | directed_multigraph_gmp base=5 |
| Directed 5-multigraphs | 81 | directed_multigraph_gmp base=6 |
| Connected directed 4-multigraphs | 80 | InvEuler(dir 4-multigraph) |
| Connected directed 5-multigraphs | 80 | InvEuler(dir 5-multigraph) |

### Newly extended OEIS sequences (opus-S50 session 2)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000612 | Hypergraphs on n nodes | 13 | 20+ | +7+ | Extending |
| A003049 | Connected Euler graphs | 88 | 110 | +22 | Complete |
| A054915 | Connected graphs, conn complement | 50 | 100 | +50 | Complete |
| A054918 | Connected digraphs, conn complement | 50 | 100 | +50 | Complete |
| A054920 | Connected binrel, conn complement | 50 | 100 | +50 | Complete |
| A054922 | Connected symrel, conn complement | 50 | 100 | +50 | Complete |
| A054933 | Digraphs up to arc reversal | 50 | 80 | +30 | Complete |

## Total impact summary

- **40+ OEIS sequences extended** with new b-file terms
- **25+ potentially new sequences** (connected variants and k-ary relations not yet in OEIS)
- **~3000+ new individual terms** across all sequences
- **Unified enumerator** handles 13 sequences in a single C file (burnside_enum_v2.c)
- **General k-ary relation enumerator** for arbitrary k (Python + C)
- **Hypergraph enumerator** via subset-orbit Burnside formula (Python)

## Key algorithmic features

1. **LCD scaling**: Work entirely in integers by scaling all contributions by LCD
2. **Bucket accumulation**: Group LCD/z contributions by t-value, shift by 2^t at end
3. **Atomic work queue**: Load-balanced multi-threading via C11 atomics
4. **Precomputed factors**: km_fac[pi][m] = k^m * m! avoids repeated GMP calls
5. **Unified framework**: Single C file handles 12 different OEIS sequences
6. **General k-uniform**: New general enumerator for any k-uniform hypergraphs
7. **Burnside orbit counting**: For k>=3, uses (1/L)*sum_{t=0..L-1} [x^k] generating function
8. **Divisor-signature Mobius**: For k=4,5 hypergraphs, avoids L iteration (64x-130x speedup)
9. **Pair orbit GF**: For triangle sequences (A008406, A052283), tracks individual orbit sizes

## Enumerator files

| File | Sequences | Notes |
|------|-----------|-------|
| burnside_enum_v2.c | 12 sequences (A000568/273/595/88/666/1174/2785/171/3086/5639/2499/2854) | Compressed partition, base-2 and base-3 |
| a000665_gmp_enum.c | A000665 only | Specialized k=3 with closed-form edges |
| a002785_gmp_enum.c | A002785 only | Self-comp oriented graphs |
| k_uniform_gmp_enum.c | Any k-uniform | General, uses Burnside orbit counting |
| k_uniform_fast_enum.c | Any k-uniform (fast) | Divisor-signature Mobius approach |
| a051240_gmp_enum.c | A051240 (k=4, fast) | Closed-form c_4, 64x speedup |
| a051249_gmp_enum.c | A051249 (k=5, fast) | Closed-form c_5, 17x speedup |
| a008406_gmp.c | A008406 triangle | Graphs by edges, pair orbit GF |
| a052283_gmp.c | A052283 triangle | Digraphs by arcs, directed pair orbit GF |
| euler_transform.py | Derived sequences | InvEuler, OGF inversion |
| derive_trivial_sequences.py | A007869, A054928 | Averages of existing sequences |
| k_ary_relations.py | A000662, A001377, A051241, k≥6 | General k-ary Burnside (Python) |
| k_ary_relations_gmp.c | Same sequences (fast) | C+GMP, orders of magnitude faster |

## Sequence-specific formulas

All core sequences use: a(n) = (1/n!) * sum_{p partition} permcount(scale*p) * base^edges(p) * extra

| Sequence | Parts of | Part filter | Scale | Base | edges formula |
|----------|----------|-------------|-------|------|---------------|
| A000568 | n | odd only | 1 | 2 | sum_{i<j} gcd + sum floor(v_i/2) |
| A000273 | n | all | 1 | 2 | 2*sum_{i<j} gcd + sum(v_i-1) |
| A000595 | n | all | 1 | 2 | 2*sum_{i<j} gcd + sum v_i |
| A000088 | n | all | 1 | 2 | sum_{i<j} gcd + sum floor(v_i/2) |
| A000666 | n | all | 1 | 2 | sum_{i<j} gcd + sum(floor(v_i/2)+1) |
| A001174 | n | all | 1 | 3 | sum_{i<j} gcd + sum floor((v_i-1)/2) |
| A002785 | n/2 | odd only | 2 | 2 | 2*sum_{i<j} gcd + sum v_i |
| A000171 | n/4 | all | 4 | 2 | 4*sum_{i<j} gcd + 2*sum v_i |
| A003086 | n/2 | all | 2 | 2 | 4*sum_{i<j} gcd + 2*sum v_i - #v |
| A005639 | n | all | 1 | 3 | gcd_even + self_even (see code) |
| A002499 | n | all | 1 | 2 | gcd + gcd_even + self_2499 |
| A002854 | n | all | 1 | 2 | gcd + self(floor(r/2)-1) + [has_odd_part] |

For k-uniform hypergraphs (k>=3): edges = orbits of k-subsets, computed via Burnside on cyclic group.

For triangle sequences:
- A008406: GF(x) = prod over pair orbits of (1 + x^{orbit_size})
- A052283: GF(x) = prod over directed pair orbits of (1 + x^{orbit_size})
