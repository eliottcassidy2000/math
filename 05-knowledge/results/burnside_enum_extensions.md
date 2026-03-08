# OEIS B-file Extensions via Burnside/Polya Partition Enumeration

## Summary of extensions computed by opus-2026-03-08-S48

### Core sequences (burnside_enum_v2.c)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000568 | Tournaments | 80 | 200+ | +120+ | Complete |
| A002785 | Self-comp oriented graphs | 14 | 300 | +286 | Complete |
| A000171 | Self-comp graphs | 100 | 370+ | +270+ | Running |
| A000273 | Digraphs | 65 | 98 | +33 | Complete |
| A000595 | Binary relations | 51 | 98 | +47 | Complete |
| A001174 | Oriented graphs | 50 | 80 | +30 | Complete |
| A000088 | Simple graphs | 88 | 95 | +7 | Complete |
| A000666 | Symmetric relations | 81 | 95 | +14 | Complete |

### k-uniform hypergraph sequences (a000665_gmp_enum.c, k_uniform_gmp_enum.c)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000665 | 3-uniform hypergraphs | 29 | 81 | +52 | Complete |
| A051240 | 4-uniform hypergraphs | 19 | 49+ | +30+ | Running |
| A051249 | 5-uniform hypergraphs | 16 | 41 | +25 | Complete |
| A309860 | 6-uniform hypergraphs | 15 | 36 | +21 | Complete |
| A309861 | 7-uniform hypergraphs | ~15 | 43 | +28 | Complete |
| A309862 | 8-uniform hypergraphs | ~15 | 32+ | +17+ | Running |
| A309863 | 9-uniform hypergraphs | ~14 | 30+ | +16+ | Running |
| A309864 | 10-uniform hypergraphs | ~14 | 30+ | +16+ | Running |
| A309865 | k-uniform triangle | 15 rows | 25 rows | +10 rows | Running |

### Derived sequences via inverse Euler transform (euler_transform.py)

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A001349 | Connected graphs | 76 | 93 | +17 | Complete |
| A003085 | Weakly connected digraphs | 64 | 95 | +31 | Complete |
| A051337 | Strongly connected tournaments | 50 | 200 | +150 | Complete |
| A086345 | Connected oriented graphs | 51 | 80 | +29 | Complete |
| A054919 | Connected binary relations | 51 | 97 | +46 | Complete |
| A054921 | Connected symmetric relations | 87 | 95 | +8 | Complete |

### New sequences (not yet in OEIS)

| Description | Terms | Derived from |
|-------------|-------|--------------|
| Connected 3-uniform hypergraphs | 80 | InvEuler(A000665) |
| Connected 4-uniform hypergraphs | 47 | InvEuler(A051240) |
| Connected 5-uniform hypergraphs | 40 | InvEuler(A051249) |

## Key algorithmic features

1. **LCD scaling**: Work entirely in integers by scaling all contributions by LCD
2. **Bucket accumulation**: Group LCD/z contributions by t-value, shift by 2^t at end
3. **Atomic work queue**: Load-balanced multi-threading via C11 atomics
4. **Precomputed factors**: km_fac[pi][m] = k^m * m! avoids repeated GMP calls
5. **Unified framework**: Single C file handles 8 different OEIS sequences
6. **General k-uniform**: New general enumerator for any k-uniform hypergraphs
7. **Burnside orbit counting**: For k≥3, uses (1/L)*sum_{t=0..L-1} [x^k] generating function

## Enumerator files

| File | Sequences | Notes |
|------|-----------|-------|
| burnside_enum_v2.c | 8 core sequences | Compressed partition, base-2 and base-3 |
| a000665_gmp_enum.c | A000665 only | Specialized k=3 with closed-form edges |
| a002785_gmp_enum.c | A002785 only | Self-comp oriented graphs |
| k_uniform_gmp_enum.c | Any k-uniform | General, uses Burnside orbit counting |
| euler_transform.py | Derived sequences | InvEuler, OGF inversion |

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

For k-uniform hypergraphs (k≥3): edges = orbits of k-subsets, computed via Burnside on cyclic group.
