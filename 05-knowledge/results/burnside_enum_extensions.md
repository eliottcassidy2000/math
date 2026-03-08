# OEIS B-file Extensions via Burnside/Polya Partition Enumeration

## Summary of extensions computed by opus-2026-03-08-S48

All computed using `burnside_enum_v2.c` (C+GMP, multi-threaded, bucket accumulation).

| Sequence | Name | OEIS had | We computed | New terms | Status |
|----------|------|----------|-------------|-----------|--------|
| A000568 | Tournaments | 80 | 200 | +120 | Complete |
| A002785 | Self-comp oriented graphs | 14 | 291 | +277 | Running |
| A000171 | Self-comp graphs | 100 | 352+ | +252+ | Running |
| A000273 | Digraphs | 65 | 88 | +23 | Complete |
| A000595 | Binary relations | 51 | 89 | +38 | Complete |
| A001174 | Oriented graphs | 50 | 80 | +30 | Complete |
| A000088 | Simple graphs | 88 | 92+ | +4+ | Running |
| A000666 | Symmetric relations | 81 | 92+ | +11+ | Running |
| A051337 | Strongly connected tournaments | 50 | 150 | +100 | Complete |

## Key algorithmic features

1. **LCD scaling**: Work entirely in integers by scaling all contributions by LCD
2. **Bucket accumulation**: Group LCD/z contributions by t-value, shift at end
3. **Atomic work queue**: Load-balanced multi-threading via C11 atomics
4. **Precomputed factors**: km_fac[pi][m] = k^m * m! avoids repeated GMP calls
5. **Unified framework**: Single C file handles 8 different OEIS sequences

## Sequence-specific formulas

All use: a(n) = (1/n!) * sum_{p partition} permcount(scale*p) * base^edges(p) * extra

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
