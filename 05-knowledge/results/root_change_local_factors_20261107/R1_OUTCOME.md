# R1 outcome: two root changes to 2^32

Pre-registered in [PREREGISTRATION.md](PREREGISTRATION.md) (0cda244632):

> **R1.** F_<=2 has no odd miss in (10^9, 2^32]: the misses below 2^32 are exactly the 140 below 3,780,967.

**Outcome: CONFIRMED.**

## Run

The enumerator is the S14 auditor's independent segmented one,
`../root_change_carmichael_chamber_20261107/audit/census_audit.c`. It builds every word letter by letter, with no
closed form and no code shared with `rc_census.c`. It was compiled unchanged and run as

    census_audit 4294967296 2 28 4 0 1000000001 512 805306368

The arguments mean:
- X = 2^32 and K = 2;
- segments of 2^28 integers, from 805,306,368, which overlaps the earlier 10^9 census;
- 4-bit saturating counters (cap 15), and no per-level bitmaps;
- miss list from 10^9 + 1.

The run took 5681.1 s on 12 threads, shared with the audit's jobs. The raw output is [census_2p32.txt](census_2p32.txt).

## Result

- **Misses.** Union misses (M_2 = 0) in [805,306,368, 2^32]: **0**. The miss list from 1,000,000,001 is empty.
- **Minimum M_2 per window.**

  | Window | Minimum M_2 |
  | --- | --- |
  | [805,306,368, 2^30) | 14 |
  | [2^30, 2^31) | at least 15 (counter cap) |
  | [2^31, 2^32) | at least 15 (counter cap) |

- **Low values.** No odd N in the counted range has M_2 <= 12. The auditor's 10^9 census found only two values with
  M_2 <= 12 in [2^28, 10^9]: 277,438,727 (8) and 301,189,927 (12).

Combined with the S14 census (P1, 10^9) and its independent audit: **the odd values missed by two root changes below 2^32
are exactly the 140 below 3,780,967.** This is FINITE-EXACT. The cutoff conjecture (HYP-9506 (E)) remains OPEN beyond
2^32.
