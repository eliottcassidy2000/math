# Spine 14: the two-apex census (FINITE-EXACT)

Session klein-2026-11-16-wheels2. Tool `scripts/two_apex_slice.c` (the census of
`chanpak_span3_growth_two_apex_seam_20261115/scripts/two_apex_small.c` with a slice parameter), run as
`two_apex_slice 14 3 8388608 8 k` for `k = 0..7` (eight slices by `a mod 8`, about one hour each in parallel), merged by
`scripts/merge_slices.py` (validated on spine 11: receipt `validate_slices_spine11.out`, every `V` line identical to the
earlier census). Data `data/tapexs_14.txt.gz` (format `V h count minpop orient`, `orient` = 1 for `A → B`, 2 for `B → A`,
3 for both); receipts `tapexs_14_slices.out` (header, line count and SHA-256 of each slice, and of the merge),
`two_apex_analysis_14.normal.out` / `.optimized.out` (byte-identical) from `scripts/spine14_analysis.py`.

## The family

Spine = the transitive tournament on `0 < 1 < … < 13`; apex `A` beats the spine vertices in a set `S_A` (any of the
`2^14` masks), apex `B` beats the spine vertices in a set `S_B` with `|S_B| ≤ 3` (470 masks), and the arc between the
apexes is taken in both orientations: `2 · 16,384 · 470 = 15,400,960` tapes. `minpop` is the least `|S_B|` realising the
value and `orient` the set of apex orientations realising it.

## Results

| spine | tapes | every odd `N < 2^k` but 7, 21 | first miss above 21 | source |
|---|---|---|---|---|
| 11 | `2^23` (all `S_B`) | `k = 16` | 77,533 | HYP-9803 (a) |
| 12 | `|S_B| ≤ 3` | `k = 16` | 117,345 | HYP-9803 (a) |
| 13 | 6,193,152 | `k = 17` (misses below `2^18`: 7, 21, 256,521) | 256,521 | HYP-9803 (a) |
| **14** | **15,400,960** | **`k = 19`** (misses below `2^20`: 7, 21 and 137 values from 729,421 on) | **729,421** | this packet |

- **Completeness to `2^19`.** Spine 14 attains **every odd `N < 524,288` except 7 and 21** (receipt): the pre-registered
  P-S1 (every odd `N < 2^18`) holds with a full extra octave. Below `2^20` the misses are 7, 21 and 137 values, the first
  729,421 = 7 · 11 · 9473, then 798,023, 816,547, 831,415, 841,509, 843,641, …; both neighbours 729,419 and 729,423 are
  attained (9 and 11 tapes, overlap 2). Coverage of the odd numbers by dyadic octave: `[2^18, 2^19)` 100 %, `[2^19, 2^20)`
  99.95 %, `[2^20, 2^21)` 85.6 %, `[2^21, 2^22)` 15.6 %; 1,138,644 distinct values in all, maximum `H = 5,889,173`.
- **The first miss grows faster than doubling.** 77,533 → 117,345 → 256,521 → 729,421 has ratios 1.51, 2.19, 2.84, so the
  pre-registered P-S2 (first miss in `[500,000, 560,000]`, extrapolating the last doubling) is **refuted**: the first miss
  is 1.39 · `2^19`. If the ratio keeps growing, spine 15 should reach every odd number below `2^21` and spine `m`
  roughly below `2^(m + 5)` or further; the one-apex family alone (`minpop = 0`) supplies only 6,245 values below `2^19`.
- **Overlap and orientation.** Below `2^18` the minimal overlap distribution is 6,245 / 83,818 / 40,937 / 70 for
  `|S_B| = 0, 1, 2, 3`: only 70 values (0.053 %, all ≥ 164,333) need overlap 3, and every value is attained in both
  orientations, so `A → B` suffices — the pre-registered P-S3 holds. Below `2^19` the distribution is
  6,245 / 115,454 / 136,209 / 4,234, and orientation starts to matter: 18 values are attained only with `B → A` (the first
  388,061) and 18 only with `A → B` (the first 400,199); the rest are attained both ways. Over the whole census the
  minimal overlaps are 6,245 / 148,165 / 482,482 / 501,752, so overlap 3 carries the top of the range.
- **What it means for HYP-9803 (a).** The conjecture (every odd `N ≠ 7, 21` is a two-apex count with overlap ≤ 3) is now
  FINITE-EXACT below `2^19` with spine 14 alone, and with spines ≤ 14 the only misses below `2^20` are the 137 values from
  729,421 on, all above `1.39 · 2^19`. Spine 15 (about 32 million tapes, four hours sliced) is the next test; the
  explicit digital formula for the count at `|S_B| ≤ 3` (HYP-9803's second test) is still not written.

## Pre-registered verdicts

P-S1 **confirmed** (and improved to `2^19`); P-S2 **refuted** (first miss 729,421, ratio 2.84 instead of 2.2);
P-S3 **confirmed** (70 values of 131,070 need overlap 3; `A → B` suffices below `2^18`, not below `2^19`).
