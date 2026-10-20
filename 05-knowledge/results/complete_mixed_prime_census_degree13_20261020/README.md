# Complete odd mixed-prime symmetry census through degree thirteen

**Status:** The permutation-subgroup classification is proved by the orbit-image and subdirect-product argument in [COMPLETENESS.md](COMPLETENESS.md), with an explicit finite certificate for the degree-thirteen transitive image. The entire 48-row census is **FINITE-EXACT**, with independent controls described in [VERIFICATION.md](VERIFICATION.md). No degree-fourteen completeness or all-degree sign assertion is made.

The universe is odd-order subgroups H ≤ S_n whose orders have at least two distinct prime divisors, modulo ambient S_n-conjugacy, with 1 ≤ n ≤ 13. The degree, fixed vertices, and permutation embedding are retained. There are **48 degree/class pairs**: 10 through degree ten, 8 at eleven, 10 at twelve, and **20 at thirteen**. Thus “18 through eleven” is cumulative, not the count at eleven.

## Principal completion and correction

For C15 acting on 5+3+1^5, the previously omitted dimension-24 cube gives

- relative tournament orbits R = 17,824;
- full tournament types tau_T = **17,564**;
- full sign-even graph types tau_E = **17,480**;
- defect tau_T − tau_E = **84**, not the preregistered 88 (nor the alternative 64).

The exact balance is 84 = F_E + Z − F_T = 282 + 62 − 260. The 62 escaping relative graph orbits lie over **60** odd graph types: 58 of multiplicity one and two of multiplicity two. The globally even multiplicity histogram is 1:17201, 2:277, 3:1, 4:1; the tournament histogram is 1:17304, 2:260. These are full-isomorphism counts, not normalizer counts.

For F21 acting on 7+1^6, the other omitted cube, the completed row is **3,040 = 3,040**, with no fusion or sign escape. The independent colored-tournament Burnside calculations give the two source counts 2 b_5(4)=17824 and b_6(2)=3040.

## Exact failure boundary

The complete mixed-prime list has exactly six failing degree/class pairs:

| Degree | Permutation embedding | tau_T | tau_E | Defect |
|---:|---|---:|---:|---:|
| 10 | C15 on 5+3+1² | 32 | 31 | 1 |
| 11 | C15 on 5+3+1³ | 172 | 168 | 4 |
| 12 | C15 on 5+3+1⁴ | 1376 | 1360 | 16 |
| 13 | C15 on 5+3+1⁵ | 17564 | 17480 | 84 |
| 13 | synchronized C15 on 5+3+3+1² | 768 | 762 | 6 |
| 13 | C5×C3² on 5+3+3+1² | 256 | 250 | 6 |

The inherited all-odd-p-group equality theorem then excludes further failures among odd subgroups in these degrees. The two double-triple deficits were already proved in the predecessor packet; this packet independently reproduces their numerical rows, not a new proof of those two deficits.

## Definitions and complete table

R denotes N_(S_n)(H)-orbits on H-fixed tournaments, equivalently relatively sign-even H-fixed graph orbits. F_T and F_E are the excess of relative orbits over full types on tournaments and globally sign-even graphs. Z counts relatively sign-even orbits of globally odd graphs. W counts their distinct full types. The identity is tau_T−tau_E=F_E+Z−F_T, with Z=W+E_wit.

“Sign-even” means that every graph automorphism has even edge-reversal character. It does not mean even vertex degrees. The precise pair-mask convention and all permutation generators are in the computation sources and compressed JSON receipts.

[TABLE_DEGREE13.md](TABLE_DEGREE13.md) gives all twenty degree-thirteen rows. [complete_census_summary.json](complete_census_summary.json) contains all 48 rows, including normalizer order, fixed-cube dimension, relative counts, full counts and all defect terms.

## Provenance, not retroactive completeness

The incoming checkpoint 1e4f40f2d5f9419ed97d1c516d0db0a83c100a38 computed **44** rows: 41 in its main receipt and three in its addenda. Its displayed table borrowed two degree-twelve rows, bringing the combined displayed table to **46**, not 48. The omitted degree-thirteen rows were exactly the two completed above.

All **18** incoming degree-thirteen rows match this audit in every reported numeric field. All **30** overlaps between the incoming receipts and previously audited packets match: 26 through degree twelve, two degree-thirteen F21-family probes, and two degree-thirteen double-triple probes. In addition, this accelerator replayed all 28 lower-degree rows and matched their audited references. See [receipt_comparison.json](receipt_comparison.json) and [all_lower_comparison.json](all_lower_comparison.json).

[CORRECTION_NOTE.md](CORRECTION_NOTE.md) records the exact source locations and two separate proof-scope cautions. It does not edit or silently overwrite the incoming theorem or hypothesis.
