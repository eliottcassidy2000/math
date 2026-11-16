# Spine 14, Zsigmondy wheels II, factor-1 shields, and the fifth root

Session klein-2026-11-16-wheels2 (2026-11-15/16), on the owner's directive: run the spine-14 two-apex census; finish the
wheel follow-up (least periods for `q ≥ 11`, the self-dual question at 120, the nested criterion); search for a factor-1
shield that still has a +2 move and for a coordinated multi-bit edit or a different seam; read 24/36/64 against
24/80/112; merge the Zsigmondy work with the Riesel and Erdős covering work; relate the `Ω(N^(1/5))` address choices to the
golden-ratio, fence and square-packing work and `{2, 3, 11}`; relate prescribing the count modulo 9 to the 8-and-9
square-sum and Catalan material. Canon: THM-5255 (new theorems), HYP-9823 and HYP-9803 updated.

## Status

- **PROVED** ([WHEELS2.md](WHEELS2.md) §1–§3). Theorem A: `W(N)` covers `Z/N` iff `a ≥ 6` or the odd wheel covers `Z/m`
  (the fibre argument on `Z/2^a × Z/m`), which settles the nested-dyadic conjecture W1 in its corrected form. Theorem B:
  for a prime-power odd part the criterion is the Kraft density `Σ q^(−j) ≥ 1`. Theorem C: for `N = 2m` a wheel is
  self-dual iff its flipper wheel covers `Z/m` (the plus sheet is a translate of the minus sheet by `(m + 1)/2`).
  Theorem D: the parity-sheet / Brier bound `Σ_keep 1/d + 2Σ_flip 1/d + ½Σ_odd 1/d ≥ 1`; it is below 1 for every period
  `≤ 150`, so **no wheel of period ≤ 150 is self-dual** (P-W2 of the previous packet proved outright), and `W(24)`,
  `W(36)`, `W(64)` cannot host Brier numbers while the period-720 Filaseta–Finch–Kozek set scores 1.476. Erdős's 1950
  covering system is, class for class, the first pattern of the Riesel wheel.
- **Least periods** (WHEELS2.md §2; exact factorisations of `Φ_(2^i q)(2)`). `a(q) ≤ 6` for every odd prime `q` (the dyadic
  part of `W(64q)`); `a(3) = 3, a(5) = a(7) = 4, a(11) = a(13) = 5, a(17) = a(19) = 6` (periods 24, 80, 112, 352, 416, 1088,
  1216); `a(23) = 6` (period 1472: `Φ_736(2)`, 106 digits, has exactly three prime factors — factordb, verified locally, receipt `nu_23_factordb_verify.out` — so `c(5, 23) = 14 < 23`). Conjecture: `a(q) = 6` for every prime `q ≥ 17`. The classical 24, 36, 64 are the cheapest
  wheels with odd parts 3, 9, 1; the least periods per prime are 24, 80, 112, 352, 416, then dyadic.
- **Spine 14** ([SPINE14.md](SPINE14.md)): spine 14 (15,400,960 tapes, eight slices merged) attains **every odd `N < 2^19` except 7 and 21**; the first miss above 21 is 729,421 (ratio 2.84 over spine 13, so the doubling guess P-S2 is refuted while P-S1 and P-S3 are confirmed); below `2^18` only 70 values need overlap 3 and `A → B` suffices, below `2^19` 18 values need `B → A`; HYP-9803 (a) is FINITE-EXACT below `2^19`.
- **Shields and moves** ([SHIELD.md](SHIELD.md)). A factor-1 shield is a transitive block, i.e. a run of ordinary sinks; it
  keeps the +2 and +4 seed swaps and multiplies the slope by `2^s`, and by BLOCKS.md no factor-1 shield can supply a new
  radix or coordinate (PROVED there). Address-level exact moves (FINITE-EXACT, `n ≤ 22`): about 2.2 % of addresses have
  an exact +2 partner address; single-bit partners exist only at CONTACT.md's critical bit and realise +2, +4 or +8
  depending on the width; 0.13–0.34 % of addresses have a +2 partner that keeps both exterior phases; every phase cell
  reaches every residue modulo 9. A third independent coordinate needs a rank-≥4 reader (two corners): OPEN.
- **Conditional theorems** ([CONDITIONAL.md](CONDITIONAL.md)). Collatz is not the famous conjecture: under `h = 2x + 1` the
  Collatz map is `h ↦ 3h` or `(h + 1)/2`, the forbidden pair is `{h(3), h(10)} = {7, 21}`, and since 10 is the only
  predecessor of 5, three quarters of all odd numbers (373,393 of the odd `h ≤ 10^6`) have orbits through 7 or 21, where any
  pullback breaks (P-C1 refuted). The famous conjecture that does carry content is Erdős–Selfridge: it forbids every
  self-dual wheel with distinct flipper parts, and a distinct-part self-dual wheel would refute it; Hough–Nielsen and BBMST
  force a part `≥ 9` divisible by 3 and a part divisible by an odd square. The triples `(5, 12, 13)`, `(7, 24, 25)` have path
  endpoints `a²` and `b²/2` (the general `8, 9`), the double `(10, 24, 26)` splits into three pieces, and 67/69 are the
  boundary of the finite exceptional set, not a Pell pair. Five bold hypotheses B1–B5 and a Vitali/determinant note.
- **Typed threads** ([THREADS2.md](THREADS2.md)). `Ω(N^(1/5))` is a two-term balance exponent (as are Friedman's fence
  bound and Erdős–Graham's `7/11`): ANALOGY of method. The golden slope: a PROPOSAL for LIVE.md's phase crossing via the
  three-distance theorem. `{2, 3, 11}`: the odd-square-bracket theorem (PROVED, cited); 132 and 352 in the wheel theory;
  otherwise NUMEROLOGY. Count modulo 9: STRUCTURAL (2 is a primitive root modulo 9; the two lowest ternary digits read run
  lengths modulo 6) against modulo 8 (Rédei plus the one- and two-cycle packings); the square-sum endpoints 8, 9 and
  Catalan's `8, 9` are NUMEROLOGY relative to the residue prescriptions.

## Files

- `scripts/`: `conjugacy_triples.py`, `collatz_tree.py` (CONDITIONAL.md), `two_apex_slice.c` (the census tool with a slice parameter), `merge_slices.py`, `spine14_analysis.py`, `nu_factor.py`, `nu23_verify.py`
  (cyclotomic factorisations), `wheel_theorems_check.py` (Theorems A–D against the exact search, `N ≤ 150`),
  `classical_covers.py` (Erdős's system, Brier bounds), `address_moves.py` (exact address-level moves and mod-9 cells).
- `receipts/`: `tapexs_14_slices.out` (header, line count and SHA-256 of each of the eight slices and of the merge; the
  70 MB of slices themselves are not stored), `validate_slices_spine11.out`, `wheel_theorems_check_150.out`,
  `nu_11/13/17/19/23.out`, `nu_23_factordb_verify.out`,
  `classical_covers.*.out` (byte-identical under `python3` and `python3 -O`), `address_moves_12..22.out`,
  `two_apex_analysis_14.*.out`.
- `data/tapexs_14.txt.gz` (the merged spine-14 census), `prediction_wheels2.txt`, [PREREGISTERED.md](PREREGISTERED.md).
