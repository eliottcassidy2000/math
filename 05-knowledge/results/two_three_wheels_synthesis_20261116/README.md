# Zsigmondy wheels, the 2–3 list doubled, and six mechanisms

Session klein-2026-11-16-two3 (2026-11-15/16), on the owner's directive: relate the Riesel wheel to `{4, 2, 5, 7, 7, 13}`
and the Leech lattice and conjecture new constructions; place the programmable-coefficient repair windows
(`Ω(X^(2/3))`, `O(√X)` bracket); take the compiled 2-and-3 problem lists of other sessions and double them with niche
ideas and questions; synthesise mechanisms and prove what they yield. Canon: THM-5254 (PROVED theorems and the exact
tables), HYP-9823 (the conjectures). Receipts byte-identical under `python3` and `python3 -O` where Python output is
deterministic (`leech_numerics`); the search scripts print timings and are replayed as single runs.

## Status

- **PROVED, new** ([WHEELS.md](WHEELS.md) §1). For `W(N)` = the primes of `2^N − 1` with one residue class modulo
  `ord_p(2)` each: (A) every dyadic level `2^i` has a prime (the Fermat numbers); (B) covering `Z/N`, `N = 2^a m`, reduces
  to covering the odd part — necessary: `Σ 1/oddpart(ord) ≥ 1` over the non-dyadic primes; sufficient: the odd parts
  cover `Z/m`; (C) for `N = 2^a q`, `q` an odd prime, `a ≤ 5`: covers iff at least `q` primes have order with odd part
  `q`; least periods 24 (`q = 3`, **because `2^6 − 1 = 63` has no primitive prime divisor**), 80 (`q = 5`), 112 (`q = 7`);
  (D) a prime serves both signs of `k·2^n ± 1` iff its order is even. The classical covering sets are the full wheels
  `W(24)`, `W(36) ∖ {109}`, `W(64)` (density exactly 1).
- **FINITE-EXACT** (WHEELS.md §2). Base 2, `N ≤ 150`: the covering periods are exactly 24, 36, 48, 60, 64, 72, 80, 84,
  96, 108, 112, 120, 128, 140, 144 (132 does not cover, by the two-arrangement argument); the nested criterion matches every decided period (P-W1). Base 3,
  `N ≤ 120`: 48, 72, 90, 96, 108, 120 (base 3 lacks modulus 2, base 2 lacks 6). No self-dual (both-sign) wheel below 150
  among the decided periods; 120 undecided after `2·10^8` nodes.
- **Typed** (WHEELS.md §4, [SYNTHESIS.md](SYNTHESIS.md) Θ3). The Riesel wheel meets the Leech/Niemeier side exactly in the
  von Staudt–Clausen primes `{2, 3, 5, 7, 13}` of `B_24` (`(p − 1) | 24 ⇒ ord_p(2) | 24`): STRUCTURAL through the divisor
  lattice of 24, nothing deeper; the literal `{4, 2, 5, 7, 7, 13}` is NUMEROLOGY; `196560 = 48·(2^12 − 1)` COINCIDENCE.
- **Audited by replay** (SYNTHESIS.md Θ5). The `Ω(X^(2/3))` disjoint three-count windows and the `O(√X)` bracket are
  `variable_reader_coverage_20261116/PROOF.md` (another session, 22:32); its checker passes here; LIMITS.md caps the
  two banks at `< 0.173·X` coverage. PROVED as scoped by its authors; Chan–Pak open.
- **The list** ([LIST.md](LIST.md)). Part A: the compiled 43 (census 20 + inventory 23). Part B: fifty further items,
  tagged classical/open/repository/question and sorted by mechanism; eight are questions posed here.
- **The mechanisms** (SYNTHESIS.md): Θ1 incommensurability, Θ2 dual primitive roots and Zsigmondy (which proves the
  wheel theorems), Θ3 the modular group `Z/2 * Z/3` (which places every 24), Θ4 "two is easy, three is hard"
  (a classifier), Θ5 additive 2 against multiplicative 3 in tournament counts (the sharp constants `2^(a+2) + 1`), Θ6 the
  dual ramification `Z[i]`/`Z[ω]`.
- **OPEN / CONJECTURED** ([HYP-9823]): W1 nested-is-without-loss; W2 least periods for `q ≥ 11` (`a(q) ≪ √q`?); W3
  self-dual wheels; W4 base-`b` analogue; plus the `[Q]` items of LIST.md.

## Files

- `scripts/zsig_wheels.py` (first survey, counts covering patterns: 48 at 24, 8352 at 36, 254016 at 48),
  `scripts/wheels2.py` (exact cover search, bases 2 and 3, Zsigmondy gaps), `scripts/wheels3.py` (exact, nested
  criterion, even-order cover), `scripts/brier_sheets.py` (two-sheet self-dual search), `scripts/exact_one.py` (single
  period, large node budget), `scripts/leech_numerics.py` (the arithmetic of §4).
- `receipts/`: `zsig_wheels_60.out`, `wheels2_base2_150.out`, `wheels2_base3_120.out`, `wheels3_base2_150.out`,
  `brier_sheets_150.out`, `exact_132.out`, `brier_96.out`, `brier_120.out`, `brier_140.out`, `leech_numerics.*.out`,
  `replay_variable_reader_verify.out`.
- `prediction_wheels.txt`, [PREREGISTERED.md](PREREGISTERED.md).
