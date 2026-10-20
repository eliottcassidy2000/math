# The 2–3 axis: every theme pairing 2 and 3 in the repository, k-valued tournaments beyond binary information, and a theory of numbers, fixed points and relations

Session mac-mini-2026-10-19-twothree (2026-10-19). Owner prompt: merge all themes relating 2 and 3 (repo and online),
continue theory-building on the structure underlying numbers, fixed points and relations; tournament enumeration is binary
relational information in a group — see what happens with numbers greater than 2. **Status: THM-4720 PROVED + FINITE-EXACT
(new sequences); HYP-9339 synthesis (typed rows, four conjectures, one of them proved). Pre-audit.**

**Sibling notes of 2026-10-19:** [Mersenne horizon](mersenne_horizon_20261019_clusters_orphans_reanchoring.md), [Farey / Mandelbrot](farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md), [Sturmian circles](duality_20261019_multiplication_addition_golden_circle_gilbreath.md), [Apollonius](apollonius_enumerative_20261019_galois_parity_cube_ford_descartes_windows.md).

## 0. Inventory of 2–3 themes in the repository (counts of files mentioning each)

`2-adic` 575, `3-adic` 236, `mod 6` 148, `log_2 3` 216, `Catalan` 159, `Gersonides` 51, `Pillai` 38, `Mahler` 112,
`Furstenberg` 14, `Stormer` 3, `Paley` 440, `cyclic triple` 49, `Redei` 222, `oriented graph` 43, `ternary` 616,
`binary` 1364, `fixed point` 619. The ones used here: THM-4484 (free cycles = Gersonides), THM-4640 (Catalan carry),
THM-4716/4717/4719 (3-smooth ladders, Sturmian circles, windows), THM-4591 (Pythagorean cycles), THM-4633 (critical
slope), THM-4556/4554/4520 (the two p-adic clocks and the level operator), THM-4522 (the `x2x3` lonely spectrum),
THM-4521 (dictionary of sixes), THM-4581 6(b) (`Z[1/6] ⋊ ⟨2,3⟩`), THM-4600/4714 (anchors and cusps), THM-001/1980/4533
(Redei and its 2-adic edge), THM-283/214/307 (Burnside formulas), the S48 Burnside enumerator.

Online: Stormer's theorem (1897; the consecutive 3-smooth pairs are 1|2, 2|3, 3|4, 8|9), Pillai's conjecture
(`|2^a − 3^b| = c` has finitely many solutions for each `c`), Furstenberg's `x2x3` theorem (1967) and conjecture (ergodic
invariant measures are Lebesgue or atomic; Rudolph 1990 for positive entropy; the zero-entropy case open, per the 2021
survey arXiv:2110.05989; a 2024 reformulation in operator algebras arXiv:2410.22701), Mahler's `3/2` problem (THM-4654),
and the 2023 "tournaments and even graphs are equinumerous" theorem (Cauchy–Frobenius counting).

## 1. Numbers greater than 2 in place of binary (THM-4720)

A tournament is a function on pairs with values in `Z/2` and the swap acting by `x ↦ 1 − x`. Replace `Z/2` by `Z/k`:
the `τ`-orbits of `Z/k` are the arc colours `{x, 1 − x}` plus, for odd `k`, the tie `1/2 mod k`. Burnside gives

    T_k(n) = Σ_(λ ⊢ n) k^(o_+(λ)) [k odd]^(o_−(λ)) / z_λ,

with `o_+ = Σ ⌊(L_i − 1)/2⌋ + Σ_(i<j) gcd(L_i, L_j)` and `o_− = #even cycles`. **The parity of `k` decides the shape:** even
`k` (no tie) sums over odd-cycle permutations only, like tournaments; odd `k` sums over all, like graphs. `k = 3` is the
oriented graphs A001174; `k = 4, …, 8` give sequences not in OEIS (`1, 2, 12, 176, 8784, 1492288, …` for 2-arc-coloured
tournaments); the rigidity defect is `2C(n,3)k^(−(2n−4))` for even `k` and `C(n,2)k^(−(n−1))` for odd `k` (matched to four
digits). The arity-3 version, 3-hypertournaments, gives `1, 1, 1, 3, 12, 1512, 6818336, 1787143624656, …` (validated by
brute force at `n = 4, 5`; not in OEIS). The Redei colour polynomial of 2-coloured tournaments has no structure beyond
parity (null result).

## 2. The synthesis (HYP-9339)

Three roles of 2 and 3 — multiplicative (`2^a 3^b`, ladders, scales, the critical slope, `x2x3`), additive (the four
Stormer/Gersonides coincidences, each a Collatz object: the two anchors, the trivial cycle, the `−5` cycle/Catalan carry),
relational (binary = tournaments, ternary = oriented graphs and hypertournaments, tie iff `2` invertible mod `k`) — with
fixed points as the common currency (Burnside, anchors, the Farey cusps, `φ² = φ + 1`, `?`'s fixed points) and `6` as the
conductor (`Z[1/6] ⋊ ⟨2,3⟩`; A073631 excludes 561 through the prime 3). Four conjectures; (2) is proved by THM-4716 +
Stormer: **balanced integer `3x+1` cycles exist exactly at the Stormer pairs `(2,1), (2,3), (4,3), (8,9)`.**

## 3. Reproduction

    python3 k_tournaments.py 14; python3 hyper_and_coloured.py 10; python3 validate_hyper.py; python3 redei_colour_polynomial.py 300

in `04-computation/experiments/two_three_20261019/`.
