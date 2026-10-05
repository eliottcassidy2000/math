# Equality in the Syracuse spectrum bound at levels 17-20 and the stitched sheets (2026-10-05): the deficits `bound - slope` fall at every level through `20` with no plateau, like `n^-1.5` to `n^-2` for `2 < q <= 4` (equality, bounded prefactor) and like `n^-0.8` to `n^-1` near `q = 2` and for `q >= 6` (equality with a slowly growing prefactor; a limit above `0.007` excluded); the centre of every odd multiplication table is `+-4^-1` and `3x-1`, `3x+1` swap the two values (the rational 2-cycle `{1/4, -1/4}` of the sign conjugacy); the A/B triples `(n, ceil(n/2), 2n+1)` carry `4y = -+1 (mod z)` and are the sheet-pair split `z = -+1 (mod 4)`; the atom operator on consecutive integers is the sibling `4I+1`; primality is blind to the sheets

**Session:** opus, `pascal-lyapunov-20261005` (fifth task), 2026-10-05.
Owner's directive: (i) decide equality in the spectrum bound of the [information-dimension bridge](collatz_information_dimension_bridge_20261005.md)
at levels `18-20`; (ii) consider how primeness versus compositeness splits the alternating members of the families
A and B, with `3x+1` and `3x-1` "stitched together when the natural numbers are conjured up"; the tournament of size
`5` under the tiling model (`6` tiles in diagonals `3, 2, 1`; the `15`-tile arrangement with an L-shaped frame
`5,0,5,0,5,4,3,2,1`; the mod-`10` multiplication table as `1/8` of it; the centres of the mod-`10`, mod-`11` and mod-`9`
tables; the three interleaved sequences and the triples `(1,1,3)`, `(2,1,5)`; the operators `p`, `q` and the atoms).
**Inherits (cited):** the bridge note (Theorem 3: `tau(q) <= min(q-1, log_3(2^q-1))`); the
[two-sheet note](collatz_two_sheet_receipts_20261005.md) (sheet pairs, the sibling automaton, cases P/M by `n mod 4`);
the [receipts four-questions note](collatz_receipts_four_questions_20261005.md) (product joins are hub coincidences,
F1); the [procgen selfie note](procgen_selfie_20261001_selfie_tournaments.md) (THM-4524, loops are the gauge; the
fixed-path tiling, `64` tilings at `n = 5`); [HYP-3244](../hypotheses/HYP-3244-tiling-half-tiling-interlocking-recursions.md)
(the fixed-path tiling model); the [Camion/Busch note](camion_busch_gaps_polyhedra_collatz_20261002.md) (`6` strong
`4`-tournaments as the octahedron; `12` classes at `n = 5`).

**Status: VERIFIED (finite numerical) for the partition functions (float64, valuations to `44`; the Syracuse law to
level `19` in full, `1.16 10^9` residues, level `20` streamed, `3.5 10^9`; mass `1` to ten digits) -- not exact
rational certificates, per the Codex audit of 2026-10-05 (MISTAKES, "atomic-prefix audit"); FINITE-EXACT for the
integer checks of sections 2-5; PROVED (Propositions 1-3, elementary); OBSERVED (the decay exponents and the verdict
on equality, which is a statement about limits read from finite levels against explicit models); typed (section 6).
Author-audited only; audit OWED. Collatz OPEN.**

**Corrections received while this note was being written (Codex, commit `162c5afe3`, MISTAKES 2026-10-05, three
entries).** (i) In the bridge note, coefficient contraction `3^k < 2^(S_k)` was identified with actual descent and the
exact valuation cylinder (modulus `2^(S_k+1)`, odd-conditional mass `2^(-S_k)`) with the coarse one (modulus
`2^(S_k)`); the carry `B` and THM-4512's exceptional representatives must be retained (`x = 1`, word `(2)`,
coefficient `3/4`, `U(1) = 1`). (ii) In the two-sheet note the sibling automaton is a sufficient common-future
detector, not a necessary one (`x = 31` fails the detector and still reaches its shadow `35`); the quarter-child
relation is reached at simultaneous time `k-3` and equality at `k-2`; an exact `-1` cell of depth `r` has `r-1`
initial valuation-one steps. (iii) The "no measure argument can cross" reading of the Vitali wall was too broad: a
computable measure with an atom at every integer makes mass-one coverage the universal statement itself (the
Codex [prefix-mass note](collatz_effective_prefix_mass_20261005.md)); what is true is only that computable points
are not random for computable atomless measures. This note uses the corrected statements; its own claims are typed
accordingly.

---

## 0. What is new, in one screen

1. **Equality, decided as far as finite levels can (section 1).** The partition functions were computed exactly to level `19` (`3^19 = 1.16 10^9` residues) and streamed at level `20` (`3.5 10^9`); the deficits `d_n(q) = bound(q) - slope_n(q)` fall at every level for every `q > 1` and nowhere plateau. For `2 < q <= 4` they decay faster than `1/n` (log-log exponents `1.6-2.2`; `d_n n^1.5` falls from `2.19` to `1.94` at `q = 2.5`; a model with a positive limit is excluded), so `tau(q) = log_3(2^q - 1)` with a bounded prefactor: **equality**. At `q = 2`, `2.1`, `1.9` and `6` the deficit is exactly `c/(n + c')` (`1/d_n` linear in `n` to `0.0-0.6%`): `d_n = 0.93/(n + 2.25)` at `q = 2`, the logarithmic correction of the linearly growing `E[rho_n^2]`; `d_n = 0.63/n` at `q = 6`: equality with a polynomial prefactor. At `q = 8` the curve `1/d_n` is concave (`0.7%` residual) and an endpoint model allows a limit up to `0.0074`: undecided, and the candidate mechanism for a strict inequality at large `q` is the coherent landing of several parity words on the residues near `-1` (the atom at `-1` carries `1.46` words' worth of mass), which the word-energy bound counts once.
2. **The stitched sheets are exact and general (section 2).** For every odd `m = 2k+1` the central block of the
   multiplication table is `[[k^2, k(k+1)], [k(k+1), (k+1)^2]] = [[d, -d], [-d, d]]` with `d = 4^-1 (mod m)`, and
   `3d - 1 = -d`, `3(-d) + 1 = d`: the `3x-1` map sends the diagonal value to the anti-diagonal one and `3x+1`
   sends it back. The owner's `25 = 36 = 3, 30 = 8` (mod `11`) and `16 = 25 = 7, 20 = 2` (mod `9`), with
   `3*3 - 1 = 8` and `3*2 + 1 = 7`, are the cases `k = 5, 4`. The reason is the rational 2-cycle `1/4 -> -1/4 -> 1/4`
   of the pair `(3x-1, 3x+1)`, i.e. the sign conjugacy `U_-(x) = -U_+(-x)`: the two sheets are one map on `Z`, and
   restricting to the positive integers ("conjuring the naturals") is what splits it in two; the centre of the
   table remembers both halves. Even moduli `2k` (`k` odd) have centre `k, k+1, k+1, k-1` (`5, 6, 6, 4` mod `10`)
   and no stitching.
3. **The 15-tile triangle (section 3).** The triangle `0 <= x <= y <= k` is a fundamental domain of the block
   `{0, +-1, ..., +-k}^2` of the mod-`m` table under the group of order `8` (swap, two sign flips), with
   `T_(k+1)` tiles and values transforming by `xy -> +-xy`: `15` tiles for mod `9` (the whole table) and for mod
   `10` without its centre cross (`19` cells, `6` more orbits: `21` in all), `21` for mod `11`. The owner's "1/8 of
   the mod-10 table" is exact for the `81` cells off the cross.
4. **The families (section 4).** `B_k = (2k-1, k, 4k-1)`, `A_k = (2k, k, 4k+1)` satisfy `4y = +1 (mod z)` on B and
   `4y = -1 (mod z)` on A, so `y = +-4^-1 (mod z)` and `3y - 1 = -y` (B), `3y + 1 = -y` (A): each triple carries
   the centre stitching of its own modulus `z` (the owner's `3(3)-1 = 8` is `B_3`, `3(2)+1 = 7` is `A_2`). The
   alternation of A and B is the sheet-pair split of the sibling automaton: `v_2(3z+1) = 1` exactly on B,
   `v_2(3z-1) = 1` exactly on A. The operator `p` sends `(A_k, B_k)` and `(B_k, A_k)` to the `x`-pair of the atom
   `y z = k(4k-1)`; the operator `q` on consecutive integers is `(4I+1, 4I+4)`, whose odd member is the Collatz
   sibling of `I` (`U(4I+1) = U(I)`).
5. **Primeness does not split the dynamics (section 4).** Primes below `10^6` sit `39,322` on B and `39,175` on A
   (Chebyshev's bias, `147`), and the valuation profiles of `3z+-1` on the primes of each family are the geometric
   law of all members to four digits: the sheets do not see primality, as the receipts notes found from the other
   side (no multiplicative fusion, product joins are hub coincidences).
6. **The tournament tiles (section 5).** On `5` vertices: `1,024` labelled tournaments, `12` isomorphism classes,
   `6` strong; a fixed Hamiltonian path leaves `6` off-path arcs in diagonals of lengths `3, 2, 1` (`2^6 = 64`
   tilings, the repo's fixed-path tiling); the selfie version has `10` arcs and `5` loops, `15 = T_5` tiles, and the
   frame around the `6` inner tiles is the `4` path arcs plus the `5` loops, `9` tiles: the L. The labels
   `5,0,5,0,5` and `4,3,2,1` are the owner's convention and are not decoded here.

---

## 1. Equality in the spectrum bound (VERIFIED, finite numerical; `collatz_syracuse_multifractal_deep_20261005.py`)

The bridge note proved `tau(q) <= min(q-1, log_3(2^q-1))` and conjectured equality for `q > 2`. The decisive
quantity is the deficit `d_n(q) = bound(q) - s_n(q)` with `s_n(q) = -log(Z_n(q)/Z_(n-1)(q))/log 3` the per-level
slope of the partition function `Z_n(q) = sum_c mu_n(c)^q`: `d_n -> 0` is equality of the exponent (a
subexponential prefactor), `d_n -> d_inf > 0` is strict inequality. The law was computed exactly to level `19`
(`3^19 = 1.16 10^9` residues, float64, a parallel kernel; mass `1` to ten digits) and level `20` streamed
(`3.5 10^9` residues, never stored); valuations to `44`.

| `q` | bound | `d_14` | `d_16` | `d_18` | `d_19` | `d_20` | `d_n n^1.5` at `13 / 20` | log-log exponent `13..20` | `1/d_n` against `n` |
|---|---|---|---|---|---|---|---|---|---|
| 1.5 | `0.5000` | `0.0145` | `0.0120` | `0.0101` | `0.0093` | `0.0085` | `0.75 / 0.76` | `1.46` | convex |
| 1.9 | `0.9000` | `0.0454` | `0.0399` | `0.0355` | `0.0335` | `0.0318` | `2.28 / 2.85` | `0.99` | `0.64/(n + 0.1)`, resid `0.4%` |
| 2.0 | `1.0000` | `0.0570` | `0.0508` | `0.0458` | `0.0436` | `0.0416` | `2.85 / 3.73` | `0.88` | `0.93/(n + 2.25)`, resid `0.0%` |
| 2.1 | `1.0832` | `0.0536` | `0.0467` | `0.0412` | `0.0388` | `0.0367` | `2.70 / 3.28` | `1.05` | `0.70/(n - 0.9)`, resid `0.6%` |
| 2.5 | `1.4003` | `0.0417` | `0.0333` | `0.0268` | `0.0241` | `0.0217` | `2.19 / 1.94` | `1.79` | convex (faster than `1/n`) |
| 3.0 | `1.7712` | `0.0356` | `0.0272` | `0.0209` | `0.0184` | `0.0161` | `1.92 / 1.44` | `2.16` | convex (faster than `1/n`) |
| 4.0 | `2.4650` | `0.0393` | `0.0318` | `0.0262` | `0.0239` | `0.0218` | `2.06 / 1.95` | `1.62` | convex (faster than `1/n`) |
| 6.0 | `3.7712` | `0.0452` | `0.0395` | `0.0351` | `0.0332` | `0.0315` | `2.28 / 2.82` | `1.01` | `0.63/(n - 0.2)`, resid `0.1%` |
| 8.0 | `5.0439` | `0.0401` | `0.0360` | `0.0329` | `0.0316` | `0.0303` | `2.00 / 2.71` | `0.79` | concave, saturating (`d_inf <= 0.0074`) |

(`q = 0.5`: deficit `-0.0010` at level `20`, i.e. the slope sits above the bound `-0.5` and approaches it from above; `q = 1`: zero exactly.)

**Verdict.** The exponent bound is attained -- `tau(q) = min(q - 1, log_3(2^q - 1))` -- for every `q` in `[1.5, 6]` at the resolution of levels `13..20`: the deficits are forced to zero (faster than `1/n` for `2 < q <= 4`, exactly `c/(n + c')` at `q = 1.9, 2, 2.1, 6`). The transition at `q = 2` is second order in the sense that `tau` is continuous with a kink (`tau'(2^-) = 1`, `tau'(2^+) = 4 ln 2/(3 ln 3) = 0.841`), and the `1/(n + 2.25)` deficit there is the logarithmic prefactor of `E[rho_n^2] ~ 0.31 n` (P7). At `q = 8` the data leave a limit between `0` and `0.0074` open. If the inequality is strict for large `q`, the mechanism is the hierarchy of residues near `-1` on which several parity words land coherently (`mu_n(-1) 2^n = 1.4622`), which the word-energy bound counts as single words; `D_inf = log_3 2` is unaffected. OBSERVED, not proved: this is the behaviour of finite levels read against two explicit models.

The maximal atom stays at `c = -1` with `mu_n(-1) 2^n -> 1.4622` (local dimension `log_3 2` with the `1/n`
correction `-log_3 1.46/n`); the atoms at `-5` and `-17` have `alpha_n` oscillating in `0.86-0.89` and
`0.92-0.95` (below their cylinder bounds `0.946`, `0.992`: coherent landings), and `1, 5, 7` descend slowly toward
`1` (`1.045, 1.043, 1.084` at level `18`).

## 2. Proposition 1: the centre of the odd multiplication table (PROVED)

**Proposition 1.** Let `m = 2k+1 >= 3` be odd and `d = 4^-1 (mod m)`. Then `k^2 = (k+1)^2 = d` and `k(k+1) = -d`
(mod `m`), and `3d - 1 = -d`, `3(-d) + 1 = d` (mod `m`). *Proof.* `4k^2 = (2k+1)(2k-1) + 1 = 1`, so `k^2 = d`;
`(k+1)^2 - k^2 = 2k + 1 = 0`; `4k(k+1) = (2k+1)^2 - 1 = -1`, so `k(k+1) = -d`; `3d - 1 + d = 4d - 1 = 0`;
`3(-d) + 1 - d = 1 - 4d = 0`. □ (Checked for all odd `m <= 401`.) In `Q` this is the 2-cycle `1/4 -> 3/4 - 1 = -1/4
-> -3/4 + 1 = 1/4` of the stitched pair, equivalently the sign conjugacy `U_-(x) = -U_+(-x)` evaluated at the
fixed point `x = 1/4` of `x -> -(3x - 1)`. The table's centre is the reduction of `+-1/4` modulo `m`.

For even `m = 2k` with `k` odd, `k^2 = k`, `(k +- 1)^2 = k + 1`, `(k-1)(k+1) = k - 1` (mod `m`) (checked, `m <=
398`): mod `10` the centre is `5` with `6, 6` on the diagonal and `4` on the anti-diagonal, and neither `3*4+1 = 3`
nor `3*6-1 = 7` closes a cycle; the stitching is a property of odd moduli, where `4` is invertible.

## 3. Proposition 2: the 15-tile triangle as a fundamental domain (PROVED + FINITE-EXACT)

Let `G` be the group of order `8` generated by `(x, y) -> (y, x)`, `(x, y) -> (-x, y)`, `(x, y) -> (x, -y)` acting
on the cells of the mod-`m` multiplication table; the table value transforms by `xy -> +-xy`.

**Proposition 2.** On the block `{0, +-1, ..., +-k}^2` every `G`-orbit has exactly one representative with
`0 <= x <= y <= k`; the block has `T_(k+1) = (k+1)(k+2)/2` orbits. For odd `m = 2k+1` the block is the whole table;
for `m = 2k` with `k` odd the cells with `x = k` or `y = k` (the centre cross, `4k - 1` cells) are not in the block.
*Proof.* Signs can be chosen nonnegative and the swap orders the pair; distinct representatives are in distinct
orbits because `|x|, |y|` are invariants. □ Counts: mod `9`: `15` orbits in all; mod `10`: `15` on the block and
`21` in all; mod `11`: `21`. The owner's `15`-tile triangle is the mod-`9` table and the mod-`10` table off its
cross, "flipped across axes by simple rules" being the sign rule.

## 4. Proposition 3: the triples, the operators and the sheets (PROVED + FINITE-EXACT)

Write `F_n = (x, y, z) = (n, ceil(n/2), 2n + 1)`: `B_k = F_(2k-1) = (2k-1, k, 4k-1)`, `A_k = F_(2k) = (2k, k, 4k+1)`.

**Proposition 3.** (a) `4y = 1 (mod z)` on B and `4y = -1 (mod z)` on A; hence `y = 4^-1` on B, `y = -4^-1` on A,
and `3y - 1 = -y (mod z)` on B, `3y + 1 = -y (mod z)` on A. (b) `x(B_k) + x(A_k) = z(B_k)`. (c) With `p(I, J) =
(x(I) + x(J), z(I) x(J))`: `p(B_k, A_k) = (4k-1, 2N)` and `p(A_k, B_k) = (4k-1, 2N-1)` with `N = k(4k-1) = y(B_k)
z(B_k)`, i.e. the second components are the `x`-pair `{2N-1, 2N}` of the atom `N`. (d) With `q(I, J) = (2(I+J) - 1,
4J)` on consecutive integers `J = I + 1`: `q = (4I + 1, 4I + 4)`, and for odd `I` the odd output is the sibling,
`U(4I+1) = U(I)`. (e) `v_2(3z+1) = 1` and `v_2(3z-1) >= 2` on B; `v_2(3z-1) = 1` and `v_2(3z+1) >= 2` on A.
*Proof.* (a) `4 ceil(n/2) = 2n + 2` or `2n`. (b)-(d) arithmetic. (e) `z = -1` or `+1 (mod 4)`. □ (All checked to
`n <= 2000`, `k <= 1000`, `I <= 10^4`.)

**Reading.** The triple `F_n` packages its own modulus `z` with the centre value `y = +-4^-1` of Proposition 1,
and the alternation of A and B is the alternation of the sign in `+-4^-1`, which is the alternation of the sheet
that takes valuation one (case P of the sibling automaton on B, case M on A). The atom `N <-> {{2N-1, 2N}, {4N-1,
4N+1}}` is the pair `(B_N, A_N)` without its middle entry, and `q` on its first pair gives `(8N - 3, 8N) =
(S(2N-1), 4 * 2N)`: the sibling ray of the odd member and the quadruple of the even one. These are exact
bookkeeping identities of the `4u + 1` ray; they contain no descent.

**Primes.** Below `10^6` there are `39,322` primes in B (`4k - 1`) and `39,175` in A (`4k + 1`): Chebyshev's bias.
The distribution of `v_2(3p+1)` over the primes `p` of A is `(0, 0.4991, 0.2503, 0.1248, 0.0630, 0.0627)` for
`a = 1..5, >= 6` against `(0, 0.5, 0.25, 0.125, 0.0625, 0.0625)` for all members of A, and likewise for
`v_2(3p-1)` on B: primality does not see which sheet is small or how small. Together with the receipts notes
(no synchronous fusion, product joins through hubs only) the split of A and B by primeness is a statement about
the primes, not about the dynamics.

## 5. The 5-tournament tiles (FINITE-EXACT)

Exhaustive: `1,024` labelled tournaments on `5` vertices, `12` isomorphism classes, `6` strong (the repo's `6`
strong at `n = 4` is a different count: the octahedral shell of score vectors). With the Hamiltonian path
`1 -> 2 -> 3 -> 4 -> 5` fixed, the `6` off-path pairs `(i, j)`, `j - i >= 2`, fall into diagonals of lengths `3, 2,
1` by `j - i`, and their orientations give the `2^6 = 64` tilings of the procgen selfie note. The selfie tournament
adds `5` optional loops: `15 = T_5` tiles, and the `9` tiles around the inner `6` are the `4` path arcs and the `5`
loops, the L. THM-4524 says the loops are a gauge (switching), so the loop arm carries no isomorphism information
while the path arm is the forced order; the owner's labels `5,0,5,0,5` and `4,3,2,1` were not decoded and are not
used.

## 6. What is structural and what is not (typed)

* Exact and general: the centre of every odd multiplication table is `+-4^-1`, swapped by the two sheets
  (Proposition 1); the triangle fundamental domain (Proposition 2); the triples carry `y = +-4^-1 (mod z)` and
  the mod-`4` sheet split (Proposition 3); `q` produces the sibling ray.
* The "stitching when the naturals are conjured": the sign conjugacy `U_-(x) = -U_+(-x)` makes `3x-1` on
  positives the restriction of `3x+1` to negatives; the natural numbers are one half of one map on `Z`, and the
  rational point `+-1/4`, visible as the centre of every odd table, is the place where the halves meet in a
  2-cycle. This is the same conjugacy that made the first trunk multiplier the sheet commutator in the two-sheet
  note.
* Not structural: primeness versus compositeness in A and B (independent of the valuation profile to four digits;
  the receipts theorems forbid multiplicative structure in the dynamics); the specific labels of the L frame; the
  identification of the `6` tiles with isomorphism classes (`12` classes, `6` strong tournaments).

## 7. Reproduction

```bash
python 04-computation/experiments/collatz_syracuse_multifractal_deep_20261005.py 19 44        # 621 s, 15.6 GB resident (level 19 in memory, level 20 streamed)
python 04-computation/experiments/collatz_syracuse_multifractal_deep_fit_20261005.py
python 04-computation/experiments/collatz_stitched_sheets_families_20261005.py                 # 2 s
```
Outputs `.out`/`.json` beside the scripts.

## 8. Verdicts

| claim | status |
|---|---|
| equality `tau(q) = min(q-1, log_3(2^q-1))` forced by the deficits for `1.5 <= q <= 6` (levels `13..20`); undecided at `q = 8` (`d_inf <= 0.0074`) | OBSERVED on VERIFIED finite numerical data (levels `13-20`) |
| centre of every odd multiplication table is `+-4^-1`, swapped by `3x-1` and `3x+1` (Proposition 1) | PROVED |
| the `T_(k+1)` triangle is a fundamental domain of the table block under the order-8 group (Proposition 2) | PROVED + FINITE-EXACT |
| triples: `4y = -+1 (mod z)`, the sheet split, `p` to the atom `yz`, `q` to the sibling (Proposition 3) | PROVED + FINITE-EXACT |
| primeness does not split the sheets (valuation profiles equal to four digits; Chebyshev bias `147` at `10^6`) | FINITE-EXACT |
| `5`-tournaments: `12` classes, `6` strong, `6` off-path tiles in diagonals `3,2,1`, `15` selfie tiles, `9`-tile L | FINITE-EXACT |
| Collatz | OPEN |
