# Six threads from the owner's directive, typed

Session klein-2026-11-15-age. Each thread: what the repository already records (file, id), the exact arithmetic
verified in this session (`scripts/clocks_check.py`, receipt `clocks_check.out`), the reasoning asked for, and a
verdict in the repository's typing (STRUCTURAL / PROVED / FINITE-EXACT / ANALOGY / COINCIDENCE / NUMEROLOGY; the
small-number rule of MISTAKE-559 applies).

## T1. The MRB constant and the logarithmic repair

**Recorded.** [`MRB.md`](../logarithmic_repairs_covers_20261116/MRB.md): with `a_n = n^(1/n) − 1` and `S_n` the
alternating tail, `a_n ≤ 2S_n ≤ a_n + (a_n − a_(n+1))` for `n ≥ 5`, hence `2nS_n = log n + O((log n)²/n)`: a logarithmic
*absolute* debit from a shrinking *relative* error. [`CORE.md`](../triangular_backbone_chanpak_20261116/CORE.md): any
induced extension that strictly increases `H` increases it by at least `2⌈log₃ H⌉` (triangle budget), and `H ≤ 3^τ`.
MRB.md compares the two scales (coefficient 1 against `2/log 3 ≈ 1.82`) and concludes that a literal MRB-type budget
cannot pay for positive induced growth at large `H`. Earlier: the Session 31 contour audit (entire perturbations preserve
branch jumps yet move the constant), the MRB pressure theorem (`O(X^(18/23))`, `O(X^(37/43))` output counts), and the
older "no logarithmic repair" line (THM-4507 for Collatz potentials, `chanpak_mixed_defects_20261115/SUPPORT_WIDTH.md`).

**Added here.** The sharp growth lemma (AUDIT.md) is the third logarithm in the picture: one added vertex multiplies `H`
by at most 39, so `log H` is a potential that rises by at most `log 39 ≈ 3.66` per vertex, and by exactly that only at the
sealed gadget (THM-5243 (2), (2′)). Combined with the triangle budget: an induced extension by `m` vertices from count `H`
to count `H′` needs `m ≥ log(H′/H)/log 39` vertices and gains at least `2⌈log₃ H⌉` at its first strict increase. Burns's
phrase for the MRB summands, "multiplicative growth, logarithmic gain" (`k^(1/k) = e^(log k/k)`,
[`MRB_LEECH_SOURCES.md`](../chanpak_mrb_leech_clock_20261114/MRB_LEECH_SOURCES.md) §1), describes gains that *decay* to
zero and alternate in sign — which is why the MRB series converges; tournament counts under extension are *monotone* with
a gain bounded by a constant that does *not* decay. The passage "almost all → every target" that the owner's Collatz
framing names is blocked on the tournament side precisely by this monotone positivity: there is no cancellation to
exploit, only the forbidden band `2⌈log₃ H⌉` and the ceiling 39.

**Verdict.** ANALOGY, with one exact new ingredient (the per-vertex log-gain ceiling `log 39` and its unique saturator)
and one exact disanalogy (decaying alternating gains versus bounded monotone gains). No theorem transfers.

## T2. The triple {7, 21, 63}, the pair {27, 63}, the chain 27 → 53 → 105

**Recorded.** THM-1370 (`{7, 21}` never attained, any order); THM-344 (63 attained at order 8 by exactly two classes);
MISTAKE-024/050 (the universal bundle is `{7, 21}`, the finite `n ≤ 7` bundle `{7, 21, 63}`);
`character_of_two_records_20261111` (63 is the only atom of the monoid `odd ∖ {7, 21}` divisible by 9, "it exists only
because 7 and 21 are missing"; elasticity 5/2 at `63·63 = 49·3⁴`); `six_tutte_collatz_20261116/COUNTS.md` (the complete
incoming boundary of `{7, 21}` under the descent map is `{13, 27, 41}`: `D(13) = 7`, `D(27) = 21`, `D(41) = 21`);
`inverse_fences_landing_20261116/INVERSE.md` (27 → 53 is one exact step; no H105 tournament contains an induced H53
tournament; one old-arc reversal is necessary and sufficient for 53 → 105; repeatable rays `1 + (h−1)·2^k` from seeds
13, 23, 41, 59, 75, 95); [`CLOCKS.md`](../logarithmic_repairs_covers_20261116/CLOCKS.md) (the ray `A_e = 13·2^e + 1`
= 27, 53, 105, 209, …; `7 | A_e` iff `3 | e`, `9 | A_e` iff `e ≡ 1 (mod 6)`, so no term is divisible by 63).

**Verified here.** `A_e` for `e = 1..12`: 27, 53, 105, 209, 417, 833, 1665, 3329, 6657, 13313, 26625, 53249; primes at
`e = 2, 8, 10`; no `A_e` with `e ≤ 1000` is divisible by 63; uncovered phases of the ray in the 24-wheel
`{2, 4, 8, 10, 16, 20}`; `27 = (4·21 − 3)/3`, `(3·27 + 3)/4 = 21`, `(13 + 1)/2 = 7`, `(41 + 1)/2 = 21`.

**Reasoning.** 27 and 63 are the two 3-adic dressings of the forbidden 21: `27` is its preimage under the ternary inverse
branch `(4d − 3)/3` and `63 = 3·21` its image under the ternary multiplier `3H + 2σ` with `σ = 0`. Both are attained
while 21 is not, so the forbidden value sits between a realised preimage and a realised image on the ternary axis. From
these two realised neighbours two rays leave: the dyadic ray `13·2^e + 1` from 27 (the even inverse `2d − 1`) and the
ternary ray `7·3^j` from 63 (ordinal sums; every `7·3^j`, `j ≥ 2`, is attained). The rays are disjoint — `63 ∤ A_e`
because `3 | e` and `e ≡ 1 (mod 6)` are incompatible — and the ternary ray's two missing initial terms are exactly
`{7, 21}`. The merge the directive asks for is therefore arithmetic, not constructive: the repository's inverse
constructions realise `27 → 53` and `53 → 105` by *different* operations (an exact one-step rewire, then one old-arc
reversal), and the all-order theorem "no H105 contains an induced H53" shows the dyadic ray cannot be walked by induced
padding at all.

**Verdict.** STRUCTURAL arithmetic (PROVED, elementary) for the dressings and the disjoint rays; the tournament content
is the cited inverse-port work, OPEN beyond it (Chan–Pak completeness, HYP-9803 (a)).

## T3. (8, 15, 17) scaled by 3, five pieces on 1…2024, and the decimal identity

**Recorded.** [`FENCE_AND_SQUARE.md`](../inverse_fences_landing_20261116/FENCE_AND_SQUARE.md): on `1, …, m−1` with
allowed sums `s, m, m+s` and `d = gcd(s, m)`, `L = m/d`, the sum graph is one path on `L−1` vertices, a second path if
`d` is even, and `⌊(d−1)/2⌋` cycles of length `2L`; for `(a, b, c) = (8, 15, 17)` scaled by `g = 3`, `s = 24² = 576`,
`m = 45² = 2025`, `d = 9`: one 224-vertex path and four 450-cycles on `1…2024`, and four seams of the new square
`46² = 2116` give one Hamiltonian path (optimal). THM-4510 (the three-square chain at 15 and the (8, 15, 17) chain on
1…224). `chanpak_collatz_mrb_fibers_20261116/README.md`: `2^u = 10^v + 3·2^w` has only the solutions `(4, 1, 1)` and
`(10, 3, 3)` (Pillai-type equation with the factor 641); `decimal_1024_extensions_20261116` and
`chanpak_decimal_rank_four_20261116/GOLDEN_DECIMAL.md` call `2^10 = 10^3 + 3·2^3` "the historical identity".

**Verified here.** The components of the `(576, 2025, 2601)`-sum graph on `1…2024`: a path on the 224 multiples of 9
(largest vertex 2016) and four 450-cycles on the residue pairs `±a (mod 9)`, `a = 1..4` (largest vertices 2024, 2023,
2022, 2021); the Pillai search reproduces `(4, 1, 1)`, `(10, 3, 3)` only; `2024 = C(24, 3) = 44·46 = 45² − 1 = 2^10 + 10^3`;
`C(n, 3) + 1` is a perfect square for `n ≤ 10^6` exactly at `n ∈ {7, 10, 24, 26, 65, 13777}`.

**Reasoning.** The "supremum equal to `3·2^3 + 10^3`" of the directive matched no quantity of the five pieces (their
sizes are 224 and 450; their largest vertices 2016 and 2021–2024). What the numbers do satisfy is a two-step decimal
chain: `24 = 3·2^3` (the scaled first leg, whose square 576 is the smallest allowed sum), `1024 = 24 + 10^3` (the unique
non-trivial solution of the Pillai equation, and the largest power of two below 2024), `2024 = 1024 + 10^3 = C(24, 3)`
(the vertex range). The five pieces themselves come from `d = g² = 9` and know nothing of 24 beyond `s = (3·8)²`.

**Verdict.** The orbit decomposition and the optimal four-seam repair are PROVED (cited). The chain `24 → 1024 → 2024`
and `2024 = C(24, 3) = 45² − 1` are COINCIDENCE: the Pillai equation has two solutions, `C(n, 3) + 1` is a square six
times below a million, and neither fact enters the component structure.

## T4. Mersenne, Brier, Riesel; 78557, 509203, {3, 5, 7, 13, 17, 241}

**Recorded.** [`CLOCKS.md`](../logarithmic_repairs_covers_20261116/CLOCKS.md): the Riesel coefficient 509203 has covering
primes `{3, 5, 7, 13, 17, 241}` = the prime divisors of `2^24 − 1 = 3²·5·7·13·17·241`, orders `(2, 4, 3, 12, 8, 24)`,
period 24; the Sierpiński coefficient 78557 uses `{3, 5, 7, 13, 19, 37, 73}`, period 36; exactly two covering patterns up
to translation, all six primes essential; **Brier obstruction for this wheel**: `ord_7(2) = 3` is odd, so 7 can never
divide both `k·2^e − 1` and `k·2^e + 1` for one `k`. No canon entry exists for Riesel/Brier numbers (the agentic sweep
found only incidental mentions and the Erdős covering reflection); Mersenne numbers have a large canon (THM-4556, 4749,
4847 and others).

**Verified here.** Both coverings hold over a full period (no uncovered `n mod 24`, resp. `mod 36`); the lcm of orders
is 24, resp. 36; `2^36 − 1 = 3³·5·7·13·19·37·73·109`, so the Sierpiński wheel omits the primitive prime 109. For
`d | 24` the primitive prime divisors of `2^d − 1` are: `d = 2: 3; 3: 7; 4: 5; 6: none; 8: 17; 12: 13; 24: 241` and
`d = 1`: none.

**Reasoning.** The Riesel wheel is exactly the set of primitive prime divisors of `2^d − 1` over the divisors `d` of 24.
Of the eight divisors two contribute nothing: `d = 1` trivially and `d = 6`, because `2^6 − 1 = 63 = 3²·7` has no
primitive prime divisor — the unique exception in Zsigmondy's theorem for base 2. So the covering set has six primes
rather than eight, and the reason is the number 63 of thread T2: the triple `{7, 21, 63} = {2^3 − 1, (2^6 − 1)/3, 2^6 − 1}`
lives inside the divisor lattice of the Zsigmondy exception, and the same 63 is why the period-24 wheel is six-pronged.
This is classical, exact, and the only place where the directive's threads T2 and T4 genuinely meet. The Brier
obstruction is a parity statement about `ord_7(2)`; a Brier number must use a wheel in which every prime has `−1` in the
subgroup generated by 2, i.e. even order for every covering prime.

**Verdict.** STRUCTURAL and classical (Zsigmondy) for the six-prong count and the 63 link; the covering facts are
PROVED (CLOCKS.md, replayed here). Nothing here constrains tournaments: the 24 is `lcm(ord_p 2)` for the chosen wheel.

## T5. The number 24: Leech lattice and Ramanujan

**Recorded.** THM-4752 (transitive tournaments as Conway's Weyl vector `(0, …, 24 | 70)` in `II_{25,1}`), THM-4758/4759
(Leech deep-hole swap theorem, umbral twining), THM-4760, THM-484 (`p² ≡ 1 mod 24`, the involution modulus behind the
Gleason orders 8 and 24), THM-486 (24 is the Pisano modulus), THM-481/487/489 (Golay, `η`-discriminant bridge), THM-4838,
HYP-9362 (`{7, 21} → {8, 24}` typed NUMEROLOGY in its numerological part), MISTAKE-559 (small-number coincidences are
NUMEROLOGY), `chanpak_mrb_leech_clock_20261114/MRB_LEECH_SOURCES.md` ("no Leech cause"), CLOCKS.md ("no Leech lattice
construction is needed to explain" the period 24).

**Reasoning.** Three different 24s: (i) the Riesel wheel's period, `lcm(2, 4, 3, 12, 8, 24)`, which is 24 because the
wheel was chosen among the primes of `2^24 − 1`; (ii) the exponent-two modulus `24 = max{n : x² ≡ 1 (mod n) for all
units x}` (THM-484), which is the arithmetic source of `η(τ) = q^(1/24)Π(1 − q^n)` and hence of the 24 in `Δ = η^24`
and of the Leech/Niemeier rank; (iii) `C(24, 3) = 2024` of thread T3. Only (ii) has a mechanism; (i) and (iii) are
choices of modulus and range. The repository has already typed every proposed bridge from `{7, 21}` or from tournament
counts to (ii) as NUMEROLOGY or COINCIDENCE, and this session found no new mechanism.

**Verdict.** NUMEROLOGY for any link between the covering wheel's 24 and the Leech/Ramanujan 24, in line with the
MISTAKE-559 policy; the classical facts (THM-484, Zsigmondy) are the only structural content.

## T6. The age-structured potential, read through the directive

The directive's first clause asked for an age-structured potential for `Z ≤ 38H`. AUDIT.md shows the lemma is proved
by a deletion injection whose data are *ages*: which of the last three sites the new vertex beats, which of the last six
the beaten head itself beats, and where the previous tail sits. CAPS.md turns the same injection into a per-profile cap
`x_s ≤ κ(s)·H` with `κ ≤ 30`, and shows (pre-registered, confirmed) that no history-refined linear cone can replace it.
The "logarithmic repair" of T1 reappears here as the fact that the potential `log H` moves by at most `log 39` per
vertex and by exactly `log 39` only at the sealed gadget, where the two beaten vertices have ages 1 and 2 and themselves
fed on three sealing sinks: age structure is what the maximiser is made of.

**Verdict.** The potential as a *proof device* for the upper bound: WITHDRAWN (superseded by the deletion proof). As a
*theorem about the state*: PROVED (CAPS.md §1–2). As a *cone certificate*: REFUTED for every history length tested.
