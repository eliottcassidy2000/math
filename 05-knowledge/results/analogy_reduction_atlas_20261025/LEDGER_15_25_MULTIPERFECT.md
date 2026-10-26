# Fifteens, twenty-fives, and multiply perfect numbers: a typed ledger

2026-10-25, klein. Each row is typed **STRUCTURE** (a map or theorem explains the numeral),
**COINCIDENCE** (same numeral, no map), or **OPEN** (a map is conceivable but not given). External
rows are CITED from the sources in `PROBLEM_CLASSIFICATION.md`; repository rows cite files.
Numerology discipline: MISTAKE-214/222/227–235 forbid promoting a shared numeral to a bridge.

## 1. Fifteens

| object | where | type | note |
|---|---|---|---|
| the 15 supersingular primes `2,3,5,7,11,13,17,19,23,29,31,41,47,59,71` | monstrous moonshine; repository: `supersingular_fus_qseries_20261022` (see §4) | STRUCTURE | primes `p` with `X_0(p)^+` of genus zero = primes dividing the Monster order |
| Sloane's ternary exponent: `2^n` has a `0` in base 3 for all `n > 15` | OEIS A004642/A102483/A346497; checked to `5.9·10^21` (Saye 2022); BB(6,4) cryptid | STRUCTURE (the 15 is the last exceptional exponent, `2^15 = 32768 = 1122221122_3`) | the conjecture is `Pi^0_1` |
| Erdős #406 ternary machine: 15 states | bbchallenge BB(15) (Stérin–Woods 2021): halts iff some `2^n`, `n > 8`, has no `2` in base 3 | COINCIDENCE with Sloane's 15 | two different ternary conjectures; the state count is a construction size |
| `15 = C(6,2)` edges of `K_6`, duads/synthemes, points of `PG(3,2)` | Petersen family (`K_6`), THM-4529 linked `K_6` | STRUCTURE (classical) | the linked-`K_6` Gale tournament is `C3[TT2,TT2,TT2]` (LCA note M2) |
| Hydra compiled as BBf(25), BB(9); Bigfoot BB(7) | cryptid sizes | not a fifteen | listed to keep sizes straight |

The supersingular fifteen and Sloane's fifteen have no map: one is a set of primes selected by the
genus of modular curves, the other the last zeroless ternary power of two. Recorded as COINCIDENCE.

## 2. Twenty-fives

| object | where | type | note |
|---|---|---|---|
| Goldbach machine, 25 states | bbchallenge BB(25) (2016, Lean-verified); now BB(21) | STRUCTURE (a construction size) | the owner's "25-state machine isomorphic to Goldbach"; size is a coordinate, not an isomorphism |
| the two regular groups of order 25 | this session's THM-4764: `Cay(Z5^2,S)[TT2]` versus `Circ(Z25)[TT2]`, depth triples `(5,21,25)` versus `(1,5,25)` | STRUCTURE | `center_transitivity_20261021c`: every vertex-transitive 25-tournament is Cayley on one of them |
| depth `25` modulo 2 | THM-4764 (`p^e` at primes dividing `p^e - 1`) | STRUCTURE | the same 25 as the number of gaps |
| `II_{25,1}` and `24 + k = 25` | THM-4752 (tournament scores in the Lorentzian lattice), THM-4758 (`(24+k)Delta`) | STRUCTURE within the Leech lane | no map to the tournament depth 25 |
| `q <= 25` good-period claim | REFUTED (THM-762/764: first good period 27) | hostile | a 25 that was a false bound |
| `27/25` ratio | THM-4239 (strict `25B < 27P`) | STRUCTURE within ports | |
| `AP(1,25,49)` centre | Session60 | STRUCTURE within squares | `25 = 5^2` |
| `c_2 = 25/8` | HYP-9367 loop law | OPEN | a derived constant |
| Antihydra counter | the walk reaches `-6`, not 25 | not a twenty-five | listed to prevent a false match |

No two rows share a map except the two tournament rows (THM-4764 and center_transitivity), which
are the same object. The Goldbach 25 and the tournament 25 are COINCIDENCE.

## 3. Multiply perfect numbers

Sieve to `2·10^6` (`multiperfect.py`, receipt `multiperfect-normal.json`): 2-perfect `6, 28, 496,
8128`; 3-perfect `120, 672, 523776`; 4-perfect `30240, 32760`. All known multiply perfect numbers are
even; whether an odd one `> 1` exists is the BB(23,3) cryptid (kg583, 2024), and the odd perfect
case is BB(37) (Reams, 2026). Both are `Sigma^0_1` existence statements.

Numerology control: are the prime factors of multiply perfect numbers supersingular? No — `8128 =
2^6·127` already uses `127`. Recorded as COINCIDENCE-FAILED; the probe cost one sieve.

Structural remark (classical, not new): `sigma(n) = kn` with `sigma` multiplicative makes
multiply perfect numbers a divisor-sum fixed-point problem, the same shape as the repository's
divisor-profile recurrences for A000568 (THM-514/516 unit-tail kernels) only in the sense that both
are sums over divisors; no map between the two is claimed.

## 4. Repository occurrences (from the bounded survey)

See the appended section below, filled from the survey of `supersingular_fus_qseries_20261022`,
THM-4752, THM-4605 and the problem ledger.

### 4a. Fifteens in the repository (typed)

| object | file | type |
|---|---|---|
| Ogg's fifteen reproduced three ways; divisor-count identities REFUTED; "the fifteen primes imply the identities" NUMEROLOGY (MISTAKE-666) | `05-knowledge/results/supersingular_divisor_shapes_qseries_20261022.md` | STRUCTURE for the list, hostile for the identities |
| the split `{29,41,59,71} / {17,19,23,31,47} / {2,…,13}` is `M \ B`, `B \ Fi22`, `Fi22`; `{29,41,59,71}` are exactly the Monster primes with no element of order `2p` (PROVED, MISTAKE-668) | `monster_levels_square_sum_7_21_20261022.md`, `session32_square_monster_20261023/SYNTHESIS.md` | STRUCTURE |
| for `p ≡ 3 (mod 4)` the class number is the Paley half-row residue imbalance; "`378 = T_27` is the sum of the fifteen" NUMEROLOGY | `ogg_triangular_chowla_20260927.md` | STRUCTURE (Dirichlet) / COINCIDENCE |
| strong-tournament Hamiltonian-path floors `f(m) = 3, 5, 9, 15, 25, 45` for `m = 3..8`; the spectrum omits `7, 21` | THM-1370 `h-spectrum-omits-7-21-all-n` | **STRUCTURE: fifteen and twenty-five are consecutive strong floors `f(6), f(7)`**, with the omitted `21` between them |
| `15 | H(T_15)`; `H = 15` lines `3·2^j - 1`; no strong self-converse order-15 tournament has black deficit six; endpoint-four repair from order 15 | THM-585, THM-4748, `nontransitive_fifteen_20261022`, `port_residual_integral_20261025/PORT_FRONTIER.md` | STRUCTURE (each its own) |
| duads/synthemes/doily `W(3,2)`, the `15 = 9 + 6` split, Mazur's 15 torsion types, 15 pentagon tilings, 15 Apéry-like sequences, with no map among the last three | `cardinality_bridges_20261018m.md`, `the-icosahedral-fifteen-s6.md`, `fifteen-kernels` | STRUCTURE (classical) / COINCIDENCE (the triple) |
| Sloane's conjecture itself is absent; Saye 2022 is CITED for `16 <= n <= 2·3^45` in the ternary-transversality lane; zeroless powers of two in base ten are THM-4580 | `collatz_procgen_20260922_transversality_foundry.md`, THM-4580, HYP-9219 | gap: the base-3 statement was not in canon before this ledger |

### 4b. Twenty-fives in the repository (typed)

| object | file | type |
|---|---|---|
| `v_T = (s | 70) in II_{25,1}` with `v·v = -2 c3(T)`, null iff transitive; score set on the level `-25` slice of the `A_24` cusp | THM-4752 | STRUCTURE |
| `(24+k)Delta` with `24 + k = 25` for `A24, D24`; `A24`'s frame shape `25·1^{-1}` and lambency `Gamma_0(25)` | THM-4758/4759 | STRUCTURE |
| the two regular groups of order 25 and the depth triples; vertex-transitive 25-tournaments; order-25 deficit-six centres vertex-transitive but none exist | THM-4764, `center_transitivity_20261021c` | STRUCTURE |
| `q <= 25` REFUTED; `25B < 27P` sharp; `AP(1,25,49)`; odd non-Parker field orders include 25; `e = 25` projectivization loss; square-sum onsets `{15, 23, 25}` | THM-764, THM-4239, Session60, Session55, Session58a, `collatz_mod6_20260922_lead_reflection3.md` | each its own; the onset pair `{15,25}` is a Pell coincidence per that file |
| `c_2 = 25/8`, the Sun mean `X^{25/24}` | HYP-9367, THM-4028 | OPEN / STRUCTURE (exponent) |

### 4c. Perfect numbers in the repository

Multiply perfect numbers, A007691 and harmonic divisor numbers have no occurrence. Perfect numbers do:
even perfect numbers are the arc-hypercube dimensions `C(2^p, 2)` at Forcade orders with a
Mersenne-prime apex (HYP-3546, part A proved) and the triangular numbers `T_{2^p - 1}`
(`collatz_aliquot_lehmer_five_20260927.md`); quasiperfect numbers must be squares or twice squares
(`ogg_triangular_chowla` Prop. 1); `C` is perfect iff `s(Cp) = C(p+2)` (`perfect_core_twin_edges_20261019`);
aliquot cycles satisfy `prod sigma(n_i)/n_i >= 2^k` with equality iff perfect (HYP-9273); THM-4746
ties Mersenne–Wagstaff pairs to Camion pseudoprimes. The odd-perfect and odd-multiply-perfect
questions enter only through the external cryptid machines above. Nothing here is a map to the
fifteens or twenty-fives.
