# How often each Hamiltonian-path count occurs: complete strong inventories below 625, Moon's chains, and which oddities line up

Session opus-2026-10-21-S50.

**Owner prompt (verbatim):**

> think how these relate to tournament metagraph structure oddities and the {7,21} gap in tournament ham path counts. besides the fact that 7 and 21 never appear, also find out how many times other values do appear and look for patterns. oddities in one aspect should always line up with oddities from another lens, just through very convoluted analogy consider prior work regarding 341=11*31 and {2,3,11}, keep exploring past threads looking for inspirational ideas to merge in as you aim for tournament structure proofs and understanding

"These" refers to the converse-Gluck exceptions `C_3`, `Z_5{1,2}`, `P_7` of S49 (THM-4743).

**Scripts.** All are in `04-computation/experiments/h_spectrum_multiplicities_20261021/`, each with its `.out`, and each prints ALL CHECKS PASSED.

| script | content | runtime |
|---|---|---|
| `h_multiplicities.py` | the class census `n ≤ 9`; `N(n,h)`; the SCC sequence law; exact multiplicity polynomials for odd `h < 125` | 11 s |
| `h_patterns.py` | strong polynomials below 125, converse pairing, late primes, residues, the repunit probe, the converse-Gluck family, Paley path orbits | 11 s |
| `strong_ears.py 225`, `strong_ears.py 375` | exact ears to orders 10–12: strong inventories below 375; the extremal classes through order 12 | 27 s / 226 s |
| `strong_ears_625.py` | exact ears to order 13: strong inventories below 625; `E(13)` | 30 min |
| `h_residue_sampling.py` | the pre-registered residue test at `n = 10, 11, 12` (the crashed first run is kept) | 6 min |
| `exact_gluings.py` | every exact gluing of two small chains | 7 min |
| `gluing_failure_modes.py` | how mixed gluings fail: (R) versus (U) | 53 s |
| `oneway_proof_checks.py` | the configurations excluded in the proof of THM-4745 (5) | seconds |
| `cross_check_macmini_census.py` | class-by-class agreement with the concurrent mac-mini census, `n ≤ 9` | 12 s |
| `audit_S50/` | the independent audit's own code | — |

**Canon.**

| item | status |
|---|---|
| THM-4744 | Moon's product law (KNOWN, with an independent proof); chain recognition (PROVED); strong inventories below 625 and the extremal classes through order 13 (FINITE-EXACT) |
| THM-4745 | exact gluings (PROVED); an independent computer-assisted reproof of the Moon–Yang theorem |
| HYP-9356 | KNOWN (Moon–Yang 2007); filed open before the literature was checked |
| HYP-9355 | OPEN (equidistribution); the 2-adic pull is refuted as a persistent law |
| MISTAKE-652 | the prior-art lapse |

## 0. What was known, what was concurrent, what is new

* **Literature.**
  * J. W. Moon, *Publ. Math. Debrecen* 19 (1972) 101–104, defined the chains `S_n`. He proved `h(S_n) = ∏(1 + 2^(p_i − 1))` and listed the minimizing patterns (p. 102, as read by the audit).
  * A. H. Busch, *EJC* 13 (2006) N3, proved that Moon's upper bound is the minimum `f(m)`.
  * J. W. Moon and L. L. M. Yang, *Utilitas Math.* 72 (2007) 193–198, "Strong tournaments with the fewest Hamiltonian paths", proved that the special chains are the **only** tournaments attaining it. This is cited via El Sahili and Ghazo Hanna, arXiv:2101.00713 (*J. Graph Theory* 2023); the primary paper was not accessed.
  * This session rediscovered the chains, the formula and the characterization from its census before checking the literature. THM-1370 had already noted that Moon "gave the chain construction".
* **Canon.**
  * THM-1370: the omissions 7, 21; Moon's monotonicity; the insertion lemma.
  * THM-4094, THM-4097, THM-4102, THM-4104: the two-prime-lane reduction; strong-ear solid intervals to 80,265.
  * The spectrum note of 2026-10-20, `spectrum_basins_20261020e/spectrum.md`:
    * Chan–Pak Conjecture 6.2;
    * the class series `F_N(z) = Σ_k (1−z)^(−k−1) Σ_(d_1⋯d_k = N) ∏ P_(d_i)(z)`;
    * `P_d` for `d ≤ 23`;
    * the logarithmic core bound.
  * THM-466, THM-478, THM-1980: the 2-adic digits of `H`.
  * S15 (`camion_busch_gaps_polyhedra_collatz_20261002.md`): `E(7) = 1` with scores `(1,1,2,3,4,5,5)`.
  * S33 and S34: the Paley path-orbit counts; 1729 at `q = 11` already typed NUMEROLOGY.
* **Concurrent lanes on the same prompt, all on 2026-10-21.**
  * **Session24.**
    * Fixed-`H` multiplicities: eventual polynomials of degree `r(H)`; self-converse and merged counts; an inventory complete for `H ≤ 73`.
    * The 7-adic degree defect `r(7·3^a) = a − 1`.
    * The tower path counts.
    * The 341 Cayley tournaments, and `{2, 3, 11}` recovered with no mechanism.
    * Metagraph intersections.
  * **mac-mini fanopsp.**
    * Camion pseudoprimes; the cyclotomic dictionary `7 = Φ_3(2)`, `21 = Φ_3(4)`, `341 = Φ_5(4)`; 1729 Carmichael; `P_19` fails Korselt.
    * HYP-9357: `P_q` is H-maximal iff `q ∈ {3, 7, 11}`.
    * Its own census for `n ≤ 9`. `cross_check_macmini_census.py` shows it agrees with this session's census class by class in `(H, strong, self-converse, |Aut|)`.
  * **Session25.** Every parity-admissible residue modulo every `M` is attained by a strong tournament.
* **New here.**
  * Complete strong inventories for every count below **625**, orders up to 12, against the concurrent 73. The audit reproduced them below 375 with independent code.
  * The exact multiplicity polynomials, determined for every odd `h < 625` and computed explicitly for `h < 375`.
  * An independent finite confirmation of Moon–Yang through order 13, with two predictions committed before computing.
  * A written proof of Moon's product law, and the recognition and rigidity of chains.
  * An independent computer-assisted reproof of the Moon–Yang uniqueness theorem (THM-4745 (5), (6)).
  * The residue data (HYP-9355), the late-seven facts, and the alignment table (§5).

## 1. How often each value occurs

* **The law (inherited).** `H` is the product over strong components, so `N(n,h)` is a finite sum of binomials `C(n − t + k, k)` over ordered `k`-tuples of nontrivial strong classes, of total size `t` and product `h`. The only input is the strong inventory `P_d`.
* **Strong inventories below 625 (THM-4744 (3)).**
  * By Busch, `f(13) = 625` and `f` is nondecreasing, so every strong tournament with `H < 625` has order at most 12.
  * Exact ears from the census give 8,769, 3,536 and 658 strong classes with `H ≤ 625` at orders 10, 11, 12.
  * **Every odd `d ∈ [3, 623]` other than 7 and 21 is the count of a strong tournament.**
  * With the solid intervals of THM-4097, THM-4102 and THM-4104 (85–2,881, 249–14,649, 429–80,265), every odd `d ≤ 80,265` other than 7 and 21 is the count of a strong tournament. This is assembled from canon.
* **Selected strong polynomials** (the full list below 375 is in `strong_ears_375.out`):

  | `d` | `P_d` | | `d` | `P_d` |
  |---|---|---|---|---|
  | 25 | `2z^6 + z^7` | | 63 | `2z^8` |
  | 27 | `2z^6 + 3z^7` | | 75 | `7z^7 + 2z^8 + 3z^9` |
  | 35 | `2z^7` | | 125 | `z^7 + 10z^8 + 2z^9 + z^10` |
  | 45 | `2z^6 + 5z^7 + 5z^8` | | 135 | `8z^7 + 14z^8 + 14z^9 + 10z^10` |
  | 49 | `7z^7 + z^8` | | 147 | `z^7 + 10z^8 + 10z^9 + 2z^10` |

* **Exact multiplicity polynomials** (`h_multiplicities.py`; valid for every `n` past the stated onset):
  * `N(n,3) = n − 2`;
  * `N(n,9) = (n−4)(n−1)/2`;
  * `N(n,15) = (n−4)(n−3)`;
  * `N(n,27) = n(n−5)(n−4)/6`;
  * `N(n,33) = n² − 5n − 1`;
  * `N(n,45) = (n³ − 12n² + 43n − 38)/2` (`n ≥ 7`);
  * `N(n,49) = 8n − 49`;
  * `N(n,63) = 2(n − 7)`;
  * `N(n,75) = n(n² − 13n + 44)/2` (`n ≥ 8`).

  The first version of the script gave quadratic laws for 45 and 75. It summed only cores of total size at most 9, and the audit caught this (MISTAKE-652).
* **Values beyond the census:**

  | `h` | `N(10,h)` | `N(11,h)` | `N(12,h)` |
  |---|---|---|---|
  | 45 | 96 | 157 | 239 |
  | 63 | 6 | 8 | 10 |
  | 75 | 70 | 121 | 192 |
  | 81 | 106 | 170 | 261 |
  | 105 | 59 | 85 | 115 |
  | 125 | 43 | 67 | 98 |
  | 135 | 146 | 278 | 490 |
  | 189 | 102 | 148 | 198 |
  | 223 | 82 | 117 | 152 |

## 2. The extremal strong tournaments of Busch's theorem

The answer is known: they are Moon's chains (Moon 1972, Moon–Yang 2007). This session reconstructs and confirms it.

* **Interval chains** `IC(p_1, …, p_k)`, Moon's `S_n`, are the transitive order `x_0 → … → x_N` with one back arc from the end to the start of each block. Consecutive blocks share an end vertex.
* **The product law (Moon).** `H(IC(p)) = ∏(2^(p_i − 1) + 1)`. THM-4744 (1) writes out a full proof:
  * glue the first block's end to the rest's start;
  * `x_0` can be entered only from the shared vertex;
  * the interior of block 2 can be entered only from below.
* **Recognition (THM-4744 (2)).**
  * All chains of order `m` have the scores `(1, 1, 2, …, m−2, m−2)` and `m − 2` cyclic triples.
  * The spine order is determined, so chains are rigid and distinct compositions are distinct classes.
  * The converse reverses the composition.
* **Optimal parts.** These are Moon's patterns. For `p ≥ 5`, `2^(p−1) + 1 > 3(2^(p−3) + 1)`, so splitting off a 2 helps. Three 2s lose to two 3s (27 > 25), and two 2s tie one 4.
* **Finite confirmation (THM-4744 (4)).**
  * The extremal classes are exactly the optimal chains for `m ≤ 13`: `E(m) = 1, 1, 2, 2, 1, 5, 3, 1, 9, 4, 1`, with self-converse counts the palindromes.
  * `E(12) = 4` and `E(13) = 1` were predicted in `fad6a40581` before they were computed. Git cannot certify the order of the order-12 run, because the 375 cap already existed in that commit; `strong_ears_625.py` did not.
* **An independent reproof of Moon–Yang (THM-4745).** Busch's Lemma 1 splits an extremal tournament into two extremal chains glued exactly at one vertex. THM-4745 proves:
  * (1) the gluing is exact iff (R) every path restricts and (U) merges are unique;
  * (2) a one-way exact gluing of chains is the end-to-start concatenation;
  * (5) for chains with parts 2, 3, 4, every exact gluing is one-way.

  The proof of (5) uses:
  * a restriction lemma, since exactness passes to the blocks containing `v`;
  * single-path constraints;
  * an explicit witness excluding junction configurations;
  * a double-merge witness excluding interior positions;
  * one finite check: blocks of 3–5 vertices at end positions, in `exact_gluings.py`.
* **The second floor** (the second-lowest strong count) at orders 5–11 is 11, 17, 27, 49, 81, 135, 243.
  * At orders 7, 9, 11 it equals `3^((m−1)/2)`, shared by 3, 5 and 8 classes, among them the all-triangle chain.
  * It is not a law: at order 13, `IC(2,2,2,3,3)` has 675 < 729, and order 15 behaves the same way (audit).
  * At order 8 the second floor is `49 = 7²`, with one class of order 8. The forbidden factor 7 returns as a square, consistent with the concurrent `r(7²) = 1`.

## 3. Patterns in the multiplicities

* **Converse pairing.**
  * `N(n,h) ≡ SC(n,h) (mod 2)` (implicit in Session24's merged counts), and every self-converse class of even order has `H ≡ 1 (mod 4)` (THM-791).
  * So a uniquely attained value belongs to a self-converse class, and is `≡ 1 (mod 4)` at even `n`.
  * At `n = 7` the uniquely attained values are 1, 61, 85, 87, 125, 127, 147, 155, 157, 171, 175, 189. The last three are the regular 7-tournaments, with `(H, |Aut|) = (171, 3), (175, 7), (189, 21)`.
* **The late seven (FINITE-EXACT).**
  * No tournament with `n ≤ 6` has `7 | H`, and at `n = 7` every tournament with `7 | H` is strong.
  * The multiples of 7 attained at `n = 7` are 35 (×2), 49 (×7), 77 (×10), 91 (×15), 105 (×4), 133 (×6), and 147, 175, 189 (×1 each).
  * The all-order structural form is the concurrent 7-adic degree defect.
  * The holes at `n = 7` in `(45, 150)` are 63, 119, 107 and 149. Two are multiples of 7 and two are primes. No product can lie in that range at `n = 7`, so these are holes of the strong order-7 spectrum, and no mechanism is claimed. A first draft "explained" them by product borrowing; the audit refuted that.
* **Residues (HYP-9355).**
  * The labelled fraction with `p | H` approaches `1/p`: within 1.64% at `n = 9` and within 0.8% at `n = 10–12`.
  * `H ≡ 1 (mod 8)` is over-represented through `n = 10` (1.179, 1.064, 1.015, 1.0030), but not at `n = 11`. There, the audit's pooled replication over 10^7 samples has residue 5 ahead (+2.6 s.e.), and residue 1 at −0.0004 ± 0.0005.
  * The monotone 2-adic "pull" of the first draft is refuted. Equidistribution stays open.
* **The repunit probe (unexplained).**
  * At `n = 8, 9` the binary repunits 31, 63, 127 and `85 = (4^4 − 1)/3` are the rarest values in their ±10% windows. 33 and 129 are the commonest there, and 65 is at the 83rd percentile.
  * A first draft blamed factorization. That is wrong: at `n = 8` the primes 131, 137 and 139 have 24, 26 and 30 classes against 10 for 127, and 29 is commoner than 31.
  * Every `2^k − 1` (`k ≥ 4`) is `≡ 15 (mod 16)` and every `2^k + 1` is `≡ 1`, and the small-`n` residue bias points the right way. But it is far too weak to account for the effect.
  * So the effect is real at these orders and **unexplained**. The cyclotomic reading is the mac-mini lane.
* **The converse-Gluck family and the maxima.**
  * `C_3` is the maximum at order 3, `Z_5{1,2}` ties at order 5, and `P_7` is the unique maximum at order 7.
  * At order 9, `C_3[C_3]` has 3159, the only class with that value, against the maximum 3357.
  * `H(Z_5{1,2}[C_3]) = 197,030,475` and `H(P_7[C_3]) = 123,156,134,945,235` (Session24's path-cover table, reproduced by DP) are 2.47 and 2.53 times the mean. No maximality claim is made.
* **Paley path orbits (NUMEROLOGY; S33/S34 and mac-mini already have them).**
  * `H(P_q)/|Aut(P_q)|` is 1, 9, 1729, 6857869865, 62293308207033 and 42745655163116955307561 for `q = 3, 7, 11, 19, 23, 31`.
  * For `q = 27`, with Frobenius (`|Aut| = 1053`), it is 397365570454381227, `≡ 3 (mod 8)`; normalized by `|AΓ_+| = 351` it is 1192096711363143681, `≡ 1`.
  * Only `q = 11` gives a pseudoprime (1729).

## 4. The pre-registered tests

| prediction (committed `fad6a40581`) | result |
|---|---|
| HYP-9356: `E(12) = 4`, the four chains, none self-converse | passed (`strong_ears_375.out`) |
| HYP-9356: no strong 12-tournament below 375; strong spectrum below 375 misses only 1, 7, 21 | passed |
| HYP-9356: `E(13) = 1`, `IC(3,3,3,3)` | passed (`strong_ears_625.out`) |
| HYP-9355 P1: `|p·freq(p | H) − 1| < 0.02` | passed (largest 0.79%) |
| HYP-9355 P2: residue 1 most frequent mod 8 | not confirmed at `n = 11` (audit replication, `audit_S50/a7*`) |
| HYP-9355 P3: `1 < 4·freq(1 mod 8) < 1.0153` at `n = 10` | passed (1.0030 ± 0.0009) |
| HYP-9355 P4: the excess decreases | not monotone |

* The order-12 and order-13 results confirm a known theorem. They are not discoveries.
* The seeds of the sampler were not pre-registered.

## 5. Which oddities line up

| oddity | where it lines up | mechanism | type |
|---|---|---|---|
| `H ∉ {7, 21}` | multiples of 7 absent for `n ≤ 6`; degree loss 2 in `N(n, 7·3^a)` (Session24); the second floor 49 at order 8 | multiplicativity over strong components plus the strong floor | PROVED (THM-1370, Session24) + FINITE-EXACT |
| four small-atom lists | path-rigid `{C_3}` (THM-2477); pure-blue `{C_3, T_5}` (THM-2454, with `T_5 = Z_5{1,2}`); signed-exceptional seeds `{C_3, T_5, P_7}`; Busch blocks `{C_3, R_4, NT_5}` | `C_3` is the only strong tournament in all four classes. Session24: exceptional ∩ path-rigid = `{C_3}` and exceptional ∩ pure-blue = `{C_3, T_5}`. Here Busch-extremal ∩ each of the three = `{C_3}`, using `H(A[C_3]) ≥ 3^|A| H(A)` for the towers | PROVED |
| `C_3`, `T_5`, `P_7` maximize `H` | converse-Gluck seeds; `tc = 1, 3, 9` (THM-2454) | none; both lists break at the next term | FINITE-EXACT coincidence |
| maximum versus minimum | order 5: `T_5` (symmetric) against `NT_5`, `IC(2,2)` (rigid); order 7: `P_7` against `IC(3,3)` | opposite ends of the symmetry scale | ANALOGY |
| Paley primes 3, 7, 11 | exceptions `P_3`, `P_7` (THM-4743); H-maximal Paley `q ∈ {3, 7, 11}` (mac-mini HYP-9357, conjectural beyond the tested `q`); `tc = 3^((q−3)/2)` breaks at 11 | not known | ANALOGY |
| 7, 21 versus `P_7` | the orders of `P_7` and of `Aut P_7` | none | NUMEROLOGY |
| 341 = 11·31, {2, 3, 11} | trunk, golden clock, Cayley tournaments on `Z_341` (Session24); Camion pseudoprimes (mac-mini); 1729 | none beyond the stated maps | NUMEROLOGY |
| repunit rarity | 31, 63, 127, 85 rare at `n = 8, 9` | unexplained (§3) | OPEN |

The principle holds where a mechanism exists: the gap at 7 organizes its multiples through multiplicativity. The other rows are coincidences of small numbers, typed so. The repunit rarity is a genuine unexplained oddity.

## 6. Next steps

1. The repunit rarity: what makes `2^k − 1` and `(4^k − 1)/3` rare at `n = 8, 9`? Test at `n = 10` by exact counts below 625.
2. HYP-9355: equidistribution, the `3 | H` excess, and the sign change of the mod-8 bias between `n = 10` and 11.
3. Read Moon–Yang (2007) to compare their proof with THM-4745 (5).
4. A general proof of Step 3 of THM-4745 (5), for blocks of any size, would extend the one-way lemma to all chains.

## Audit

The independent audit of 2026-10-21 (`audit_S50/`, own code) confirmed:
* the proofs of THM-4744 (1), (2) and THM-4745 (1), (2), (4);
* every finite-exact computation: the inventories below 375, `E(m)` through 12, all `P_d`, and the gluing census.

It found:
* the prior art (Moon 1972, Moon–Yang 2007);
* the truncation bug in the multiplicity polynomials;
* the false second-floor law;
* the failed explanations of the order-7 holes and of the repunits;
* the non-replication of P2;
* the contradicted monotone 2-adic statement;
* several credits and wording errors.

All are corrected above (MISTAKE-652). THM-4745 (5), (6) were added after the audit and have had a further audit pass.
