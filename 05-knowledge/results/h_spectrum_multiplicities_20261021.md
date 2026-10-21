# How often each Hamiltonian-path count occurs: Busch's extremal tournaments, complete strong inventories below 375, and which oddities line up

Session opus-2026-10-21-S50.

**Owner prompt (verbatim):**

> think how these relate to tournament metagraph structure oddities and the {7,21} gap in tournament ham path counts. besides the fact that 7 and 21 never appear, also find out how many times other values do appear and look for patterns. oddities in one aspect should always line up with oddities from another lens, just through very convoluted analogy consider prior work regarding 341=11*31 and {2,3,11}, keep exploring past threads looking for inspirational ideas to merge in as you aim for tournament structure proofs and understanding

"These" refers to the converse-Gluck exceptions `C_3`, `Z_5{1,2}`, `P_7` of S49 (THM-4743).

**Scripts.** All are in `04-computation/experiments/h_spectrum_multiplicities_20261021/`, each with its `.out`, and each prints ALL CHECKS PASSED.

| script | content | runtime |
|---|---|---|
| `h_multiplicities.py` | the class census `n ≤ 9`: `N(n,h)`, the SCC sequence law, multiplicity polynomials, patterns | 19 s |
| `h_patterns.py` | strong polynomials below 125, converse pairing, late primes, residues, the repunit probe, the converse-Gluck family, Paley `tc` | 11 s |
| `strong_ears.py 225`, `strong_ears.py 375` | exact ears to orders 10–12; complete strong inventories below 375; Busch-extremal classes through order 12; the chain product law | 27 s / 226 s |
| `h_residue_sampling.py` | the pre-registered residue test, `n = 10, 11, 12` | about 6 min |
| `exact_gluings.py` | every exact gluing of two chains of at most 6 vertices with `(|A|−1)(|B|−1) ≤ 16`, and of `C_3`, `R_4` with `IC(3,3)`: one-way, and a concatenation | about 7 min |
| `gluing_failure_modes.py` | how mixed gluings of `C_3`, `R_4` fail: (R)-violations versus extra merges | 53 s |
| `cross_check_macmini_census.py` | class-by-class agreement with the concurrent mac-mini census, `n ≤ 9` | 30 s |

**Canon.**
* THM-4744: PROVED, plus FINITE-EXACT through order 12.
* HYP-9355: OPEN; its predictions were committed in `fad6a40581`.
* HYP-9356: OPEN; the same commit predicted `E(12) = 4`, which then passed.
* THM-4745: PROVED. It shows that one-way exact gluings of chains are concatenations, and that HYP-9356 reduces to a one-way lemma, which is FINITE-EXACT for small chains.

## 0. Concurrent work, and what is new here

* **Already canon (credited).**
  * THM-1370: the omissions 7, 21; Moon's monotonicity of the strong floor; the insertion lemma.
  * THM-4094, THM-4097, THM-4102, THM-4104: the two-prime-lane reduction, and selected strong-ear solid intervals up to 80,265.
  * The spectrum note of 2026-10-20, `spectrum_basins_20261020e/spectrum.md`:
    * Chan–Pak Conjecture 6.2 as the open completeness question (CITED);
    * Busch's minimum (CITED);
    * the logarithmic core bound;
    * the class series `F_N(z) = Σ_k (1−z)^(−k−1) Σ_(d_1⋯d_k = N) ∏ P_(d_i)(z)`, with `P_d` for `d ≤ 23`.
  * THM-466, THM-478, THM-1980: the 2-adic digits of `H`.
* **Concurrent Session24 lane (today, before this note was drafted; `session24_research_synthesis_20261021.md`).**
  * Fixed-`H` multiplicities: eventual polynomials of exact degree `r(H)`, self-converse and merged counts, a strong inventory complete for `H ≤ 73` (census to order 8), and the **7-adic degree defect** `r(7·3^a) = a − 1`. That is the all-order structural effect of the gap.
  * The exact path counts of the exceptional towers.
  * The 341 Cayley tournaments and the recovered `{2, 3, 11}` (no mechanism).
  * Metagraph corollaries.
  * Session24 noted this session's order-9 checkpoint and left the THM-4744/4745 and HYP-9355/9356 identifiers alone.
* **Concurrent mac-mini fanopsp lane** (checkpoint `f002fc4253`, 16:50). It covers the pseudoprime and cyclotomic lens, `7 = Φ_3(2)`, `21 = Φ_3(4)`, `341 = Φ_5(4)`:
  * THM-4746: Camion pseudoprimes;
  * THM-4747: Hamiltonian-cycle congruences;
  * HYP-9357: `P_q` is H-maximal iff `q ∈ {3, 7, 11}`;
  * its own gentourng census for `n ≤ 9` of `H`, `C`, `tc`.

  `cross_check_macmini_census.py` compares that census with this one. The two agree class by class in `(H, strong, self-converse, |Aut|)` for every `n ≤ 9`, so each independently reproduces the other. The pseudoprime material is not duplicated here.
* **New here.**
  * Exact ears to orders 10–12 on top of the order-9 census, which mac-mini reproduced independently. Together they give **complete strong inventories for every count below 375**, against the concurrent 73.
  * The **Busch-extremal tournaments**:
    * the interval-chain product law (PROVED, all orders);
    * recognition and rigidity of chains (PROVED);
    * the extremal classes are exactly the optimal chains through order 12, with `E(12) = 4` pre-registered;
    * the conjectured all-order law (HYP-9356).
  * The residue distribution of `H` with its 2-adic pull (HYP-9355, with a pre-registered sampling test).
  * The finite late-seven facts.
  * Two numerology probes with negative verdicts.
  * An alignment table (§5).

## 1. How often each value occurs

* **The law (inherited).**
  * A class is its ordered sequence of strong components, and `H` is the product over components. So `N(n, h)` is a finite sum of binomials `C(n − t + k, k)` over ordered `k`-tuples of nontrivial strong classes, of total size `t` and product `h`.
  * For fixed `h`, `N(n, h)` is eventually a polynomial in `n`.
  * The only input is the strong inventory `P_d`.
* **Strong inventories below 375 (THM-4744 (3)).**
  * Busch: `f(10) = 125`, `f(11) = 225`, `f(12) = 375`, and `f` is nondecreasing.
  * So the census settles every strong class with `H < 125`. The exact ears settle `H < 375`: 1739 strong 10-classes and 411 strong 11-classes with `H ≤ 375`.
  * Every odd `d ∈ [3, 373]` other than 7 and 21 is the count of a strong tournament.
  * Selected strong polynomials:

    | `d` | `P_d` | | `d` | `P_d` |
    |---|---|---|---|---|
    | 25 | `2z^6 + z^7` | | 63 | `2z^8` |
    | 27 | `2z^6 + 3z^7` | | 75 | `7z^7 + 2z^8 + 3z^9` |
    | 35 | `2z^7` | | 125 | `z^7 + 10z^8 + 2z^9 + z^10` |
    | 45 | `2z^6 + 5z^7 + 5z^8` | | 135 | `8z^7 + 14z^8 + 14z^9 + 10z^10` |
    | 49 | `7z^7 + z^8` | | 147 | `z^7 + 10z^8 + 10z^9 + 2z^10` |

    The complete list to 373 is in `strong_ears_375.out`.
* **Exact multiplicities beyond the census**, for every odd `h < 375` and every `n`:

  | `h` | `N(10,h)` | `N(11,h)` | `N(12,h)` |
  |---|---|---|---|
  | 63 | 6 | 8 | 10 |
  | 75 | 70 | 121 | 192 |
  | 105 | 59 | 85 | 115 |
  | 125 | 43 | 67 | 98 |
  | 135 | 146 | 278 | 490 |
  | 147 | 70 | 123 | 192 |
  | 189 | 102 | 148 | 198 |
  | 195 | 126 | 238 | 396 |
  | 223 | 82 | 117 | 152 |

  These come from the inherited series. They agree with the census for every `n ≤ 9` and every `h < 375`.
* **Small polynomials.** These agree with the concurrent table where they overlap:
  * `N(n,3) = n − 2`;
  * `N(n,9) = (n−4)(n−1)/2`;
  * `N(n,27) = n(n−5)(n−4)/6`;
  * `N(n,33) = n² − 5n − 1`;
  * `N(n,45) = 6n² − 74n + 233`;
  * `N(n,49) = 8n − 49`;
  * `N(n,63) = 2(n − 7)`;
  * `N(n,75) = 4(n − 6)²`.

## 2. The extremal strong tournaments of Busch's theorem (THM-4744, HYP-9356)

Busch proves that a strong `m`-tournament has at least `f(m) = min 3^a 5^b` Hamiltonian paths, with `2a + 3b = m − 1`, and that Moon's construction attains this. He does not determine which tournaments attain it. This session does, through order 12, and conjectures the answer for all orders.

* **Interval chains.**
  * `IC(p_1, …, p_k)` is the transitive order `x_0 → … → x_N` with one back arc from the end to the start of each block.
  * Consecutive blocks share an end vertex, and the blocks are nearly transitive tournaments on `p_i + 1` vertices.
* **Product law (PROVED, all orders).** `H(IC(p)) = ∏ (2^(p_i − 1) + 1)`.
  * Glue the first block to the rest at the shared vertex `v`. Every cross arc points forward.
  * Each pair of block paths has exactly one merge.
  * Every Hamiltonian path restricts to a pair, for two reasons:
    * the first vertex `x_0` can only be entered from `v`;
    * the interior of block 2 can only be entered from below.
* **Recognition (PROVED).**
  * All chains of order `m` have the scores `(1, 1, 2, …, m−2, m−2)` and `m − 2` cyclic triples, Moon's minimum.
  * The spine order is determined by the scores and the two tie-breaking arcs, and the back arcs spell out the composition.
  * So chains are rigid, distinct compositions are distinct classes, and the converse reverses the composition.
* **Optimal parts.** Per added vertex, the block costs are `log(2^(p−1)+1)/p = 0.549, 0.536, 0.549, 0.567, …` for `p = 2, 3, 4, 5`.
  * So the optimal chains use 3s, at most two 2s, or one 4 in place of two 2s.
  * Their number is `E(3q) = q`, `E(3q+1) = 1`, `E(3q+2) = q(q+3)/2`.
* **Finite verification (FINITE-EXACT, `m = 3..12`).**
  * The extremal classes are exactly the optimal chains: `E(m) = 1, 1, 2, 2, 1, 5, 3, 1, 9, 4`.
  * The self-converse counts are `1, 1, 2, 0, 1, 1, 1, 1, 3, 0`, the palindromes.
  * The formula was read off `m ≤ 9` and confirmed at `m = 10, 11` in the same session.
  * `E(12) = 4` was committed as a prediction (`fad6a40581`) before the order-12 computation, and passed.
* **Reduction (THM-4745).**
  * Busch's Lemma 1 splits an extremal `T`, `m ≥ 6`, into two strong pieces sharing one vertex. Equality forces both pieces to be extremal, and the gluing to be exact.
  * If the cross arcs are one-way, exactness puts the shared vertex at the end of one spine and the start of the other.
    * At every other position there is a Hamiltonian path of `T` whose restriction to the first piece is not a path.
    * The constructions use the terminal strong component of `A − v`, or explicit two-path covers of `A − v`.
  * So the classification at every order follows from the one-way lemma. That lemma is exhaustive for small chains (§6).
* **The second floor.** The second-lowest strong count at orders `5, …, 11` is `11, 17, 27, 49, 81, 135, 243`.
  * At odd `m ≥ 7` it is the all-triangle chain, `3^((m−1)/2)`.
  * At `m = 8` it is `49 = 7²`, which is not a chain value. The forbidden factor 7 returns as a square one floor up, in line with the concurrent `r(7²) = 1`.

## 3. Patterns in the multiplicities

* **Converse pairing.**
  * Non-self-converse classes pair with their converses, so `N(n,h) ≡ SC(n,h) (mod 2)`. Every self-converse class of even order has `H ≡ 1 (mod 4)` (THM-791 with a fixed-point-free anti-involution).
  * **Hence a uniquely attained value always belongs to a self-converse class, and at even `n` it is `≡ 1 (mod 4)`.**
  * At `n = 7` the uniquely attained values are 1, 61, 85, 87, 125, 127, 147, 155, 157, 171, 175, 189. The last three are the three regular 7-tournaments: `|Aut| = 3, 7, 21`, with `P_7` at 189.
  * Checked for every `n ≤ 9`, and on the order-10 and order-11 inventories.
* **The late seven.**
  * No tournament with `n ≤ 6` has `7 | H`. So `7 | H(T)` forces a strong component of order at least 7 with `7 | H`.
  * At `n = 7` every tournament with `7 | H` is strong. The multiples of 7 attained there are 35 (×2), 49 (×7), 77 (×10), 91 (×15), 105 (×4), 133 (×6), 147, 175 and 189 (×1 each). Missing are 7, 21, 63, 119 and 161.
  * The odd holes of `n = 7` in `(45, 150)` are 63, 119 (multiples of 7) and 107, 149 (primes). These are exactly the values that cannot borrow a product factor at that order.
  * The all-order form of "the gap shadows its multiples" is the concurrent 7-adic degree defect. 7 is also the only small prime that is late: odd multiples of 7 are in range at order 6, but 7 first divides a strong count at order 7.
* **Residues (HYP-9355).**
  * For `p = 3, 5, 7, 11, 13` the labelled fraction with `p | H` approaches `1/p`. At `n = 9` it is within 1.6%.
  * The 2-adic distribution is pulled toward 1. At `n = 7, 8, 9` the relative mass of `H ≡ 1 (mod 8)` is 1.179, 1.064 and 1.015, and at `n = 7` mod 16 it falls from residues 1 and 9 to residues 7 and 15.
  * The pre-registered sampling test at `n = 10, 11, 12` is in §4.
* **The repunit probe (NO MAP).**
  * At `n = 8` and 9 the binary repunits 31, 63, 127 and `85 = (4^4 − 1)/3` are the rarest values in their ±10% windows, while 33, 65 and 129 are the commonest.
  * The binary digits are a red herring:
    * 31 and 127 are prime;
    * 63 needs the forbidden factor 7;
    * `85 = 5·17` needs `4 + 6 = 10` vertices to factor;
    * `2^k + 1` is divisible by 3 for odd `k` (33, 129) and by 5 for `k ≡ 2 (mod 4)` (65), and so borrows a cheap factor.
  * 255, 511, 1023, 341 and 1365 are typical.
  * Only small values are affected, through the factorization law. The cyclotomic reading `7 = Φ_3(2)`, `21 = Φ_3(4)` belongs to the mac-mini lane.
* **The converse-Gluck family and the maxima.**
  * `C_3` is the maximum at order 3, `Z_5{1,2}` ties at order 5 (two classes with 15), and `P_7` is the unique maximum at order 7.
  * At order 9, `C_3[C_3]` has 3159, below the maximum 3357 of the circulant `Z_9{1,2,3,5}`, though it is the only class with 3159.
  * The concurrent path-cover table gives `H(Z_5{1,2}[C_3]) = 197,030,475` and `H(P_7[C_3]) = 123,156,134,945,235`, reproduced here by DP. These are 2.47 and 2.53 times the mean `n!/2^(n−1)`, the same order as Paley tournaments (2.53 at `q = 19`, 2.56 at `q = 23`). No maximality claim is made at 15 or 21.
* **Paley presentation counts (NUMEROLOGY, negative).**
  * `tc(P_q) = 2H/(q(q−1))` is 1, 9, 1729, 6857869865, 62293308207033 and 42745655163116955307561 for `q = 3, 7, 11, 19, 23, 31`. For `q = 27`, with Frobenius, it is 397365570454381227.
  * Only `q = 11` gives a base-2 pseudoprime: 1729, a Carmichael number.
  * Every value is `≡ 1 (mod 8)` with the affine normalization, as HYP-9273 predicts.
  * The 341 thread ends here for this note.

## 4. The pre-registered tests

| prediction (committed `fad6a40581`) | result |
|---|---|
| HYP-9356: `E(12) = 4`, the four chains (2,3,3,3), (3,2,3,3), (3,3,2,3), (3,3,3,2), none self-converse | PASSED (`strong_ears_375.out`) |
| HYP-9356: no strong 12-tournament below 375; strong spectrum below 375 misses only 1, 7, 21 | PASSED |
| HYP-9355 P1: `|p·freq(p | H) − 1| < 0.02`, `p ≤ 13`, `n = 10, 11, 12` | PASSED (largest deviation 0.8%) |
| HYP-9355 P2: residue 1 is the most frequent mod 8 | PASSED within the stated 3-s.e. tolerance; strictly largest at `n = 10` and 12, but at `n = 11` residue 5 is ahead by 1.6 s.e. |
| HYP-9355 P3: `1 < 4·freq(1 mod 8) < 1.0153` at `n = 10` | PASSED (1.0030 ± 0.0009) |
| HYP-9355 P4 (informative): the excess decreases from `n = 10` to `n = 12` | NOT monotone (+0.0030, −0.0018, +0.0032). The pull is resolved only through `n = 10` |

* The sampler first drew 64-bit arc codes, which overflow at `n = 12`, where there are 66 arcs.
  * That run crashed after `n = 11`. Its tables are kept in `h_residue_sampling_run1_crashed_at_n12.out`, and its `n = 10, 11` verdicts were the same.
  * The committed version draws one byte per arc.
* The side observation that `3 | H` is over-represented decays: 1.016 at `n = 9`, 1.0034 ± 0.0007 at `n = 10`, 1.0010 ± 0.0010 at `n = 11`.

## 5. Which oddities line up

This answers the owner's principle that oddities in one lens should line up with oddities in another. Each row says whether a mechanism is known.

| oddity | where it lines up | mechanism | type |
|---|---|---|---|
| `H ∉ {7, 21}` | multiples of 7 absent for `n ≤ 6`; holes 63, 119, 161 at `n = 7`; degree loss 2 in `N(n, 7·3^a)`; second floor 49 at order 8 | multiplicativity over strong components plus the strong floor | PROVED (THM-1370, Session24 defect) + FINITE-EXACT |
| four small-atom lists | path-rigid strong `{C_3}` (THM-2477); pure-blue strong `{C_3, R_5}` (THM-2454); signed-exceptional seeds `{C_3, R_5, P_7}`; Busch blocks `{C_3, R_4, NT_5}` | Session24 proves exceptional ∩ path-rigid = `{C_3}` and exceptional ∩ pure-blue = `{C_3, R_5}`. Here Busch-extremal ∩ path-rigid = Busch-extremal ∩ pure-blue = Busch-extremal ∩ exceptional = `{C_3}`.
* The only strong path-rigid tournament is `C_3`, and `H(R_5) = 15 > 9`.
* For the towers, choosing one of the 3 triangle paths in each block of a lifted path gives `H(A[C_3]) ≥ 3^|A| H(A)`. So `H(X[C_3^[k]]) ≥ 3^((n−d)/2) H(X)`, and this exceeds `f(n) = 3·5^((n−3)/3)` for every `k ≥ 1` (ratio `≥ (√3/5^(1/3))^(n−3)` for `X = C_3`; 1.94 at `R_5[C_3]`; 8.8 at `P_7[C_3]`). The seeds are also above the floor: `R_5` has 15 > 9 and `P_7` has 189 > 25.

So `C_3` is the one strong tournament in all four classes. At orders 5 and 7 the maximum (`R_5`, `P_7`: symmetric, exceptional) and the minima (`IC(2,2)` and `NT_5`; `IC(3,3)`: rigid, palindromic) sit at opposite ends | PROVED |
| `C_3, C_5, P_7` maximize `H` | converse-Gluck seeds; `tc = 1, 3, 9 = 3^((p−3)/2)` (THM-2454) | none; both lists break at the next term (order 9; `tc(P_11) = 1729`) | FINITE-EXACT coincidence |
| 7 and 21 | orders of `P_7` and of `|Aut P_7|`; exceptional orders `7·3^k` | none (Session24 states this too) | NUMEROLOGY |
| 341 = 11·31 and {2, 3, 11} | Collatz trunk, golden clock, Cayley tournaments on `Z_341` (Session24); Carmichael `tc(P_11) = 1729`; Camion pseudoprimes (mac-mini THM-4746: 21 passes, 341 fails at 31) | none beyond the stated maps | NUMEROLOGY |
| Paley primes 3, 7, 11 | converse-Gluck exceptions `P_3, P_7` (THM-4743); H-maximal Paley tournaments `q ∈ {3, 7, 11}` (mac-mini HYP-9357, concurrent); `tc = 3^((q−3)/2)` breaks at 11 | not known; 11 is the first Paley order that is not exceptional, and the last H-maximal one | ANALOGY |
| uniquely attained values | self-converse classes | converse pairing | PROVED |
| the extremal and the exceptional | Busch-extremal chains are rigid, with one back arc per block; exceptions are maximally symmetric | the two ends of the symmetry scale; nothing transfers | ANALOGY |

The principle holds where a mechanism exists: the gap at 7 really does organize the multiples of 7, through multiplicativity. The other rows are coincidences of small numbers, and are typed so.

## 6. Next steps

1. **HYP-9356 for all orders: one lemma remains.**
   * THM-4745 proves:
     * Busch's merge inequality with its equality case: exact iff (R) every path restricts and (U) merges are unique;
     * that a **one-way** exact gluing of two chains is the end-to-start concatenation;
     * that the full classification follows from the one-way lemma by induction.
   * The **one-way lemma** says that an exact gluing of two extremal chains has all cross arcs in one direction.
     * It is exhaustive for chain pairs of at most 6 vertices with `(|A|−1)(|B|−1) ≤ 16`, and for `C_3`, `R_4` with the optimal 7-chain. This covers every induction step up to order 10.
     * Its proof must use both (R) and (U): for `C_3 ∪ C_3`, four of the fourteen mixed patterns violate only (U).
   * Proving it would close the classification of Busch's equality case at all orders.
2. Orders 13 and 14 by ears (prediction `E(13) = 1`, `E(14) = 14`); this needs a strong 12-inventory below 625.
3. HYP-9355: is the 3-excess (`3 | H` slightly over-represented through `n = 11`) a symmetry effect? `|Aut|` divides `H`, and order-3 automorphisms are the commonest.
