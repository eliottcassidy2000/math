# Tournament pseudoprimes: the Camion test, Hamiltonian-cycle congruences, Paley extremality, and the 7 / 21 / 341 dictionary

Session mac-mini-2026-10-21-fanopsp. **Canon:** THM-4746, THM-4747. **Hypotheses:** HYP-9357, HYP-9358. **Independently audited** ([report](tournament_pseudoprimes_audit_20261021.md)); corrections logged as MISTAKE-651.

**Owner prompt.**
* The owner sent three pasted blocks: the seven-point sphere with the Fano analogy and arXiv:2610.26686; the converse-symmetry exceptions `C_3, C_5{1,2}, P_7`; and the `{7, 21}` Hamiltonian-path gaps with their multiplicities, `341 = 11·31` and `{2, 3, 11}`.
* The added instruction: **"think pseudoprimes"**.
* The owner's standing thesis: "oddities in one aspect should always line up with oddities from another lens, just through very convoluted analogy".

**Concurrent lanes on the same prompt (credited, not redone).**

| lane | work |
|---|---|
| opus-2026-10-21-S49 | THM-4742 (the seven-point sphere is two spheres; Fano's seventh line); THM-4743 (Gluck with the converse) |
| Session24 | HYP-9353 PROVED: the converse exceptions are exactly `X[C_3^[k]]`, `X ∈ {C_3, Z_5{1,2}, P_7}` |
| session23b | the tetrahedral Fano bridge |
| opus-2026-10-21-S50 | THM-4744/4745, HYP-9355/9356: the multiplicity laws of the H spectrum (interval chains `H = Π(2^(p_i−1)+1)`; Busch-extremal classes = optimal chains; H equidistributes mod every integer with a 2-adic pull toward 1). S50's census agrees with this session's class by class (commit 29d15e4ae1). |

* I sent S50 a lane message (`agents/opus/inbox/MSG-2395…`) and took only the pseudoprime lens.
* The link arXiv:2610.26686 is Bucić–He–Lee's proof of the Erdős–Szemerédi sunflower conjecture (§7.6), not the Gluck paper S49 cites.

**Scripts.** All are in `04-computation/experiments/tournament_pseudoprimes_20261021/`, each with an `.out`; the independent audit's scripts are in `audit/`.

## 0. Results

1. **Camion pseudoprimes (THM-4746, PROVED + FINITE-EXACT, audited).**
   * At prime order `p`, rotations force #strong labelled tournaments ≡ #circulants `= 2^((p−1)/2) ≡ (2/p) (mod p)`.
   * The strong count is quasi-periodic: `s(m+p−1) ≡ (2/p) s(m) (mod p)`.
   * The audit's Lemma B adds that `s(m+p(p−1)) ≡ (2/p) s(m) (mod p²)` for `m ≥ 3`.
   * This gives a **Korselt-type criterion** whose free residues `1, 2, 3` are K1, the missing strong K2, and the **triangle**.
   * Consequences:
     * every Carmichael number, every `2p`, every `3p`, and every `3b(3b−2)` with both factors prime is a pseudoprime;
     * this includes `4^q − 1` whenever `2^q − 1` and `(2^q+1)/3` are both prime;
     * `21` is one; `341` is not, because its residue at 31 is the Wagstaff factor 11; `3·341 = 1023 = 2^10 − 1` is one;
     * `p²` passes iff `p² | D(p)`, which below 10000 happens only for `p = 3, 13`.
   * The quasi-periodicity itself is shared with connected graphs (Mani–Stones 2016). What is tournament-specific is the circulant target and its free residues.
2. **Hamiltonian-cycle congruences (THM-4747, PROVED + FINITE-EXACT, audited).**
   * An automorphism of prime order with a fixed vertex divides `C(T)`.
   * Circulants on `Z_p`: `C ≡ (p−1)/2 (mod p·μ_S)`; Paley: `C(P_p) ≡ (p−1)/2 (mod p(p−1)/2)`.
   * Reflection-symmetric paths: `ρ_0 ≡ μ_S (mod 2μ_S)`.
   * Summing over tournaments gives Wilson's theorem (the classical rotation proof, weighted).
   * Composite orders 9, 15, 21 pass no circulant version of the naive test.
3. **Paley extremality (HYP-9357, FINITE-EXACT + CITED, audited).**
   * For prime `q`, `P_q` is the circulant H-maximum exactly at `q = 3, 7, 11`; the rotational circulant beats it at 19 and 23. This is in line with Alon's and Wormald's conjecture that the rotational tournament is the maximiser for large odd `n`.
   * `{3, 7, 11}` are exactly the `q ≡ 3 (mod 4)` whose Paley design is 2-transitive (Kantor), i.e. Galois' exceptional primes.
   * **New:** at `n = 13` a rigid regular tournament (`|Aut| = 1`) has `H = 3719893` and `C = 168361`, beating every circulant on both counts.
4. **Census `n ≤ 9` (FINITE-EXACT, lens-specific).**
   * The presentation counts `tc = H/|Aut|` miss exactly `7` and `21` below 60 at `n = 7, 8, 9`.
   * The Hamiltonian-cycle spectrum has no floor gaps, only ceiling gaps (§6).
5. **Dictionary (§7).** Typed bridges:
   * the cyclotomic/Cipolla reading of `7, 21, 341` (mechanism PROVED for the Camion split);
   * `{2, 3, 11}`;
   * the shared atoms `K1, C3` across five canon lists;
   * the Fano/tetrahedron group statement;
   * sunflowers and cap sets.

## 1. Inheritance pass

* **Closest proved mechanisms.**
  * THM-1370 (1): `H` is multiplicative over the ordered strong components (Wright's recurrence is the same decomposition).
  * LEM-003: automorphisms act freely on Hamiltonian paths.
  * THM-2477: the path-rigid strong atoms are `K1, C3`.
  * THM-4743 / Session24.
  * In the literature: Specker–Blatter (ultimate periodicity of MSO-definable labelled counts modulo every `m`) and Mani–Stones (connected graphs modulo prime powers).
* **Canonical hostile examples.**
  * `341`: the first base-2 pseudoprime, which fails the strong test.
  * The rooted count `m·s(m)`, which is not `(p−1)`-quasi-periodic (period `p(p−1)`).
  * *Corrected by the audit:* connected graphs are **not** a hostile control. They obey `cg(m+p−1) ≡ (2/p)cg(m)` for `m ≥ 2`; the pre-audit test demanded it from `m = 1` (MISTAKE-651).
* **Corrected near misses.**
  * (i) I first read the `n = 7` interior holes of the H spectrum (`63, 107, 119, 149`) as "prime-like" (realisable only by strong tournaments). That is vacuous: at `n = 7` every `H ≥ 47` comes only from strong tournaments, by vertex count. Withdrawn before writing.
  * (ii) `tc(P_11) = 1729` is a Carmichael number, but `tc(P_19)` fails Korselt at 7: NUMEROLOGY.
  * (iii) A program bug (zero classes allocated only for `p ≤ N/3`) truncated a first list run; it was fixed and rerun.
  * (iv) The audit's three major corrections:
    * the connected-graph control (above);
    * the Mersenne table (squareful `k` silently skipped, `k = 59` miscalled undecided);
    * "`D` is not periodic mod `p²`", a test begun at `m = 1` where Lemma B starts at `m = 3`.
* **Least-used sidecar.** The Hurwitz ring of binomial convolution over `F_p`, and the group algebra `K[(Z/p)]` inside it.

## 2. Camion pseudoprimes (THM-4746)

**Definitions.**
* `s(m)` counts strong labelled tournaments: A054946 `1, 0, 2, 24, 544, 22320, 1677488, …`. By Camion's theorem these are exactly the tournaments with a Hamiltonian cycle.
* `c(m)` counts circulants on `Z_m`: `2^((m−1)/2)` for odd `m`, `0` for even `m`.
* `D = s − c`. A composite `n` is a **Camion pseudoprime** if `s(n) ≡ c(n) (mod n)`.

**Why primes pass.** The rotation `x ↦ x+1` of `Z_p` acts on labelled strong tournaments. Its fixed points are the circulants, all strong. So `s(p) ≡ c(p) = 2^((p−1)/2) ≡ (2/p)` by Euler's criterion.

**Quasi-periodicity (proof in THM-4746).**
* `2^C(m,2)` is quasi-periodic mod `p` with multiplier `(2/p)`, so it is a combination of exponentials `E_g = (g^m)`, with `g` running over an additive group `G ≅ Z/p` in `F_p` or `F_(p²)`.
* Wright's recurrence is `(e − S) ∗ T = e` in the Hurwitz ring. `K[G]` is a local ring, so `e − S` is again a combination of the same exponentials.
* Hence `s(m+p−1) ≡ (2/p)s(m)`.
* The same algebra is closed under the shift, so connected graphs obey the law too, from `m = 2` (audit Lemma C; Mani–Stones). Strong semicomplete digraphs obey it with `(3/p)`. Rooted counts are the genuine negative (`camion_controls.out`, corrected).
* **Mod `p²` (Lemma B).** `s(m+p(p−1)) ≡ (2/p)s(m)` for `m ≥ 3`. The failures are exactly at `m = 1, 2` (`camion_even_sq.out`, corrected).

**Korselt-type criterion.** For `n` with squarefree odd part:

    s(n) ≡ c(n) (mod n)  ⟺  p | D(n mod (p−1)) for every odd prime p | n.

* With Lemma B it extends to cube-free odd parts (residues mod `p(p−1)`).
* It agrees with direct computation on all 1,321 composites `n ≤ 2000` with squarefree odd part, and on all 346 cube-free cases `≤ 2500` (audit).

**Zero classes.**
* `Z_p = {r : p | D(r)} ⊇ {1, 2, 3}` for `p ≥ 5`, because `D(1) = D(2) = D(3) = 0`.
* For `r ≥ 4`, `D(r) > 0`, and the extra (sporadic) classes are Poisson with mean `(p−4)/p` (HYP-9358(b)). Examples: `43: 39`; `47: 17`; `59: 21`; `103: 37, 43`; `139: 77, 103`; ….
* For odd `n`, `3` and `5` are blind primes: the odd residues mod `p−1 ∈ {2, 4}` are only `1` and `3`. For `4 | n` the prime 5 obstructs (`D(4) = 24`): `20, 40, 60` fail.
* `7` is the first prime with a nonzero odd class, `r = 5`: `D(5) = 540 = 2²·3³·5 ≢ 0 (mod 7)`.

**Families (all PROVED by the criterion).**

| family | residues | examples |
|---|---|---|
| Carmichael | all `r = 1` | 561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841 (all ≤ 20000 pass) |
| `2p` | `r_p = 2` | 6, 10, 14, 22, 26, … (1228 below 20000) |
| `3p` (`p > 3`) | `r_p = 3`, `r_3 = 1` | 15, 21, 33, 39, 51, 57, 69, … (857 below 20000); `9` also passes |
| `3b(3b−2)`, `b ≥ 5`, both prime | `3` at both large primes, `1` at 3 | 195, 399, 1023, 1443, 4623, 12099, 16383, 19599 |
| `4^q − 1` (Mersenne–Wagstaff double primes) | as above | `q = 5, 7, 13, 17, 19, 31, …`; among `2^k − 1`, `k ≤ 64`, exactly `k = 4, 6, 10, 12, 14, 26, 34, 38, 62` pass (audit) |
| sporadic | uses `Z_p \ {1,2,3}` | odd: 1239, 3765, 8401, 12019, 14317; even: 186, 530, 1898, 3290, 7202, 8282, 17066 |

**The 60 exotic squarefree members up to 20000** (not of the form `2p` or `3p`).

* **Odd (37):** 105 = 3·5·7, 195 = 3·5·13, 231 = 3·7·11, 399 = 3·7·19, 561 = 3·11·17, 595 = 5·7·17, 1023 = 3·11·31, 1105 = 5·13·17, 1239 = 3·7·59, 1443 = 3·13·37, 1729 = 7·13·19, 2465 = 5·17·29, 2821 = 7·13·31, 3003 = 3·7·11·13, 3315 = 3·5·13·17, 3765 = 3·5·251, 4465 = 5·19·47, 4623 = 3·23·67, 5187 = 3·7·13·19, 5565 = 3·5·7·53, 6601 = 7·23·41, 7395 = 3·5·17·29, 8401 = 31·271, 8463 = 3·7·13·31, 8911 = 7·19·67, 10153 = 11·13·71, 10585 = 5·29·73, 10815 = 3·5·7·103, 12019 = 7·17·101, 12099 = 3·37·109, 14317 = 103·139, 14763 = 3·7·19·37, 15841 = 7·31·73, 16383 = 3·43·127, 18915 = 3·5·13·97, 19599 = 3·47·139, 19803 = 3·7·23·41.
* **Even (23):** 30, 182, 186, 530, 870, 1122, 1406, 1898, 2210, 3290, 3458, 3782, 4930, 5402, 5642, 6734, 7202, 8282, 11102, 13202, 17066, 17390, 17822.

**Prime squares (PROVED by Lemma B).**
* `p²` passes iff `p² | D(p)`. The Camion–Wieferich primes below 10000 are `3` and `13` only, so `9` and `169` pass.
* Mod 9: `D(m+6) ≡ −D(m)` for `m ≥ 3`. So the factor 9 never obstructs, and `63 = 4³ − 1` and `4095 = 2^12 − 1` pass.

**The tournament-count test** `2^C(n,2) ≡ c(n) (mod n)` (odd `n`) is the Fermat test to the base `c(n) = 2^((n−1)/2)`.
* It contains all odd base-2 pseudoprimes, plus `85 = (4⁴−1)/3, 133, 481, 511 = 2⁹−1, 949, …`. This sequence is not in the OEIS.
* Its only universal residue is `1`, while the strong test has `1, 3`. **The universal residue sets differ exactly by the triangle class.**
* The tests also disagree for non-universal reasons (audit).

| | count test | strong (Camion) test |
|---|---|---|
| 21 = 3·7 = (4³−1)/3 | fails | **passes** (triangle class at 7) |
| 85 = 5·17 = (4⁴−1)/3 | passes | fails |
| 341 = 11·31 = (4⁵−1)/3 | **passes** (base-2 pseudoprime) | fails at 31 (residue 11) |
| 1023 = 3·11·31 = 4⁵−1 | fails | **passes** (triangle class at 11 and 31) |
| 561, 1105, 1729 | pass | pass |

**Burnside reading.** The prime-order formulas `(T(m) + (m−1)c(m))/m` and `(s(m) + (m−1)c(m))/m` count rotation classes of all, resp. strong, tournaments. The two pseudoprime notions are exactly the composite `m` where these formulas still give integers. The first is kind-pasteur's `|Tour(m)/Z_m|`.

## 3. Hamiltonian cycles under automorphisms (THM-4747)

**Statements.**
* An automorphism fixing a directed Hamiltonian cycle rotates it, so it has no fixed vertex.
* For `Cay(Z_p, S)` with `μ = |{a : aS = S}|` (odd: `M_S` misses `−1`):
  * the `(p−1)/2` arithmetic cycles form orbits of size `μ`;
  * every other cycle has a free orbit of size `pμ`;
  * hence `C ≡ (p−1)/2 (mod pμ)`.
* The signed group of order `2pμ` acts on paths, and only the reflections `x ↦ 2c−x` fix paths. Burnside gives `ρ_0 ≡ μ (mod 2μ)`.
* An involutory anti-automorphism always reverses an odd number of Hamiltonian paths (Rédei).

**Checks.** All 262 circulant classes `3 ≤ n ≤ 23` pass (the audit re-derived every value, with brute force over all permutations for every connection set `n ≤ 11`). At `n = 9, 15, 21`, no class satisfies `C ≡ #(S ∩ Z_n^*) (mod n)`.

**Paley table.**

| `q` | `H` | `C` | `tc` | `ρ_0` |
|---|---|---|---|---|
| 3 | 3 | 1 | 1 | 1 |
| 7 | 189 | 24 | 9 | 9 |
| 11 | 95095 | 5505 | 1729 | 185 |
| 19 | 1172695746915 | 34358763933 | 6857869865 | 573057 |
| 23 | 15760206976379349 | 374127973957716 | 62293308207033 | 63871533 |

## 4. Paley extremality (HYP-9357)

* **Data.** The circulant H-maximum is:
  * `P_q` at `q = 7, 11`;
  * the rotational `Z_n{1..m}` at `n = 13, 15, 17, 19, 21, 23`, beating `P_19` by 0.98% and `P_23` by 1.6%;
  * `Z_9{2,3,4,8}` (`H = 3357`, the global maximum) at `n = 9`.
* **Ratios.** `H·2^(n−1)/n!` rises from 2.40 (`P_7`) to 2.60 (rotational 23). For orientation, Adler–Alon–Ross give `e − o(1)` as a lower bound for the maximum, and Wormald's blow-ups of rotational tournaments about 2.8559.
* **`n = 13` (new).** Steepest ascent over cyclic-triangle reversals, 160 restarts, finds a regular tournament with `|Aut| = 1`, not self-converse, with `H = 3719893` and `C = 168361`. The best circulant values are `3711175` and `166848`. So the global maximum at 13 is not a circulant, and `max H(13) ≥ 3719893` (OEIS A038375 stops at 11).
  * The maximisers' symmetry orders run 3, 5/3, 21, 9, 55 for `n = 3..11`, then 1 at `n = 13`.
  * At `n = 15`, 16 restarts found nothing above the rotational circulant (a weak search).
* **Why the list `{3, 7, 11}` is not arbitrary.**
  * It is Galois' list of exceptional primes for `PSL(2,p)` acting on `p` points, restricted to `q ≡ 3 (mod 4)`.
  * It is also Kantor's list for 2-transitive Paley designs: the trivial design, the Fano plane, the biplane. The deduction was repaired after the audit, via Sylow conjugacy of regular cyclic subgroups of `PGL(m,2)`.
* **Nesting.** Converse-Gluck `{3, 7}` ⊂ H-maximal `{3, 7, 11}` = 2-transitive design.
* **Credit.** Part 3 is a special case of Alon's and Wormald's conjecture that the rotational tournament maximises for large odd `n`. The `n = 13` record shows that "large" exceeds 13.
* **Mechanism.** A HEURISTIC: design symmetry balances the 5- and 6-vertex statistics while it exists, i.e. through 11.

## 5. What the pseudoprime lens says about the 7/21 gaps (honest scope)

* **Not explained here.** The gaps are THM-1370: multiplicativity plus the strong floor `3, 5 | 9`. Nothing here gives a new proof.
* **What is shared, after the audit.** The ordered strong decomposition makes `H` and `tc = H/|Aut|` multiplicative over components. That is why the presentation-count spectrum inherits exactly the gaps `7, 21` (`analyze_census.out`: no other odd value below 60 is missing at `n = 7, 8, 9`).
  * The quasi-periodicity of the strong count is **not** a consequence of that decomposition. It is a property of the exponential family generated by `2^C(m,2)` mod `p`, and connected graphs share it.
  * The pre-audit claim of a single shared mechanism is withdrawn (MISTAKE-651).
* **Shared atoms.** The orders where every strong tournament is a circulant are `1` and `3`: K1 and C3. They recur across the canon:
  * the universal Camion residues;
  * the path-rigid atoms (`tc = 1`, THM-2477/4741);
  * the first two members of the pure-blue alphabet `{K1, C3, T5}` (THM-4730);
  * the substituted towers of every converse-Gluck exception (Session24);
  * the first two H-values `1, 3`.
  * S50's checkpoint 93356203d1 records a proved four-way intersection `{C3}` of its own atom lists.
* **The 7.** Order 5 is where the shared structure breaks.
  * The strong 5-tournaments jump to `H ≥ 9`, so `7` is skipped.
  * There are 544 strong labelled ones against 4 circulants, so `7` is the first prime that sees them in the Camion test.
  * `T5 = Z_5{1,2}` is pure-blue and converse-exceptional, but not path-rigid (`tc = 3`).
  * This is ANALOGY: both statements are about strong 5-tournaments, but there is no map between them.

## 6. Census `n ≤ 9` (FINITE-EXACT; reproduced by the audit and by S50)

* **Method.** `tour_census.c` reads every gentourng class. Class counts are 2, 4, 12, 56, 456, 6880, 191536 (A000568 ✓). For each class it computes `H`, `C`, strongness, `|Aut|`, the anti-automorphism count, and `c3`.
* **H.** Below the ceiling, the missing odd values are exactly `7, 21` at `n = 8, 9`.
  * At `n = 6` the missing values are `7, 21, 35, 39`; at `n = 7` they are `7, 21, 63, 107, 119, 149` plus the ceiling run.
  * The strong spectra have extra holes just above the Moon–Busch floor: `47, 55, 67` at `n = 8`; `77, 79, 83` at `n = 9`.
  * Class multiplicities for `n ≤ 6` are in `analyze_census.out`; the multiplicity laws are S50's lane.
* **C (directed Hamiltonian cycles).**
  * Every value from 0 to the maximum occurs at `n = 3..6` and `n = 8`.
  * At `n = 7` the values `16, 18–23` are missing (maximum 24, at `P_7` only).
  * At `n = 9`, 21 values in `[215, 248]` are missing (maximum 249).
  * There are no floor gaps: the cycle spectrum has no analogue of 7 and 21.
* **tc = H/|Aut|.** The missing odd values below 60 are exactly `{7, 21}` at `n = 7, 8, 9` (THM-4741, HYP-9352).
* **Labelled multiplicities.** `N_n(h) = Σ n!/|Aut|` is divisible by `n!/gcd(n!, h)` (forced by `|Aut| | h`). There are no violations at `n ≤ 9`.
* **Divisibility frequencies.** The frequencies of `p | H` match `1/p` for `p = 3, 5, 7, 11, 13`, so 7 has no special bias (also S50's HYP-9355).
* **Strong tilings.** Tournaments containing the path `1→…→n` that are strong: `σ(n) = 1, 0, 1, 5, 50, 903, 30773, …` (not in the OEIS). They compose as an ordinary generating function, and `σ(p) mod p` shows no Fermat-type law (`strong_tilings.out`). The labelled count with its rotation action is what carries the congruence.

## 7. Dictionary (typed)

Each line gives source → target, the map, what is preserved, and the status.

**7.1 Cyclotomic / Cipolla.**
* `7 = Φ_3(2) = |PG(2,2)|`, `21 = Φ_3(4) = Φ_3(2)Φ_3(−2) = 7·3 = |PG(2,4)|`, `341 = Φ_5(4) = 31·11 = |PG(4,4)|`.
* The Cipolla numbers are `N_q = (4^q−1)/3 = M_q·W_q` (Mersenne times Wagstaff); they are base-2 pseudoprimes for prime `q ≥ 5`.
* **PROVED mechanism (THM-4746).** At `M_q` the residue of `N_q` is `W_q`. So `N_q` passes the strong test at `M_q` iff `M_q | D(W_q)`. That holds at `q = 3` because `W_3 = 3` is the triangle class, and fails at `q = 5` (`31 ∤ D(11)`).
* The PG-cardinality reading is DICTIONARY: there is no tournament map to projective planes.

**7.2 `{2, 3, 11}`.**
* `1023 = 2^10 − 1 = 3·11·31` passes with the triangle residue at 11 and 31; `341 = 1023/3` fails with residue `11`. So `{2, 3, 11}` = (base, triangle, Wagstaff residue) is DICTIONARY with a PROVED core.
* The repo's other `{2, 3, 11}` (Collatz integer-cycle periods, THM-4591; odd-square brackets; the `(2,3,11)` triangle group) have no map to this. The exception is `PSL(2,11)`, a `(2,3,11)`-generated group that acts on the biplane of `P_11` (HYP-9357): STRUCTURAL, but only as a group action.

**7.3 Fano and the seven-point sphere.**
* The seven points of the Fano plane are the cosets `PSL(2,7)/S_4`: Galois' exceptional action.
* The stabiliser `S_4` is the full symmetry group of the tetrahedron, and Fano minus a point is the edge matroid of `K_4` (session23b).
* So the Fano plane is "the tetrahedron plus one point" in a precise group sense. It is also the out-neighbourhood design of `P_7`, the third primitive converse exception.
* S49's THM-4742 shows that the pasted seven points of a tetrahedron are not cospherical: six labels lie on one sphere, and the vertex is separated by a sum of squares. Both are seven-from-a-tetrahedron constructions in which one incidence is special.
* ANALOGY, with the group statement PROVED (classical).

**7.4 Sign involution.**
* The Camion test sees `(2/p)`: the quadratic character of 2, i.e. whether doubling preserves or reverses Paley orientation.
* `×2` is an automorphism of `P_7, P_23, P_31` (`2 ∈ QR`) and an anti-automorphism of `P_3, P_11, P_19` (`2 ∈ NQR`).
* In the Cipolla pair, the Mersenne factor carries the doubling-automorphic Paley tournament and the Wagstaff factor the doubling-anti-automorphic one: `(P_7, P_3)`, `(P_31, P_11)`, `(P_127, P_43)`.
* `P_7[C_3]`, of order 21, is the converse exception built on the `q = 3` pair. `P_31[P_11]`, of order 341, is not exceptional.
* STRUCTURAL for the doubling statements (`(2/q) = 1` iff `q ≡ ±1 mod 8`; Mersenne primes are `≡ 7`, Wagstaff primes `≡ 3 mod 8`). ANALOGY for the rest.

**7.5 Pseudoprimes of the stacking monoid.**
* A composite `h` all of whose factorisations use `7` or `21` (`35, 49, 63, 77, 91, 119, 133, 161, …, 343`) is realised only by tournaments with a single non-singleton strong component (corollary of THM-1370).
* This is the stacking monoid's analogue of a pseudoprime. It is left to S50's strong-atom lane.

**7.6 Sunflowers (arXiv:2610.26686).** Bucić–He–Lee prove that families in `2^[n]` with no `k`-sunflower have size `≤ 3(2 − 1/(16k³))^n`.
* The OCF `H = I(Ω, 2)` counts families of pairwise disjoint odd cycles with weight `2^k`, i.e. empty-kernel sunflowers among odd cycles.
* The converse towers `C_3^[k] = Cay(F_3^k, lexicographic)` live on the cap-set space. Every line is a cyclic triangle, though not conversely. Its cyclic-triangle-free sets are its transitive subtournaments, of maximum size exactly `2^k` (two of three blocks at each level). Cap sets (no line) can be larger, `≥ 2.2^k`; the `k = 3` sunflower bound comes from them (Alon–Shpilka–Umans).
* ANALOGY; cheapest test: whether a sunflower-type bound controls the number of odd-cycle packings, i.e. `I(Ω, x)` for large `x`.

## 8. Next obligations

1. HYP-9358(d′): the mod-`p^k` law for `k ≥ 3` (derivatives up to order `2(k−1)`); a faster `p²` search for Camion–Wieferich primes beyond `10^4`.
2. HYP-9357: `H(P_27)`, `H(P_31)` by an orbit-reduced or meet-in-the-middle count; a stronger search at `n = 15, 17`; certify or improve `max H(13) ≥ 3719893`.
3. HYP-9358(e): the exact thresholds for shifted and rooted decompositions.

## 9. Reproduction

```bash
cd 04-computation/experiments/tournament_pseudoprimes_20261021
clang -O3 -o circ_counts circ_counts.c && for n in 3 5 7 9 11; do ./circ_counts $n; done > circ_small.out
./circ_counts 19 paley > paley19.out; ./circ_counts 23 paley > paley23.out; ./circ_counts 23 > circ_23.out   # ~4 min
clang -O3 -o tour_census tour_census.c && for k in 3 4 5 6 7 8 9; do gentourng -q $k | ./tour_census $k > census_n$k.txt; done   # n = 9 stored gzipped
clang -O3 -o camion_psp camion_psp.c && ./camion_psp 1200 > camion_psp_1200.out
clang -O3 -o camion_verify camion_verify.c && ./camion_verify 199 2500 1500 > camion_verify.out
clang -O3 -o camion_list camion_list.c && ./camion_list 20000 3000 > camion_list_20000.out   # ~6 min
clang -O3 -o camion_wieferich camion_wieferich.c && ./camion_wieferich 3 10000 > camion_wieferich_10000.out
clang -O3 -o hmax_search hmax_search.c && ./hmax_search 13 40 7 > hmax_search_13.out; ./hmax_search 13 120 23 > hmax_search_13b.out; ./hmax_search 15 16 11 > hmax_search_15.out
for s in camion_period camion_split camion_even_sq camion_D_factor camion_controls camion_sporadic check_congruences analyze_census analyze2 analyze3 mersenne_camion strong_tilings; do python3 $s.py > $s.out; done
cd audit && sh run_audit.sh   # the independent audit, ~5.5 min
```

## 10. Status

| item | status |
|---|---|
| quasi-periodicity, Euler criterion, rotation lemma, mod 3 / mod 5, Korselt criterion, families, count-test comparison, semicomplete and connected-graph analogues, mod-`p²` lemma (THM-4746) | PROVED (audited; Lemmas B, C by the auditor) |
| criterion vs direct, lists, zero classes, Camion–Wieferich ≤ 10000, Mersenne `k ≤ 64`, controls | FINITE-EXACT |
| fixed-vertex lemma, circulant/Paley cycle congruences, `ρ_0` law, Wilson (THM-4747) | PROVED (audited) |
| circulant and Paley data `n ≤ 23`; `n = 13` rigid record | FINITE-EXACT |
| Paley H-max iff `{3, 7, 11}` iff 2-transitive design (HYP-9357) | OPEN (data + CITED Kantor; design part PROVED given Kantor) |
| Camion–Wieferich infinitude, Poisson classes, counting, mod `p^k` for `k ≥ 3`, thresholds (HYP-9358) | OPEN |
| §7 bridges | as typed there |
| novelty of THM-4746 | no prior statement found of the explicit law, the criterion, or Camion pseudoprimes; prior art for periodicity in principle (Specker–Blatter) and for connected graphs (Mani–Stones) credited |

## 11. Audit

The independent adversarial audit reproduced every number and checked every proof. Its report is [here](tournament_pseudoprimes_audit_20261021.md); its scripts are in `audit/`.

**Corrections applied (MISTAKE-651).**
* **Three MAJOR.**
  * Connected graphs obey the law (Lemma C; Mani–Stones), so the "linear condensation is the mechanism" thesis is withdrawn.
  * The Mersenne table had omitted `k = 6, 12` and miscalled `k = 59`.
  * "`D` is not periodic mod `p²`" is replaced by Lemma B. Lemma B proves the mod-9 law, the `p²` equivalence and the cube-free criterion.
* **Ten MINOR.**
  * The class count is 262.
  * The Kantor citation and the Sylow step of HYP-9357 were fixed.
  * The equivalence "`p²` passes ⟺ `p² | D(p)`" is now proved.
  * The pointer to the exotic list now points to §2, and the even sporadic members were added.
  * The 142/163 even-case count was corrected.
  * "3 and 5 test nothing" now reads "for odd `n`".
  * The missing hypotheses `p ≥ 5`, `odd n` and `μ` odd were added.
  * The "climbs towards e" remark was replaced.
* **Credits added.** Specker–Blatter; Mani–Stones; Alon (1990) and Wormald; Adler–Alon–Ross; Peterson/Sagan and Anderson–Benjamin–Rouse for the rotation proofs; Gessel.
* **Overclaim fixed.** "The triangle exactly separates the two tests" now reads "the universal residue sets differ by the triangle class".

## 12. Continuation: every odd prime power (owner request, 2026-10-21)

The owner asked to prove and test HYP-9358(d′): `s(m + φ(p^k)) ≡ (2/p)·s(m) (mod p^k)` for `m ≥ 2k − 1`.

**Proof.** The concurrent **Session25** lane published the proof first ([additive proof](session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md), [independent review](session25_spectrum_wiring_20261021/session25-additive-review/REVIEW.md)). It is now PROVED in THM-4746 §7.

* This session found the same proof independently, by the same method:
  * write `2^C(m,2) = ω^C(m,2)·(1+pu)^C(m,2)`, with `ω` the Teichmüller lift;
  * expand in divided derivatives `C(m,d)·x^(m−d)` over Teichmüller nodes, with the coefficient of degree `d` divisible by `p^⌈d/2⌉`;
  * show that this filtration is closed under products and that `T` is invertible in it (Neumann series);
  * conclude that the delta terms vanish mod `p^k` once `m ≥ 2k − 1`.
* It is credited to Session25 and not filed as a separate theorem.
* Session25 also proved the general-base, shift (connected graphs, `n ≥ 2r`) and unrestricted CRT corollaries.

**New here (THM-4746 §7b, statement 9).**

* **Prime-power towers (PROVED).** `D(p^a) ≡ (2/p)D(p^(a−1)) (mod p^a)`. So the exponents with `p^a` a Camion pseudoprime form an initial segment `{1, …, h(p)}`, and `h(p)` equals the frozen eventual value of `v_p(D(p^a))`.
  * `h(p) ≥ 2` iff `p` is Camion–Wieferich.
  * `h(3) = h(13) = 2`, and `h(p) = 1` for every other odd `p ≤ 10^4`.
  * So for `p ≤ 10^4` the odd prime-power Camion pseudoprimes are exactly the primes, `9` and `169`.
  * This mirrors the classical fact that `p^a` is a base-2 Fermat pseudoprime iff `a ≤ v_p(2^(p−1) − 1)`.
* **Independent replication (FINITE-EXACT).**
  * 48 pairs `(p, k)`, every index up to `φ(p^k) + 4k + 20`, with no failure at `m ≥ 2k − 1`.
  * The unrestricted CRT criterion against direct `s(n) mod n` for all 2569 composites `n ≤ 3000`.
  * Towers to `3^9, 5^6, 7^5, 11^4, 13^3, 17^3, 19^3, 23^3`.
  * Connected graphs, 15 pairs.
  * Valuations and ten threshold entries were recomputed by an independent Python implementation.
* **Thresholds (data for HYP-9358(e)).** `2k − 1` is sharp for generic `p`. It is lower for `p < 2k − 1`; for example `τ(3,10) = 13` against the sufficient 19. The table and the conjectures (e1), (e2) are in HYP-9358.
* **Literature.** The connected-graph congruence of Mani–Stones, as quoted by Kotek–Makowsky (Prop. 1.11; the original was not accessible here), is stated for `n ≥ p^r`. The Hurwitz argument extends it to `n ≥ 2r`.

**Scripts.** `04-computation/experiments/tournament_pseudoprimes_20261021/modpk/`:
* C: `pk_period.c`, `run_grid.sh`, `korselt_full.c`, `tower_deep.c`, `cg_pk.c`;
* independent Python: `freeze_check.py`, `tau_check.py`;
* each with its `.out`.

## 13. Continuation 2: the Camion function (owner request: hunt h(p) ≥ 3, HYP-9357 at 27 and 31, long hypothesis/test loops)

The canon statement is THM-4750 (pre-audit); the open laws are HYP-9361. This section records how the loops went.

**Loop 1: h(p) ≥ 3 hunt.**
* The fast Camion–Wieferich test `s(p) ≡ ε + p·u_p (mod p²)` (`cw_fast.c`) took the search to `10^5`.
* It found one new prime, **22409**, confirmed by a second code that computes `D(p) mod p³` directly.
* Measuring its tower height required the new idea of the session. `D(p^a)` cannot be reached directly, but the higher Kummer congruences make `l ↦ ε^(−l)D(p + l(p−1))` a polynomial mod `p^K`. Evaluating it at the `p`-adic point `l = −p/(p−1)` (that is, `n = 0`) gives the **Camion limit `Λ_p`**, and `h = v_p(Λ_p)`.
* `h(22409) = 2`. No prime `≤ 10^5` has `h ≥ 3`. The heuristic says there are finitely many such primes (`Σ 1/p² ≈ 0.09`).

**Loop 2: what the congruences represent.**
* The Mahler coefficients decay at half the Kummer rate: one power of `p` per two steps, because each power of the Fermat drift `u` arrives with the degree-2 polynomial `C(n,2)`.
  * This is sharp at every tested non-Wieferich prime `p ≥ 7`.
  * At `p = 3, 5, 7` there is an extra `v_p((2k−1)!!)` (`camion_speed.out`; audit finding 28).
* That suggested the curvature formula `a_2/p ≡ −q_2(p)X(m−2)`, proved and then checked in every class of every prime `≤ 10^4`.
* Its prediction is that at Wieferich primes the counts become linear mod `p²`, so each zero class has one second-order root. At other primes the roots of a quadratic come in pairs (0 or 2).
  * The prediction was made before the test and confirmed exactly: 4 of 4 classes at 1093 and 4 of 4 at 3511 have one root.
  * Over 1226 primes the split:inert ratio is 1:1. The triangle class is always split, by its exact zero `D(3) = 0`.
* Two further predictions also passed:
  * **base 3** (strong semicomplete digraphs) is at full speed exactly at the base-3 Wieferich prime 11;
  * **`r = 3`** (connected 3-uniform hypergraphs) gives cubic class sequences with random-cubic root statistics.
* Recorded but not claimed: a deficit of double roots, and a slight underdispersion of the sporadic counts.

**Loop 3: sparks from the random wide-net survey** (122 files sampled across canon, results, reflections and write-ups).
* *Spark: "(E^(p−1) − ε)^r s ≡ 0 mod p^r".* REFUTED except at Wieferich primes. Half speed gives only `(E^(p−1) − ε)^(2k−1) s ≡ 0 mod p^k`.
* *Spark: "a 2-adic Camion face via Moon's odd-order theorem".* Confirmed and sharpened into THM-4750 (6):
  * `T ≡ Φ_3 (mod 2)` forces the parity law for strong isomorphism classes;
  * the groupoid count `f(n) = Σ 1/|Aut|` is 2-adically quasi-periodic with period `3·2^j`, since `F_4^×` replaces `F_p^×`;
  * on `n ≡ 2 (mod 3)` it has a single simple 2-adic zero `n* ≡ 3250 (mod 4096)`.
  * The spark's guessed formula `3 + v_2((n−2)/6)` fails from `n = 50`. The truth is `2 + v_2(n − n*)`.
* Repunit-side sparks (one-apex stratum, Fermat-type monoid, free-vertex doubling, carry adjacency, and others) were forwarded to opus S51, whose lane that is. Session25/26 already have the one-apex recurrence.
* The HYP-9357 spark "Paley + sink switching class" was not run.

**Loop 4: repunit archaeology.**
* About 45 threads were catalogued and forwarded to opus S51 as raw material for THM-4748/4749 (`agents/opus/inbox/MSG-2399`).
* The unifying reading: a base-`b` repunit truncates the fixed point `−1/(b−1)` of `x ↦ bx + 1`.
* Three unrecorded errata were found and fixed (MISTAKE-653).

**Loop 5: HYP-9357.** The rotational circulant beats Paley at `q = 27` and `q = 31` (two independent exact codes; ratios 1.0205 and 1.0239). The gap keeps growing with `q`.

**What it all represents.**
* Mod `p`, a labelled count in the `2^C(n,2)` family is a twisted exponential sum over the Teichmüller group. That is Camion pseudoprimality.
* Mod `p^k` it is a `p`-adic analytic function of the class index. Its speed is `1/(degree of the edge count)`, and its curvature is the Fermat quotient of the number of edge states.
* The tower heights are the valuations of one regularized value, `Λ_p`, the branch through `K1` evaluated at "zero vertices".
* At `p = 2` the same picture is controlled by Moon's theorem and by `Φ_3`, the mod-2 shadow of the tournament EGF.

**Audit (2026-10-21; [report](camion_function_audit_20261021.md); MISTAKE-654).**
* The core proofs are sound, and every FINITE-EXACT number reproduced with independent code (audit b1–b11).
* Two MAJOR fixes were applied:
  * **M1.** A double root mod `p` can be a ramified pair of zeros: the triangle class at `p = 5`, class 2 at `p = 7`.
  * **M2.** The `n*` rational-reconstruction evidence was replaced by a computation mod `2^140`: no rational of height `≤ 2^69.5`.
* Sixteen MINOR fixes and four upgrades were applied: the congruence mod `p^(a+1)`, the cheaper `Λ_p` formula, the exact double-factorial excess, and an entire-function proof sketch for the 2-adic statements.
* Post-audit continuation, not covered by that audit:
  * Teichmüller rigidity, `Λ_p = σ_p − ⟨2⟩^(−1/2)`, and the **Camion lift** `κ(2) = ω(2)σ_p^(−2)` with `h(p) = v_p(2 − κ(2))` (THM-4750 §5b);
  * the `r = 4` `S_4` statistics (§6b);
  * the Camion–Wieferich lists for base 3 (17, 31) and connected graphs (5, 197);
  * the FLINT search to `10^6`.

