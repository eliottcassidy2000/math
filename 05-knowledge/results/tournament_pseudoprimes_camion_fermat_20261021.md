# Tournament pseudoprimes: the Camion test, Hamiltonian-cycle congruences, Paley extremality, and the 7 / 21 / 341 dictionary

Session mac-mini-2026-10-21-fanopsp. **Canon:** THM-4746, THM-4747. **Hypotheses:** HYP-9357, HYP-9358.

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
| opus-2026-10-21-S50 | THM-4744/4745 and HYP-9355/9356 reserved: the multiplicity laws of the H spectrum |

* I sent S50 a lane message (`agents/opus/inbox/MSG-2395…`) and took only the pseudoprime lens.
* The link arXiv:2610.26686 is Bucić–He–Lee's proof of the Erdős–Szemerédi sunflower conjecture (§7.6), not the Gluck paper S49 cites.

**Scripts.** All are in `04-computation/experiments/tournament_pseudoprimes_20261021/`, each with an `.out`.

* C: `circ_counts.c`, `tour_census.c`, `camion_psp.c`, `camion_verify.c`, `camion_list.c`.
* Python: `camion_*.py`, `check_congruences.py`, `analyze*.py`, `mersenne_camion.py`.
* Reproduction commands are in §9.

## 0. Results

1. **Camion pseudoprimes (THM-4746, PROVED + FINITE-EXACT).**
   * At prime order `p`, rotations force #strong labelled tournaments ≡ #circulants `= 2^((p−1)/2) ≡ (2/p) (mod p)`.
   * The strong count is quasi-periodic: `s(m+p−1) ≡ (2/p) s(m) (mod p)`.
   * This gives a **Korselt-type criterion** whose free residues `1, 2, 3` are K1, the missing strong K2, and the **triangle**.
   * Consequences:
     * every Carmichael number, every `2p`, every `3p`, and every `3b(3b−2)` with both factors prime is a pseudoprime;
     * this includes `4^q − 1` whenever `2^q − 1` and `(2^q+1)/3` are both prime;
     * `21` is one; `341` is not, because its residue at 31 is the Wagstaff factor 11; `3·341 = 1023 = 2^10 − 1` is one.
2. **Hamiltonian-cycle congruences (THM-4747, PROVED + FINITE-EXACT).**
   * An automorphism of prime order with a fixed vertex divides `C(T)`.
   * Circulants on `Z_p`: `C ≡ (p−1)/2 (mod p·μ_S)`; Paley: `C(P_p) ≡ (p−1)/2 (mod p(p−1)/2)`.
   * Reflection-symmetric paths: `ρ_0 ≡ μ_S (mod 2μ_S)`.
   * Summing over tournaments gives Wilson's theorem.
   * Composite orders 9, 15, 21 pass no circulant version of the naive test.
3. **Paley extremality (HYP-9357, FINITE-EXACT + CITED).**
   * `P_q` is the circulant H-maximum exactly at `q = 3, 7, 11`; the rotational circulant beats it at 19 and 23.
   * These are exactly the `q ≡ 3 (mod 4)` whose Paley design is 2-transitive (Kantor), i.e. Galois' exceptional primes.
   * The converse-Gluck Paley exceptions `{3, 7}` sit inside this list.
4. **Census `n ≤ 9` (FINITE-EXACT, lens-specific; S50 owns the multiplicity laws).**
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
  * THM-1370 (1): `H` is multiplicative over the ordered strong components. The same ordered decomposition is Wright's recurrence, which drives THM-4746.
  * LEM-003: automorphisms act freely on Hamiltonian paths.
  * THM-2477: the path-rigid strong atoms are `K1, C3`.
  * THM-4743 / Session24.
* **Canonical hostile examples.**
  * `341`: the first base-2 pseudoprime, which fails the strong test.
  * Connected graphs: same total `2^C(m,2)`, unordered components, no quasi-period.
* **Corrected near misses (this session).**
  * (i) I first read the `n = 7` interior holes of the H spectrum (`63, 107, 119, 149`) as "prime-like" (realisable only by strong tournaments). That is vacuous: at `n = 7` every `H ≥ 47` comes only from strong tournaments, by vertex count. Withdrawn before writing.
  * (ii) `tc(P_11) = 1729` is a Carmichael number, but `tc(P_19)` fails Korselt at 7: NUMEROLOGY.
  * (iii) A program bug (zero classes allocated only for `p ≤ N/3`) truncated a first list run; it was fixed and rerun.
* **Least-used sidecar.** The Hurwitz ring of binomial convolution over `F_p`, and the group algebra `K[(Z/p)]` inside it.

## 2. Camion pseudoprimes (THM-4746)

**Definitions.**
* `s(m)` counts strong labelled tournaments: A054946 `1, 0, 2, 24, 544, 22320, 1677488, …`. By Camion's theorem these are exactly the tournaments with a Hamiltonian cycle.
* `c(m)` counts circulants on `Z_m`: `2^((m−1)/2)` for odd `m`, `0` for even `m`.
* `D = s − c`. A composite `n` is a **Camion pseudoprime** if `s(n) ≡ c(n) (mod n)`.

**Why primes pass.** The rotation `x ↦ x+1` of `Z_p` acts on labelled strong tournaments. Its fixed points are the circulants, all strong. So `s(p) ≡ c(p) = 2^((p−1)/2) ≡ (2/p)` by Euler's criterion.

**Quasi-periodicity (proof in THM-4746).**
* `2^C(m,2)` is quasi-periodic mod `p` with multiplier `(2/p)`, so it is a combination of exponentials `E_g = (g^m)`, with `g` running over an additive group `G ≅ Z/p` in `F_p` or `F_(p²)`.
* Wright's recurrence is `(e − S) ∗ T = e` in the Hurwitz ring.
* `K[G]` is a local ring, so `e − S` is again a combination of the same exponentials.
* Hence `s(m+p−1) ≡ (2/p)s(m)`.
* Linearity of the condensation is essential:
  * strong semicomplete digraphs obey the same law with `(3/p)`;
  * connected graphs obey none (`camion_controls.out`).

**Korselt-type criterion.** For `n` with squarefree odd part:

    s(n) ≡ c(n) (mod n)  ⟺  p | D(n mod (p−1)) for every odd prime p | n.

* It agrees with direct computation for all 472 odd squarefree composites `≤ 2500` and all 142 even `n ≤ 400` (`camion_verify.out`, `camion_even_sq.out`).

**Zero classes.**
* `Z_p = {r : p | D(r)} ⊇ {1, 2, 3}`, because `D(1) = D(2) = D(3) = 0`.
* For `r ≥ 4`, `D(r) > 0`, and the extra (sporadic) classes are Poisson with mean `(p−4)/p` (HYP-9358(b)). Examples: `43: 39`; `47: 17`; `59: 21`; `103: 37, 43`; `139: 77, 103`; ….
* `3` and `5` are blind primes: the odd residues mod `p−1 ∈ {2, 4}` are only `1` and `3`.
* `7` is the first prime with a nonzero odd class, `r = 5`: `D(5) = 540 = 2²·3³·5 ≢ 0 (mod 7)`.

**Families (all PROVED by the criterion).**

| family | residues | examples |
|---|---|---|
| Carmichael | all `r = 1` | 561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841 (all ≤ 20000 pass) |
| `2p` | `r_p = 2` | 6, 10, 14, 22, 26, … (1228 below 20000) |
| `3p` | `r_p = 3`, `r_3 = 1` | 9, 15, 21, 33, 39, 51, 57, 69, … (857 below 20000) |
| `3b(3b−2)`, `b ≥ 5`, both prime | both `3` | 195, 399, 1023, 1443, 4623, 12099, 16383, 19599 |
| `4^q − 1`, `2^q−1` and `(2^q+1)/3` prime | both `3` | `q = 5, 7, 13, 17, 19, 31, …` (`k = 10, 14, 26, 34, 38, 62` among `2^k − 1`, `k ≤ 64`) |
| sporadic | uses `Z_p \ {1,2,3}` | 1239 (59: 21), 3765 (251: 15), 8401 (271: 31), 12019 (101: 19), 14317 (103: 37; 139: 103) |

**Prime squares.**
* `9` and `169` pass; `25, 27, 49, 81, 121, 125` fail.
* The Camion–Wieferich primes (`D(p) ≡ 0 mod p²`) below 3000 are only `3` and `13`.
* Modulo 9, the bad odd `m` are exactly `m ≡ 1 (mod 6)`, `m ≥ 7` (to 700).
* `D` is not periodic modulo `p²` with any period dividing `2p(p−1)`.

**The tournament-count test** `2^C(n,2) ≡ c(n) (mod n)` is the Fermat test to the base `c(n) = 2^((n−1)/2)`.
* It contains all odd base-2 pseudoprimes, plus `85 = (4⁴−1)/3, 133, 481, 511 = 2⁹−1, 949, …`. This sequence is not in the OEIS.
* Its only universal residue is `1`. **The triangle is exactly the difference between the two tests.**

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
* For `Cay(Z_p, S)` with `μ = |{a : aS = S}|`:
  * the `(p−1)/2` arithmetic cycles form orbits of size `μ`;
  * every other cycle has a free orbit of size `pμ`;
  * hence `C ≡ (p−1)/2 (mod pμ)`.
* The signed group of order `2pμ` acts on paths, and only the reflections `x ↦ 2c−x` fix paths. Burnside gives `ρ_0 ≡ μ (mod 2μ)`.
* An involutory anti-automorphism always reverses an odd number of Hamiltonian paths (Rédei).

**Checks.** All 264 circulant classes `3 ≤ n ≤ 23` pass. At `n = 9, 15, 21`, no class satisfies `C ≡ #(S ∩ Z_n^*) (mod n)`. Summing over labelled tournaments yields Wilson's theorem, so the summed cycle count has no pseudoprimes.

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
* **Ratios.** `H·2^(n−1)/n!` climbs towards `e`: from 2.40 (`P_7`) to 2.60 (rotational 23).
* **Why the list `{3, 7, 11}` is not arbitrary.**
  * It is Galois' list of exceptional primes for `PSL(2,p)` acting on `p` points, restricted to `q ≡ 3 (mod 4)`.
  * It is also Kantor's list for 2-transitive Paley designs: the trivial design, the Fano plane, the biplane.
  * By the multiplier-group argument, `QR_q` is a Singer set only for `q = 3, 7`.
* **Nesting.** Converse-Gluck `{3, 7}` ⊂ H-maximal `{3, 7, 11}` = 2-transitive design.
* **Mechanism.** A HEURISTIC: design symmetry balances the 5- and 6-vertex statistics.

## 5. What the pseudoprime lens says about the 7/21 gaps (honest scope)

* **Not explained here.** The gaps are THM-1370: multiplicativity plus the strong floor `3, 5 | 9`. Nothing here gives a new proof.
* **Shared mechanism.** The same object, the **ordered strong decomposition**, produces three things:
  1. the multiplicativity of `H`, and so the gaps;
  2. the quasi-periodicity of `s`, and so the Camion pseudoprimes;
  3. the multiplicativity of `tc = H/|Aut|` over components. That is why the presentation-count spectrum inherits exactly the same gaps `7, 21` (`analyze_census.out`: no other odd value below 60 is missing at `n = 7, 8, 9`).
* **Shared atoms.** The orders where every strong tournament is a circulant are `1` and `3`: K1 and C3. They recur across the canon:
  * the universal Camion residues;
  * the path-rigid atoms (`tc = 1`, THM-2477/4741);
  * the first two members of the pure-blue alphabet `{K1, C3, T5}` (THM-4730);
  * the substituted towers of every converse-Gluck exception (Session24);
  * the first two H-values `1, 3`.
* **The 7.** Order 5 is where the shared structure breaks.
  * The strong 5-tournaments jump to `H ≥ 9`, so `7` is skipped.
  * There are 544 strong labelled ones against 4 circulants, so `7` is the first prime that sees them in the Camion test.
  * `T5 = Z_5{1,2}` is pure-blue and converse-exceptional, but not path-rigid (`tc = 3`).
  * This is ANALOGY: both statements are about strong 5-tournaments, but there is no map between them.

## 6. Census `n ≤ 9` (FINITE-EXACT; lens-specific)

* **Method.** `tour_census.c` reads every gentourng class. Class counts are 2, 4, 12, 56, 456, 6880, 191536 (A000568 ✓). For each class it computes `H`, `C`, strongness, `|Aut|`, the anti-automorphism count, and `c3`.
* **H.** Below the ceiling, the missing odd values are exactly `7, 21` at `n = 8, 9`.
  * At `n = 6` the missing values are `7, 21, 35, 39`; at `n = 7` they are `7, 21, 63, 107, 119, 149` plus the ceiling run.
  * The strong spectra have extra holes just above the Moon–Busch floor: `47, 55, 67` at `n = 8`; `77, 79, 83` at `n = 9`.
  * The class multiplicities for `n ≤ 6` are in `analyze_census.out`; the laws are S50's lane.
* **C (directed Hamiltonian cycles).**
  * Every value from 0 to the maximum occurs at `n = 3..6` and `n = 8`.
  * At `n = 7` the values `16, 18–23` are missing (maximum 24, at `P_7`).
  * At `n = 9`, 21 values in `[215, 248]` are missing (maximum 249).
  * There are no floor gaps: the cycle spectrum has no analogue of 7 and 21.
* **tc = H/|Aut|.** The missing odd values below 60 are exactly `{7, 21}` at `n = 7, 8, 9` (THM-4741, HYP-9352).
* **Labelled multiplicities.** `N_n(h) = Σ n!/|Aut|` is divisible by `n!/gcd(n!, h)` (forced by `|Aut| | h`). There are no violations at `n ≤ 9`.
* **Divisibility frequencies.** The frequencies of `p | H` match `1/p` for `p = 3, 5, 7, 11, 13`, so 7 has no special bias.

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
* The converse towers `C_3^[k] = Cay(F_3^k, lexicographic)` live on the cap-set space. Its lines are cyclic triangles. Its cyclic-triangle-free sets are its transitive subtournaments, of maximum size exactly `2^k` (two of three blocks at each level). Cap sets (no line) can be larger, `≥ 2.2^k`; the `k = 3` sunflower bound comes from them (Alon–Shpilka–Umans).
* ANALOGY; cheapest test: whether a sunflower-type bound controls the number of odd-cycle packings, i.e. `I(Ω, x)` for large `x`.

## 8. Next obligations

1. HYP-9358: a mod-`p²` refinement (prime squares; the mod-9 law); search for Camion–Wieferich primes beyond 3000 with a faster `p²` recurrence.
2. HYP-9357: `H(P_27)`, `H(P_31)` by an orbit-reduced or meet-in-the-middle count; a local search for non-circulant maxima at `n = 13, 15`.
3. Quasi-periodicity for other linearly decomposed classes (HYP-9358(e)); a Camion analogue for presentation counts (rotation classes of closed tilings).
4. Cross-check S50's multiplicity tables against `census_n*.txt` when they land.

## 9. Reproduction

```bash
cd 04-computation/experiments/tournament_pseudoprimes_20261021
clang -O3 -o circ_counts circ_counts.c && for n in 3 5 7 9 11; do ./circ_counts $n; done > circ_small.out
./circ_counts 19 paley > paley19.out; ./circ_counts 23 paley > paley23.out; ./circ_counts 23 > circ_23.out   # ~4 min
clang -O3 -o tour_census tour_census.c && for k in 3 4 5 6 7 8 9; do gentourng -q $k | ./tour_census $k > census_n$k.txt; done
clang -O3 -o camion_psp camion_psp.c && ./camion_psp 1200 > camion_psp_1200.out
clang -O3 -o camion_verify camion_verify.c && ./camion_verify 199 2500 1500 > camion_verify.out
clang -O3 -o camion_list camion_list.c && ./camion_list 20000 3000 > camion_list_20000.out   # ~6 min
for s in camion_period camion_split camion_even_sq camion_D_factor camion_controls camion_sporadic check_congruences analyze_census analyze2 analyze3 mersenne_camion; do python3 $s.py > $s.out; done
```

## 10. Status

| item | status |
|---|---|
| quasi-periodicity, Euler criterion, rotation lemma, mod 3 / mod 5, Korselt criterion, families, count-test comparison, semicomplete analogue (THM-4746) | PROVED |
| criterion vs direct, lists, zero classes, Camion–Wieferich ≤ 3000, Mersenne `k ≤ 64`, controls | FINITE-EXACT |
| fixed-vertex lemma, circulant/Paley cycle congruences, `ρ_0` law, Wilson (THM-4747) | PROVED |
| circulant and Paley data `n ≤ 23` | FINITE-EXACT |
| Paley H-max iff `{3, 7, 11}` iff 2-transitive design (HYP-9357) | OPEN (data + CITED Kantor) |
| Camion–Wieferich infinitude, Poisson classes, counting, mod 9, general principle (HYP-9358) | OPEN |
| §7 bridges | as typed there |
| novelty of THM-4746 | NOT CHECKED beyond OEIS (no congruence recorded on A054946) |
