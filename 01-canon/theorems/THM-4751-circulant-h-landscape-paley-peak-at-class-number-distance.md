---
id: THM-4751
title: "The circulant H-landscape: a Paley peak at class-number distance from the rotational peak. (1) For a prime q = 3 mod 4, q > 3, the Paley tournament P_q and the rotational tournament R_q = Cay(Z_q, {1,...,(q-1)/2}) are B_q = ((q-1)/2 - (2 - (2/q)) h(-q))/2 difference-class flips apart in the circulant class graph (minimum over multipliers), and agree on a fraction 1/2 + (2 - (2/q)) h(-q)/(q-1) of their arcs (Dirichlet; the half-row imbalance of the Ogg note, the circulant companion of THM-4666). (2) B_q <= 1 iff q in {3, 7, 11}; B_q >= 2 for every q >= 19 (elementary). (3) FINITE-EXACT, odd n <= 25: the H-maximal circulant is unique and is R_n or a one-flip neighbour of R_n (P_7, Z_9{1,2,3,5}, P_11); the Paley class is a strict local maximum of H at q = 7, 11, 19, 23 and 31 (all its flip-neighbours are isomorphic); R_n is a strict local maximum at n = 13..31. (4) If every H-maximal circulant is within one flip of R_n (HYP-9362(a)), then P_q is never the circulant H-maximum for q >= 19, the circulant 'only if' of HYP-9357. (5) KNOWN facts, alignment NUMEROLOGY: among primes q = 3 mod 4 the spectral ring Z[(1+sqrt(-q))/2] of P_q is norm-Euclidean (covering radius^2 = (q+1)^2/(16q) < 1) exactly for q = 3, 7, 11, and a PID exactly for the Heegner primes 3, 7, 11, 19, 43, 67, 163. (6) The two peaks are the two extremes of the additive energy E(S) of the connection set: c_4 = n(E - m)/4 on circulants; R_n is the unique maximiser for every odd n (Lev; RLT c_4 theorem), P_q the unique minimiser for prime q = 3 mod 4 (flat spectrum, Galois); so the class number measures the distance between the energy extremes (prior art THM-140, THM-162). (7) FINITE-EXACT mechanism: in H = sum 2^k alpha_k the energy minimum wins alpha_1 (and alpha_2 from n = 13), R_n wins every alpha_k with k >= 3; at n = 13 R_n overtakes only through alpha_4."
status: >
  AUDITED (first pass, 2026-10-22: 05-knowledge/results/lattice_bridges_audit_20261021.md; no false statement, statuses and
  prior art corrected; MISTAKE-659, MISTAKE-660). Second pass 2026-10-22 (lattice_bridges_audit_pass2_20261022.md):
  statement 7, the HYP-474 refutation and the R_27..R_31 local maxima CONFIRMED by independent code.
  PROVED: statements 1 (Dirichlet's class number formula; the overlap A_q - B_q = (2 - (2/q))h(-q) is classical and already
  in the repo as the "half-row imbalance" of P_q in 05-knowledge/results/ogg_triangular_chowla_20260927.md; new here is the
  multiplier-minimised flip distance), 2 (elementary, least quadratic non-residue bound), 4 (formal consequence), 6
  (identities and the extremes; uniqueness of the maximiser uses the CITED RLT_n theorem).
  KNOWN: statement 5 (each side classical); its alignment with (2) is NUMEROLOGY (MISTAKE-559).
  FINITE-EXACT: statement 3 (complete multiplier-class tables of circ_counts.c for odd n <= 23 and of ie_ham.c for n = 25,
  reproduced class by class by the audit's independent a3_hcount.c; two independent H computations at q = 31 (ie_ham.c;
  audit a3_hcount.c); R_27, R_29, R_31 against all their flip-neighbours), statement 7 (alpha_decomp.c, which reproduces
  every H in the tables), and the check of statement 1 for all 1135 primes 7 <= q < 20000.
  PRIOR ART IN THE REPO:
  - the March 2026 "orientation cube" program on the same landscape: THM-134-paley-local-max-parseval.md (Paley the unique local maximum on the Parseval
    simplex, p = 7, 11; its continuous-extension claim is OPEN since Session28's 2026-10-22 scope repair), THM-135 (P_19 is not the circulant maximum; the interval wins), THM-137-paley-crossover-mechanism.md
    (Paley vs interval crossover, two-mechanism reading), THM-139 ("Paley is #2 (local max, no single flip improves)" at
    p = 19), THM-140 (interval maximum and Paley minimum of the additive energy E, E-H sign reversal), HYP-464, HYP-480
    (interval universal H-maximiser, exhaustive for primes p <= 19 and p = 23);
  - THM-4666 (the natural order of P_p has ranking excess p h(-p)) and the Ogg note (h(-p) as the half-row imbalance);
  - THM-481 (Paley tournaments generate the extended QR codes and the Gleason ring; Kim-Sole 2008), the proved form of the
    Paley -> E_8 / Golay bridge used in the results note;
  - Savchenko, J. Graph Theory 83 (2016) 44-77: c_5 + 2c_4 = n(n^2-1)(n^2-9)/160 for every regular tournament; statement
    6's c_4 = n(E - m)/4 is its circulant translation.
  CAUTION: every condition that is a threshold in (11, 19) on primes q = 3 mod 4 yields {3, 7, 11}, because no such prime
  lies between 11 and 19; this includes statement 5 and the alpha_4 threshold of statement 7.
session: mac-mini-2026-10-21-fanopsp (owner prompt: 1729, 24 dimensions, Heegner numbers, Niemeier/Leech, {7,21} vs {8,24})
source: 05-knowledge/results/lattices_codes_tournaments_20261021.md
scripts:
  - 04-computation/experiments/lattice_bridges_20261021/class_number_distance.gp (+ .out)
  - 04-computation/experiments/lattice_bridges_20261021/circ_flip_distance.py, circ_landscape.py, paley_flip_margin.py (+ .out; data ../tournament_pseudoprimes_20261021/circ_*.out)
  - 04-computation/experiments/lattice_bridges_20261021/circ_landscape_big.py (+ circ25_H.txt, circ25_landscape.out), rot_local_max.py (+ rot27/29/31_*), P31_oneflip.out
  - 04-computation/experiments/lattice_bridges_20261021/spectral_peaks.py, energy25.py, energy_peaks_table.py (+ .out)
  - 04-computation/experiments/lattice_bridges_20261021/alpha_decomp.c, alpha_all.py (+ alpha_*.out)
  - audit: 04-computation/experiments/lattice_bridges_20261021/audit/
depends_on:
  - Dirichlet's class number formula for Q(sqrt(-q)), q = 3 mod 4, q > 3: (2 - (2/q)) h(-q) = sum_(a < q/2) (a/q) (CITED)
  - class number one: Heegner, Math. Z. 56 (1952) 227-253; Baker, Mathematika 13 (1966) 204-216; Stark, Michigan Math. J. 14 (1967) 1-27 (CITED)
  - Lev, Linear equations over F_p and moments of exponential sums, Duke Math. J. 107 (2001) 239-263 (CITED, checked at abstract level)
  - c_4 is maximised among n-tournaments, odd n >= 5, exactly by the regular locally transitive tournament RLT_n = R_n (Beineke-Harary, Canad. Math. Bull. 8 (1965); Kotzig 1968; as stated in Savchenko, arXiv:2403.19555, Section 1) (CITED)
  - THM-118 (tr A^k = k c_k for k = 3, 4, 5)
related:
  - 05-knowledge/hypotheses/HYP-9362-one-flip-conjecture-paley-peak-and-the-7-21-to-8-24-dictionary.md
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
  - 01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md
  - 01-canon/theorems/THM-481-paley-gauge-generates-gleason.md
  - 01-canon/theorems/THM-137-paley-crossover-mechanism.md
  - 01-canon/theorems/THM-140-energy-H-correlation-reversal.md
---

# THM-4751 — The circulant H-landscape: Paley peak and rotational peak

**One line.** Among circulant tournaments on `Z_q`, the number of Hamiltonian paths has a rotational peak and a separate Paley peak. They are the two extremes of the additive energy of the connection set: an arithmetic progression against a flat-spectrum set. The two peaks lie `B_q` difference-class flips apart, and `B_q` is the class number of `Q(√−q)` up to an affine change. The Paley peak is the global one only when it is adjacent to the rotational peak, which happens exactly for `q = 3, 7, 11`. In the OCF the flat set wins on single cycles and the progression on packings of three or more disjoint cycles.

## Notation

* `q` is a prime with `q ≡ 3 (mod 4)`, and `m = (q − 1)/2`. For general odd `n`, `m = (n − 1)/2`.
* `QR`, `NQR` are the nonzero squares and non-squares mod `q`.
* `P_q = Cay(Z_q, QR)` is the Paley tournament and `R_n = Cay(Z_n, {1, …, m})` the rotational tournament.
* A circulant tournament `Cay(Z_n, S)` contains exactly one of `k`, `n − k` for each `k ≤ m`. The **circulant class graph** has as vertices the classes of `S` up to multipliers `a ∈ Z_n^×`. Two vertices are adjacent when the sets differ in one difference class (one "flip").
* `A_q = |QR ∩ [1, m]|`, `B_q = m − A_q`.

## 1. Overlap and flip distance (PROVED)

By Dirichlet's class number formula, `A_q − B_q = (2 − (2/q)) h(−q)` for `q > 3`. This half-row imbalance is already in the repo (Ogg note). Hence

    A_q = (m + (2 − (2/q))h(−q))/2,   B_q = (m − (2 − (2/q))h(−q))/2.

**Arcs.** `P_q` and `R_q` share the arc `i → i + k` (`k ≤ m`) iff `k ∈ QR`. So they agree on a fraction `1/2 + (2 − (2/q))h(−q)/(q − 1)` of their arcs.

**Flip distance.** The minimum over multipliers `a` of `#{k ≤ m : −k ∈ aQR}` is `B_q`.
* `a ∈ QR` fixes `QR`, and since `−1 ∈ NQR` the count is then `B_q`.
* `a ∈ NQR` gives `A_q > B_q`, by `h > 0`.
* Multipliers act on connection sets by automorphisms of the flip hypercube. So a path of `k` flips between two classes lifts to `k` flips from `R_q` to some `aP_q`, and this minimum is the class-graph distance. The audit confirmed it by BFS for every class, odd `n ≤ 25`.

**Checked.** All 1135 primes `7 ≤ q < 20000`, with PARI `qfbclassno` (`class_number_distance.out`). The audit confirmed them by counting reduced forms.

**How far apart.** For `q > 11` the excess `A_q − B_q = (2 − (2/q))h(−q)` is at least 3, so `B_q ≤ (m − 3)/2`. Equality holds exactly when `(2 − (2/q))h(−q) = 3`:
* at the Heegner primes 19, 43, 67, 163 (`B_q = 3, 9, 15, 39`);
* at `q = 23, 31` (`h = 3`, `(2/q) = 1`; the only such primes ≡ 7 (mod 8), by Watkins).

**Companion.** THM-4666 measures `h(−p)` with the linear order on `P_p` (ranking excess `p·h(−p)`). Here the class number measures the distance to the cyclic order's tournament `R_q`.

## 2. Adjacency only for q = 3, 7, 11 (PROVED, elementary)

`B_3 = 0` and `B_7 = B_11 = 1`. `B_q ≥ 2` for every `q ≥ 19`:

* **`q ≡ 3 (mod 8)`.** 2 and 8 are non-residues, and both are `≤ m` once `q ≥ 17`.
* **`q ≡ 7 (mod 8)`.** With `n_0` the least non-residue, `n_0` and `2n_0` are non-residues.
  * The elementary bound is `n_0 < 1/2 + √(q + 1/4)`. With `m_0 = ⌈q/n_0⌉`, the number `m_0n_0 − q ∈ (0, n_0)` is a residue, so `m_0` is a non-residue and `m_0 ≥ n_0`.
  * This gives `2n_0 ≤ m` for `q ≥ 23`.

The audit checked every step for all primes `q ≡ 3 (mod 4)` below `2·10^6`.

## 3. The landscape (FINITE-EXACT)

Sources:
* the complete multiplier-class tables `circ_counts.c` for odd `n ≤ 23` and `ie_ham.c` for `n = 25` (208 classes);
* the audit reproduced these class by class with an independent inclusion–exclusion code.

**Global maximum (unique in every case).**
* It is `R_n` for `n = 13, …, 25`.
* It is a one-flip neighbour of `R_n` for `n = 7` (`P_7`), `n = 9` (`Z_9{1,2,3,5}`, the class of `Z_9{2,3,4,8}`) and `n = 11` (`P_11`).

**Strict local maxima.** These are counted in the hypercube of connection sets, where a flip that lands in the same multiplier class counts as a tie.
* `n = 7, 9, 11`: the global maximum only.
* `n = 13, 15`: `R_n` only.
* `n = 17`: `R_17` and one class at distance 2.
* `n = 19`: `R_19`, `P_19` (at `B_19 = 3`) and one class at distance 2.
* `n = 21`: 5 classes.
* `n = 23`: 6 classes, among them `P_23` (at `B_23 = 4`).
* `n = 25`: 6 classes.

If only distinct neighbour classes count (the class graph without loops), five more classes are strict. Each has a self-loop flip:
* `Z_11{1,2,3,4,6}`;
* `Z_15{1,2,5,7,9,11,12}`;
* `Z_19{1,2,5,9,11,12,13,15,16}`;
* `Z_21{1,2,3,7,8,9,10,15,16,17}`;
* `Z_23{1,2,5,6,7,8,9,11,13,19,20}`.

**The rotational peak.** `R_n` is a strict local maximum for `13 ≤ n ≤ 31`. For `n ≥ 27` this was tested against all `m` flip-neighbours by inclusion–exclusion (`rot*_local_max.out`). From `n = 25` on, the best neighbour flips `{m, m + 1}`, and its ratio is `1 − 0.104/n` to three digits.

**The Paley peak.** Residue multipliers are automorphisms of `P_q` and permute its `m` difference classes simply transitively, since `−1 ∉ QR`. So all flip-neighbours of `P_q` are isomorphic, and there is no self-loop. THM-139 already noted the one-flip peak at `p = 19`. The margins `H(neighbour)/H(P_q)` (`paley_flip_margin.out`) are:

| `q` | margin |
|---|---|
| 7 | 0.92593 |
| 11 | 0.97825 |
| 19 | 0.99246 |
| 23 | 0.99542 |
| 31 | 0.99822 |

The `q = 31` value comes from `H = 19841263005616787646549771` against `H(P_31) = 19876729650849384218015865`. Both are reproduced by the audit; `H(P_31)` also has the audited S34 value. All margins are below 1, so `P_q` is a strict local maximum.

## 4. Consequence for HYP-9357 (PROVED, conditional)

If HYP-9362 (a) holds in the form "every `H`-maximal circulant class is within one flip of `[R_n]`", then by §2 `P_q` is not the circulant `H`-maximum for any prime `q ≥ 19`. For `q = 19, 23, 31` this is already FINITE-EXACT.

## 5. Euclidean and Heegner primes (KNOWN; alignment NUMEROLOGY)

The nontrivial eigenvalues of `P_q`, `(−1 ± √−q)/2`, generate `O_K` with `K = Q(√−q)`.

* **Euclidean.** `O_K` is norm-Euclidean exactly for `q = 3, 7, 11`. The circumcentre `z* = 1/2 + ((q − 1)/(4q))√−q` of the triangle `0, 1, (1 + √−q)/2` lies in `K`. It is the deep hole, at squared distance `(q + 1)²/(16q)`, and this is `< 1` iff `q < 7 + √48`. So the computation proves the statement directly (context: Lemmermeyer, Expo. Math. 13 (1995) 385–416).
* **Heegner.** `O_K` is a PID exactly for `3, 7, 11, 19, 43, 67, 163` (Heegner, Baker, Stark).

So `{3, 7, 11}` is cut out by covering radius, while the Heegner primes 19, 43, 67, 163 are PIDs without the Euclidean property, and their Paley tournaments are far from rotational (§1). The coincidence of this 3-element set with (2) is NUMEROLOGY under MISTAKE-559: see the caution in the status.

## 6. The two peaks are the two extremes of additive energy (PROVED; prior art THM-140, THM-162)

Let `E(S) = #{(a, b, c, d) ∈ S⁴ : a + b ≡ c + d (mod n)}`, and `c_k` be the number of directed `k`-cycles.

**(6a) Identities.** On circulant tournaments of odd order `n`,

    c_4 = n(E − m)/4,   2c_5 + nE = (2m⁵ + 5m⁴ + 5m² + 3m)/5.

*Proof.*
* `A` is circulant, with eigenvalues `λ_j = Σ_(s∈S) ω^(js)` and `λ_0 = m`.
* Since `S ⊔ −S = Z_n \ {0}`, `λ_j = −1/2 + iy_j` for `j ≠ 0`.
* Parseval gives `Σ_j y_j² = m(m + 1/2)`, and `Σ_j |λ_j|⁴ = nE(S)`.
* THM-118 gives `tr A^k = k c_k` for `k = 4, 5`. Expanding `Re λ_j⁴` and `Re λ_j⁵` in `y_j` gives both formulas. ∎

The identities were checked on every class, `5 ≤ n ≤ 25`, with `c_4`, `c_5` counted by DFS in the audit. The second is the circulant case of Savchenko's identity `c_5 + 2c_4 = n(n² − 1)(n² − 9)/160`, which holds for every regular tournament. THM-162 already observed `corr(Σ|λ|⁴, c_5) = −1` at `p = 11`.

**(6b) The extremes.**
* **Maximum.** `R_n` is the unique maximiser of `E`, up to multipliers, for every odd `n ≥ 5`.
  * For prime `n` the inequality is Lev's theorem: arithmetic progressions maximise the number of solutions of `a_1 + a_2 − a_3 − a_4 = 0` with the variables in sets of given sizes.
  * For every odd `n`, uniqueness follows from the CITED theorem that `RLT_n ≅ R_n` uniquely maximises `c_4` among `n`-tournaments, together with (6a).
  * `Aut(R_n) = Z_n`, since the stabiliser of 0 fixes the transitive out- and in-neighbourhoods pointwise. So an isomorphism of circulants onto `R_n` is affine, i.e. a multiplier.
* **Minimum.** For prime `q ≡ 3 (mod 4)`, `P_q` is the unique minimiser.
  * Cauchy–Schwarz on `Σ_(j≠0)|λ_j|⁴`, with `Σ_(j≠0)|λ_j|²` fixed, gives the bound, with equality iff the spectrum is flat.
  * A flat spectrum forces `λ_j ∈ Q(√−q)`, so `S` is `QR`-invariant, i.e. `S ∈ {QR, NQR}`.
* **FINITE-EXACT, odd `7 ≤ n ≤ 25`.** The energy minimum is unique except at `n = 15` (four classes) and `n = 21` (two).
* THM-140 (kind-pasteur S58) already had "interval always has maximum E, Paley always has minimum E", with `E(Int) = m(2m² + 1)/3` and `E(QR) = m(m² + 1)/2`.

So statement (1) reads: **the two energy extremes are `B_q` flips apart.**

**(6c) Which extreme wins** (FINITE-EXACT; complete tables for odd `n ≤ 25`; `energy_peaks_table.out`).
* The global `H`-maximum is the energy minimum for `n ≤ 11` and `R_n` for `13 ≤ n ≤ 25`.
* For `17 ≤ n ≤ 25` every energy minimiser is also a strict local maximum of `H`, a second peak.
* This is THM-140's `E`–`H` sign reversal, now located on the whole tables.

## 7. Mechanism by disjoint packings (FINITE-EXACT, new)

Write `H = Σ_k 2^k α_k`, where `α_k` is the number of sets of `k` pairwise vertex-disjoint odd cycles (the OCF). The code `alpha_decomp.c` reproduces every `H` in the tables.

* **Per-`k` champions, over all classes** (`alpha_all.out`, `alpha_all_21.out`; every `H` checked).
  * For `n = 17, 19, 21` an energy minimiser is the unique maximiser of `α_1` and `α_2` (at `n = 21` it is `Z_21{1,2,3,4,6,9,10,13,14,16}`).
  * `R_n` is the unique maximiser of every `α_k` with `3 ≤ k ≤ ⌊n/3⌋`, except the perfect triangle factors `α_(n/3)` when `3 | n` (`n = 15, 21`), which are maximised elsewhere.
  * For `n = 13, 15`, `R_n` uniquely maximises `α_3` and `α_4`.
  * The ratio `α_k(R_n)/α_k(flat)` grows with `k`: at `n = 19` it is 1.043, 1.126, 1.317, 1.846 for `k = 3..6`. The `H`-contributions `2^kΔα_k` themselves shrink after `k = 3`.
* **The split of `ΔH = H(R_n) − H(flat)`**, flat = the energy minimum (at `n = 21`: `Z_21{1,2,3,7,10,12,13,15,16,17}`, one of the two energy minima, the other giving `ΔH = 1,834,373,228,646`; `n = 15`, with four energy minima, is omitted) (`alpha_decomp_*.out`):

| `n` | `2Δα_1 + 4Δα_2` | `Σ_(k≥3) 2^kΔα_k` | `ΔH` |
|---|---|---|---|
| 9 | −162 | +72 | −90 |
| 11 | −4,620 | +2,552 | −2,068 |
| 13 | −90,636 | +98,800 | +8,164 |
| 17 | −189,570,298 | +318,784,408 | +129,214,110 |
| 19 | −15,235,138,144 | +26,752,215,992 | +11,517,077,848 |
| 21 | −831,973,445,620 | +2,686,870,223,200 | +1,854,896,777,580 |

* **The switch at 13.** Without `α_4`, `R_13` would lose (`+81,120 < 90,636`). Four disjoint odd cycles need 12 vertices, so the switch from the flat class to `R_n` comes exactly when they fit.
* **Corrections to the March reading** (MISTAKE-660).
  * THM-137's "the interval's cycles are more independent (higher `α_2`)" holds only for `n ≤ 11`. From `n = 13` the interval has fewer pairs and wins through `α_k`, `k ≥ 3`.
  * HYP-474 ("Paley maximises the `k`-cycle count at every odd `k`") fails at `k = p = 19, 23`: `c_19(R_19) = 34,841,356,485 > c_19(P_19) = 34,358,763,933`, and `R_23` has the most Hamiltonian cycles of any circulant class at `n = 23`. For every odd `5 ≤ k ≤ p − 2` Paley is still the unique maximiser.
* Independently verified in the second audit pass: `c_k` and the Hamiltonian-cycle count by Held–Karp and by inclusion–exclusion (all classes, `n = 19, 23`); `α_k` by explicit enumeration (`n ≤ 13`) and by ranked zeta transforms (`n ≤ 21`) (`05-knowledge/results/lattice_bridges_audit_pass2_20261022.md`).
