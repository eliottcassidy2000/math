# Lattices, codes, class numbers and tournaments (2026-10-21, mac-mini-2026-10-21-fanopsp)

**Owner prompt.** Merge the following into the tournament program:
* earlier work on 1729, the 24 of string theory and `q^(1/24)`, and Smale's h-cobordism;
* the Heegner numbers 1, 2, 3, 7, 11, 19, 43, 67, 163;
* the 24 Niemeier lattices and sphere packing, the deep holes of Leech, and the analogue in 8 dimensions;
* `II_{25,1}` (25 = 24 + 1, Lorentzian against Euclidean);
* `E_8`: 248 = 8 + 240, and its Langlands self-duality;
* the modular forms `Δ = (E_4³ − E_6²)/1728`, `E_8 = E_4²`, `E_10 = E_4E_6`, `E_14 = E_4²E_6`;
* `{7, 21}` against `{8, 24}`;
* again, why repunits are rare, and the earlier repunit work;
* four attached preprints (Niemeier genus 5, ADE moment rigidity, two on dimension 96).

**Three sessions received this prompt.** This note is the mac-mini lane and also an index to the other two.

| session | lane | where |
|---|---|---|
| mac-mini fanopsp (this note) | class number = Paley–rotational overlap; circulant `H`-landscape; Paley → QR codes → `E_8`/Golay; `{7,21} → {8,24}`; Heegner/Euclidean; `E_8` via `A_8`; an independent check of the four PDFs (§11) | THM-4751, HYP-9362 (pre-audit); scripts in `04-computation/experiments/lattice_bridges_20261021/`; errata MISTAKE-659 |
| opus S52 (MSG-1913; lane split agreed in MSG-2403) | `II_{25,1}` score lifts; Kervaire–Milnor Mersenne factors; 1728/1729 and eta tests; a different set of three PDFs | THM-4752, THM-4753, HYP-9363 (reserved) |
| Session27 | repunit recovery; typed modular/lattice synthesis; the Lorentzian score bridge; 1729 lanes; Smale source check; full review of the same four PDFs | `05-knowledge/results/session27_repunits_lattices_20261022/` (audited; MISTAKE-658) |

**Status words** follow the repo: PROVED, KNOWN (classical), FINITE-EXACT, DICTIONARY (exact facts on both sides, no map), NUMEROLOGY (no mechanism; recorded so it is not rediscovered). The repo policy is MISTAKE-559: a coincidence of small numbers is NUMEROLOGY unless a map carries the mechanism.

## 0. Summary

| item | status | where |
|---|---|---|
| the circulant `H`-landscape has a rotational peak and a Paley peak; their flip distance is `B_q = ((q−1)/2 − (2 − (2/q))h(−q))/2` | PROVED (distance) + FINITE-EXACT (peaks, `n ≤ 25`, `q ≤ 31`) | THM-4751, §1 |
| the `H`-maximal circulant is within one flip of rotational, and this would give HYP-9357's circulant "only if" | FINITE-EXACT (`n ≤ 25`; `R_n` a strict local maximum to `n = 31`) + conjecture; refines HYP-480 | HYP-9362 (a), THM-4751 (4) |
| the two peaks are the two extremes of additive energy; the flat class wins single cycles and pairs, the progression every packing of `≥ 3` disjoint odd cycles; HYP-474 fails at `p = 19, 23` | PROVED (extremes; prior art THM-140, THM-162, Savchenko) + FINITE-EXACT (packing law) | §13, THM-4751 (6)–(7), MISTAKE-660 |
| `P_7`, `P_23`, `P_47` generate the extended QR codes `[8,4,4]`, `[24,12,8]`, `[48,24,12]`, hence `E_8` (theta `E_4`) and `N(A_1^24)` | KNOWN; proved in the repo as THM-481 (Kim–Solé 2008); re-checked here | §2 |
| `{7, 21}` are the point counts of the only extendable projective planes; one extension gives 8, three give 24, and the extensions are the Steiner systems behind `E_8` and Leech | NUMEROLOGY (MISTAKE-559); the chains are old (S80, S88, S100, the 2026-05-02 reflection), new is only the Hughes/Lam selection; the "Baer tower" reading HYP-1348 was refuted | HYP-9362 (d), §3 |
| `{3, 7, 11}` = norm-Euclidean `Q(√−q)` = `H`-maximal Paley = Kantor = Galois | NUMEROLOGY, except the one-flip lens (PROVED) and the packing threshold (FINITE-EXACT, §13) | §4 |
| scores and `II_{25,1}` | PROVED identity (Session27, credited to S52) + my cross-checks | §5 |
| `E_8` roots via `A_8`: `240 = 72 + 168` | numerical DICTIONARY | §6 |
| modular forms; 1729; Smale | Session27 and S52 lanes; pointers and cross-checks only | §7–§9 |
| repunit rarity | PROVED by opus S51 (THM-4748/4749, audited) | §10 |
| the four attached preprints | Session27's review, plus my additions | §11 |

## 1. The circulant H-landscape and the class number (THM-4751, HYP-9362)

**Peaks.** On circulant tournaments `Cay(Z_n, S)`, one move reverses a difference class `{k, n − k}`. With this move the number `H` of Hamiltonian paths has:
* a rotational peak `R_n = Cay(Z_n, {1, …, m})`, where `m = (n − 1)/2`;
* for prime `q ≡ 3 (mod 4)`, a Paley peak `P_q`.

**The Paley peak is strict.**
* All one-flip neighbours of `P_q` are isomorphic, because the residue multipliers act simply transitively on the difference classes.
* The margins `H(neighbour)/H(P_q)` are 0.9259, 0.9783, 0.9925, 0.9954, 0.9982 at `q = 7, 11, 19, 23, 31` (`paley_flip_margin.out`).

**Distance between the peaks.**
* It is `B_q = |NQR ∩ [1, m]|`, minimised over multipliers.
* By Dirichlet's class number formula, `B_q = (m − (2 − (2/q))h(−q))/2`. Checked for all 1135 primes `q ≡ 3 (mod 4)` with `7 ≤ q < 20000` (`class_number_distance.out`; the audit confirms by counting reduced forms).
* The overlap of the two tournaments is `1/2 + (2 − (2/q))h(−q)/(q − 1)` of all arcs.
* Prior art: the half-row imbalance (`ogg_triangular_chowla_20260927.md`) and its linear-order companion THM-4666 (ranking excess `p·h(−p)`); the orientation-cube program THM-134/135/137/139/140 and HYP-480 (Paley against interval, March 2026; THM-139 already notes that no single flip improves `P_19`).

**Which peak wins.**
* `B_q ≤ 1` exactly for `q = 3, 7, 11`. For `q ≤ 31` these are exactly the primes at which `P_q` is the circulant `H`-maximum (FINITE-EXACT; conjecturally all `q`, HYP-9357).
* The Heegner primes 19, 43, 67, 163 have `h = 1`, so for their size they are as far from `R_q` as possible: `B_q = 3, 9, 15, 39` (so are `q = 23, 31`, where `h = 3` and `(2/q) = 1`).
* On odd `n ≤ 25` the global maximum is always `R_n` or a one-flip neighbour of it. If that holds for every `n` (HYP-9362 (a)), Paley is never the circulant maximum for `q ≥ 19`.

**Reading.** The class number measures how far the Paley peak sits from the rotational peak. Paley wins only when the two peaks touch.

## 2. Paley tournaments, codes and lattices (KNOWN; THM-481)

**When it works.** Take `q ≡ 7 (mod 8)`. Among `q ≡ 3 (mod 4)` these are exactly the primes with `(2/q) = +1`, the Camion multiplier of THM-4746.
* The out-neighbourhoods of `P_q` span the binary QR code. For `q = 7` they are the 7 Fano lines.
* Its parity extension is doubly even and self-dual.
* THM-481 (June 2026, claudebox and kind-pasteur) proved this for all such `q`, rediscovering Kim–Solé 2008.

Re-checked here (`paley_codes.out`):

| `q` | extended code | weight enumerator | lattice |
|---|---|---|---|
| 7 | `[8,4,4]` | `1 + 14y⁴ + y⁸` | Construction A gives `E_8`, theta `E_4 = 1 + 240q + 2160q² + 6720q³ + 17520q⁴` |
| 23 | `[24,12,8]` | `1 + 759y⁸ + 2576y¹² + 759y¹⁶ + y²⁴` | Construction A gives the Niemeier lattice `N(A_1^24)`, theta `E_4³ − 672Δ = 1 + 48q + 195408q² + …` |
| 31 | `[32,16,8]` | 620 weight-8 words | — |
| 47 | `[48,24,12]`, extremal | 17296 weight-12 words | — |

**Leech.**
* The Leech lattice comes from the Golay code by Conway's twisted construction, not by plain Construction A (Session27 makes the same point).
* Its 196560 minimal vectors are `4·C(24,2) + 759·2^7 + 24·2^12`. Its theta series is `E_4³ − 720Δ`, and its character is `j − 720`.
* The repo's code–modular dictionary THM-489 has `W_g24 = W_ê8³ − 42P_24`, with `P_24 ↦ 16Δ`.

**Higher dimensions.**
* Nebe's extremal 72-dimensional lattice is built over `Z[(1 + √−7)/2]` (HYP-2433..2435; THM-487 "the 72 frame"). That ring is the spectral ring of `P_7`.
* The two attached 96-dimensional preprints also work over it (§11).
* The `[72,36,16]` code remains the open "72 problem". The `P_71` gauge gives `d = 12`.

## 3. {7, 21} → {8, 24}: extendable projective planes (NUMEROLOGY; HYP-9362 (d))

**Facts.**
* A projective plane of order `n` extends to a 3-design only if `n ∈ {2, 4, 10}` (Hughes 1965). There is no plane of order 10 (Lam–Thiel–Swiercz 1989).
* So the extendable planes are `PG(2,2)`, with 7 points, and `PG(2,4)`, with 21 points.
* `PG(2,2)` extends once: `PG(2,2) → S(3,4,8)` (`7 + 1 = 8`). `S(3,4,8)` is the set of 14 weight-4 words of `[8,4,4]`, whose Construction A is `E_8`. It extends no further, since there is no `S(4,5,9)`.
* `PG(2,4)` extends three times: `PG(2,4) → S(3,6,22) → S(4,7,23) → S(5,8,24)` (`21 + 3 = 24`). The 759 octads of `S(5,8,24)` are the weight-8 Golay words from which Leech is built (Witt; `M_21 = PSL(3,4)` inside `M_24`).

**The tournament side.** `H = 7` and `H = 21` are never taken by `H` (THM-029, THM-079, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, PROVED). That no other odd value is missed is a CONJECTURE (FINITE-EXACT below 80,265, THM-4749).

**Status.**
* No map is known, so by MISTAKE-559 this is NUMEROLOGY.
* The chains themselves are old: opus S80 (`fano_hamming_forbidden.out`), opus S88 (`golay_baer_bridge_88.py`, "7 → 21 → 24"), kind-pasteur S100 (`projective_geometry_deep.out`, "PG(2,4) + 3 = 24"), the 2026-05-02 reflection, and concurrently opus S52 (also NUMEROLOGY). New here is only the Hughes/Lam selection.
* Session27 noted that the proposal "still owes a map", and that the QR length rule `q → q + 1` pairs 7 with 8 but 23, not 21, with 24. The extendable-plane reading is the pairing that does send 21 to 24. It is selected by a classification theorem (Hughes, Lam), not by fitting.
* Unlike the refuted "Baer tower" `|PG(2, 2^(2^k))|` of HYP-1348 (273 is attained), this reading predicts nothing at `q = 16`.
* The unaudited reflection `sphere-packing-eight-twentyfour.md` had both chains with several errors. They are listed in its errata banner (MISTAKE-659).

**What would make it mathematics.** A map from the OCF obstruction behind `H ∉ {7, 21}` to the extendability of `PG(2,2)` and `PG(2,4)`. Both point counts are `Φ_3` values: `7 = Φ_3(2)` and `21 = Φ_3(4)`. But `13 = Φ_3(3)` and `31 = Φ_3(5)` are attained, and `43 = Φ_3(6)` too (HYP-1317). So the selection must use extendability or the prime-power order, not `Φ_3` alone.

## 4. Heegner and Euclidean (KNOWN facts; alignment NUMEROLOGY; THM-4751 (5))

**The spectral ring.** The nontrivial eigenvalues of `P_q`, `(−1 ± √−q)/2`, generate `O_K` with `K = Q(√−q)`.
* It is a PID exactly for the Heegner primes `3, 7, 11, 19, 43, 67, 163`.
* It is norm-Euclidean exactly for `3, 7, 11`. The deep hole of `Z[(1 + √−q)/2] ⊂ C` has squared distance `(q + 1)²/(16q)`, which is `< 1` exactly then: a covering-radius condition in two dimensions.

**Characterizations of `{3, 7, 11}` among primes `q ≡ 3 (mod 4)`:**
1. `P_q` is the circulant `H`-maximum (FINITE-EXACT for `q ≤ 31`; among all tournaments known only for `q ≤ 11`, A038375; conjecturally exactly `{3, 7, 11}`, HYP-9357);
2. the Paley design is 2-transitive (Kantor);
3. Galois' exceptional primes;
4. `P_q` is within one flip of `R_q` (PROVED; THM-4751 (2));
5. `Q(√−q)` is norm-Euclidean (classical);
6. `H(P_q)/|Aut P_q| − 1` is a cube: 0, 8, 1728. It is not a cube at 19, 23, 27, 31 (`numerology_checks.out`; S52's 1729 lane has the same test).

Only (4) and the packing threshold of §13 come with a tournament mechanism. (4) is conditional on HYP-9362 (a). The threshold says the flat class wins while the OCF stops at `α_3`, i.e. `n ≤ 11`, and is FINITE-EXACT.

**Caution.** No prime `≡ 3 (mod 4)` lies between 11 and 19, so any threshold in `(11, 19)` gives this set. The repo has rejected Heegner-based LRC and chromatic-number claims (MISTAKE-229, MISTAKE-230, MISTAKE-568, MISTAKE-576); the historical entries HYP-1202/HYP-1220 ("Heegner numbers are tournament constants") remain CONFIRMED numerology. Session27 confirms that the owner's nine-value Heegner list, including 1, is the correct radicand list.

## 5. Score vectors, permutohedra and II_{25,1} (Session27 + S52 lane; my cross-checks)

**Primary.** Session27, `LORENTZIAN_SCORE_BRIDGE.md` (credited to S52, MSG-1913); S52, THM-4752 (reserved).
* For a 25-vertex tournament, `v(T) = (s | 70) ∈ II_{25,1}` has norm `−2c_3(T)`.
* It is null exactly for the transitive tournaments.
* Every transitive lift is a coordinate permutation of the Weyl vector `w = (0, 1, …, 24 | 70)`, whose null quotient `w^⊥/Zw` is the Leech lattice.

**Answer to THM-872's open question.**
* THM-872 asks which Niemeier lattice, if any, the scores select at the 24-rung. Its Euclidean ladder `E_8 = D_8⁺`, `D_16⁺`, … only gives a sum-zero slice, which is not unimodular.
* The Lorentzian lift at `n = 25` does select one: the only null score lifts are the transitive ones, and they all give Leech.
* This is a reading of Session27's proved identity, not a new theorem.

**My cross-checks** (`score_vectors_forests.out`, and the verifier's `f_checks.out`).
* Labelled score vectors are the lattice points of the permutohedron `conv{σ(0, …, n−1)}`, and their number is the number of labelled forests on `n` vertices: 1, 2, 7, 38, 291, 2932, 36961, 561948 for `n ≤ 8`. This is Landau's theorem plus Stanley's count.
* The permutohedron is the Voronoi cell of `A_{n−1}^*`, checked exactly for `n ≤ 6`.
* Reversing an arc `u → v` moves the score vector by the `A_{n−1}` root `e_v − e_u`.
* For every `n`, `c_3(T) = (n³ − n)/24 − |s − s̄|²/2`. So `c_3` is a paraboloid on the permutohedron, zero exactly at its vertices. At `n = 25` this is the Lorentzian norm identity, since `−2c_3 = Σs_i² − 4900`.
* By the cannonball problem, `Σ_(k≤N) k²` is a square only for `N = 1, 24`. Checked for `N < 10^6`; Watson proved it.

## 6. E_8 through A_8 (numerical DICTIONARY)

`E_8 = A_8 ∪ (A_8 + [3]) ∪ (A_8 + [6])`. Its roots are:
* the 72 roots `e_i − e_j` of `A_8`: the ordered pairs of 9 points, i.e. the arcs of `K_9` in both directions;
* the `2·C(9,3) = 168` vectors `±(e_i + e_j + e_k − (1/3)Σe)`, two for each 3-subset.

So `240 = 72 + 168` and `248 = 80 + 84 + 84` (`su(9) ⊕ Λ³ ⊕ Λ⁶`).

**Tournament reading.** The 72 are the arcs of the complete digraph on 9 vertices, and 168 equals the number of cyclic triangles counted with both orientations. The bijection is only numerical: `S_9` preserves the two glue cosets but permutes triangle orientations transitively. `168 = |PSL(2,7)|` is a further coincidence.

**Prior repo results.**
* THM-868 (`THM-868-e8-bridge-score-lattice.md`): the score half-vectors of 8-tournaments lie in `E_8`, with the 70 sum-zero roots as the near-regular shell.
* Session27's `SCORE_ROOT_COLLISION.md`: reversing a directed triangle of `P_7` keeps the scores but moves `H` from 189 to 171. So the `E_8` point does not determine `H`.

**The 8-dimensional analogue of the deep holes.** Session27, `PRIOR_1729_AND_E8_HOLES.md` §B:
* the null quotient of `II_{9,1}` is always `E_8`, with one cusp;
* the `E_8` deep hole is a 16-vertex cross-polytope;
* in rank 24, by contrast, there are 24 Niemeier cusps and 23 deep-hole types.

**Self-duality.** `E_8` is its own Langlands dual (simply laced, trivial centre) and its root lattice is unimodular. Session27 warns that Langlands duality and lattice duality are different operations. The closest tournament notion is self-converse (e.g. `P_q` via `x ↦ −x`). That is ANALOGY only.

## 7. Modular forms (pointers)

**Identities.**
* `Δ = (E_4³ − E_6²)/1728 = q − 24q² + 252q³ − 1472q⁴ + …`. `M_8`, `M_10`, `M_14` are one-dimensional, so `E_8 = E_4²`, `E_10 = E_4E_6`, `E_14 = E_4²E_6`.
* `θ_(E_8) = E_4`. `θ_(E_8⊕E_8) = θ_(D_16⁺) = E_4²`: Milnor's isospectral tori (THM-480, and THM-872 "the metagraph hears one of the two drums").
* Niemeier: `θ = E_4³ + (24h − 720)Δ`. Leech: `E_4³ − 720Δ`.
* `Δ = η²⁴`, with `η = q^(1/24)∏(1 − qⁿ)`.
* The normalizations are in Session27's `session27-modular/report.md`.

**Proved tournament–modular links in the repo:**
* THM-4652/HYP-9274: the Paley wheel `W_4` sum is governed by `η(2z)⁴η(4z)⁴` (weight 4, level 8). Its exponent `(2·4 + 4·4)/24 = 1` is the `q^(1/24)` bookkeeping. The same eta product is Mason's `M_24` class-4A form (THM-4669).
* HYP-9275: the Paley run sums through 12 steps lie in `Z[q] + Z[q]a_q`.
* THM-4653 (PROVED reduction) and HYP-9276 (FINITE-EXACT): the `K_6 − P_6` count is the adjoint of the weight-3 level-32 newform.

## 8. 1729 (Session27 + S52 lanes; pointers)

* **Session27** (`session27-modular/PRIOR_1729_AND_E8_HOLES.md`) separates three earlier lanes:
  * Ramanujan's near-miss recurrence (an exact all-`n` identity);
  * the Ono–Trebat-Leder K3 surface (cited);
  * the tournament ratio `H(P_11)/|Aut P_11| = 95095/55 = 1729` (HYP-9028). The quotient counts free orbits of Hamiltonian paths.

  It also notes that 1729 is a Carmichael number and hence a Camion pseudoprime (THM-4746 (5)).
* **Earlier refutations.** "1729 = j(i) + 1" as a modular law was refuted as HYP-2306. The taxicab property is a theorem about Eisenstein norms (THM-463).
* **Cross-check.** The cube pattern stops after 11 (§4 (6)). `|M_12| = 1728·55` gives `H(P_11) = |M_12| + |Aut P_11|`, but the analogue fails at `q = 3` and `q = 23`: NUMEROLOGY.

## 9. Smale, h-cobordism, exotic spheres (Session27 + S52 lanes; corrected reading)

**Primary.** Session27, `session27-smale/ADDENDUM.md`, checked against Milnor 2011; MISTAKE-563 corrected the old `bP` script.
* The odd parts of `|bP_(4m)|` contain `2^(2m−1) − 1 = 7, 31, 127, 511, 2047, …`.

**Against THM-4749.** Every one of these with exponent `≤ 40` is forbidden (7) or has degree at most that of its Fermat partner `2^(2m−1) + 1`:
* by degree;
* at `k = 11` by leading coefficient;
* at `k = 23, 29` the degrees tie and the coefficients are not yet decided.
* For odd `k ≥ 17` other than 21 and 25 the degrees are conditional on `V`-membership of large prime factors (THM-4749 status).

**Why.** The tournament reason is the triangle: `3 = Φ_2(2) = H(C_3)` divides `2^k + 1` for odd `k`.
* The even exponent `63 = 2^6 − 1` is rare as well (`Φ_3(2) = 7` divides it), though it is not a `bP` factor.
* `15` and `255` are commoner than `17` and `257`.

So the `bP` list is the odd-exponent half of a pattern governed by the cyclotomic factors of `2^k ± 1`. Topology adds no mechanism (NUMEROLOGY).

**Correction.** An earlier draft of this note said 511 was typical. That is wrong: THM-4749 makes 511 rarer than 513.

## 10. Why repunits are rare (opus S51, audited; and the catalogs)

* **THM-4748.** For `k ≥ 5`, every tournament with `H = 2^k − 1` needs at least three reversed arcs from transitive, while `2^k + 1` needs one. The bound is attained by the Fibonacci fan.
* **THM-4749.** The number of tournament classes with `H = h` grows like `n^w(h)`, where `w(h)` is the longest factorization of `h` into strong values. Strong values are all odd numbers `≥ 3` except 7 and 21, FINITE-EXACT below 80,265. Read through `2^k ± 1 = ∏ Φ_d(2)`:
  * `2^k − 1` is the rarer value for `k = 5, 6, 7, 9, 13, 15`;
  * it is the commoner value for `k = 4, 8, 16` (Fermat products) and `k = 12` (`4095 = 3²·5·7·13`, the 7 absorbed into 35);
  * the degrees tie at `k = 10, 11, 14`; THM-4749 decides `k = 10` (1023 commoner) and `k = 11` (2047 rarer); `k = 14` is not computed.
* **HYP-9359.** The arithmetic bias comes from divisibility by 3 (the triangle) and from order-3 automorphisms.
* **Catalogs of earlier repunit work:**
  * mine, about 45 threads (MSG-2399, previous session);
  * Session27's measure-aware recovery (`session27-repunits/REPORT.md`, with bounded remote search).

  The unifying reading is that a base-`b` repunit truncates the fixed point `−1/(b−1)` of `x ↦ bx + 1`.
* **Lattice side (ANALOGY only).** `2^k − 1 = |PG(k−1, 2)|` is the length of the Hamming codes, the nontrivial perfect binary codes besides the Golay code of length 23 (up to parameters; Tietäväinen, van Lint 1973).

## 11. The four attached preprints ("SRFP311T1 collaboration", 2026)

**Primary referee report:** Session27, `session27-papers-review/REVIEW.md`, with independent lattice-controls and Niemeier audits. Its decisive points:
* **Niemeier genus 5.** The hyperelliptic nonvanishing claims contradict Codogni (Doc. Math. 21 (2016), Thm 6.1). Equal root counts force the theta differences to vanish on the hyperelliptic locus in every genus, for exactly the five equal-`h` pairs.
* **Extremal 96.** The `q⁶` coefficient is wrong. The Lagrangian count uses the symplectic formula, though the form is orthogonal in characteristic 7. The no-short-vector certificate was not reproduced, and the later paper retreats to "minimum in {2, …, 10}".
* **ADE compression.** The rank-43 four-moment threshold and its unique collision were reproduced by an independent enumeration, with an explicit finite-difference mechanism.

**This session's independent verification** (scripts and outputs in `04-computation/experiments/lattice_bridges_20261021/pdf_verify/`) agrees with all of these. Items marked **[new]** are not in Session27's review.

### Extremal 96 ("An Extremal Even Unimodular Lattice in Dimension 96")

* **The `q⁶` coefficient** of the unique weight-48 form with `q¹..q⁴` zero is `2,972,108,280,960,000`, computed in Python and with PARI `mfinit`. This paper prints `184,139,089,776,640`; the dual-horizon paper is right. `q⁵ = 565,866,362,880` in both. (`a_extremal_theta*.out`)
* **[new] The base lattice is not of minimum 28** (FINITE-EXACT, inside a reconstruction).
  * The paper's twisting unit `u = −u_4u_8u_16 ∈ O_F^×` is undefined, since `O_F^× = {±1}`. The verifier therefore built `J` from all 128 totally positive twists, and all 128 give `7·`Leech. The four examined in full (candidates 0–3) satisfy every property the paper states for `J` (Prop. 2.2).
  * In two of these reconstructions, the base lattice `E_8 ⊗_(O_F) J` (96-dimensional, even, positive definite, `det = 7^48`, minimum diagonal 28) contains the explicit vector `x_1⊗y_1 + x_2⊗y_2` of norm **24**. Here `x_i` are roots and `y_i` are minimal vectors, with `h_E(x_1,x_2) = (1 − 3√−7/7)/2` and `h_J(y_1,y_2) = −7 − 3√−7`, so the norm is `56 + 4·Re(bc) = 24`.
  * I re-verified the integer witness myself (`w24_independent_check.out`).
  * This is an explicit witness for Session27's logical objection: Prop. 3.2 bounds only pure tensors. It does not touch the extremality of the glued lattice, since `24 > 10`.
* **[new] The union bound is mis-scaled** (HEURISTIC).
  * The "≈ 10^12 forbidden vectors" matches the Gaussian count at norm `≤ 4`, not at norm `≤ 8`.
  * At norm `≤ 8` a random isotropic gluing would keep about `10^6` short vectors, not `5·10^−9`.
  * For scale: by Siegel–Weil, `E_48` gives a mass-weighted average of about 21 norm-6 and `1.6·10^7` norm-8 vectors per even unimodular 96-dimensional lattice. An extremal lattice needs both to be 0. (`e_checks.out`, `a_extremal_theta.out`)
* **Lagrangian count.** The form is symmetric in characteristic 7, so the hyperbolic count is `∏_(i=0)^(23)(7^i + 1)`, not the symplectic `∏_(i=1)^(24)(7^i + 1)`. The verifier brute-forced it over `F_3, F_5, F_7`. The probability that a line lies in a random Lagrangian is `1/(7^23 + 1)`.
* **[new] Smaller errors.**
  * The Galois description `(Z/35)^×/{±1, ±2, …}` is garbled. The subgroup fixing `Q(√−7)` is `{k : (k/7) = 1}`, cyclic of order 12, generated by 2. It does not contain −1.
  * Both 96-dimensional papers cite Nebe's `Γ_72` as Ann. of Math.; it is J. reine angew. Math. 673 (2012) 237–247.
  * The minimum bound is Mallows–Odlyzko–Sloane, J. Algebra 36 (1975).
  * The `E_8 ⊗ J` design mirrors Sakamoto, arXiv:2609.15633: two extremal rank-88 lattices from a Hermitian tensor product over `Q(√−23)`. It is not cited.

### Dual horizon 96

* The explicit `H_E` is `E_8` as a rank-4 Hermitian `O_F`-lattice: the trace form is even, unimodular, with 240 roots and 2160 norm-4 vectors.
* The Schur identity `7A − RR^T = I_2` holds.
* The list of 12 Galois elements is right.
* Session27 notes that this paper leaves the minimum open in `{2, …, 10}`, and that its matrix is not the old one. Dimension 96 remains open.

### Niemeier genus 5

* Table 2.1, Theorem 2.8 and all four values of Theorem 4.11 are reproduced by direct frame counts. So are all the `A_5`, `D_5`, `D_4` multiplicities, including `E_7(A_5) = 1344 = 336 + 1008`.
* **[new] Genus 3 and 4.** On root-type coefficients:
  * at genus 3, equal-`h` pairs agree;
  * at genus 4, `F_h = (m_h/896)·F_30` on all 11 rank-4 classes, with `m = (1, 9, 20, 105, 896) = Δm(D_4)/5`;
  * `a(I_4)`, the number of ordered orthogonal root 4-frames, takes 24 distinct values, so it separates all Niemeier lattices as claimed.

  Non-root coefficients were not tested.
* **Normalization conflict.** `F_30^(4) = E_4·(Θ_(E_8²) − Θ_(D_16⁺))` exactly, so the companion paper's `J_8` normalization is off by the factor 3360. Session27 found the same.
* **Decisive.** Session27's Codogni contradiction: the hyperelliptic claims fail.

### ADE compression

* **Census.** It matches: 5,565,832 root systems of rank `≤ 42` with injective signatures. At rank 43 there are 1,556,185 systems and exactly one collision, `R_1 = 2A_1 + 2A_4 + A_5 + 2A_7 + A_8 + D_6` and `R_2 = 2A_2 + A_3 + 2A_6 + A_9 + 3D_5`. This makes three independent enumerations: the authors', Session27's and this session's. (`b_paper2_dp.out`)
* **[new] What the four moments are** (PROVED).
  * The paper defines `c_m` only by printed polynomials (p. 3).
  * The verifier proved that `c_1..c_4` are the first four coefficients of `log Σ_j O_j z^j/j!`, where `O_j` counts ordered `j`-tuples of mutually orthogonal roots. This holds for all ranks, by polynomial degree; `O_m` was brute-forced to `m = 6`.
  * So the theorem says: the numbers of orthogonal root frames of sizes 1–4 determine an ADE root system up to rank 42, and first fail at rank 43.
* **[new] The printed `c_5` is not the fifth such coefficient** for any of the 30 types tested.
  * Example: for `A_1`, `log(1 + 2z)` gives `32/5`, but `64/5` is printed.
  * True formulas: `l_5(A_k) = (16/5)k(k+1)(42k⁴ − 176k³ + 326k² − 296k + 105)` and `l_5(D_k) = (32/5)k(k−1)(672k⁴ − 5784k³ + 19301k² − 29489k + 17295)`.
  * For the rank-43 pair, `Δl_5 = 1,082,880` and `ΔO_5 = 129,945,600`, not the printed `158,592`. The paper's `ΔN_5 = 480·Δc_5` has no basis. The pair still separates at order 5. (`b_paper2_formulas.out`, `b_O5_bruteforce.out`)
* **[new] The weighted 15-design claim is correct** (PROVED). A degree-12 harmonic theta series of Leech is a multiple of `Δ²`, which forces `w_3/w_2 = (6/4)^6/48 = 243/1024`. (`d_leech.out`)

### The tournament bridge in the compression paper (DICTIONARY with exact structure)

The frame series of a root system is an independence polynomial.
* A set of pairwise orthogonal roots is an independent set of the non-orthogonality graph `N_R` (`r ~ s` iff `⟨r, s⟩ ≠ 0`).
* So `Σ_j O_j z^j/j! = I(N_R, z)`, and the paper's cumulants are the coefficients of `log I(N_R, z)`.
* For an orthogonal direct sum, `N_(R ⊕ R')` is the disjoint union of `N_R` and `N_R'`. So `I` multiplies and `log I` adds; this is why the `c_m` are additive.

The tournament program has the same structure, with:
* odd cycles in place of roots;
* vertex-disjointness in place of orthogonality;
* the chain of strong components in place of the orthogonal decomposition.

Then the OCF gives `H(T) = I(Ω(T), 2)`, and `Ω(T)` is the disjoint union of the `Ω(T_i)` over the strong components, so `H` is multiplicative along the chain.

**What the parallel says.**
* The ADE theorem uses four Taylor coefficients of `log I`.
* `H` keeps one evaluation of `I`, at `x = 2`. The multiplicity counts of THM-4749 measure exactly how much that one evaluation forgets.

**Question.** How many coefficients of `log I(Ω(T), x)` are needed to recover the multiset of strong components, within a parametrized family of strong tournaments? This would be the tournament version of the rank-43 threshold.

## 12. What would turn the dictionaries into mathematics

* Prove the packing law HYP-9362 (e), or the one-flip conjecture HYP-9362 (a), or at least that `R_n` is always a strict local maximum. That alone does not exclude Paley, but it is the first step. The complete class table at `n = 25` is in §13.
* Find a map from the OCF obstruction behind `H ∉ {7, 21}` to the extendability of `PG(2,2)` and `PG(2,4)` (§3).
* Find a tournament statistic whose value at `P_q` is the covering radius of `O_(Q(√−q))` (§4).
* Compute the tournament version of the rank-43 threshold (§11, bridge paragraph): within a parametrized family of strong tournaments, how many coefficients of `log I(Ω, x)` determine the strong-component multiset?

## 13. The two peaks are the two extremes of additive energy (extremes KNOWN/PROVED, THM-140; landscape and packings FINITE-EXACT to n = 25)

**Setting.** Let `E(S) = #{(a, b, c, d) ∈ S⁴ : a + b ≡ c + d (mod n)}` be the additive energy of the connection set.
* **Exact relations.** On circulant tournaments, `c_4 = n(E − m)/4` and `2c_5 + nE` is constant. They follow from `tr A^k = k c_k` for `k = 4, 5` (THM-118), the normality of `A`, and `Σ_j |λ_j|⁴ = nE(S)`. Checked over every class table (`spectral_peaks.out`).
* So maximising `E` is the same as maximising the 4-cycle count `c_4` and the fourth spectral moment, and the same as minimising `c_5`.
* **Prior art.** The relations are the circulant case of Savchenko's identity `c_5 + 2c_4 = n(n²−1)(n²−9)/160` for every regular tournament (J. Graph Theory 83 (2016) 44–77; Prop. A.1 of arXiv:2403.19555); new is `c_4 = n(E − m)/4` on circulants. The energy extremes and the two-mechanism reading are THM-140 (kind-pasteur S58: interval maximum and Paley minimum of `E`, `E`–`H` sign reversal) and THM-137; the Paley end is THM-162 and HYP-464.

**The two extremes.**
* **Maximum.** For prime `p`, Lev's rearrangement theorem (Duke Math. J. 107 (2001)) says arithmetic progressions maximise the number of solutions of linear equations with variables in sets of given size. So `R_p` maximises `E` among all half-sets. It is the unique maximiser up to multipliers for every odd `n ≥ 5`: `RLT_n = R_n` uniquely maximises `c_4` among `n`-tournaments (arXiv:2403.19555 §1), `c_4 = n(E − m)/4`, and `Aut(R_n) = Z_n` makes isomorphic circulants multiplier-equivalent (FINITE-EXACT check `n ≤ 25`).
* **Minimum.** For prime `q ≡ 3 (mod 4)`, `P_q` minimises `E` uniquely, by Cauchy–Schwarz with equality exactly for a flat spectrum, and a flat spectrum forces `S` to be `QR`-invariant.
* **Consequence.** The class-number formula of §1 says that **the two energy extremes are `B_q` flips apart**.

**The landscape against energy** (complete class tables; `energy_peaks_table.out`; `n = 25` new here, `circ25_landscape.out`). Local maxima are counted in the hypercube convention, where a flip landing in the same multiplier class is a tie; THM-4751 (3) lists the five extra classes of the loopless class-graph convention.

| `n` | global `H`-maximum | energy minimum |
|---|---|---|
| 7, 9, 11 | the energy minimum (`P_7`, `Z_9{1,2,3,5}`, `P_11`) | the only strict local maximum (at `n = 11` the loopless class graph also has `Z_11{1,2,3,4,6}`) |
| 13, 15 | `R_n` = the energy maximum | not a local maximum (at `n = 15` one of the four energy minimisers, `Z_15{1,2,5,7,9,11,12}`, is strict in the loopless class graph) |
| 17, 19, 21, 23, 25 | `R_n` | a strict local maximum of `H`: a second peak, `P_19` and `P_23` at prime `q ≡ 3 (mod 4)`, `H`-rank 3, 3, 11–12, 8, 23 |

At `n = 25` (208 classes) the global maximum is `R_25`, `H = 2,418,453,569,285,650,675`. There are 6 strict local maxima, and the energy minimum `Z_25{1,2,3,4,5,7,9,10,13,14,17,19}` is the highest of the others.

**Rotational local-maximum tests** (`rot_local_max.py`; inclusion–exclusion, exact):
* `R_27`, `R_29` and `R_31` are strict local maxima. This completes the partial test pre-registered in HYP-9362 (a).
* In every case, including `n = 25`, the best neighbour flips the last class `{m, m + 1}`. Its ratio is 0.99585, 0.99613, 0.99641, 0.99667 at `n = 25, 27, 29, 31`, so `n(1 − ratio) ≈ 0.104` (`rot27/29/31_local_max.out`).
* The rotational peak's margin therefore decays like `1/n`. The Paley peak's margin decays at least quadratically: `(1 − margin)q²` falls from 3.63 to 1.72 (HYP-9362 (b)).

**Reading.** `H` on circulants trades off two kinds of structure:
* the pseudorandom end: flat spectrum, the most 5-cycles, minimum energy, Paley;
* the additive end: an arithmetic progression, the most 4-cycles, maximum energy, rotational.

Short odd cycles decide small `n`, and the progression wins from `n = 13` on. The class number measures the distance between the two ends.

**Mechanism: disjoint packings** (FINITE-EXACT; `alpha_decomp.c`, which reproduces every `H` in the class tables).
* Split `H = Σ_k 2^k α_k` (OCF), where `α_k` counts sets of `k` pairwise vertex-disjoint odd cycles.
* Against the energy minimum ("flat"), `R_n` has fewer odd cycles (`α_1`) and, from `n = 13` on, fewer disjoint pairs (`α_2`).
* But `R_n` has more packings of every size `k ≥ 3`, and the gap grows with `k`. At `n = 19`, for instance, `α_6` is 4,100,656 for `R_19` against 2,221,081 for `P_19`.

Split of `ΔH = H(R_n) − H(flat)` (`alpha_decomp_9_13.out`, `alpha_decomp_17_19.out`, `alpha_decomp_21.out`):

| `n` | flat class | `2Δα_1 + 4Δα_2` | `Σ_(k≥3) 2^k Δα_k` | `ΔH` |
|---|---|---|---|---|
| 9 | `Z_9{1,2,3,5}` | −162 | +72 | −90 |
| 11 | `P_11` | −4,620 | +2,552 | −2,068 |
| 13 | `Z_13{1,2,3,5,6,9}` | −90,636 | +98,800 (of which `α_4`: +17,680) | +8,164 |
| 17 | `Z_17{1,2,3,5,7,8,11,13}` | −189,570,298 | +318,784,408 | +129,214,110 |
| 19 | `P_19` | −15,235,138,144 | +26,752,215,992 | +11,517,077,848 |
| 21 | `Z_21{1,2,3,7,10,12,13,15,16,17}` | −831,973,445,620 | +2,686,870,223,200 | +1,854,896,777,580 |

**Champions of each packing number** (all classes, `n ≤ 19`; `alpha_all.out`):
* For `n = 17, 19` the energy minimum is the unique maximiser of `α_1` and `α_2`, and `R_n` is the unique maximiser of every `α_k` with `3 ≤ k ≤ ⌊n/3⌋`.
* For `n = 13, 15`, `R_n` is the unique maximiser of `α_3` and `α_4`, and an energy minimiser maximises `α_1`.
* At `n = 15` the triangle factors (`α_5`) are maximised elsewhere.

So the two peaks are the two champions: of few-cycle counts (pseudorandom) and of many-cycle packings (progression).

**The crossover.**
* At `n = 13` the rotational class wins only through `α_4`, since `+81,120 < 90,636` without it.
* So `R_n` overtakes the flat class exactly when four disjoint triangles fit, at `n ≥ 12`.
* Within these tables, the flat (pseudorandom) class is `H`-maximal exactly when the OCF stops at `α_3`, i.e. `n ≤ 11`.

**Consequences.**
* This gives characterization (1) of `{3, 7, 11}` in §4 a tournament mechanism: they are the primes `q ≡ 3 (mod 4)` with `deg I(Ω(P_q), x) = ⌊q/3⌋ ≤ 3`. It is FINITE-EXACT, and the small-numbers caveat still applies, since 12 lies in `(11, 19)`.
* **Corrections to the March reading (MISTAKE-660).**
  * HYP-474 ("Paley maximises the number of `k`-cycles at every odd `k`", CONFIRMED at `p = 7, 11`) fails at `k = p`: `c_19(R_19) = 34,841,356,485 > c_19(P_19) = 34,358,763,933`, and at `n = 23` `R_23` has the most Hamiltonian cycles of any circulant class (the `C` column of `circ_19.out`, `circ_23.out`; `alpha_decomp_17_19.out`).
  * THM-137's "the interval's cycles are more independent (higher `α_2`)" holds only for `n ≤ 11`. From `n = 13` the interval has fewer pairs and wins through `α_k`, `k ≥ 3`.
  * THM-161's "Paley maximizes every `c_k` (Savchenko's theorem for DRTs)" misreads Savchenko, whose results compare doubly regular tournaments with `RLT_n`.

**Reading.** A pseudorandom connection set maximises the number of cycles and of pairs. An arithmetic progression lays its odd cycles side by side around the circle, so it maximises the large disjoint packings. The OCF weights a packing of size `k` by `2^k`, so for `n ≥ 13` the progression wins.
