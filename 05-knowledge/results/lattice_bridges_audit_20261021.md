# Audit of THM-4751, HYP-9362, the lattices/codes note (§0–10, §13) and MISTAKE-659 (2026-10-21/22)

**Auditor.** Independent, refutation-first subagent (Claude), for the mac-mini-2026-10-21-fanopsp session. Work done only in `/Users/e/Documents/GitHub/math-wt-fanopsp`. No existing file was edited; no git operation was run.

**Audited versions** (sha1 at audit start):

| file | sha1 |
|---|---|
| `01-canon/theorems/THM-4751-circulant-h-landscape-paley-peak-at-class-number-distance.md` | `650f4b82` (= commit `f804c3e24d`) |
| `05-knowledge/hypotheses/HYP-9362-one-flip-conjecture-paley-peak-and-the-7-21-to-8-24-dictionary.md` | `985f5ba5` (= `f804c3e24d`) |
| `05-knowledge/results/lattices_codes_tournaments_20261021.md` | `d0f6429e` (§0–10 unchanged in the later `47145651`; §11–13 kept growing during the audit) |
| `07-reflections/sphere-packing-eight-twentyfour.md` (errata banner) | `b053ef11` |
| `01-canon/MISTAKES.md` (MISTAKE-659 entry) | `2c5bc60c` |

**Scripts and outputs.** All in `04-computation/experiments/lattice_bridges_20261021/audit/`; none reuses session code. `a3_hcount.c` and `a7_hspectrum.c` build with `cc -O3 -o a3_hcount a3_hcount.c` (binaries removed after use).

| script | output | content |
|---|---|---|
| `a1_class_number.py`, `a1_class_number.gp` | `a1_class_number.out`, `a1_class_number_gp.out` | `h(−q)` by counting reduced forms, all `q ≡ 3 (4)`, `7 ≤ q < 20000`; PARI `qfbclassno`, `quadclassunit` and the Dirichlet sum on every one |
| `a1b_classgraph.py` | `a1b_classgraph.out` | BFS in the circulant class graph, `[R_q]` to `[P_q]`, `q ≤ 31` |
| `a2_nonresidue.py` | `a2_nonresidue.out` | every step of the THM-4751 (2) proof, all `q ≡ 3 (4)` below `2·10^6` |
| `a3_hcount.c`, `a3_landscape.py` | `a3_landscape_17.out` (n = 3..17), `a3_landscape_19_23.out`, `a3_landscape_25_25.out`, `a3_classes_n.txt` | complete class tables, odd `n ≤ 25`: `H` by Karp inclusion–exclusion over necklaces (all `n`) and Held–Karp (`n ≤ 21`); class graph, BFS distances, local maxima, `E`, `c_4`, `c_5` (DFS) |
| `a3b_selfloops.py` | `a3b_selfloops.out` | classes whose strictness depends on the self-loop convention |
| `a4_q31_sets.txt` | `a4_q31.out` | `H(P_31)` and two different one-flip neighbours |
| `a5_misc.py` | `a5_misc.out` | QR codes, theta series, THM-489, Leech count, `E_8 ⊃ A_8`, forests, `c_3`, cannonball, cubes, `bP`, THM-4749 degrees |
| `a6_errata.py` | `a6_errata.out` | MISTAKE-659 items |
| `a7_hspectrum.c` | `a7_hspectrum.out` | exhaustive `H`-spectrum of all labelled tournaments, `n ≤ 7`, with witnesses |
| `a8_designs.py` | `a8_designs.out` | Witt chain from the `P_7`, `P_23` codes; Hughes divisibility |
| `a9_regular_c4c5.py` | `a9_regular_c4c5.out` | `c_5 + 2c_4` on all regular tournaments; symbolic proof |
| `a10_rot_sets.txt`, `a10_ratios.py` | `a10_rot.out`, `a10_ratios.out` | `H(R_27)`, `H(R_29)`, `H(R_31)`; margins; ratios; normalised values |
| `a11_covering.py` | `a11_covering.out` | covering radius of `Z[(1+√−q)/2]` |

Severity: **high** = a false statement typed PROVED or FINITE-EXACT; **medium** = a status or label that makes a conjecture look proved, a mislabelled number, missing prior art that pre-states the result, or a conflict with an earlier conjecture; **low** = imprecision, convention or citation detail.

## Summary verdicts

| # | item | verdict | severity |
|---|---|---|---|
| 1 | THM-4751 (1): Dirichlet `A_q − B_q = (2 − (2/q))h(−q)`, 1135 primes `7 ≤ q < 20000` | CONFIRMED | — |
| 2 | (1): flip distance `= B_q`; class-graph distance `=` min over multipliers | CONFIRMED (the lifting step is unstated) | low |
| 3 | (1): arc-overlap fraction `1/2 + (2 − (2/q))h/(q − 1)` | CONFIRMED | — |
| 4 | (1): "Heegner primes … least correlated with `R_q`" | CORRECTED (shared with `q = 23, 31`) | low |
| 5 | (2): `B_q ≤ 1` iff `q ∈ {3,7,11}`; least-non-residue proof | CONFIRMED (the bound alone works from `q = 23`) | — |
| 6 | (3): `H` of every circulant class, odd `n ≤ 23` (and 25) | CONFIRMED | — |
| 7 | (3): global maxima `P_7`, `Z_9{1,2,3,5}`, `P_11`, `R_13..R_25` | CONFIRMED | — |
| 8 | (3): strict-local-maximum lists | CORRECTED (true only if a self-loop flip counts as a tie; "class graph" convention gives 5 more) | medium |
| 9 | (3)/HYP-9362 (b): Paley margins `q = 7, 11, 19, 23` | CONFIRMED (two algorithms) | — |
| 10 | (3): all one-flip neighbours of `P_q` isomorphic | CONFIRMED (proof below; `q = 31` checked) | — |
| 11 | `q = 31` values `H(P_31)`, `H(P_31` one flip`)` | CONFIRMED (no longer single-source) | — |
| 12 | (4): logic | CONFIRMED, given "every maximiser" in HYP-9362 (a) | low |
| 13 | (5): deep hole `(q+1)²/(16q)`; norm-Euclidean iff `3, 7, 11`; PID iff Heegner | CONFIRMED (math); label CORRECTED (KNOWN facts + NUMEROLOGY alignment); citations missing | low |
| 14 | HYP-9362 (a) evidence (`n ≤ 23`) | CONFIRMED (also `n = 25`) | — |
| 15 | HYP-9362 (b): "`1 − margin` ≈ `2.5/q²`" | CORRECTED (`(1 − margin)q²` falls from 3.63 to 1.72) | low |
| 16 | HYP-9362 (c): ratios 1.0098 … 1.0239 | CONFIRMED | — |
| 17 | HYP-9362 (c): "Paley normalised value (2.55, 2.60, 2.63, 2.66)" | REFUTED (these are the rotational values) | medium |
| 18 | HYP-9362 (c) vs HYP-9028 / HYP-9274 (A) | CORRECTED (prior conjectures predict ratio → 1; not cited) | medium |
| 19 | HYP-9362 (d): Hughes, Lam–Thiel–Swiercz, Witt chains, no `S(4,5,9)` | CONFIRMED (computed from the Paley codes) | — |
| 20 | (d) and note §3: "7 and 21 are the only values never taken (THM-1370)" | CORRECTED (completeness is a CONJECTURE; repeats MISTAKE-559 (1)) | medium |
| 21 | (d) status DICTIONARY | CORRECTED (NUMEROLOGY per MISTAKE-559; (d)'s own refutation test is already met) | medium |
| 22 | (d): `Φ_3` values 13, 31, 43 attained | CONFIRMED (witnesses at `n = 5, 6, 6`) | — |
| 23 | (d): codes cited as "THM-4751 (5)" | CORRECTED (THM-481) | low |
| 24 | note §2: QR code table, 620, 17296, kernels | CONFIRMED | — |
| 25 | note §2, §7: `E_4`, `Δ`, `E_4³ − 672Δ`, Leech theta, `j − 720`, THM-489 | CONFIRMED | — |
| 26 | note §2, §6: 196560; `240 = 72 + 168`; `248 = 80 + 84 + 84` | CONFIRMED | — |
| 27 | note §5: forests 1 … 561948; `c_3` formula; cannonball | CONFIRMED | — |
| 28 | note §1: "exactly the `H`-maximal Paley primes"; "1135 primes … below 20000" | CORRECTED | low |
| 29 | note §4: Heegner MISTAKE citations | CORRECTED (225, 233 unrelated; HYP-1202/1220 not rejected) | low |
| 30 | note §9: THM-4749 quotes | CORRECTED (511/513, `k = 11`, ties 23/29 right; sentence self-contradictory, conditionality dropped) | low |
| 31 | note §10: THM-4749 quotes; perfect codes | CORRECTED (`k = 12` not a Fermat product; `k = 14` undecided; trivial perfect codes) | low |
| 32 | note: other cited files (THM-872/868/481/489/480/463/487, Session27 files, …) | CONFIRMED (number collisions THM-868, THM-1370, THM-137) | low |
| 33 | MISTAKE-659: `h(E_6) = 12`; Golay `[23,12,7]`; `P_7`: 80 cycles, `α_2 = 7`, `H = 189`; two Fano planes; dodecahedron; regular 5-tournaments `H = 15` | CONFIRMED | — |
| 34 | MISTAKE-659 / banner: `Z_E8(s) = 240ζ(s)ζ(s−3)` | CORRECTED (factor `2^(−s)` for the reflection's own `Σ|v|^(−2s)`) | low |
| 35 | banner completeness | CORRECTED (four further errors in the reflection) | low |
| 36 | prior art | CORRECTED (March 2026 orientation-cube program THM-134/135/137/139/140, HYP-464/480; earlier 7 → 21 → 24 chains; HYP-9028) | medium |
| 37 | addendum (6)(i): `c_4 − nE/4`, `c_5 + nE/2` constant | CONFIRMED; KNOWN in general (Savchenko 2016: `c_5 + 2c_4 = n(n²−1)(n²−9)/160` for every regular tournament) | medium (prior art) |
| 38 | (6)(ii): `R_p` maximises `E`; Lev; uniqueness | CONFIRMED; Lev covers it; uniqueness is PROVED for every odd `n` via the known `c_4` extremal theorem | low |
| 39 | (6)(ii): `P_q` minimises `E` | CONFIRMED (KNOWN: THM-162, THM-140, HYP-464); uniqueness provable | — |
| 40 | (6)(iii): finite table `7 ≤ n ≤ 23` | CONFIRMED (also `n = 25`); prior art THM-140 | — |
| 41 | (6)(iv) | CONFIRMED (restatement of (1)) | — |

No high-severity error: every computed number in THM-4751, HYP-9362 and §0–10 that I could test is right. The problems are statuses, conventions, one mislabelled row, and prior art.

## A. THM-4751 (1)

**1. Dirichlet formula — CONFIRMED.**
* `h(−q)` was obtained by counting reduced forms `(a,b,c)`, `b² − 4ac = −q`, with no L-function. `B_q` was obtained by Euler's criterion.
* All 1135 primes `q ≡ 3 (mod 4)`, `7 ≤ q < 20000`: 0 mismatches.
* All 1135 values also equal PARI `qfbclassno`, `quadclassunit(-q).no` and `Σ_(a<q/2)(a/q)/(2 − (2/q))` (`a1_class_number_gp.out`).
* The formula is false at `q = 3` (`A − B = 1 ≠ 3`), as the theorem's `q > 3` says.
* There are 1136 such primes below 20000 if `q = 3` is counted (`a1_class_number.out`). THM-4751 says "7 ≤ q < 20000", which is right. The note says "all 1135 primes `q ≡ 3 (mod 4)` below 20000", which needs the same qualifier (#28).

**2. Flip distance — CONFIRMED.**
* For `a ∈ QR` the count `#{k ≤ m : −k ∈ aQR}` is `B_q`; for `a ∈ NQR` it is `A_q > B_q` (`a1_class_number.out`, all 1135 primes).
* Class-graph distance equals the multiplier minimum. Multipliers act on connection sets by automorphisms of the flip hypercube: if `S, S′` differ in the class `{k, −k}`, then `aS, aS′` differ exactly in `{ak, −ak}`. So a class-graph path of length `k` lifts to `k` flips from `R` to some `aP`, and conversely.
* BFS confirms this:
  * for `[R_q]` to `[P_q]`, `q = 7, 11, 19, 23, 31` (distances 1, 1, 3, 4, 6; `a1b_classgraph.out`);
  * for every class to `[R_n]` and `[P_n]`, odd `n ≤ 25` (`a3_landscape_*.out`, "BFS distance == min-over-multipliers: True").
* THM-4751 §1 never states this lifting step. One sentence fixes it (edit T3).
* Aside, not claimed in the files: `[P_q]` is at maximal class distance from `[R_q]` for `q ≤ 31` (eccentricity of `[R_q]` = `B_q`).

**3. Arc overlap — CONFIRMED.** It is `A_q/m` by definition. Checked at `q = 23`: 161 of 253 arcs, which equals `1/2 + 3/22` (`a1_class_number.out`).

**4. "Least correlated" — CORRECTED (low).**
* For `q > 11`, `A_q − B_q = (2 − (2/q))h ≥ 3`, so `B_q ≤ (m − 3)/2`.
* Equality holds for `q ∈ {19, 23, 31, 43, 67, 163}` below 20000 (`a1_class_number.out`): the four Heegner primes `≡ 3 (mod 8)`, and also `q = 23, 31`, which have `h = 3` and `(2/q) = 1`.
* The list is complete: by Watkins/Oesterlé, the only `h = 3` primes `≡ 7 (mod 8)` are 23 and 31.
* So "least correlated" is not a property of the Heegner primes alone. The note's §1 line "as far from `R_q` as possible" has the same issue.

## B. THM-4751 (2) — CONFIRMED

Each step was checked for every prime `q ≡ 3 (mod 4)` below `2·10^6` (`a2_nonresidue.out`, 0 failures):

* For `q ≡ 3 (mod 8)`, `q ≥ 19`: 2 and 8 are non-residues and `≤ m`.
* For `q ≡ 7 (mod 8)`:
  * `m_0 = ⌈q/n_0⌉` gives `0 < m_0n_0 − q < n_0`;
  * `m_0n_0 − q` is a residue, so `m_0` is a non-residue and `m_0 ≥ n_0`;
  * hence `n_0(n_0 − 1) < q`, i.e. `n_0 < 1/2 + √(q + 1/4)`;
  * `2n_0` is a non-residue (since `(2/q) = 1`) and `2n_0 ≤ m`.

The direct counts give `B_q ≤ 1` only at `q = 3, 7, 11` (`q < 20000`).

Remark: the bound alone gives `2n_0 ≤ m` from `q = 23` on, since `2⌊1/2 + √23.25⌋ = 10 ≤ 11`. THM-4751's "for `q ≥ 31`" plus the separate `q = 23` check is correct, just weaker than needed (edit T2, optional).

## C. THM-4751 (3) — the landscape

**6. Class tables — CONFIRMED.**
* My `a3_hcount` computes `H` by inclusion–exclusion over binary necklaces (iterative FKM, weight = period, arithmetic mod `2^128`). This is a different algorithm from the session's `circ_counts.c` (Held–Karp from vertex 0).
* It agrees with my own Held–Karp on every class, odd `n ≤ 21`.
* It agrees class by class with the session data:
  * `circ_*.out`: `n = 3..23`, with 1, 1, 2, 4, 4, 6, 16, 16, 30, 88, 94 classes;
  * `circ25_H.txt`: 208 classes at `n = 25`, 208 agreements.
* Spot values: `H(P_23) = 15760206976379349`, `H(R_25) = 2418453569285650675`.

**7. Global maxima — CONFIRMED.**
* `P_7` (189), `Z_9{1,2,3,5}` (3357, the class of `{2,3,4,8}`, since `7·{2,3,4,8} = {1,2,3,5}`) and `P_11` (95095).
* `R_n` for `n = 13, …, 25`. The maximiser is unique in every case.
* Distance to `[R_n]` is 1 at `n = 7, 9, 11`.

**8. Strict local maxima — CORRECTED (medium; convention).**
* THM-4751 §3 and HYP-9362 (b) say "strict local maximum of `H` on the circulant class graph".
* The stated lists (7/9/11: global max only; 13, 15: `R_n` only; 17: `R` + one class at `d = 2`; 19: `R`, `P_19`, one class at `d = 2`; 21: 5 classes; 23: 6 classes including `P_23`) are exactly right only in the **hypercube convention**. There a flip that lands in the same multiplier class (a self-loop of the class graph) counts as a tie, so a class with a self-loop is never strict. This is what `circ_landscape.py` computes.
* In the class graph (loops dropped, the usual convention for "neighbour"), five more classes are strict local maxima (`a3b_selfloops.out`):

| `n` | class | `H` | `d_R` | self-loop flip `{k, n−k}` → multiplier |
|---|---|---|---|---|
| 11 | `Z_11{1,2,3,4,6}` | 93467 | 1 | `k = 3 → a = 2`, `k = 4 → a = 6` |
| 15 | `Z_15{1,2,5,7,9,11,12}` (an energy minimiser) | 197068455 | 2 | `k = 5 → a = 11` |
| 19 | `Z_19{1,2,5,9,11,12,13,15,16}` | 1167653827609 | 3 | `k = 1 → 15`, `k = 5 → 14` |
| 21 | `Z_21{1,2,3,7,8,9,10,15,16,17}` (2nd-highest class) | 125056114727355 | 3 | `k = 7 → 8` |
| 23 | `Z_23{1,2,5,6,7,8,9,11,13,19,20}` | 15720696951606199 | 3 | `k = 1 → 11`, `k = 2 → 21` |

* `n = 25` has no self-loop classes, and both conventions give the session's 6 maxima.
* The Paley statements are unaffected: `P_q` has no self-loop and a single neighbour class.
* The §13 table rows "7, 9, 11 … the only strict local maximum" and "13, 15 … not a local maximum" inherit the same convention. At `n = 11` and `n = 15` they are false in the class graph.

**9. Paley margins — CONFIRMED.**
* `H(P_q)` and `H(P_q` with `{1, −1}` flipped`)` are 189/175, 95095/93027, 1172695746915/1163848601583, 15760206976379349/15688038105929693.
* Margins: 0.925926, 0.978253, 0.992456, 0.995421 (`a10_ratios.out`).
* `q = 23` is not single-source: I computed it by inclusion–exclusion, not by a `2^23` DP. The session's value came from Held–Karp.

**10. Isomorphic flip-neighbours — CONFIRMED (proof).**
* Difference classes are the cosets `{±k}` of `{±1}` in `Z_q^×`.
* Because `q ≡ 3 (mod 4)`, `−1 ∉ QR`. So `QR ∩ {±1} = {1}`, and `QR → Z_q^×/{±1}` is an isomorphism. Hence `QR` acts simply transitively on the `m` classes.
* Each `a ∈ QR` is an automorphism of `P_q`. So `x ↦ ax` maps `P_q` with `{±k}` flipped onto `P_q` with `{±ak}` flipped.
* Every flip of `P_q` therefore lies in one multiplier class: one neighbour class, and no self-loop, since `aQR = QR` for `a ∈ QR`.
* Checked at `q = 31`: flipping `{1, 30}` and flipping `{3, 28}` both give `19841263005616787646549771` (`a4_q31.out`).

**11. `q = 31` — CONFIRMED.**
* My `a3_hcount` gives `H(P_31) = 19876729650849384218015865` and one flip `19841263005616787646549771`, so the margin is 0.998216 (69,273,668 necklaces, about 3 min per tournament).
* Independent earlier sources for `H(P_31)`:
  * opus S34 (2026-10-14) and its audit (`paley_h31_audit_20261014/REPORT.md`: Karp inclusion–exclusion modulo three 58-bit primes, 397 s);
  * mac-mini `paley_extremal_20261021/ie_ham.c`.
* S51's `hyp9357_rotational_vs_paley.py` does not recompute `H(P_31)`; it imports the audited S34 value.
* My run also reproduces `H(R_27) = 426990256001863046739`, `H(R_29) = 87133267870433699742519` (= session `rot29_H.txt`) and `H(R_31) = 20351875726012821643802839` (`a10_rot.out`).

## D. THM-4751 (4), (5)

**12. (4) — CONFIRMED with a wording condition (low).**
* If every `H`-maximal circulant class of order `q` is within one flip of `[R_q]`, then a Paley maximiser would force `B_q ≤ 1`, i.e. `q ≤ 11`, by (2).
* HYP-9362 (a) says "the circulant tournament maximising `H` … is `R_n` or differs from it in one difference class". This presupposes a unique maximiser, which holds for `n ≤ 25`. It should say "every `H`-maximal circulant". Otherwise a tied Paley maximiser far from `R` is not excluded.
* For the record, (4)'s conclusion is already FINITE-EXACT at `q = 19, 23, 31`.

**13. (5) — mathematics CONFIRMED; label and citations CORRECTED (low).**
* The triangle `0, 1, ω = (1 + √−q)/2` is isosceles and acute. Its circumcentre `z* = 1/2 + ((q − 1)/(4q))√−q` lies in `K`, with `|z*|² = (q + 1)²/(16q)`.
* A grid search confirms that `z*` is the deep hole for `q = 3 … 163` (`a11_covering.out`).
* `(q + 1)² < 16q` iff `q < 7 + √48 = 13.93`. Among primes `q ≡ 3 (mod 4)` below `10^4`, that is `q ∈ {3, 7, 11}`.
* Since `z* ∈ K`, this proves "norm-Euclidean exactly for 3, 7, 11" directly; no citation is needed. Context: Lemmermeyer, Expo. Math. 13 (1995) 385–416; Motzkin, Bull. AMS 55 (1949) 1142–1146 (19, 43, 67, 163 are not Euclidean for any function).
* "PID exactly for the Heegner primes" needs Heegner (Math. Z. 56 (1952) 227–253), Baker (Mathematika 13 (1966) 204–216) and Stark (Michigan Math. J. 14 (1967) 1–27). THM-4751 lists them only as "(CITED)" without references. Numerically, `h = 1` occurs below 20000 only at 3, 7, 11, 19, 43, 67, 163.
* Label: both facts are KNOWN. Their alignment with (2) is a coincidence of a 3-element set, and the theorem itself notes that any threshold in `(11, 19)` gives it. MISTAKE-559's rule types this NUMEROLOGY, not DICTIONARY.

## E. HYP-9362

**14. (a) — CONFIRMED** for `n ≤ 23`, and the session's `n = 25` table is reproduced (`R_25`). Prior art is HYP-480 (#36).

**15. (b) — margins CONFIRMED; the decay claim CORRECTED (low).**
* `(1 − margin)·q² = 3.630, 2.631, 2.723, 2.422, 1.715` at `q = 7, 11, 19, 23, 31` (`a10_ratios.out`).
* There is no stable constant 2.5. At `q = 31` the constant is 1.7, so the decay is at least quadratic over this range.

**16–18. (c).**
* The ratios are right: 1.009821, 1.015947, 1.020468, 1.023905.
* **REFUTED (medium):** "the Paley normalised value `H·2^(q−1)/q!` (2.55, 2.60, 2.63, 2.66)". Those are the **rotational** values (2.552, 2.598, 2.632, 2.658). The Paley values are 2.527, 2.557, 2.579, 2.596 (`a10_ratios.out`; same as S51's `hyp9357_rotational_vs_paley.out` and HYP-9357's table).
* **Conflict (medium):** two earlier conjectures together predict `H(R_q)/H(P_q) → 1`:
  * HYP-9028 (kind-pasteur S134): both the rotational and the Paley excess tend to `e`;
  * HYP-9274 (A): `R(q) → e` for Paley.

  So, given HYP-9274 (A), (c)'s "limit > 1" is equivalent to the rotational excess tending to a limit above `e`, i.e. to HYP-9028 failing for `R_n`.
* The data lean that way. The rotational values 2.632, 2.645, 2.658 at `n = 27, 29, 31` rise by 0.0138 and then 0.0122, and HYP-9028's fitted `α` has fallen to about 0.69 at `n = 31`. Neither HYP-9028 nor HYP-9274 is cited.
* `q = 27` is a prime power, and `P_27` is not a circulant on `Z_27`. Its ratio is not a circulant statement (low).
* "Wormald's blow-up constant for the asymptotic maximum (about 2.8559)": this is the best known lower-bound constant, not the maximum (HYP-9357 says this correctly).

**19. (d) facts — CONFIRMED (`a8_designs.out`).**
* An extension of a plane of order `n` is a `3-(n²+n+2, n+2, 1)` design. Its block count is an integer iff `(n + 2) | 12`, i.e. `n ∈ {1, 2, 4, 10}`. This is the counting condition of Hughes, Amer. J. Math. 87 (1965) 761–778.
* Lam–Thiel–Swiercz, Canad. J. Math. 41 (1989) 1117–1123.
* From the code of `P_7`: the 14 weight-4 words form `S(3,4,8)`, and the derived design at the parity point is the Fano plane.
* From the code of `P_23`: the 759 octads form `S(5,8,24)`. Derived designs are `S(4,7,23)` (253 blocks) and `S(3,6,22)` (77 blocks). At three points one gets a 2-(21,5,1) design with 168 hyperovals, the projective plane of order 4, unique and hence `PG(2,4)`.
* `S(4,5,9)` and `S(6,9,25)` fail divisibility. Witt, Abh. Math. Sem. Hamburg 12 (1938).

**20. "Only values never taken" — CORRECTED (medium).**
* THM-1370 (h-spectrum file) proves that 7 and 21 are never `H` values.
* Its status says "Completeness of the range remains OPEN", and it lists "CONJECTURE (spectrum completeness)".
* THM-4749 has the strong values FINITE-EXACT below 80,265.
* HYP-9362 (d) ("THM-1370: 7 and 21 are the only values never taken by `H`") and note §3 ("are the only values never taken by `H` (THM-029, THM-079, THM-1370)") repeat MISTAKE-559 (1) word for word in substance.

**21. DICTIONARY → NUMEROLOGY (medium).**
* MISTAKE-559 states the rule ("a coincidence of small numbers is NUMEROLOGY unless a map carries the mechanism"). It explicitly kept "the NUMEROLOGY verdict on {7, 21}".
* opus S52 typed the projective-plane reading NUMEROLOGY (`lattices_modular_repunits_20261021.md` §5).
* HYP-9362's own refutation test ("an independent characterization of the forbidden values … that makes no reference to the projective plane structure") is met by THM-1370's proof. There 7 and 21 sit in the strong-floor crevasses `(5, 9)` and `(15, 25)`, and the hypothesis admits this.
* So (d) is NUMEROLOGY. The same applies to the note's §0 rows for `{7,21}` and `{3,7,11}`, and to THM-4751 (5)'s alignment.

**22. `Φ_3` remarks — CONFIRMED** by exhaustive enumeration of all labelled tournaments on `n ≤ 7` (`a7_hspectrum.out`, explicit witness arcs printed):
* `H = 13 = Φ_3(3)` is first attained at `n = 5`;
* `H = 31 = Φ_3(5)` and `H = 43 = Φ_3(6)` are first attained at `n = 6`;
* 7 and 21 never occur for `n ≤ 7`;
* 63 is absent at `n = 7`, consistent with THM-1370.

**23.** (d)'s "Two Paley routes … (THM-4751 (5))" should cite THM-481. (5) is the Euclidean/Heegner statement (low).

## F. Results note §1–10

**24. Codes — CONFIRMED (`a5_misc.out`).**
* Ranks of `A(P_q)` over GF(2) are `(q + 1)/2` for `q = 7, 23, 31, 47, 71`.
* Extended codes, all doubly even and self-dual:
  * `[8,4,4]`: `1 + 14y⁴ + y⁸`;
  * `[24,12,8]`: `1 + 759y⁸ + 2576y¹² + 759y¹⁶ + y²⁴`;
  * `[32,16,8]`: `A_8 = 620`, `A_12 = 13888`, `A_16 = 36518`;
  * `[48,24,12]`: `A_12 = 17296` (extremal).
* Kernels: `ker A(P_7) = [7,3,4]` and `ker A(P_23) = [23,11,8]`, each equal to the even-weight subcode of the row space. The row spaces `[7,4,3]` and `[23,12,7]` are perfect.

**25. Modular forms — CONFIRMED.**
* `E_4 = 1, 240, 2160, 6720, 17520`; `Δ = q − 24q² + 252q³ − 1472q⁴`.
* `E_8 = E_4²`, `E_10 = E_4E_6`, `E_14 = E_4²E_6`.
* `E_4³ − 672Δ = 1 + 48q + 195408q² + …`, which equals the Construction A theta of `[24,12,8]` computed from the weight enumerator.
* `E_4³ − 720Δ = 1 + 196560q² + …`; `θ_Λ/Δ = q⁻¹ + 24 + 196884q`.
* THM-489: `W_g24 = W_e8³ − 42P_24` holds, and `P_24 ↦ 16Δ` under Construction A.

**26. Lattices — CONFIRMED.**
* `4·276 + 759·128 + 24·4096 = 196560`.
* `E_8 = A_8 + Z[3]` was constructed explicitly. It has 72 + 84 + 84 = 240 norm-2 vectors (84 in each glue class, found by brute force), and the Gram determinant of an HNF basis is 1.
* `248 = 80 + 84 + 84`.

**27. §5 — CONFIRMED.**
* Distinct labelled score vectors, computed by DP over arcs (not by Landau's criterion), are 1, 2, 7, 38, 291, 2932, 36961, 561948 for `n ≤ 8`. These equal the labelled forest counts from the tree-recursion.
* `c_3 = (n³ − n)/24 − |s − s̄|²/2`: proved symbolically and checked on 560 random tournaments, `n ≤ 30`.
* `Σ_(k≤N) k²` is a square, for `N < 10^6`, only at `N = 1, 24` (Watson, Messenger of Math. 48 (1918)).
* `P_7` with the cyclic triangle `{0,1,3}` reversed has unchanged scores and `H = 171`, matching Session27's `SCORE_ROOT_COLLISION.md`.

**28. §1 overclaims (low).**
* "These are exactly the `H`-maximal Paley primes": circulant `H`-maximality is FINITE-EXACT only for `q ≤ 31`; global maximality is known only for `q ≤ 11` (A038375); beyond that it is HYP-9357's conjecture.
* §4 item 1 has the same scope problem.
* "1135 primes … below 20000" needs "`7 ≤ q`" (#1).

**29. §4 Heegner citations (low).**
* MISTAKE-225 (skew parity and tournament games) and MISTAKE-233 (theta and modular curves for LRC) do not mention Heegner numbers.
* MISTAKE-229, -230 and -576 are on point; MISTAKE-568 and -569 are Heegner-related and uncited.
* The repo has not rejected "Heegner numbers as tournament constants": HYP-1202 and HYP-1220 are still CONFIRMED (as numerology) in `INDEX-HISTORICAL-THROUGH-2026-07-21.md`.

**30. §9 THM-4749 quotes (low).**
* Recomputed `w(2^k ± 1)` for `k ≤ 40` with `V` = odd `≥ 3` minus `{7, 21}` (`a5_misc.out`):
  * the THM-4749 table for `k ≤ 16` is reproduced exactly;
  * odd-`k` ties occur only at 11, 23, 29;
  * the Mersenne side is never larger for odd `k`;
  * the certified `k` are all `k ≤ 16` plus 18, 20, 21, 22, 24, 25, 26, 30, 36, exactly as in THM-4749's status.
* So "511 rarer than 513" (`w = 1` vs 4), "`k = 11` by leading coefficient" (THM-4749: 2049 commoner) and "ties at `k = 23, 29`" are accurate.
* Two problems with the §9 sentence:
  * "Every one of these … is either forbidden (7) or rarer than its Fermat partner" contradicts its own third bullet (`k = 23, 29` undecided);
  * it drops THM-4749's caveat that the odd `k ≥ 17` other than 21 and 25 are conditional on `V`-membership of large factors.

**31. §10 (low).**
* "commoner value for `k = 4, 8, 12, 16` (Fermat products)": `2^12 − 1 = 4095` is not a Fermat product. Its commonness comes from absorbing 7 into 35 (THM-4749 statement 3).
* "ties at `k = 10, 11, 14` are decided by the leading coefficients": THM-4749 decides only `k = 10` (1023 commoner) and `k = 11` (2047 rarer). The `k = 14` values (16383, 16385) exceed its `h < 2100` polynomials.
* "the perfect binary codes besides the Golay code": this omits the trivial perfect codes (repetition codes, the whole space). Up to parameters the statement is Tietäväinen–van Lint (1973); non-linear codes with Hamming parameters also exist.
* The THM-4748 summary is accurate: `fa ≥ 3` for `H = 2^k − 1`, `k ≥ 5`; the Fibonacci fan; `2^k + 1 = H(NT_(k+2))`.

**32. Cited files — CONFIRMED.**
* Each of these exists and says what is claimed:
  * THM-872: "OPEN which of the 24 Niemeier lattices, if any, is score-selected"; `D_24⁺`-coset, not unimodular; "the metagraph hears one of the two drums";
  * `THM-868-e8-bridge-score-lattice.md`: half-vectors in `E_8`'s half-integer coset, 70 = C(8,4) roots as the near-regular shell;
  * THM-481: Kim–Solé 2008 rediscovery, claudebox + kind-pasteur, 2026-06-11;
  * THM-489: `P_24 ↦ 16Δ`;
  * THM-480 (Milnor pair); THM-463 (Eisenstein norms); THM-487 (Nebe over `Z[(1+√−7)/2]`; Paley gauge stalls at `d = 12`);
  * THM-4652/4669 (`η(2z)⁴η(4z)⁴`, Mason class 4A, cycle shape `2⁴4⁴`); HYP-9275/9276; HYP-9028 (1729); THM-4746 (5) (Carmichael numbers pass); HYP-1348, HYP-1317, HYP-2306;
  * the Session27 files `LORENTZIAN_SCORE_BRIDGE.md` (S52 concurrent, MSG-1913), `SCORE_ROOT_COLLISION.md`, `EXPLORATION.md` ("still owes a map … QR-code length rule q→q+1 instead pairs 7 with 8 and 23 with 24"), `session27-modular/report.md` and `PRIOR_1729_AND_E8_HOLES.md` (16-vertex cross-polytope; 23 Leech deep-hole types), `SYNTHESIS.md` (nine Heegner radicands including 1), `session27-smale/ADDENDUM.md`, `session27-repunits/REPORT.md`, `session27-papers-review/REVIEW.md`.
* `|bP_4m|` for `m = 2..6` (28, 992, 8128, 261632, 1448424448) is divisible by `2^(2m−1) − 1`.
* Number collisions mean citations should give the path:
  * two THM-868 files;
  * two THM-1370 files (the other is the elliptic Jacobian theorem);
  * two THM-137 files (crossover mechanism; eigenvector interaction).

## G. MISTAKE-659 and the errata banner

**33. CONFIRMED (`a6_errata.out`, `a5_misc.out`).** Each listed error is present in the reflection text, and each correction is right:
* "Coxeter number of E_6" is at reflection line 75. `h = |Φ|/rank = 72/6 = 12`.
* The perfect Golay code is `[23,12,7]`; the kernel code is `[23,11,8]`, the even subcode; `[24,12,8]` is not perfect (lines 89–90, 243).
* `P_7` has `c_3 = 14`, `c_5 = 42`, `c_7 = 24`, so `α_1 = 80`; `α_2 = 7`; and `1 + 160 + 28 = 189` equals the DP value (line 110).
* The cyclic triples are `dev{0,1,3} ∪ dev{0,1,5}`: disjoint, each a Fano plane, 14 distinct triples, each with one orientation (lines 52–54).
* The dodecahedron has `(V, E, F) = (20, 30, 12)` (line 158).
* The mass-formula paragraph is garbled (line 175). The `ζ_M` paragraph rests on the k-nacci numerology, which is MISTAKE-028.
* All 24 labelled regular 5-tournaments have `H = 15`. HYP-1363's "13" is wrong, and the index line is quoted correctly.
* The `one_twenty_eight_atlas` edit is right, and the original "likewise" wording was wrong.
* `golay_baer_bridge_88.py` (lines 155–157) is correct as stated.

**34. `Z_E8` normalisation — CORRECTED (low).**
* The reflection defines `Z_E8(s) = Σ_(v≠0)|v|^(−2s)`. With roots of norm 2 the shell `|v|² = 2m` has `240σ_3(m)` vectors, so `Z_E8(s) = 2^(−s)·240ζ(s)ζ(s−3)`.
* `240ζ(s)ζ(s−3)` is the series `Σ(|v|²/2)^(−s)`.
* Checked at `s = 6`: the partial sum is 4.58589, against 4.58589 for the `2^(−s)` form and 293.497 without it.
* The only pole is at `s = 4` (residue 16.235). At `s = 1` the product is finite, because `ζ(−2) = 0`.

**35. Further errors in the reflection, not in the banner (low):**
* (a) §1.1, lines 63–64: "The automorphism group of the oriented Fano plane is G_2". `G_2` is `Aut(O)`, a 14-dimensional Lie group. The collineations of the Fano plane preserving `P_7`'s line orientations form a group of order 21 (`a6_errata.out`).
* (b) §3.1, line 221: "Bott periodicity (Atiyah-Hirzebruch, 1961)". It is due to Bott (1957/59).
* (c) §3.3, lines 281–282: "At s=2: Z_E8(2) = π^4/90 × (mass formula constant)". `s = 2` is outside the half-plane of convergence; the continued value is `60ζ(2)ζ(−1) = −5π²/6`.
* (d) §1.2, line 68: "regular (= self-complementary)". These are not equivalent; the transitive `TT_5` is self-converse.
* The §2.2 sentence "21 … is also the complete conflict graph `Ω(T_7)`" (line 139) is a consequence of the `K_7` error and should be flagged with it.

## H. Prior art (medium)

THM-4751, HYP-9362 and §1/§13 of the note cite only THM-4666, the Ogg note, THM-481 and HYP-9357. The repo already has a March 2026 "orientation cube" program on exactly this landscape. The orientation cube `{±1}^m` is the flip hypercube.

* **THM-135** (opus S58): `P_19` is not the circulant maximum; the interval `{1..9}` wins. It records `H(C_19) = 1184212824763` and `H(T_19) = 1172695746915`. HYP-464 ("Paley maximises `H` among circulants") is refuted at 19.
* **THM-137** (`THM-137-paley-crossover-mechanism.md`, kind-pasteur S56c): Paley maximises `H` among circulants at `p = 3, 7, 11` but not at 19 or 23. It proposes a two-mechanism reading (cycle counts favour Paley, independence favours the interval). Open question 2 is "Does the interval maximize `H` among ALL circulant tournaments for all `p ≥ 19`? (HYP-480)".
* **THM-139**, table footnote: "At p=19: … Paley is #2 (local max, no single flip improves)". This is a prior statement of the Paley one-flip peak at `p = 19`.
* **THM-134**: Paley is the unique local maximum of `H` on the Parseval simplex, `p = 7, 11`.
* **THM-140** (kind-pasteur S58): additive energy against `H`.
  * "Interval always has maximum E", "Paley always has minimum E".
  * `corr(E, H) = −1.000, −0.088, +0.509` at `p = 7, 11, 13`.
  * `E(Int) = m(2m² + 1)/3` and `E(QR) = m(m² + 1)/2`; both check against my tables, e.g. 891 and 671 at `n = 23`.
  * The same pseudorandom-versus-additive reading as §13.
* **HYP-480**: "INTERVAL UNIVERSAL H-MAXIMIZER", exhaustive for primes `p ≤ 19` and `p = 23`. HYP-9362 (a) refines it to all odd `n` and to "within one flip".
* **HYP-464** (Paley maximises `c_5` among circulants) and **HYP-485** (Paley orientation is the top Walsh-2 eigenvector "by QR multiplication equivariance"). HYP-485 uses the same mechanism as #10.
* **HYP-9028** and **HYP-9274 (A)** for HYP-9362 (c) (#18).
* **For (d) and §3:**
  * opus S80 (`05-knowledge/results/fano_hamming_forbidden.out`): a Steiner chain, garbled as `S(2,3,7) → S(3,6,22) → …`;
  * opus S88 (`04-computation/golay_baer_bridge_88.py`): "Fano → PG(2,4) → Golay → Leech" and "7 → 21 → 24";
  * kind-pasteur S100 (`05-knowledge/results/projective_geometry_deep.out`): "PG(2,4) + 3 = 24 → GOLAY CODE … S(5,8,24) can be built from PG(2,4) by adding 3 points";
  * the 2026-05-02 reflection's Chains A/B;
  * concurrently, opus S52 (`lattices_modular_repunits_20261021.md` §5, NUMEROLOGY), which also attributes the reading to mac-mini (commit `8a3a6bbd2e`).

  New in HYP-9362 (d) is only the use of Hughes and Lam as the selection principle. The note's "new reading" in §0 should say so.
* **Searched with no earlier hit:**
  * class number as a distance between circulant tournaments;
  * "flip distance" or "one-flip" in this sense (the hits are unrelated: strategy cubes, FKN, deletion decks);
  * Hughes's theorem (the only other hits are a different Hughes in THM-4529).

## Addendum: proposed THM-4751 (6) (energy extremes)

**37. (i) — CONFIRMED, and KNOWN in greater generality (medium, prior art).**
* On circulants, with `λ_j = −1/2 + iy_j` and `Σ_j y_j² = m(m + 1/2)` (Parseval), I derive exactly:
  * `c_4 = n(E − m)/4`;
  * `2c_5 + nE = (2m⁵ + 5m⁴ + 5m² + 3m)/5`.
* Both were checked on every class `5 ≤ n ≤ 25` with `c_4`, `c_5` counted by DFS, not by zero-sum counting (`a3_landscape_*.out`; constants `−nm` and, e.g., 79189 at `n = 23`).
* Normality is not needed. For **every** regular tournament, write `A = (J − I + B)/2` with `B = A − Aᵀ` skew and `BJ = JB = 0`. Then `(J + B − I)^k = (B − I)^k + α_kJ`, so `tr A^4` and `tr A^5` are both affine in `tr B⁴`. Eliminating it:
  * `c_5 + 2c_4 = n(n² − 1)(n² − 9)/160`;
  * `c_4 = (n⁴ − 4n³ + 3n + tr B⁴)/64`.
* Checked symbolically, exhaustively on all 2640 regular 7-tournaments, and on random regular tournaments for `n = 9, 11, 13` (`a9_regular_c4c5.out`).
* The identity is **Savchenko, J. Graph Theory 83 (2016) 44–77**, restated with a spectral proof as Prop. A.1 of Savchenko, arXiv:2403.19555. The repo's INV-053 cites that paper ("c5, c6 formulas").
* THM-162 already observed `corr(Σ|λ|⁴, c_5) = −1` at `p = 11`.
* New in (i) is only the circulant translation `c_4 = n(E − m)/4`.

**38. (ii) maximum — CONFIRMED; uniqueness is PROVED, not only FINITE-EXACT.**
* Lev, Duke Math. J. 107 (2001) 239–263 (abstract): for a linear equation with nonzero coefficients over `F_p` and variables in arbitrary sets of given sizes, the number of solutions is at most the number when the sets are replaced by arithmetic progressions centred at 0.
* For `a_1 + a_2 − a_3 − a_4 = 0` with all `a_i ∈ S`, put `A_1 = A_2 = S` and `A_3 = A_4 = −S` (same size). The bound is then `E(P)`, `P` a centred progression of length `|S|`, and `E(P) = E({1..m})`. So Lev covers it.
* I checked Lev only at abstract level (Project Euclid). The full text was not accessible. The same reduction is used for Lev's `Z` version in Aaronson, arXiv:1801.07135.
* A stronger, uniqueness-giving route exists for every odd `n`, composite included:
  * `c_4 = s_4` is maximised over all `n`-tournaments, odd `n ≥ 5`, exactly by the regular locally transitive tournament `RLT_n ≅ R_n`. This is stated in arXiv:2403.19555 §1, which builds on Beineke–Harary (Canad. Math. Bull. 8 (1965)) and Kotzig (1968) for the `s_4` bounds. I read the statement there and did not check the original proofs.
  * With `c_4 = n(E − m)/4` this makes `R_n` the unique `E`-maximiser up to isomorphism.
  * Isomorphism implies multiplier-equivalence here. `Aut(R_n) = Z_n`, because the stabiliser of 0 fixes the transitive out- and in-neighbourhoods `{1..m}`, `{m+1..2m}` pointwise. So an isomorphism `Cay(Z_n, S) → R_n` normalises `Z_n` and is affine.
* My tables confirm a unique `E`-maximum at `R_n` for every odd `7 ≤ n ≤ 25`.

**39. (ii) minimum — CONFIRMED (KNOWN: THM-162, THM-140, HYP-464).**
* Cauchy–Schwarz on `Σ_(j≠0)|λ_j|⁴` with `Σ_(j≠0)|λ_j|²` fixed gives the bound, with equality iff the spectrum is flat.
* For prime `q`, flatness forces `λ_j = (−1 ± √−q)/2 ∈ Q(√−q)`. Then `σ_a(λ_j) = λ_{aj} = λ_j` for `a ∈ QR`, so `S` is `QR`-invariant, hence `S ∈ {QR, NQR}`. So `P_q` is the unique `E`-minimiser up to multipliers. The note does not state this uniqueness; it can be added as PROVED.

**40. (iii) — CONFIRMED** (`a3_landscape_*.out`), on the complete class tables:
* `E` has a unique maximum at `R_n` for `7 ≤ n ≤ 25`.
* The `E`-minimum is unique except at `n = 15` (4 classes, `E = 183`) and `n = 21` (2 classes, `E = 514`).
* It is Paley at `n = 7, 11, 19, 23`.
* The `H`-maximal class is the `E`-minimiser at `n = 7, 9, 11` and `R_n` (`E`-maximiser) for `13 ≤ n ≤ 25`.
* The §13 table's `H`-ranks of the energy minimum (3, 3, 11–12, 8, 23 for `n = 17..25`) are right. "`Z_25{…}` is the highest of the others" is right.
* For the `n = 11` and `n = 15` rows, see #8.
* `n = 7, 11, 13` are already in THM-140.

**41. (iv) — CONFIRMED.** It is a restatement of (1) once (ii) identifies the extremes. "The two energy extremes `R_q`, `P_q` are `B_q` flips apart" holds for prime `q ≡ 3 (mod 4)`, with both extremes unique (#38, #39).

**§13, other.**
* The snapshot's sentence "Paley maximises every `c_k` among doubly regular tournaments (THM-161)" misreads Savchenko. Savchenko's results compare `DR_n` with `RLT_n` (arXiv:2403.07629); `c_k` for small `k` is constant across DRTs; THM-161's own wording overstates this.
* The current version replaces the sentence with a "Caveat to THM-161", which is the right direction. The later `α`-decomposition text was not audited.

## Required edits (old → new)

**THM-4751**

* **T1** (§1, line 61). `These are the Paley tournaments least correlated with `R_q`.` → `For q > 11 the excess A_q − B_q = (2 − (2/q))h(−q) is at least 3, so B_q ≤ (m − 3)/2; equality holds exactly when (2 − (2/q))h(−q) = 3, i.e. at the Heegner primes 19, 43, 67, 163 and at q = 23, 31 (h = 3, (2/q) = 1).`
* **T2** (§2, line 72, optional). `This gives `2n_0 ≤ m` for `q ≥ 31`.` → `This gives `2n_0 ≤ m` for `q ≥ 23`.`
* **T3** (§1, after "The minimum over multipliers … is `B_q`."). Add: `Multipliers act on connection sets by automorphisms of the flip hypercube, so a path of k flips between two classes lifts to k flips from R_q to some aP_q; hence this minimum is the class-graph distance.`
* **T4** (§3). `**Strict local maxima.**` → `**Strict local maxima** (in the hypercube of connection sets: a flip that lands in the same multiplier class counts as a tie).` Then append: `Against distinct neighbour classes only (the class graph without loops), five more classes are strict: Z_11{1,2,3,4,6}, Z_15{1,2,5,7,9,11,12}, Z_19{1,2,5,9,11,12,13,15,16}, Z_21{1,2,3,7,8,9,10,15,16,17} and Z_23{1,2,5,6,7,8,9,11,13,19,20}; each has a self-loop flip.` Use the same convention in HYP-9362 (b) and §13's table.
* **T5** (§4). `If the one-flip conjecture HYP-9362 (a) holds` → `If HYP-9362 (a) holds in the form "every H-maximal circulant class is within one flip of [R_n]"`.
* **T6** (title and §5 heading). `(5) DICTIONARY:` → `(5) KNOWN facts; their alignment with (2) is NUMEROLOGY:`. Also `## 5. Euclidean and Heegner primes (DICTIONARY)` → `(KNOWN; alignment NUMEROLOGY)`, and status `CAUTION: statement 5 is a DICTIONARY.` → `CAUTION: the alignment in statement 5 is NUMEROLOGY (MISTAKE-559).`
* **T7** (depends_on). Replace "(CITED)" with references: Heegner, Math. Z. 56 (1952) 227–253; Baker, Mathematika 13 (1966) 204–216; Stark, Michigan Math. J. 14 (1967) 1–27. Note that the Euclidean statement is proved by §5's covering-radius computation, with Lemmermeyer, Expo. Math. 13 (1995) 385–416 for context.
* **T8** (status, PRIOR ART IN THE REPO). Add THM-134, THM-135, `THM-137-paley-crossover-mechanism.md`, THM-139 ("Paley is #2 (local max, no single flip improves)" at `p = 19`), THM-140 (interval max `E`, Paley min `E`, `E`–`H` sign reversal), HYP-464 and HYP-480 (interval universal `H`-maximiser, exhaustive `p ≤ 19`, `p = 23`).
* **T9** (status). `one H computation at q = 31` → `two independent H computations at q = 31 (ie_ham.c; audit a3_hcount.c)`.

**HYP-9362**

* **Y1** ((b), line 46). `` `1 − margin` decays roughly like `2.5/q²` over this range. `` → `` `(1 − margin)·q²` = 3.63, 2.63, 2.72, 2.42, 1.72 at `q = 7, 11, 19, 23, 31`: at least quadratic decay, with no stable constant. ``
* **Y2** ((c), line 52). ``Wormald's blow-up constant for the asymptotic maximum (about 2.8559) and the Paley normalised value `H·2^(q−1)/q!` (2.55, 2.60, 2.63, 2.66) suggest a finite limit.`` → ``Wormald's blow-up constant (about 2.8559, the best known lower-bound constant for max H) and the normalised values `H·2^(q−1)/q!` (rotational 2.552, 2.598, 2.632, 2.658; Paley 2.527, 2.557, 2.579, 2.596) suggest a finite limit. HYP-9274 (A) conjectures that the Paley value tends to e, and HYP-9028 that the rotational value does too; together they predict ratio → 1. So (c) amounts to the rotational limit exceeding e, i.e. to HYP-9028 failing for R_n. q = 27 is a prime power (P_27 is not a circulant).``
* **Y3** ((d), line 61). `` * THM-1370: `7` and `21` are the only values never taken by `H`, at any `n`. `` → `` * THM-1370: `7` and `21` are never taken by `H`, at any `n` (PROVED). That they are the only odd values missed is THM-1370's completeness CONJECTURE (FINITE-EXACT below 80,265: THM-4749); MISTAKE-559 (1). ``
* **Y4** ((d), line 62). `(THM-4751 (5))` → `(THM-481)`.
* **Y5** (status, line 12; and the last paragraph). `it is recorded as DICTIONARY.` → `it is recorded as NUMEROLOGY (MISTAKE-559 kept the NUMEROLOGY verdict on {7, 21}).` And `so (d) is at best a dictionary unless a map is found.` → `so by this criterion (d) is NUMEROLOGY unless a map is found.`
* **Y6** ((d)). Add prior art: opus S80 `fano_hamming_forbidden.out`, opus S88 `golay_baer_bridge_88.py` ("7 → 21 → 24"), kind-pasteur S100 `projective_geometry_deep.out` ("PG(2,4) + 3 = 24"), the 2026-05-02 reflection, and S52 `lattices_modular_repunits_20261021.md` §5. New here: the Hughes/Lam selection.
* **Y7** ((a)). `the circulant tournament maximising the number H of Hamiltonian paths is` → `every circulant tournament maximising the number H of Hamiltonian paths is`. Add: `Refines HYP-480 (interval universal maximizer, primes p ≤ 19, and p = 23) and THM-137 open question 2.`

**Results note**

* **N1** (§0 row 4). `DICTIONARY (new reading; the "Baer tower" reading HYP-1348 was refuted)` → `NUMEROLOGY (MISTAKE-559); the chains are old (S80, S88, S100, the 2026-05-02 reflection), the Hughes/Lam selection is new; the "Baer tower" reading HYP-1348 was refuted`.
* **N2** (§0 row 5). `DICTIONARY with a small-numbers caveat` → `NUMEROLOGY, except (4) (PROVED)`.
* **N3** (§1, line 51). `Checked for all 1135 primes `q ≡ 3 (mod 4)` below 20000` → `Checked for all 1135 primes `q ≡ 3 (mod 4)` with `7 ≤ q < 20000``.
* **N4** (§1, line 56). `These are exactly the `H`-maximal Paley primes.` → `For q ≤ 31 these are exactly the primes at which P_q is the circulant H-maximum (FINITE-EXACT; conjecturally all q, HYP-9357).`
* **N5** (§1, line 57). After `as far from `R_q` as possible`, add `(so are q = 23, 31, where h = 3 and (2/q) = 1)`.
* **N6** (§1, line 53). Append to the prior-art line: `; the orientation-cube program THM-134/135/137/139/140 and HYP-480 (Paley vs interval, March 2026; THM-139 already notes that no single flip improves P_19)`.
* **N7** (§3, line 96). `are the only values never taken by `H` (THM-029, THM-079, THM-1370).` → `are never taken by `H` (THM-029, THM-079, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, PROVED); that no other odd value is missed is a CONJECTURE (FINITE-EXACT below 80,265, THM-4749).`
* **N8** (§3, line 99). `* No map is known, so this is a DICTIONARY.` → `* No map is known, so by MISTAKE-559 this is NUMEROLOGY.`
* **N9** (§4, line 113). `1. `P_q` is `H`-maximal (FINITE-EXACT; HYP-9357);` → `1. `P_q` is the circulant `H`-maximum (FINITE-EXACT for q ≤ 31; among all tournaments known only for q ≤ 11, A038375; conjecturally exactly {3,7,11}, HYP-9357);`
* **N10** (§4, line 122). `The repo has already rejected "Heegner numbers as tournament constants" and the LRC Heegner claims (MISTAKE-225/229/230/233, MISTAKE-576).` → `The repo has rejected Heegner-based LRC and chromatic-number claims (MISTAKE-229, MISTAKE-230, MISTAKE-568, MISTAKE-576); the historical entries HYP-1202/HYP-1220 ("Heegner numbers are tournament constants") remain CONFIRMED numerology.`
* **N11** (§9, line 194 ff.). `Every one of these with exponent `≤ 40` is either forbidden (7) or rarer than its Fermat partner `2^(2m−1) + 1`:` → `Every one of these with exponent ≤ 40 is forbidden (7) or has degree at most that of its Fermat partner 2^(2m−1) + 1:` Then keep the three bullets, and add: `For odd k ≥ 17 other than 21 and 25 the degrees are conditional on V-membership of large prime factors (THM-4749 status).`
* **N12** (§10, lines 212–213). `it is the commoner value for `k = 4, 8, 12, 16` (Fermat products);` → `it is the commoner value for k = 4, 8, 16 (Fermat products) and k = 12 (4095 = 3²·5·7·13, the 7 absorbed into 35);` And `ties at `k = 10, 11, 14` are decided by the leading coefficients.` → `the degrees tie at k = 10, 11, 14; THM-4749 decides k = 10 (1023 commoner) and k = 11 (2047 rarer); k = 14 is not computed.`
* **N13** (§10, line 220). `the perfect binary codes besides the Golay code of length 23` → `the nontrivial perfect binary codes besides the Golay code of length 23 (up to parameters; Tietäväinen, van Lint 1973)`.
* **N14** (§13 prior art, line 314). Replace with: `* **Prior art.** (i) is the circulant case of Savchenko's identity c_5 + 2c_4 = n(n²−1)(n²−9)/160 for every regular tournament (J. Graph Theory 83 (2016) 44–77; Prop. A.1 of arXiv:2403.19555); new is c_4 = n(E − m)/4 on circulants. The energy extremes and the two-mechanism reading are THM-140 (kind-pasteur S58) and THM-137; the Paley end is THM-162 and HYP-464.`
* **N15** (§13, line 317). `Its uniqueness up to multipliers is FINITE-EXACT: `n ≤ 25`, composite `n` included.` → `It is the unique maximiser up to multipliers for every odd n ≥ 5: RLT_n = R_n uniquely maximises c_4 among n-tournaments (arXiv:2403.19555 §1), c_4 = n(E − m)/4, and Aut(R_n) = Z_n makes isomorphic circulants multiplier-equivalent (FINITE-EXACT check n ≤ 25).` Under **Minimum**, add `uniquely (flat spectrum forces S to be QR-invariant)`.

**MISTAKE-659 and the banner**

* **M1** (banner). `The Epstein zeta function of `E_8` is `Z_(E8)(s) = 240 ζ(s) ζ(s−3)`, a Dirichlet series with its pole at `s = 4`.` → `The Epstein zeta function Z_(E8)(s) = Σ_(v≠0)|v|^(−2s) is 2^(−s)·240 ζ(s) ζ(s−3) (roots of norm 2; 240ζ(s)ζ(s−3) is Σ(|v|²/2)^(−s)), a Dirichlet series with its only pole at s = 4; the "s = 2" value is the continuation −5π²/6.` Make the same correction in MISTAKE-659 (`It is `240ζ(s)ζ(s−3)`, with its pole at 4.`).
* **M2** (banner and MISTAKE-659, add):
  * §1.1: the oriented Fano plane's automorphism group has order 21; `G_2` is `Aut(O)`;
  * §3.1: Bott periodicity is Bott (1957/59), not Atiyah–Hirzebruch;
  * §1.2: "regular (= self-complementary)" is false as an equivalence;
  * §2.2's "complete conflict graph `Ω(T_7)`" falls with the `K_7` error.

## Not done

* §11 (the four preprints) was not audited, as instructed.
* The `α`-decomposition text added to §13 after the audit started (`alpha_decomp*.out`, `alpha_all.out`) was not audited.
* Lev (2001) was checked from its abstract only.
* Savchenko's statements were checked through arXiv:2403.19555 (HTML) and its abstract, not through the 2016 journal text.
