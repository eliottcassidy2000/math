# Independent audit of mac-mini-2026-10-22-sierp: THM-4754, THM-4755, THM-4756, HYP-9364 and the Sierpinski note (2026-10-22)

* **Auditor:** independent refutation-first audit (subagent, Claude Opus 5.5), 2026-10-22. No edits to existing files; no git operations.
* **Scope.** The committed texts of `26967fe4cb` (THM-4754, THM-4755, HYP-9364, note) and `a7d2e23c06` (THM-4756), plus what was added while the audit ran: commit `a6c47712da` and the uncommitted edits seen at about 02:10 (THM-4756 second level, THM-4755 addendum on Conjecture 1, HYP-9364 addendum (f), note section 10).
* **Method.** Only the auditor's own code. The session's scripts were read for file formats but not run or reused. Scripts and outputs: `04-computation/experiments/sierpinski_ocf_20261022/audit/`.
* **Severity.** HIGH = wrong mathematics in a PROVED statement. MEDIUM = wrong priority or credit, or evidence overstated in canon or a hypothesis. LOW = wording, digits, citations.

**Bottom line.**
* No mathematical error was found in THM-4754, THM-4755 or THM-4756.
* THM-4756 is a correct and complete proof of THM-086: Theorem A for every `n`, the Eulerian corollary, and Theorem B qualitatively. Its mod-2 paragraph also completes THM-094, which it does not cite. No earlier complete proof exists in the repo.
* What needs correcting:
  * two novelty claims that are already in Bhattacharjee's release: the nine-vertex circulant observation, and Conjecture 1 for `d ≤ 9`;
  * a mislabelled runner-up at `d = 9`;
  * the digits of `κ_S`;
  * extrapolation precision in HYP-9364 (b)–(d), which depends on the choice of points;
  * the outright refutation of HYP-9028's universality clause, which none of the files states;
  * stale or misattributed sentences in the note.

## Summary verdicts

| # | item | verdict | severity |
|---|---|---|---|
| **A** | **THM-4756** | | |
| 1 | Step 1: `F = Σ c_m (x−1)^(n−m)`, `c_m` ordered path covers | CONFIRMED | — |
| 2 | Step 2: exact identity `Σ_X H(X)H(X^c) = Σ_R 4^|R| 2^(n−|V(R)|)`, bijection and multiplicities | CONFIRMED | — |
| 3 | Step 3: `c_2 ≡ ((−1)^n+1)H (mod 3)` | CONFIRMED | — |
| 4 | `F ≡ H·A_n (mod 3)`, the `A_n` pattern, `F ≡ (1+x)^(n−1) (mod 2)` | CONFIRMED | — |
| 5 | Does it prove all of THM-086? | CONFIRMED for Theorem A, the corollary and Theorem B (qualitative); THM-086's "30–35%" and its deletion–contraction sketch are faulty, and the two files use `c_j` in different senses | LOW |
| 6 | Step 2 identity is prior art (procgen selfie note, THM-4524 package) and is not cited | CORRECTED | LOW |
| 7 | Completes THM-094 (mod 2), which is not cited | CORRECTED | MEDIUM |
| 8 | Priority: no earlier complete proof of THM-086 in the repo | CONFIRMED | — |
| 9 | Second level `u_3 ≡ H + hc` (odd `n`), `0` (even `n`) mod 3 (added later) | CONFIRMED | — |
| 10 | "the failure for p = 5 is explained" | CORRECTED (heuristic) | LOW |
| **B** | **THM-4755** | | |
| 11 | Section 1: walk `(J+S)/d`, covariance on `1^⊥`, `ρ = Π(d²+σ²)/(d²−σ²)` | CONFIRMED (and by exact finite-`n` blow-up expectations) | — |
| 12 | Agreement with Wormald's table | CONFIRMED | — |
| 13 | "six other 9-vertex circulants": Wormald lists seven; `{1,2,3,5}` is a tournament and matches; the whole observation is Bhattacharjee's | CORRECTED | MEDIUM |
| 14 | "New relative to the release" list (note §5) | REFUTED for the circulant item; the rest needs qualification | MEDIUM |
| 15 | Lexicographic law | CONFIRMED | — |
| 16 | Tower product equals the determinant | CONFIRMED (`k ≤ 4`; the file claims `k ≤ 3`) | — |
| 17 | `κ_S = 2.79150988819834634766` | CORRECTED → `2.79150988819834776750` | LOW |
| 18 | `H(S_k)`, `v_3`, normalised values, digit counts | CONFIRMED for `k ≤ 5`; `k = 6` not recomputed | LOW |
| 19 | "So the towers' own normalised counts converge …" stated as fact | CORRECTED (conjectural) | LOW |
| 20 | Conjecture-1 addendum: `d = 9` runner-up is not `C_3[C_3]`; the note says the release "does not treat Conjecture 1" | REFUTED | MEDIUM |
| 21 | `|Aut|` divides `H` used without citing LEM-003 | CORRECTED | LOW |
| 22 | `((7+21)/(25+1))^3` in a canon theorem without a NUMEROLOGY label | CORRECTED | LOW |
| **C** | **THM-4754** | | |
| 23 | Statements 1–6 | CONFIRMED | — |
| 24 | Statement 7 no longer conditional | CONFIRMED (already fixed in `a6c47712da`) | — |
| 25 | Prior-art citations | CONFIRMED, minor caveats | LOW |
| **D** | **HYP-9364 and note §§4–6** | | |
| 26 | Normalised values (every `H` recomputed) | CONFIRMED | — |
| 27 | (b) rotational fits | Reproduced, but unstable; `e` is excluded, the identification is good only to about `10^(−2)` | MEDIUM |
| 28 | (c) Paley fit `2.718332` | Reproduced, but depends on omitting `P_27` | MEDIUM |
| 29 | (d) ratio fit `1.04953` | Reproduced; spread about `3·10^(−3)` | LOW |
| 30 | Circulant `κ/ρ` bands, rank correlation, `argmax ρ = R_n` | CONFIRMED | LOW (wording) |
| 31 | Tower deviations `0.0159, 0.0055, 0.0019` | CORRECTED (last value is 0.0018) | LOW |
| 32 | Conflict with HYP-9028: its universality clause is already refuted by Wormald's theorems, which is not said | CORRECTED | MEDIUM |
| 33 | Relation to HYP-9274 (A) and HYP-9362 (c) | CONFIRMED; "settles in principle" overstates | LOW |
| 34 | (b) labelled "Wormald's Conjecture 5" | CORRECTED | LOW |
| 35 | Addendum (f), the second-order law: the data | CONFIRMED (fits are in-sample) | LOW |
| 36 | Note §10: all regular `n = 9, 11` | CONFIRMED | — |
| **E** | **Note §7 and the rest** | | |
| 37 | Arithmetic of the number board | CONFIRMED | — |
| 38 | Numerology presented as more than numerology | CORRECTED (2 places) | LOW |
| 39 | Citations: Session25, THM-343, THM-079/115, THM-4502 "41-tail" | CORRECTED | LOW |
| 40 | Note §9 still says "THM-086 … still a proof sketch" | CORRECTED (stale) | LOW |
| 41 | Interval sequence `Σ_N`, OEIS absence, the 27-vertex table | CONFIRMED | — |

## Details

### A. THM-4756

**1–4. CONFIRMED.**

*Step 1.* Cut each ordering at its unchosen steps. A pair (ordering, `j` chosen forward steps) then corresponds bijectively to an ordered `(n−j)`-tuple of vertex-disjoint directed paths; the inverse is concatenation.

*Step 2.* The OCF is applied to `T[X]` and to `T[V∖X]`, with `H(∅) = 1` and `H = 1` on one or two vertices.
* A pair `(R_1 ⊆ X, R_2 ⊆ X^c)` is the same as an independent `R = R_1 ∪ R_2` together with a side for each cycle. That gives `2^|R|` sides times the OCF weight `2^|R|`.
* Each vertex outside `V(R)` gets a free side: `2^(n−|V(R)|)`.

*Step 3.* Modulo 3, `4 ≡ 1` and `2 ≡ −1`, and `|V(R)| ≡ |R| (mod 2)`. Removing the terms `X = ∅` and `X = V` gives `c_2 ≡ ((−1)^n − 2)H ≡ ((−1)^n + 1)H`. For `m ≥ 3`, `m! | c_m` because the paths are distinct.

*Conclusion.* For even `n`, `x − 2 ≡ x + 1`.

Script `a_thm4756_check.c`, output `a_thm4756_check.out`. Zero failures in each of seven checks:
1. the exact integer identity of Step 1;
2. the exact identity of Step 2, with directed cycles counted by an independent Held–Karp DP;
3. the OCF;
4. `c_2 (mod 3)`;
5. `F ≡ H·A_n (mod 3)`;
6. `F ≡ (1+x)^(n−1) (mod 2)`;
7. `u_3`.

The checks were run on:
* every labelled tournament with `n ≤ 6`, brute force over all `n!` orderings;
* all classes at `n = 7` (456) and `n = 8` (6880), brute force;
* all 191,536 classes at `n = 9`, by subset DP (53 s);
* random tournaments: `n = 9` (20, brute force), `n = 10` (3000), `n = 11` (1000), `n = 12` (15) and `n = 13` (10), by DP.

`A_n (mod 3)` and `A_n (mod 2)` follow the stated pattern for `n ≤ 80` (`a_eulerian_mod3.py`, `.out`).

**5. What THM-086 states.**
* **Theorem A:** `3 | c_j` for `j < val(n)`, where THM-086's `c_j` is the coefficient of `(x−1)^j`. The proof covers every `n ≥ 1`, not only `n ≥ 3`.
* **Corollary:** `3 | A(n,k) ⇒ 3 | F_k`. Proved.
* **Theorem B:** follows. The coefficient of `(x−1)^val(n)` is `≡ H` (odd `n`) or `≡ −H` (even `n`), and `H = 1` for the transitive tournament.
* **Theorem B's "approximately 30–35%"** is exactly the frequency of `3 | H`. Exact labelled values are 25% (`n = 3, 4`), 41.4% (`n = 5`) and 29.7% (`n = 6`), so `n = 3, 4, 5` fall outside the stated range. Over classes it is 37.3%, 35.1% and 34.6% at `n = 7, 8, 9`.
* **The deletion–contraction sketch is invalid as written.** It assumes "T/e is an (n−1)-vertex tournament", but under THM-083's contraction (`w` takes the in-arcs of `u` and the out-arcs of `v`) the cyclic triangle `u→v→z→u` already gives arcs both ways between `w` and `z`. THM-4756 does not use this route.
* **Notation clash.** THM-086's `c_j` is THM-4756's `c_(n−j)`. Both files call them `c`.
* **THM-085 is consistent.** Its `c_0, c_1, c_2` are the path-cover numbers `c_n, c_(n−1), c_(n−2)`, divisible by `n!`, `(n−1)!` and `(n−2)!`.

**6. Prior art for Step 2.** `05-knowledge/results/procgen_selfie_20261001_selfie_tournaments.md` (THM-4524 package, 2026-10-01), §2, already states: "Check: U_T(1,1) = sum over ordered 2-block covers of H(T[B1]) H(T[B2]) = sum_S 4^|S| 2^(N−|V(S)|) (VERIFIED N <= 7)". It is also the principal specialisation `U_T(1^2)` of the Grinberg–Stanley formula. The mod-3 use and the link to `c_2` are new.

**7. THM-094.** THM-094 (kind-pasteur S38) states `F ≡ (1+x)^(n−1) (mod 2)` with status "PROOF SKETCH". THM-4756's mod-2 paragraph (`m! | c_m`, Rédei) proves it in two lines, but THM-4756 does not cite THM-094.

**8. Priority.** Searched for "THM-086", "Taylor zero", "almost-tournament claim", "Eulerian conjecture", "forward-edge polynomial", "U_T(1,1)", and "path cover" together with "mod 3".
* THM-086 is still a sketch everywhere: INV-123 says "PROOF SKETCH COMPLETE", and the 2026-03-10 synthesis lists the almost-tournament claim as open.
* There is no earlier complete proof.
* A brief web search found no literature statement of `F ≡ H·A_n (mod 3)`. This search is not exhaustive.

**9. Second level (in `a6c47712da`). CONFIRMED.**
* Inclusion–exclusion gives `c_3 = U(1,1,1) − 3U(1,1) + 3H`.
* Modulo 9, the only terms `6^|R| 3^(n−|V(R)|)` of `U(1,1,1)` that survive are those with `|R| = 1` and `V(R) = V`, i.e. Hamiltonian cycles.
* Then `6u_3 ≡ 6(H + hc)` (odd `n`) or `0` (even `n`) mod 9.
* Check 7 of `a_thm4756_check` passes on every case listed above.

**10.** The `p = 5` sentence is a heuristic, not a proof.

### B. THM-4755

**11. CONFIRMED.**
* `P = (I + 2A)/d = (J + S)/d`. The identity `(I−Q)^(−1) + (I−Qᵀ)^(−1) − I = (I−Q)^(−1)(I−QQᵀ)(I−Qᵀ)^(−1)` holds, and with `QQᵀ = −S²/d²` the stated ratio follows.
* Independent test (`b_blowup_exact.py`, `.out`). The exact finite-`n` value of `E[H'(T(T*,n))]/g'(dn)`, computed as P(walk balanced)/P(i.i.d. balanced) and extrapolated in `1/n`:

  | `T*` | extrapolated | `ρ` |
  |---|---|---|
  | `C_3` | 2.000000 | 2 |
  | `R_5` | 2.315794 | 44/19 = 2.315789 |
  | `R_7` | 2.463026 | 2.463038 |
  | `P_7` | 2.370318 | 64/27 = 2.370370 |

* Remark: Wormald's Theorem 4 and this computation evaluate the same limit. So once the local-limit step is written rigorously, `d/𝓗 = ρ²` follows wherever `𝓗 ≠ 0`. The non-vanishing still needs a direct Hessian argument (Bhattacharjee's Proposition 3.1).

**12. CONFIRMED** (`b_rho_check.gp`, `.out`).
* Exact equality for `R_3, R_5, R_7, R_9, R_11, R_13` and `R'_9 = C_3[C_3]` (5488/2197).
* Wormald's 7-vertex example: its completion to a regular tournament is unique (checked by brute force), with `ρ = 3088/1287`.
* `ρ_T(R_21), ρ_T(R_41), ρ_T(R_81), ρ_T(R_101)` = 2.721518859…, 2.786703868…, 2.820801899…, 2.827747088…, matching every printed digit.
* Wormald's `ρ_R(R_81) = 2.855842467` should read …468. Bhattacharjee already notes this.

**13.** Wormald (p. 6) lists **seven** sets `D`.
* The six sets `{1,3,4,5}`, `{1,2,4,5}`, `{2,3,4,5}`, `{1,3,5,6}`, `{1,4,5,6}`, `{2,3,5,6}` each contain both `r` and `9 − r`, so they are not tournaments; their values lie below `(10/8)^4 = 2.4414`.
* `D = {1,2,3,5}` is a tournament with `ρ = 451696/183937 = 2.4557104`, matching Wormald's 2.4557. THM-4755 omits it.
* All of this, including the `{1,2,3,5}` value and the Jensen-bound argument against Wormald's "1.63 g(n)", is the Remark on Wormald's circulants in Bhattacharjee's preprint (`preprint/main.tex` of the cited repo, release v1.0.0).

**14. The "new relative to the release" list** (release v1.0.0, published 2026-10-14T09:27Z, MIT; confirmed via the GitHub API).
* **Circulant observation:** not new (item 13).
* **Spectral moment expansion:** a term-by-term Taylor expansion of the release's eq. (trace) `log ρ = Σ log((d²+s_k)/(d²−s_k))`, so new only in form.
* **Lexicographic law, tower product, self-blow-up principle:** grep of the preprint finds no lexicographic, composition, substitution or tower content, so these are new.
* **The note's description of the release** omits its upper bound `U_d` (limit 2.8574012, via Grzesik–Král'–Lovász–Volec) and its exhaustive `d ≤ 9` table with `R_d` maximal.
* The note's "settling Conjecture 6 modulo his Proposition 3.1" understates the release: Proposition 3.1 is proved on paper (implicit function theorem plus Schur complement) and is only unformalised in Lean. THM-4755's own status line states this correctly.

**15. CONFIRMED.**
* `S(Q[B]) = S_Q⊗J + I⊗S_B` acts as `b·S_Q` on `u⊗1` and as `S_B` on `u⊗w` with `w ⊥ 1`.
* Numerically, `ρ(det) = law` to 20 digits for `R5[C3]`, `C3[R5]`, `P7[C3]`, `C3[P7]` and `R5[P7]`.
* Cosmetic: the displayed formula has an unbalanced parenthesis, and its product runs over the pairs `±iσ`.

**16. CONFIRMED.** `ρ(S_k)` from the determinant equals the level product exactly for `k ≤ 4`, both for `C_3^[k]` and for `Cay(Z/3^k, sign 1)` (`b_rho_check.out`).

**17. CORRECTED.**
* `κ_S = 2.7915098881983477674999548…` (`b_kappaS_mpmath.py` with `artanh`, and PARI at 60 digits).
* THM-4755 prints `2.79150988819834634766`, which is wrong from the 15th decimal. The note's `2.79150988819834…` is right.
* The `log κ_S` series and the excess `4.20·10^(−7)` over `(5488/2197)e^(1/9)` are right. That excess is `≈ 6·3^(−15)`: `e^(1/9)` is the leading term of the tail product.

**18. CONFIRMED for `k ≤ 5`.**

*Evidence.*
* `H(S_3) = 417382500592116023859` three ways:
  * inclusion–exclusion on `Z_27` over rotation orbits (`b_ie_circulant.c`, `b_ie_S3.out`, 10.9 s);
  * the auditor's top-level formula `H(C_3[B,B,B]) = 3 Σ_a c_a (c_a² + c_a c_(a+1) + c_(a+1)²)`, with `c(S_2)` from a direct DP;
  * the auditor's shape-sum law.
* The shape-sum law was validated against brute force on 40 random `C_3[B_0, B_1, B_2]`. The `c_m(S_3)` it produces give a palindromic, nonnegative `F(S_3, x)` summing to `27!`.
* `H(S_4)` (98 digits) comes from that law (`b_tower_H.py`, `.out`).
* `H(S_5)` (404 digits) comes from the same law modulo 64 primes plus CRT (`b_tower_k5_modp.py`, `.out`).
* All values equal Session24's and the session's.

*Values.*
* `v_3(H(S_k)) = 1, 5, 13, 40, 121`.
* `κ(S_k)` = 2, 2.2285714, 2.5723575, 2.7135298, 2.7648911.
* Digit counts 1, 4, 21, 98, 404 (and 1554 for `k = 6`, read from the session's value).

*Not independently checked.* `H(S_6)` was not recomputed. Its `v_3 = 364` is consistent with `3^364 = |Aut|` dividing `H` (LEM-003).

**19. Overstatement in canon.**
* Richardson through `k = 3..6` reproduces 2.791509575761 (`−3.12·10^(−7)`).
* Through `k = 4..6` it gives `−5.2·10^(−6)`; through `k = 2..6`, `+2.3·10^(−7)` (`d_tower_fit.py`, `.out`).
* Strong evidence, but not a proof. Section 4's last bullet states convergence as fact.

**20. Conjecture-1 addendum.**
* My enumeration (`d_regular_rho.py`, `d_regular_rho11.py`) reproduces `ρ` for `d = 5, 7, 9, 11`, including the runner-ups 2.3993783994 (`d = 7`) and 2.5682251025 (`d = 11`). `d = 13` was not checked.
* At `d = 9` the runner-up `2.4981158008 = 443488/177529` is **not** `C_3[C_3]`. `C_3[C_3]` is `5488/2197 = 2.4979517524`, third of 15; I confirmed it is the module-decomposable class. Bhattacharjee's table agrees.
* "Bhattacharjee's release does not treat Conjecture 1" is false. The preprint discusses Wormald's Conjecture 1, proves `((d+1)/(d−1))^((d−1)/2) ≤ ρ_T ≤ U_d` over regular tournaments, and states that an exhaustive computation shows `R_d` maximal among regular tournaments for `d ≤ 9`. Only `d = 11, 13` are new.

**21.** The `v_3` table (§4) and the note's `Aut(Σ_N)` paragraphs use "`|Aut|` divides `H`", which is LEM-003 (freeness on Hamiltonian paths). LEM-003 is not cited.

**22.** Section 3, bullet 2 presents `((7+21)/(25+1))^3` without a label. The note's §7 correctly says this reading is NUMEROLOGY.

### C. THM-4754

**23. CONFIRMED** (`c_thm4754_check.py`, `.out`).
1. `Cay(Z/3^k, sign 1)` is a tournament equal to `C_3^[k]` under the base-3 digit map, with the lowest digit as the top level (`k ≤ 5`).
2. `x ↦ 3x+b` is an isomorphism onto class `b`, and class `b` dominates class `b+1` (`k ≤ 6`). `−1/2 = …111` and `(3x+1)^k(x) = 3^k x + (3^k−1)/2` both hold.
3. For every unit `λ` (`k ≤ 5`): `λ ≡ 1` gives an automorphism and `λ ≡ 2` reverses every arc.
4. Arc-flip law: 0 violations in 20,000 random pairs `x ≡ y (mod 2^t)`, `t ≤ 40`, for both `C` and the shortcut `T`. As a control, 964/2000 pairs without the congruence violate it. The proof is correct: equal parity for `t` steps, then the difference is `3^o (y−x)/2^e` and `1/2 ≡ −1 (mod 3)`.
5. Steiner reduction: equals `(−1)^(n−1) Σ C(n−1,i) x_i` for `n < 60`. It depends only on the endpoints exactly when `n − 1 ∈ {0, 1, 3, 9, 27}`.
6. `2^k = 3J_k + (−1)^k` holds. The odd `n < 10^6` with `3n+1` a power of 2 are exactly `J_(2m)`.

**24.** After THM-4756, statement 7 is unconditional. Commit `a6c47712da` already rewrote it ("By THM-4756"). No further edit needed.

**25. Citations.** Checked: line 253 of the pseudoprimes note, THM-4743, Session24 PATHCOVERS (`|Aut| = 3^((3^k−1)/2)`; values to `k = 5` in `receipt.json`), HYP-1028..1030 (rows in INDEX-HISTORICAL), THM-467/T274, THM-1425 (Jacobsthal), THM-4517 and the procgen trunk notes. Minor caveats:
* HYP-9353 names `C_3[C_3]` but never gives `Z_9{1,3,4,7}`.
* THM-1425 does not contain `I(P_4, 2) = 21`, and its own path indexing is off by one. That is pre-existing.

### D. HYP-9364 and note §§4–6

**26. CONFIRMED.**
* Every `H` used was recomputed by inclusion–exclusion (`d_ie_small.out`, `d_ie_large.out`): `R_n` for every odd `n ≤ 31`, and `P_7`, `P_11`, `P_19`, `P_23`, `P_31`. All match the repo files. Not recomputed: `R_33` (session's `R33.out`) and `P_27` (on `F_27`, from `P27.out`).
* `κ(R_n)` = 2.551967, 2.576702, 2.597757, 2.615845, 2.631571, 2.645362, 2.657557 for `n = 19..31`, and 2.668417 at `n = 33` (session's `R33.out`).
* `κ(P_q)` = 2.527148, 2.556980, 2.595512 for `q = 19, 23, 31` (`d_fits.out`).

**27. (b) rotational fits** (`d_fits.py`, `.out`).
* The quoted fits reproduce exactly: 2.8535955 (`n = 27..31`) and 2.8571531 (`n = 25..31`).
* They are not stable:
  * through `n = 23..31` the fit gives 2.8707, and through `21..31` it gives 2.9613;
  * with `H(R_33)` the fits through `29..33`, `27..33` and `25..33` give 2.85386, 2.85503 and 2.84838.
* So "bracketing" is not evidence at `10^(−3)`.
* The data do exclude `e`. `n²·Δκ/2` is 4.95, 5.10, 5.21, 5.32, 5.41, 5.49, 5.56 for `n = 20, …, 32` and still rising. Reaching `e` would need a total further increase of `e − κ(R_33) = 0.050`, while increments of this size sum to at least about 0.17.
* Conclusion: the contradiction with HYP-9028 is well supported, but the value `cosh 1/cos 1` is supported only to about `10^(−2)`.

**28. (c) Paley fit.**
* The fit through `q = 19, 23, 31` reproduces 2.7183322.
* The repo also has the Paley tournament on `F_27` (`paley_extremal_20261021/P27.out`, `H = 418425945688463432031`). HYP-9274 (A) covers prime powers, and HYP-9362 (c) itself quotes `q = 27`.
* Including it gives 2.72566 (`q = 19, 23, 27, 31`) and 2.72050 (`q = 23, 27, 31`). The `5·10^(−5)` agreement is therefore an artifact of the chosen points.
* `(1 − κ/e)·q` = 1.336, 1.365, 1.386, 1.400 rises toward HYP-9274's 3/2.

**29. (d) ratio fit.** 1.0495313 reproduces. With `q = 27` the fits give 1.04755 and 1.04894, so the spread is about `3·10^(−3)`.

**30. CONFIRMED** (`d_circ_band.py`, `.out`).

| `n` | classes | `κ/ρ` range | Spearman |
|---|---|---|---|
| 13 | 6 | [0.9245, 0.9643] | +0.600 |
| 15 | 16 | [0.9269, 0.9680] | −0.035 |
| 17 | 16 | [0.9375, 0.9722] | +0.118 |
| 19 | 30 | [0.9425, 0.9791] | −0.009 |
| 21 | 88 | [0.9454, 0.9781] | −0.025 |
| 23 | 94 | [0.9505, 0.9819] | −0.044 |
| 25 | 208 | [0.9537, 0.9817] | +0.092 |

* `argmax ρ = R_n` among circulants for every `n = 7..25`.
* `n(1 − κ/ρ)` stays between 0.98 and 1.16 at the lower edge and about 0.46 at the upper edge, so "roughly like `1/n`" is fair.
* Wording: the "208 classes" are multiplier classes of connection sets, not isomorphism classes. At `n = 9` the data list both `Z_9{1,3,4,7}` and `Z_9{1,4,6,7}`, and both are `≅ C_3[C_3]`.

**31.** The deviations `1 − κ/ρ` are 0.01586, 0.005451, 0.001836 at `n = 81, 243, 729`, so the last should read 0.0018. `n·(1 − κ/ρ)` = 1.285, 1.325, 1.339, consistent with "about 1.35/n".

**32. HYP-9028.**
* HYP-9028 conjectures `exc → e` for `R_n` and Paley. Its "universality upgrade" adds that "every sequence of regular tournaments has exc → e".
* That clause is already **refuted**, independently of HYP-9364, by Wormald's Theorems 1, 2 and 5. `R(C_3, n)` (regular blocks) is supported on regular tournaments on `3n` vertices with `E[H'] ∼ 2e^(1/3) g'(3n)`, where `2e^(1/3) = 2.7912 > e`. So for infinitely many `n` some regular tournament has `exc ≥ 2.79 − o(1)`. The exact tower value `exc(S_6) = 2.7826 > e` also exceeds `e`.
* HYP-9364 and the note withdraw only "the rotational half".
* HYP-9364 also predicts that HYP-9028's control family `C_3[R_k,R_k,R_k]` tends to `2(cosh(1/3)/cos(1/3))^3 ≈ 2.7917`, not to `e`.

**33.**
* HYP-9274 (A), Paley `→ e`: consistent.
* HYP-9362 (c) records the conflict "ratio → 1" and makes no prediction; the value 1.0506 is new.
* "This settles in principle the conflict" should be conditional.

**34.** "(Wormald's Conjecture 5)" is imprecise. Wormald's Conjecture 5 is `lim ρ_R(R_d) = lim h(n)/g(n) ≈ 2.855958`. The statement about `H(R_n)/g(n)` is the unnumbered conjecture after it, with path analogues; the value `cosh 1/cos 1` is Bhattacharjee's.

**35. Addendum (f)** (`d_second_order.py`, `.out`).
* `corr(y, M4)` = 0.982, 0.983, 0.958, 0.987, 0.992 for `n = 17..25`. The per-`n` fits give `a = 0.315–0.349` and `b = 4.58–5.33`; the file's "b ≈ 4.3–5.3" is slightly wide.
* Spearman with the corrected predictor is 0.735–0.884. The maximum relative error is `2.4–2.6·10^(−3)` at `n = 19, 23, 25` and `1.1·10^(−2)` at `n = 21`.
* The maximiser is found every time, but with two fitted parameters per `n` (in-sample).
* "`R_13` first wins, by 0.22%" is the margin over the flattest class `{2,4,5,6,10,12}`. The runner-up circulant `{3,4,5,6,11,12}` trails by only 0.10%.

**36. CONFIRMED** (`hcount.c`, `d_regular_principle.py`, `.out`).
* `n = 11` (1223 classes): `κ/ρ ∈ [0.9166, 0.9804]`, Spearman `−0.743`.
* `n = 9`: Spearman `−0.710` against the note's `−0.72` (tie handling).

### E. Note §7 and the rest

**37. CONFIRMED** (`e_numberboard.py`, `.out`).
* **Collatz:** 27 takes 111 steps, 41 of them odd and 70 halvings; `U(9) = 7`, `U(21) = 1`, `C(9) = 28`.
* **Number identities:** `3^k = k^3` only at `k = 3`; `k³+1 = (k+1)Φ_6(k)`, with `Φ_6(k) = 1 + 2C(k,2)`, `7 = Φ_6(3) = 28/4` and `21 = Φ_6(5) = 126/6`; `{111_2, 111_3, 111_4} = {7, 13, 21}` and `(3x+1)^3(0) = 13`.
* **Tournament counts:** `c_3(regular, 2N+1) = Σ_(i≤N) i²`; Watson's squares occur only at `N ∈ {0, 1, 24}` (checked to `2·10^6`) and `Σ_(i≤24) i² = 70²`; `H(R_7) = 175` and `H(P_7) = 189`.
* **H-values:** `H = 13` occurs at `n = 5`; the profiles for `H = 7` and `21` are exactly as listed.

**38.** Two places present numerology as more than that:
* THM-4755 §3, bullet 2 (item 22);
* note §4: "The factors are exactly the numbers of the prompt … `(7 + 21)/(25 + 1)`".

The appearance of 28/26 is forced. Its decomposition into the prompt's pairs is not.

**39. Citations.**
* "Session25 shows … and on strong classes of orders 5–8, `I(Ω, −1)` is never 0": the second half is from `repunit_rarity_hyp9357_20261021.md` §8 (opus-2026-10-21-S51), not Session25.
* THM-343 is "Lean-formalised modulo 7 cited axioms", which include the OCF.
* THM-079 and THM-115 treat `H = 21` without Busch or crevasses. That mechanism is in THM-1370 and the camion_busch note only.
* THM-4502's "41-tail" is the orbit value `27 → 82 → 41`, not the 41 odd steps.
* HYP-1043 is accurate. "The F_21 note" is `05-knowledge/results/mersenne_switch_parity_f21_compression_20261006.md`.
* The camion_busch list, Theorem P, the least-resistance note, the-two-staircases and Session29 are accurate.

**40.** Note §9, last bullet: "Decide THM-086 in full; it is still a proof sketch." This contradicts §1 and THM-4756.

**41. CONFIRMED.**
* The first 20 terms of `H(Σ_N)` match the note's first 14 (`hcount.c interval`, `e_interval_seq.out`).
* OEIS has no match for `1,3,5,11,45,159,621,3159,14877`, for `3,3159,417382500592116023859`, or for `3,15,175,3267,93027,3711175`.
* 27-vertex table: `ρ(R_27) = 2.751119`; `(14/13)^13 = 2.620601`; `κ = 2.631571, 2.572358, 2.578788`.
* The `Aut(Σ_N)` recursion argument is sound: the three blocks are the maximal modules of a strong `C_3[A,B,C]`. It was not recomputed with nauty.

## Required edits (exact old → new)

**THM-4756** (`01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md`)

1. Status, last line:
   * old: `  The failure for p = 5 noted in THM-086 is explained: modulo 5 the weights 2^(-|C|) depend on the cycle length.`
   * new: `  The failure for p = 5 noted in THM-086 is consistent with this (heuristic, not a proof): modulo 5 the weights 2^(-|C|) depend on the cycle length.` followed by these new status lines:
     * `  PRIOR ART: the exact identity of step 2 is already stated and checked (N <= 7) as "U_T(1,1) = sum over ordered 2-block covers of H(T[B1]) H(T[B2]) = sum_S 4^|S| 2^(N-|V(S)|)" in 05-knowledge/results/procgen_selfie_20261001_selfie_tournaments.md (THM-4524 package, section 2); its mod-3 use is new here.`
     * `  The mod-2 paragraph completes THM-094 (F = (1+x)^(n-1) mod 2, previously a proof sketch). Theorem B of THM-086 also follows (the coefficient of (x-1)^val(n) is +-H mod 3, and H = 1 for the transitive tournament); its "30-35%" is the frequency of 3 | H (labelled: 25%, 25%, 41.4%, 29.7% at n = 3..6).`
     * `  NOTATION: THM-086's c_j is the coefficient of (x-1)^j, i.e. c_(n-j) here. THM-086's deletion-contraction sketch assumed that T/e is a tournament, which fails for THM-083's contraction; this proof does not use it.`
   * Audit trail: the "independently audited" line is for the session to add.
2. `related:` block:
   * old: `  - 01-canon/theorems/THM-085-f-omega-mod9-universal.md`
   * new: `  - 01-canon/theorems/THM-085-f-omega-mod9-universal.md` followed by `  - 01-canon/theorems/THM-094-fk-mod2-tournament-independent.md` and `  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md`
3. Body:
   * old: `every tournament carries the binary Sierpinski row.`
   * new: `every tournament carries the binary Sierpinski row. This proves THM-094.`

**THM-4755** (`01-canon/theorems/THM-4755-blow-up-gains-skew-spectrum-lexicographic-law-ternary-tower.md`)

4. Title:
   * old: `the limit is kappa_S = 2.79150988819834634766...`
   * new: `the limit is kappa_S = 2.79150988819834776750...`
5. Section 3 display:
   * old: `= 2.79150988819834634766…`
   * new: `= 2.79150988819834776750…`
6. Status:
   * old: `Wormald's six other 9-vertex` / `'circulants' are not tournaments as printed (each connection set contains both r and 9 - r), and their values lie below the` / `Jensen bound 2.44, so they cannot be regular tournaments.`
   * new: `Wormald lists seven 9-vertex` / `'circulants' (p. 6): six are not tournaments as printed (each connection set contains both r and 9 - r) and lie below the` / `Jensen bound 2.44; the seventh, D = {1,2,3,5}, is a tournament with rho = 451696/183937 = 2.45571, matching his 2.4557.` / `This observation is Bhattacharjee's (Remark on Wormald's circulants in the release preprint), not new here.`
7. Section 3:
   * old: ``* The second is `((3³ + 1)/(3³ − 1))^3 = ((7 + 21)/(25 + 1))^3`.``
   * new: ``* The second is `((3³ + 1)/(3³ − 1))^3`. (Reading 28 and 26 as `7 + 21` and `25 + 1` is NUMEROLOGY.)``
8. Section 4:
   * old: ``The values for `k ≤ 3` also come from inclusion–exclusion (`ternary_H.out`).``
   * new: ``The values for `k ≤ 3` also come from inclusion–exclusion (`ternary_H.out`). The values for `k ≤ 5` were recomputed independently in the audit (`sierpinski_audit_20261022.md`). Since `Aut` acts freely on Hamiltonian paths (LEM-003), `v_3(H) ≥ v_3|Aut| = (3^k − 1)/2`.``
9. Section 4:
   * old: `* So the towers' own normalised counts converge to the limit of their blow-up gains (HYP-9364 (a)).`
   * new: `* This supports the conjecture HYP-9364 (a) that the towers' own normalised counts converge to the limit of their blow-up gains. Richardson through 3 or 5 levels gives −5.2·10^(−6) and +2.3·10^(−7), so the agreement is robust to about 5–6 digits; it is not a proof.`
10. Addendum, table row:
    * old: ``| 9 | 15 | 2.5475143943 | 2.4981158008 (`C_3[C_3]`) |``
    * new: ``| 9 | 15 | 2.5475143943 | 2.4981158008 (443488/177529; `C_3[C_3]` = 5488/2197 = 2.4979517524 is third) |``
11. Addendum, first sentence:
    * old: ``has a unique maximum at the rotational `R_d` (`wormald_conj1.out`, `wormald_conj1_13.out`):``
    * new: ``has a unique maximum at the rotational `R_d` (`wormald_conj1.out`, `wormald_conj1_13.out`). The cases `d ≤ 9` are already in Bhattacharjee's release (its table of the regular tournaments on `d ≤ 9` vertices, with `R_d` maximal); `d = 11, 13` are new:``

**HYP-9364** (`05-knowledge/hypotheses/HYP-9364-self-blow-up-principle-tower-limit-rotational-cosh-over-cos.md`)

12. Title:
    * old: `(Wormald's Conjecture 5);`
    * new: `(Wormald's conjecture following his Conjecture 5, path form; the value is Bhattacharjee's);`
13. Status:
    * old: `0.0159, 0.0055, 0.0019 at n = 81, 243, 729`
    * new: `0.0159, 0.0055, 0.0018 at n = 81, 243, 729`
14. Status:
    * old: `at n = 25 (208 classes)`
    * new: `at n = 25 (208 multiplier classes)`
15. Status:
    * old: `  - (b) exact H(R_n) for odd n <= 31. Polynomial fits in 1/n give 2.85360 (3 points, n = 27..31) and 2.85715 (4 points,` / `    n = 25..31), bracketing cosh(1)/cos(1) = 2.85596. This contradicts the rotational half of HYP-9028, which predicts e.`
    * new: `  - (b) exact H(R_n) for odd n <= 33. Polynomial fits in 1/n give 2.85360 (3 points, n = 27..31) and 2.85715 (4 points,` / `    n = 25..31), but are not stable: 5 and 6 points give 2.8707 and 2.9613, and with H(R_33) the 3-, 4- and 5-point fits give` / `    2.8539, 2.8550, 2.8484. The limit is fixed only to about 1e-2. The data do exclude e: n^2 (kappa(R_(n+2)) - kappa(R_n))/2` / `    is 5.56 at n = 31 and still rising, while reaching e needs only 0.050 more. This contradicts HYP-9028.`
16. Status:
    * old: `  - (c) exact H(P_q), q = 19, 23, 31: the 3-point fit gives 2.718332, against e = 2.718282.`
    * new: `  - (c) exact H(P_q), q = 19, 23, 31: the 3-point fit gives 2.718332, against e = 2.718282. With the Paley tournament on F_27` / `    (P27.out) the fits give 2.72566 (q = 19..31) and 2.72050 (q = 23, 27, 31): the 5e-5 agreement depends on the points chosen.`
17. Status:
    * old: `gives 1.04953, against 1.050648.`
    * new: `gives 1.04953, against 1.050648 (with q = 27: 1.04755 and 1.04894).`
18. Body:
    * old: `and agree to 7 digits.`
    * new: `and agree to about 6 digits (Richardson through 3, 4, 5 levels: −5.2e-6, −3.1e-7, +2.3e-7).`
19. Body:
    * old: `* This settles in principle the conflict recorded in HYP-9362 (c).`
    * new: `* If true, it resolves the conflict recorded in HYP-9362 (c).`
20. Body:
    * old: ``* The rotational half of HYP-9028 (limit `e`) should then be withdrawn, and HYP-9274 (A) (Paley limit `e`) kept.``
    * new: ``* The rotational half of HYP-9028 (limit `e`) should then be withdrawn, and HYP-9274 (A) (Paley limit `e`) kept. HYP-9028's universality clause ("every sequence of regular tournaments has exc → e") is already refuted, independently of this hypothesis: by Wormald's Theorems 1, 2 and 5, the regular blow-ups `R(C_3, n)` have `E[exc] → 2e^(1/3) = 2.7912`, so regular tournaments with `exc ≥ 2.79 − o(1)` exist for infinitely many `n`; also `exc(S_6) = 2.7826 > e`. This hypothesis predicts HYP-9028's control family `C_3[R_k, R_k, R_k]` → `2(cosh(1/3)/cos(1/3))^3 ≈ 2.7917`.``

**Note** (`05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md`)

21. §0 table:
    * old: `| data agree to 1e-3, 5e-5, 1e-3 |`
    * new: `| selected 3–4-point fits land within 3e-3, 5e-5, 1.1e-3; the fit-to-fit spread is about 1e-2 (audit); the data exclude e for R_n |`
22. §1:
    * old: ``Session25 shows that every odd residue modulo every `M` is a strong `H`-value, and on strong classes of orders 5–8, `I(Ω, −1)` is never 0.``
    * new: ``Session25 shows that every odd residue modulo every `M` is a strong `H`-value, and `repunit_rarity_hyp9357_20261021.md` (opus-2026-10-21-S51) that on strong classes of orders 5–8, `I(Ω, −1)` is never 0.``
23. §2:
    * old: ``  * THM-343 (`H = 7`, Lean-formalised);``
    * new: ``  * THM-343 (`H = 7`, Lean-formalised modulo cited axioms, including the OCF);``
24. §2:
    * old: `  * THM-079, THM-115 and THM-1370 (all `n`, via Busch's strong floor`
    * new: `  * THM-079 and THM-115 (`H = 21`), and THM-1370 (both values, all `n`, via Busch's strong floor`
25. §4:
    * old: ``The second factor is `(27 + 1)/(27 − 1) = (7 + 21)/(25 + 1)`, cubed.``
    * new: ``The second factor is `(27 + 1)/(27 − 1)`, cubed; reading 28 and 26 as `7 + 21` and `25 + 1` is NUMEROLOGY.``
26. §5:
    * old: `settling Conjecture 6 modulo his Proposition 3.1;`
    * new: `settling Conjecture 6; its Proposition 3.1 (Wormald's analytic H equals a determinant) is proved on paper and checked by exact computation, but not formalised in Lean;`
27. §5:
    * old: ``* the spectral moment expansion `log ρ = −Σ tr(S^(4m+2))/((2m+1) d^(4m+2))`;``
    * new: ``* the spectral moment expansion `log ρ = −Σ tr(S^(4m+2))/((2m+1) d^(4m+2))` (a Taylor expansion of the release's `log ρ = Σ log((d²+s_k)/(d²−s_k))`);``
28. §5:
    * old: `* the observation that Wormald's six other 9-vertex "circulants" are not tournaments as printed and lie below the Jensen bound;`
    * new: delete this line, and add to the list of the release's claims: `* that six of Wormald's seven 9-vertex "circulants" are not tournaments, while the seventh, {1,2,3,5}, gives 451696/183937, matching his value; an upper bound U_d (limit 2.8574012); and an exhaustive table for d ≤ 9 with R_d maximal;`
29. §6:
    * old: `Fits in `1/n` over `n ≤ 31` give 2.8536 and 2.8572, bracketing 2.85596.`
    * new: `Fits in `1/n` over `n ≤ 33` give 2.848–2.871 depending on the points (2.8536 and 2.8572 for the top 3 and 4 points up to 31): consistent with 2.85596 at the 1e-2 level, while the increments exclude `e`.`
30. §6:
    * old: `The 3-point fit gives 2.718332, against 2.718282 (HYP-9274 (A)).`
    * new: `The 3-point fit (19, 23, 31) gives 2.718332, against 2.718282; adding `P_27` moves the fits to 2.7205–2.7257 (HYP-9274 (A)).`
31. §6:
    * old: `The 3-point fit gives 1.04953 (HYP-9362 (c), now with a conjectured value).`
    * new: `The 3-point fit gives 1.04953, or 1.0476–1.0489 with `q = 27` (HYP-9362 (c), now with a conjectured value).`
32. §7:
    * old: `Prior analysis is THM-4504 and THM-4502 (the 41-tail family).`
    * new: `Prior analysis is THM-4504 and THM-4502 (whose "41-tail" is the orbit value 27 → 82 → 41, not the 41 odd steps).`
33. §9:
    * old: `* Decide THM-086 in full; it is still a proof sketch.`
    * new: `* Point THM-086, THM-094 and INV-123/INV-124 to THM-4756.`
34. §10, table row:
    * old: ``| 9 | 15 | 2.5475143943 | 2.4981158008 (`C_3[C_3]`) |``
    * new: ``| 9 | 15 | 2.5475143943 | 2.4981158008 (443488/177529; `C_3[C_3]` = 2.4979517524 is third) |``
35. §10:
    * old: `Bhattacharjee's release does not treat Conjecture 1.`
    * new: `Bhattacharjee's release treats Conjecture 1 for regular bases: it proves ((d+1)/(d−1))^((d−1)/2) ≤ ρ ≤ U_d and checks exhaustively that R_d is maximal for d ≤ 9; d = 11, 13 are new here.`
36. §10:
    * old: ``*Consequence.* Since `|Aut|` divides `H`,``
    * new: ``*Consequence.* Since `|Aut|` divides `H` (LEM-003),``

**Recommended outside the audited files** (the auditor may not edit them):
* **THM-086:** status "PROVED by THM-4756". Its deletion–contraction sketch has a gap (`T/e` is not a tournament).
* **THM-094:** status "PROVED (THM-4756, mod-2 paragraph)".
* **HYP-9028:** mark the universality upgrade REFUTED (Wormald Theorems 1, 2 and 5; item 32), and the rotational half as contradicted by HYP-9364 and the increment data.
* **INVESTIGATION-BACKLOG:** close INV-123 and INV-124.

## Files (all under `04-computation/experiments/sierpinski_ocf_20261022/audit/`)

**THM-4756**
* `a_thm4756_check.c` / `.out`: steps 1–3, mod 2, `u_3`, `3 | H` frequencies.
* `a_eulerian_mod3.py` / `.out`

**THM-4755, statement 1 and Wormald**
* `b_rho_check.gp` / `.out`: Wormald's table, the 9-vertex circulants, towers, lexicographic law.
* `b_kappaS_mpmath.py` / `.out`
* `b_blowup_exact.py` / `.out`: finite-`n` blow-up expectations.

**Tower counts**
* `b_ie_circulant.c`, `d27.txt`, `b_ie_S3.out`: `H(S_3)` by inclusion–exclusion.
* `b_tower_H.py` / `.out` (`k ≤ 4`) and `b_tower_k5_modp.py` / `.out` (`k = 5`).

**THM-4754**
* `c_thm4754_check.py` / `.out`

**HYP-9364 and the note**
* `d_ie_small.out`, `d_ie_large.out`: `R_n` and `P_q` recomputed.
* `d_fits.py` / `.out`, `d_tower_fit.py` / `.out`: extrapolations.
* `d_circ_band.py` / `.out`: bands and rank correlations.
* `d_regular_rho.py`, `d_regular_rho11.py` / `.out`: Conjecture 1 for `d = 7, 9, 11`.
* `d_regular_principle.py` / `.out`, `d_second_order.py` / `.out`

**Number board**
* `e_numberboard.py` / `.out`
* `hcount.c`, `e_interval_seq.out`
