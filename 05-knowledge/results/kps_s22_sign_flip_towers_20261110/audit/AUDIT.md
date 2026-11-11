# Independent audit of kind-pasteur-2026-11-10-S22 (THM-5130, HYP-9662, packet `kps_s22_sign_flip_towers_20261110`)

Audited state: commit `b5912dbd67` (checkpoint 2, "pre-audit"). Auditor: independent adversarial agent, 2026-11-11.
Nothing outside this `audit/` folder was created or changed by the audit. (During the audit the working tree acquired
an uncommitted edit of `05-knowledge/hypotheses/HYP-9658-...md` that is not mine; see minor finding m24.)

## Verdict: PASS WITH CORRECTIONS

- **Mathematics.** Every proof of THM-5130 (A)-(H) was checked line by line and is correct; two short gaps (H4
  equality case, the remark on 15) are filled below. No false theorem, no sign error, no wrong edge case was found.
- **Computation.** Every claim was recomputed with programmes written for this audit (plain permutation enumeration,
  an independent subset recursion, literal path generation). 0 mismatches: gadget law (5,876 cases incl. an
  exhaustive small universe), Theorem B, all 4,098 certificates, D4 (including two literal counts at order 24),
  the 9-grid and 43-grid data, the domination ray, the abundancy bounds, the neutral-gadget classification on all
  labelled tournaments of order <= 7, and the pre-registered runs E1, E2, E3 (re-derived from scratch).
- **Corrections needed.** Three major items, all about credit / the E4 pre-registration, none about a theorem:
  M1 (a "not found elsewhere" claim that is false), M2 (two concurrent notes that cover parts of Theorem A are not
  cited), M3 (two of the four E4 predictions were decidable from data already on disk; one outcome sentence is false
  as written). 26 minor items.

## 1. Proofs (task 1)

| item | verdict | note |
|---|---|---|
| Theorem A, bijection with walk systems of `L'` | CORRECT | blocks entered only from `B`, left only towards `v`; (o1), (o2) right |
| A, one walk: `r + 2M` (one visit), `M` (two), none with three | CORRECT | `M + h` and `M - W`; ending at `z`: `M - h + r`, `M` |
| A, two walks glued at `z`: `r + 2M`, `2M` | CORRECT | marked-visit bijection; both marked positions give `w3` after relabelling the free block |
| A, Case I / Case II: `2M + h + r - E`, `4M + h + r`, `2M` | CORRECT | recomputed term by term: Case I one visit `S + (M - S) + h`; Case II one visit `S + (h - E) + (M - S - W)`; Case I two visits `S + (M - S) + h + (M - S) + S`; Case II two visits `S + h + (M - S) + S + (M - S - W)`; three visits `M + M`; four impossible |
| A, the four rows | CORRECT | `u - sigma`, `nu - eps` count nonempty ordered pairs with the marked block in either slot |
| Remarks (1)-(3) | CORRECT | empty-base state `(1,0,1,0)` checked by hand and for all labelled `X` of order <= 4 |
| Theorem B (all nine constants, `w`, `w3`, `u3(j)`) | CORRECT | including `i = 1`, `m = 0`, `a = delta + 1`, `j = a` |
| Corollary B | CORRECT | "only if" is terse (m7) |
| Theorem C | CORRECT | (ii) needs `j >= 1` as stated |
| D1, D2, D3 | CORRECT | D1 = THM-4748 shared-endpoint formula with `s = L - 1`, `L - 2` (checked algebraically, also for general `delta`) |
| D4 | CORRECT | state recounted three ways; towers by hand (`M_1 = 646`, `M_2 = 52448`, `M_3 = 2^23` for every `n`) and literally at order 24 |
| F1, F2 (a)-(c), F3, F4 | CORRECT | F2 (c) uses "strong and `H = 1` implies one vertex" (standard) |
| H1 | CORRECT | |
| H2 (domination ray) | CORRECT | conditional on Session261 EXTENSION4 (1), which was spot-checked on all 5,405 rooted labelled tournaments of order <= 5 |
| H3 (iii) port lemma | CORRECT | port unselected, in-degree 1 in `B`: at most one initial visit and one visit after the in-neighbour, so every count is `c_1 + c_2 (2^a - 2)` |
| H4 merge inequality | CORRECT | the merge-sort rule keeps both orders; the pair is recovered by restriction, so the map is injective |
| H4 "`a(X) = 2^n` iff `X` transitive" | CORRECT WITH A GAP | asserted without proof; fix in m5 |
| "15 lies on no ray of offset +1 or -1" | CORRECT WITH A GAP | uses `H >= 5` for strong tournaments of order >= 4 without saying so, and the case `H(B) = 3` is only summarised; full argument in m6 |

Edge cases examined by hand and by exhaustion: `O` empty, `O = V(B)`, `P` empty, `P = V(X)`, `|B| = 1`, `|X| = 1`,
`B` empty (Remark 1), `i = 1`, `m = 0`, `a = delta + i`, `[i = 1][p = n]`, `[i_2 = 1][m_1 = 0]`.

## 2. Independent recomputation (task 2)

All programmes are in this folder; none imports the session's code. Outputs are the `.out` files.

| task | programme | result |
|---|---|---|
| counters agree | `aud_a_gadget_law.py` T0, `aud_c_selftest.py` | subset recursion = permutation enumeration on 1,230 + 435 pairs (orders 1..8); transitive controls at orders 12..20 |
| (a) gadget law | `aud_a_gadget_law.py` T1, T1e | 400 random `(B, O, X, P)`, orders <= 4 each (377 literal states by permutations, offsets -17..5), and all 5,476 cases with `B`, `X` labelled of order <= 3 and every `O`, `P`: 0 failures |
| (b) Theorem B | T2, T2b | 79 gadgets (`delta <= 5`, `i <= 4`, `m <= 3`), by ordered set partitions and by permutations; level map on 235 literal states, `delta = 1..4` |
| (c) certificates | `aud_c_certs.py` + `aud_pairstate.c` | all 2,049 + 2,049 lines rebuilt and recounted, 0 mismatches; 157 + 155 sampled lines (150 random each + the named `c`) recounted by literal path generation and closed (`T_1`, `T_2`): 0 mismatches |
| (d) D4 | `aud_defg.py` D0-D7, `aud_d4_tower24.py` | state `(-1, 1024, 217, 64)`; `2047, 4095, 8191`; out-lists = code; towers `n = 1, 7` at order 24: `H = 16,777,215` by literal generation of every path |
| (e) H3 | `aud_defg.py` E1-E8 | five 9-grid pairs as tabulated, `N(9, r, s)` for `r = 1, 2`; carrier: port in-degree 1 (in-neighbour 7), not beaten by the root, states `(1, 360, 135, 92)`, `(1, 704, 263, 178)`, `(1, 1392, 519, 350)`, `(1, 2768, 1031, 694)` for `a = 1..4`; four low sources; `H = 175`; bookkeeping: nothing left out |
| (f) H2 | `aud_defg.py` F1-F3, `aud_f_structure.py` S3 | `d(G_t) = rho_dom + 2^t delta_dom`, `t = 1, 2, 3`, with `rho_dom = d(G_0) - iota - delta_dom + 1`, on all 1,023 labelled pairs of order <= 4; `M >= delta_dom`, equality iff `B` transitive, on all labelled pairs of order <= 6 (2,064,384 at order 6) |
| (g) H4 | `aud_defg.py` G1, G2 | all 1,099 labelled tournaments of order <= 5: merge inequality for every `S`, merge map injective, `(n+1) H <= Pi_2 <= 2^n H`, equality on the right iff transitive |
| (G) neutral gadgets | `aud_neutral.c` | all labelled tournaments of order <= 7, every nonempty `P` (266,338,304 gadgets at order 7): `kappa = 0` only for `|P| = 2`, always `W = 2H`; with `H = 1` exactly distance 1 or 2; labelled counts agree with `neutral.out` class counts |
| (F), remark on 15 | `aud_f_structure.py` | F1-F4 on all labelled tournaments of order <= 5 (2,166 unit sets); seed slopes below 10 are 1, 2, 4, 6, 8; anti-seed slopes below 9 are 3, 5, 6 (order <= 6) |
| session scripts | `rerun_*.out` | `verify_s22.py` (18 PASS) and `sources_check.py` (8 PASS) pass and reproduce the committed `.out` files; so do `certs_verify.py` (with a build of their `pairstate.c`), `ray_cover.py`, `mersenne_towers.py`, `exploratory/ozj.py` |

43-grid, logical completeness: carrier `1 + 2^t (16 + 172 * 2^a)` gives `s = t + 4 >= 5`, `r = a - 2 >= 1` (module
sizes `a >= 3`); sources: `s = 2` (`t = r`), `s = 3` (`t = r + 2`), `s = 4` (`t = r + 3`) for all `r >= 1`; `s = 1`
for `r >= 2` (`t = r - 1`); head `(1, 1) = 175`. Every `(r, s)` with `r, s >= 1` is covered once the unit offset of
the carrier is known for all `a`; that is the port lemma, which is correct.

## 3. Numbers (task 3)

Cross-checked against `.out` / `.err` / certificate files, and where possible recomputed (`aud_numbers.py`,
`aud_reach.c`, `aud_reach_e2.py`, `aud_unit_states.c`): certificate counts (2,049 each; 1,969 / 63 / 17 and
1,959 / 88 / 2; `v2 = 0`: 1,225 and 1,187; max order 14; max `v2` 7 and 10; the 33 exceptional `c` with orders
10 (2), 11 (4), 12 (16), 13 (9), 14 (2)), E1 (805, 1,612, 80, 1,330, 15, uncovered list), E2 (339 start states, 89 /
1,028 / 1,872 and 70 / 1,124 / 1,782, 32 of 33, 3773 missing), E3 (614, 1,357, 1,592; 54, 1,128, 1,554; 5, 73, 85 at
depth 6 with `v2` = 32, 31, 30), 1,649, 3,722 of 4,095, 62,422, 596, 702, 77, 2,166, 1,540 / 980, 32,767, the 9-grid
table, the eight-numbers table. All match except the items listed as minor findings m1 and m8-m14.

## 4. Pre-registration (task 4)

| block | predictions committed | outcomes committed | order |
|---|---|---|---|
| E1 / E2 (R1-R5) | `86a1fc6e5a`, 2026-11-10 16:20:30 | `33114ad576`, 2026-11-11 07:37:17 | predictions first |
| E3 (R6a-c) | `24076b49e2`, 2026-11-10 16:42:39 | `33114ad576` | predictions first |
| E4 (G1-G4) | `33114ad576`, 07:37:17 | `b5912dbd67`, 07:54:02 | predictions first, but see M3 |

Both early commits sit in the linear history before Session250 (16:51) and every other concurrent packet, which is
consistent with "pushed" (push times themselves cannot be read from the repository). Nothing above an outcome line
was altered by a later commit. Outcomes R1-R6 and G1-G4 match the data files, and R1-R6 were re-derived by my own
programmes. Caveat for E4: M3.

## 5. Credit / novelty (task 5)

Answers to the five questions.

(i) No other file states the four-row transport law. `ONE_EXIT.md` Theorem 1 (`H = hH + FU - sE - dW`) is row `h`
only, but for an arbitrary one-exit attachment and with exactly the base coordinates `(H, U, E, W) = (h, M, E, h - r)`;
together with the closing law it implies that some matrix exists (M2). `MARKED_RESELECTION.md` / `MODULE_REBASE.md`
give the module half as a principle (dependence on the marked cover tensor only), with no coefficients.

(ii) Not contradicted. Lines `c 2^s - 1` already covered for all large `s` by other families: `c = 1, 3, 5`
(THM-4748; THM-5085), `c = 9` (Session261 MODULE_REBASE section 2, `H_A(m) = 9 * 2^m - 1`), `c = 15` (multiplicand
carrier 21, `H_m = 15 * 2^m - 1`), `c = 3` again (Session260, `24 * 2^k - 1`). See m19.

(iii) Yes, absent. Session261 NINE_OFFSET_EXTENSION has levels 1-6 and states levels `s >= 7` OPEN; no 43-grid,
slope 172 family or level-7 source `(129, 288)` occurs in the multiplicand, seed_recursion, Session260-262 packets.

(iv) The merge inequality and `Pi_2 <= 2^n H` were not found in canon. Session197 has both lower bounds
(`B_T(2) >= (n+1) H(T)` and `B_T(2) >= 2^n`); multiplicativity of `Pi_2` over ordinal sums is Session216
three-core/PROOF.md (1) (`q' = q1 q2`), see m20.

(v) Overstated sentences: M1, M2, m15, m16, m17, m18.

## 6. Findings

### MAJOR

**M1. A "not found elsewhere" claim is false (F2 (a), (b)).**
`01-canon/theorems/THM-5130-...md:29-31`: "F1 in the exact form, F2 (a), (b) and F4 were not found elsewhere."
- F2 (a) is `session258_endpoint_seed_20261111/seed_decision/CORE_CLASSIFICATION.md` equation (1),
  `K_ba(B) = E * U_a(A) * F_b(C) * prod Q_2(S_l)`, with the same bounds (`>= H(B)`, and (3)); commit `63b86513f8`,
  2026-11-10 22:40. The same file is cited by the session for F2 (c).
- F2 (b), offset part, is in the proof of that theorem (`W_B(O) = H(exterior) W_S(O)`) and in
  `session250_sixteen_20261110/constructive/PROOF.md` Theorem 5 (`c_B = c_C product H(D)`), 2026-11-10 16:53.
- Both, including the slope formula of F2 (b), are instances of the ordinal-sum law (1) of
  `session216_overnight_constructions_20261107/three-core/PROOF.md` (`A' = A1 q2 + h1 A2`,
  `T' = A1 G2 + T1 h2 + h1 T2 + E1 S02`).
- `M >= |O| H(B)` is Session250 Lemma 4 (16:53), before Session258 (THM line 30, SOURCES line 35).
*Correction:* "KNOWN: F2 (a) (Session258 CORE_CLASSIFICATION (1)), F2 (b) (Session250 Theorem 5, Session258; both
are the join law of Session216 three-core (1)), F2 (c) (Session250, Session258), `M >= |O| H(B)` (Session250
Lemma 4, Session258 (5)), F3 (Session260). Not found elsewhere: F1 as an identity (the defect of the cut injection)
and F4."

**M2. Credit for Theorem A is incomplete: two concurrent notes and one special case are not cited.**
`THM-5130:9-10` ("New here ... the whole state is transported linearly through a new closer"), `README.md:34-35`
("Not found elsewhere: the transport law for the whole state"), `SOURCES_AND_DOMINATION.md:47-50` ("Specific to this
session ... and the items of this note").
- (a) `spectrum_multiplicity_signed_tables_20261111/returns/ONE_EXIT.md` Theorem 1 (commit `b3b3591701`, 2026-11-10 23:33):
  for any tournament `X` with one exit vertex `x` attached to `(G, A)`, `H = hH + FU - sE - dW`. The base enters
  through the same four coordinates. With `X = X' => v`, `x = v` one gets `(h, F, s, d) = (H(X'), Pi_2(X'), 0, H(X'))`,
  i.e. row `h`. Since `T_t(B+, P)` is itself a one-exit attachment onto `(B, O)`, Theorem 1 with the closing law
  (and its (10) for `t = 0`) already implies that `r+`, `M+`, `E+` are linear in `(h, M, E, W)`. The entries are not
  there. Checked: `aud_credit_checks.py` C2.
- (b) `session261_creative_chanpak_20261111/constructors/MARKED_RESELECTION.md`, marked-context theorem (commit
  `40eba89fd5`, 2026-11-11 02:11), with `MODULE_REBASE.md` (`79ec1b74a5`, 01:49): the host's `(H, W, U)` depend on a marked module only through the
  tensor `(C_r, S_r, E_r, K_r)`. The nine gadget constants are entries of it: `H_X = C_1`, `Pi = C_2 + 2 C_1`,
  `w = K_2`, `w3 = K_3`, `sigma = S_1`, `u = S_1 + S_2`, `u3 = S_3`, `eps = E_1`, `nu = E_1 + E_2` (C3).
- (c) H1 at the point base is `NINE_OFFSET_EXTENSION.md` section 3 (1), `U_x = 2^(p+1) 3^q + 2^(q+1) - 0^p` (C4).
What survives: the explicit matrix (rows `r`, `M`, `E`) and the closure of the four coordinates under "close, then
re-select inside `X`" are stated nowhere else (Session216 star/PROOF.md section 6 says the closed output is not
re-inserted as a source), and the statement was pre-registered at 16:20 on 2026-11-10, before (a) and (b).
*Correction:* add to the status of (A) "CONCURRENT: ONE_EXIT Theorem 1 (row `h` for every one-exit attachment, same
four base coordinates; with the closing law it gives the existence of the matrix); Session261 MARKED_RESELECTION
(the module enters through its marked cover tensor). New: the explicit rows `r`, `M`, `E`." Add (c) to the status of
(H). Replace "Not found elsewhere: the transport law for the whole state" by "the explicit transport matrix", and
"the items of this note" by a list (Lemma 1, Proposition 4, levels 7-11, the 43-grid, the upper bound of the
abundancy index).

**M3. E4 pre-registration: two of the four predictions were not blind, and one outcome sentence is false as written.**
- `PREREGISTRATION.md:123-124`: "Nothing below has been computed: no state `E` of the pairs with
  `(h, M) = (15, 2^k)` ... no pair with `H(B) = 31` and a power of two as slope." The files
  `runs/e2_states_order_le7.txt` (the start states of E2, run on 2026-11-10) and `runs/e2_states_order8.txt`, both
  committed in the same commit as this block, contain the state `1 16 64 6`, i.e. `E = 6`, and the unit states of
  order <= 8 (base counts of seeds with a power of two as slope: 1, 5, 9, 15, 19, 27, 51, 351; no 31). The order-9 census of unit
  pairs (`pt_9.txt`, used for the certificates of the same commit) also existed. So G1 and G3 were decidable from
  data on disk. G1 failed, so the data were evidently not consulted, but G3 "CONFIRMED" is a look-up, not a blind
  prediction. G2 and G4 were blind.
- `SOURCES_AND_DOMINATION.md:222-223`: "G3 CONFIRMED (no pair of order at most 9 with `H(B) = 31` and slope 64 or
  128)". False without `r = 1`: the pair (order 6, code 20612, mask 34), quoted in `runs/rays_through.out` and in
  verify V10, has state `(15, 128, 31, 8)`.
*Correction:* add after the E4 outcome block: "Audit note: the quantities of G1 and G3 were already present in
`runs/e2_states_order_le7.txt`, `runs/e2_states_order8.txt` and the order-9 census; they count as retrospective
checks. G2 and G4 were blind." In SOURCES section 8 write "no pair of order at most 9 with `r = 1`, `H(B) = 31`
and slope 64 or 128".

### MINOR

Proof gaps and statements.

- **m1.** `THM-5130:249-250` and `README.md:134-137`: "1,969 seeds and 1,959 anti-seeds are pairs of order at most
  9; the others are towers". 47 of the 80 seed towers and 54 of the 90 anti-seed towers also have order <= 9
  (`make_certs.py` keeps a tower on ties). Write: "2,016 seed and 2,013 anti-seed certificates have order <= 9
  (1,969 and 1,959 quoted from the census, 47 and 54 written as towers); the other 33 and 36 are towers of order
  10 to 14."
- **m2.** `THM-5130:316`, status line 32, `SOURCES:16-17`: "`M = delta_dom` exactly when `B` is transitive" is
  FINITE-EXACT only in the direction "only if". "If" is provable: for `O = {o_1 < ... < o_k}` in `TT_n` the subsets
  with least element `o_j` contribute `2^(k-j) 2^(c_j)` and `k - j + c_j = n - o_j`. Row 17 of the SOURCES table
  ("zero-defect interfaces are the one-apex ones, PROVED (S7)") mixes the proved state formula with the
  finite-exact characterisation.
- **m3.** `THM-5130:316-323` and status line 32: H3 is labelled FINITE-EXACT, but (ii) and (iii) are statements for
  all `r` (and all `r, s`): "PROVED from FINITE-EXACT states (closing law; port lemma)", as in SOURCES line 21.
- **m4.** The certificates of the four low 43-grid sources and of the head appear only inside `sources_check.py`.
  Add to SOURCES section 4: level 1 (7, 532488, {1,3}) state `(3, 172, 51, 26)`; level 2 (7, 313346, {1,3})
  `(5, 172, 49, 22)`; level 3 (6, 16392, {0,5}) `(9, 86, 17, 2)`; level 4 (6, 17034, {0,3}) `(17, 86, 33, 14)`;
  head: order 7, code 17976, `H = 175`.
- **m5.** `THM-5130:326-327`: add the proof of the equality case: if every term equals `H(X)` the merge map is a
  bijection for every `S`; for a Hamiltonian path `v_1 ... v_n` and `S = {v_i, v_j}`, `i < j`, the restriction is a
  path, so `v_i -> v_j`; hence `X` is transitive.
- **m6.** `THM-5130:348-352`: add (a) "a strong tournament other than `C_3` and a point has `H >= 5`"; (b) for
  `H(B) = 3`, `M <= 8` forces `|O| = 2`; both in `C_3` gives `r = 1`; both outside gives `W` in {6, 12, ...};
  one in `C_3` and one outside at distance `k` in the component order gives `W = 3 * 2^k + [k = 0]`, so `k = 0`,
  slopes `9 * 2^j` (after) and `15 * 2^j` (before).
- **m7.** `THM-5130:176-177`: for "only if" say that sink padding gives two bases with the same `r` and different
  `M`, so `kappa_X = 0` is forced (and then `r_X = -1`).

Numbers and reproducibility.

- **m8.** `PREREGISTRATION.md:84`: "`2^N - 1` paths on `N >= 24` vertices" for all 15 states. The three states of
  slope `2^25` give `2^(N+1) - 1` on `N = 25, 26` vertices (`mersenne_towers.out`, last three lines).
- **m9.** `PREREGISTRATION.md:15`: the displayed coefficient `(1 - 2^(delta-1))` of `r` is 0 at `delta = 1`; the
  true value is -1 (THM-5130 has it right). The block cannot be edited; add an erratum line under the outcomes.
- **m10.** `SOURCES:16`, `PREREGISTRATION.md:149`: no script or output in the packet produces "3,978 pairs of order
  <= 6"; `sources_check.py` S5 covers labelled order <= 5, and tests the ray form only for orders <= 3 and only
  through `d(G_2) - d(G_1)`. The audit closes both (section 2 (f)); cite it or add the script.
- **m11.** `scripts/src10.c` skips `r < 5`, but `runs/src43grid.txt` contains 15 lines with offset 3: the file was
  not produced by the committed programme.
- **m12.** `THM-5130:339-341` (same control in `README.md:221-222`): "pairstate.c ... is correct to order 25". The programme declares
  `unsigned out[24], inn[24]`, so order 25 is out of bounds (undefined behaviour), and its header says "orders up
  to 20". Say "valid to order 24". No effect: certificates have order <= 14.
- **m13.** `README.md:227`: "about two minutes"; `verify_s22.py` runs in about 3 seconds here.
- **m14.** `SOURCES:110-112`: at `s = 5` the least-order (order 6) one-vertex sources have both the first type
  (slope 36, excess 3) and the second (slope 72, excess 39; `runs/src9grid.txt`, line `6 568 16 33 33 72`).

Wording, scope, status.

- **m15.** `README.md:37-39`: "the only universal operation of this kind (up to order 7)". The scan supports: "the
  only tail gadgets of order <= 7 that keep `|r| = 1` for every base are the transitive pairs at distance 1 or 2"
  (sink and source padding, and compositions, are also universal).
- **m16.** `README.md:82-84`: "there is no gadget of order at most 7 that keeps the sign" needs "for every base"
  (the one-vertex gadget keeps the sign of a positive offset).
- **m17.** `SOURCES:126`: "No uniform family:" should read "No uniform family was found:" (scope: ports of unit
  roots of order <= 9).
- **m18.** `SOURCES:145`: "The grids 3, 21, 35 ... are finished exactly this way." Grids 3 and 21 were completed
  through stable ports and head grids (seed_recursion R3), grid 35 through offset-transport carriers and a product
  replacement (multiplicand); "in the same pattern (one carrier for large `s`, finitely many sources and heads)".
- **m19.** `THM-5130:27-28, 253`, `README.md:143`: "the anti-seeds are new" should carry "as a certified list; the
  lines of `c = 1, 3, 5, 9, 15` were already covered for large `s`" (answer (ii)).
- **m20.** `THM-5130:32-33, 327`: add "multiplicativity of `Pi_2`: KNOWN, Session216 three-core (1)"; Session197
  also has `Pi_2 >= 2^n`.
- **m21.** `README.md:147-148`: "the first tournament of a ray has at least `4 H(B)` paths" holds for `|O| >= 2`
  only (for `|O| = 1` the bound is `3 H(B)`, attained when `o` is a sink).
- **m22.** `HYP-9662:54-55`: "11, 23, 91, 95, 125, 127 do not". In the E3 box the odd `c <= 129` not reached are
  11, 23, 29, 91, 95, 125, 127; either add 29 or restrict to the eight numbers of R6a (11, 23, 91, 95, 127).
  `HYP-9662:35`: "of order `X^(2/3)`" is `X^(2/3) / (13 ln X + 102)`. `HYP-9662:40-41`: "a tower over a transitive
  pair" is used for "a tower over a transitive base" (in THM-5130 a transitive pair has two selected vertices).
- **m23.** `SOURCES:90-91`: "`w` = 10 does not occur for at most four selected vertices": it never occurs (four or
  more selected vertices give `w >= 14`). `SOURCES:196-197`: "integers only at 2, 3, 5, 7 and 19" has no range;
  true for all primes up to 30,000 (audit).
- **m24.** Working tree, `05-knowledge/hypotheses/HYP-9658-...md` (uncommitted, not by the audit): the new
  2026-11-11 paragraph is inserted in the middle of the 2026-11-10 continuation, so the sentences "All 2,049 lines
  ... are now complete ..." now read as part of the S22 paragraph. Move the new paragraph after "... bounded
  universes."
- **m25.** `PREREGISTRATION.md:53, 90, 130`: the label G3 is used for the decision rule and for an E4 prediction.
  The message of commit `b5912dbd67` says "outcomes of E3 and E4"; the E3 outcome was already in `33114ad576`.
- **m26.** `THM-5130:37`: "Independent audit: PENDING" to be replaced by the result of this audit.

### Observations (no action required)

- Session260's five-vertex carrier `(24,1,3,6,14)` is the one-closer tournament `C(TT_2, {1,2}; TT_2)`, and its
  marked families for `v = 1, 2` are exactly the towers of D3 over `(TT_2, {1,2})` with gap `c = delta - 1`
  (`aud_credit_checks.py` C1: 108 marked modules, Session260 (4)-(5) = THM-5130 (1)). The session's label
  CONCURRENT for the bracket and for D3 is therefore accurate, and the two results agree on the same tournaments.
- The pre-registration of 16:20 already lists the families `15 * 3^m` and `30 * 3^m - 1` (as hand computations),
  before endpoint_coverage_compiler (`70429f4c72`, 2026-11-10 22:42) and Session260 coverage/PROOF.md
  (`4315a74fe0`, 2026-11-11 00:55). First-commit times used in this report: CORE_CLASSIFICATION `63b86513f8` 2026-11-10 22:40;
  Session250 constructive/PROOF.md `42969527be` 16:53; ONE_EXIT `b3b3591701` 23:33; MODULE_REBASE `79ec1b74a5`
  2026-11-11 01:49; MARKED_RESELECTION `40eba89fd5` 02:11; NINE_OFFSET_EXTENSION `f01a711377` 02:36.
- `rays_through.out` (15, 17, 110, 156 rays) was not recomputed (it needs class lists); 62,422 is the right number
  of class-set pairs of order <= 7.

## 7. The audit's programmes

Build: `gcc -O2 -o aud_pairstate aud_pairstate.c`, likewise `aud_neutral.c`, `aud_reach.c`, `aud_unit_states.c`
(the Python drivers expect `aud_pairstate.exe` in this folder). Python 3, no third-party package. Binaries were
deleted after the runs. Runs: `aud_neutral N` for `N = 1..7`; `aud_reach 10 4 10 4 4097` (E1),
`aud_reach 5 2 10 6 4097` (E3), `aud_reach 10 4 10 3 4097 pairs FILE` (transitive-pair starts);
`aud_unit_states N | sort -u` for `N = 1..7`.

| file | what it does |
|---|---|
| `aud_core.py` | permutation enumeration of `H`, states, `Q`, gadget constants; a pull-form subset recursion; `close`, `transitive`, codes |
| `aud_pairstate.c` | state by subset recursion; state by literal depth-first generation of every path (no table); `H` only; literal `H` for orders up to 30 |
| `aud_a_gadget_law.py` | tasks 2 (a), (b), remarks (1)-(3) |
| `aud_c_selftest.py`, `aud_c_certs.py` | the C counter against permutations; all certificates; file statistics |
| `aud_defg.py`, `aud_d4_tower24.py` | tasks 2 (d)-(g) |
| `aud_f_structure.py` | (F), slopes of unit pairs of order <= 6, H3 (i) at order 6, quoted states |
| `aud_neutral.c` | (G) on all labelled tournaments of order <= 7 |
| `aud_numbers.py` | quoted run statistics from the packet's data files |
| `aud_reach.c`, `aud_unit_states.c`, `aud_reach_e2.py` | E1 and E3 from scratch (1.77e9 and 2.69e10 states); the exploratory run with transitive-pair starts (37 / 620, union 1,649: `aud_reach_e1x.out`); the 339 unit states of order <= 7 from all labelled tournaments; E2 |
| `aud_credit_checks.py` | computations behind M2 and the observations |
| `rerun_*.out`, `aud_g1_gadget_sources.out` | reruns of the session's scripts; G1 from the unit states of order <= 8 |
