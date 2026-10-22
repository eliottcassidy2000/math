# Audit pass 2: Hamiltonian cycles (HYP-474), the packing split (THM-4751 §7, note §13, MISTAKE-660), THM-137/161, rotational local maxima (2026-10-22)

**Auditor.** The independent refutation-first subagent of pass 1 (`lattice_bridges_audit_20261021.md`), on the coordinator's second-pass request. No existing file was edited; no git operation was run.

**Audited versions** (sha1 at the start of pass 2):

| file | sha1 |
|---|---|
| `01-canon/theorems/THM-4751-circulant-h-landscape-paley-peak-at-class-number-distance.md` | `2cc088a2` (commit `678577d631`, §6–7 added) |
| `05-knowledge/results/lattices_codes_tournaments_20261021.md` | `28a61be5` |
| MISTAKE-660 | in `01-canon/MISTAKES.md`, line 19637 |

**Scripts and outputs.** All in `04-computation/experiments/lattice_bridges_20261021/audit/`; nothing reuses `circ_counts.c` or `alpha_decomp.c`. Build with `cc -O3 -o a12_cycles a12_cycles.c` and `cc -O3 -o a3_hcount a3_hcount.c`; the binaries are removed after use.

| script | output | content |
|---|---|---|
| `a12_cycles.c` (mode `hk`) | `a12_RP19_hk.out`, `a12_RP23_hk.out`, `a12_all19_hk.out`, `a12_all23_hk.out` | Held–Karp from vertex 0: `H`, and every `c_k` via `c_k = (n/k)·#(k-cycles through 0)`, for all classes at `n = 19, 23` |
| `a12_cycles.c` (mode `iec`) | `a12_RP19_iec.out`, `a12_RP23_iec.out` | Hamiltonian cycles by inclusion–exclusion over closed walks from 0, a second method |
| `a12b_analysis.py` | `a12b_analysis.out` | HYP-474 per odd `k`, against all classes; trace checks for `c_3`, `c_4`, `c_5`; the session's `C=` column |
| `a12c_order.py` | `a12c_order.out` | `H`-order against `c_n`-order (THM-161, fact 5) |
| `a13b_alpha_explicit.py` (method E) | `a13b_alpha_explicit_{7,9,11,13}.out` | explicit DFS of every directed odd cycle; explicit enumeration of disjoint collections |
| `a13_alpha_ie.py` (method Z) | `a13_alpha_ie_{7,…,21}.out`, `.json` | `f(X)` from the Held–Karp dump plus rotation; ranked zeta transform; inclusion–exclusion; three primes and CRT |
| `a13c_alpha_claims.py` | `a13c_alpha_claims.out` | the split table, champions, THM-137 check, comparison with `alpha_decomp_*.out` |
| `a14_rotflip_sets.txt`, `a14_rot31flip_sets.txt` | `a14_rotflip.out`, `a14_rot31flip.out` | best and first flips of `R_27`, `R_29`, `R_31` (`a3_hcount ie`) |

Severity as in pass 1.

## Summary verdicts

| # | item | verdict | severity |
|---|---|---|---|
| P1 | HYP-474 fails at `k = p`: `c_19(R_19) = 34,841,356,485 > c_19(P_19) = 34,358,763,933` | CONFIRMED (two methods) | — |
| P2 | `n = 23`: `C(R_23) = 382,137,819,329,585`, the unique maximum over all 94 classes; `C(P_23) = 374,127,973,957,716` | CONFIRMED (two methods for `R`, `P`; Held–Karp for all classes) | — |
| P3 | HYP-474 for odd `k < p` | CONFIRMED: `P_q` is the unique maximiser of `c_k` for every odd `5 ≤ k ≤ p − 2` over all circulant classes, at `p = 19` and at `p = 23`; `c_3` is constant | — |
| P4 | session `C=` column, `circ_19.out` and `circ_23.out` | CONFIRMED (30/30, 94/94 classes) | — |
| P5 | repo paraphrase of Savchenko (TANGENTS T134, INV-053: "c_m(DR_n) > c_m(RLT_n) for ALL m ≡ 1,2,3 mod 4") | REFUTED at `m = n = 19, 23`; the cited arXiv:2403.07629 treats only `m = 8` | medium |
| P6 | `α_k` at `n = 9, 11, 13` (all classes) | CONFIRMED (method E = method Z = `alpha_decomp`; OCF reproduces `H`) | — |
| P7 | `α_k` at `n = 15, 17, 19` (all classes) and `n = 21` (`R`, both energy minima) | CONFIRMED (method Z; OCF reproduces `H` for every class) | — |
| P8 | split table (`n = 9, 11, 13, 17, 19, 21`) | CONFIRMED, every entry exact | — |
| P9 | "at `n = 13` `R` wins only through `α_4`" | CONFIRMED (without `α_4`, `ΔH = −9,516`) | — |
| P10 | champions: `n = 17, 19`; `n = 13, 15`; `α_5` at `n = 15` | CONFIRMED | — |
| P11 | THM-137 refinement: `α_2(R) > α_2(flat)` only for `n ≤ 11` | CONFIRMED (`n = 7, 9, 11` greater; `13, 15, 17, 19, 21` smaller) | — |
| P12 | what THM-137 and THM-161 say | CONFIRMED as quoted in MISTAKE-660. THM-161's fact 5 ("H-class ordering follows c_11") is nearly but not exactly true at `n = 19, 23` | low |
| P13 | THM-4751 §7 / note §13 wording | CORRECTED (`n = 21` "energy minimum" is two classes; "gap grows with k" holds for ratios only; the `{3,7,11}` "mechanism" is a threshold restatement) | low |
| P14 | `R_27`, `R_29` best flips (`{13,14}`, `{14,15}`) and first flips | CONFIRMED (4 of 4 values reproduced exactly) | — |
| P15 | `R_31` best flip `{15,16}`: `H = 20,284,103,532,402,676,770,852,133`, ratio 0.996670 | CONFIRMED | — |

## 1. Hamiltonian cycles and HYP-474

**P1, P2 — the counts.** Two independent methods agree:
* Held–Karp from vertex 0, `C = Σ_(v→0) dp[full][v]`;
* inclusion–exclusion `C = Σ_(X∋0) (−1)^(n−|X|) (A_X^n)_(00)`, modulo `2^128`.

| `n` | class | `C`, Held–Karp | `C`, inclusion–exclusion |
|---|---|---|---|
| 19 | `R_19` | 34,841,356,485 | 34,841,356,485 |
| 19 | `P_19` | 34,358,763,933 | 34,358,763,933 |
| 23 | `R_23` | 382,137,819,329,585 | 382,137,819,329,585 |
| 23 | `P_23` | 374,127,973,957,716 | 374,127,973,957,716 |

Over all 94 classes at `n = 23`, `R_23` is the unique maximiser of `c_23`. The Held–Karp `H` agrees with the pass-1 tables on every class, and the trace identities `tr A^k = k c_k` (`k = 3, 4, 5`, THM-118) hold on every class (`a12b_analysis.out`).

**P3 — only `k = p` fails.** At `p = 19` and at `p = 23`, `P_q` is the unique maximiser of `c_k` over all circulant classes for every odd `5 ≤ k ≤ p − 2`; `c_3 = n(n²−1)/24` is constant. Margins at `p = 19`:

| `k` | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 |
|---|---|---|---|---|---|---|---|---|
| `P_19` | 11628 | 424080 | 12156390 | 249208902 | 3280900392 | 23662379790 | 69401425077 | 34358763933 |
| `R_19` | 10488 | 391362 | 10807884 | 224515210 | 2961329208 | 21901889133 | 66503305202 | **34841356485** |

At `p = 23` the same holds through `k = 21` (`753,133,060,365,334 > 734,346,446,683,604`), and `k = 23` reverses. MISTAKE-660 states the `p = 19` part; the `p = 23` part ("odd `5 ≤ k ≤ 21` still favour Paley") is worth adding.

**P4.** My `c_n` equals the session's `C=` column on 30/30 classes at `n = 19` and 94/94 at `n = 23`.

**P5 — a second casualty (medium).**
* TANGENTS T134 and INV-053 (kind-pasteur S21) record "Savchenko (2024) proves `c_m(DR_n) > c_m(RLT_n)` for ALL `m = 1,2,3 mod 4` (including all odd m)".
* `P_19` and `P_23` are doubly regular, and `R_n` is `RLT_n`, the regular locally transitive tournament.
* So `c_19(R_19) > c_19(P_19)` and `c_23(R_23) > c_23(P_23)` refute that paraphrase at `m = n`.
* The cited paper, arXiv:2403.07629, is about `m = 8` only: `c_8(DR_n) > c_8(RLT_n)` for `11 ≤ n ≤ 35`, reversed for `n ≥ 39` (abstract).
* MISTAKE-660 says that THM-161 "overstates Savchenko, whose results compare doubly regular tournaments with `RLT_n`". That is right, but it leaves T134/INV-053 standing. A March-2026 snapshot of §13 said "R_19 is not [a DRT], so there is no conflict". That reasoning was wrong, because `RLT_19 = R_19` is exactly Savchenko's comparison object.

## 2. The packing split

**Methods.** Both are independent of `alpha_decomp.c`, which uses a min-vertex path DP over all `2^n` masks and the min-element recursion `F[W]`.
* **Method E** (`n ≤ 13`):
  * every directed odd cycle is enumerated as a vertex sequence by DFS, then grouped by vertex set;
  * all collections of pairwise disjoint vertex sets are enumerated explicitly, weighted by products of cycle counts.
* **Method Z** (`n ≤ 21`):
  * `g(Z)` for `Z ∋ 0` is obtained by Held–Karp from vertex 0;
  * `f(X) = g(X − min X)` by rotation;
  * ranked zeta transform `F[r](W) = Σ_(X⊆W, |X|=r) f(X)`, and `P_k = F^k` as rank polynomials;
  * the ordered disjoint `k`-tuples are `Σ_Y Σ_j (−1)^j C(n−|Y|, j) P_k[|Y|+j](Y)`, i.e. inclusion–exclusion in place of a Möbius transform, divided by `k!`;
  * arithmetic modulo three primes, combined by CRT, with the third prime as a check.

**P6, P7 — values.**
* Methods E and Z agree on every class at `n = 7, 9, 11, 13`.
* Method Z agrees with `alpha_decomp_*.out` on every class printed there (`n = 9, 11, 13, 17, 19, 21`).
* For every class computed (`n ≤ 19` all classes; `n = 21`: `R_21` and both energy minima), `Σ 2^k α_k` equals the pass-1 `H`, and `α_1 = Σ_X f(X)`.

Selected rows (method Z):

| `n` | class | `α_1, …, α_⌊n/3⌋` |
|---|---|---|
| 13 | `R_13` | 606027, 436748, 87568, 3224 |
| 13 | `Z_13{1,2,3,5,6,9}` | 645937, 439452, 77428, 2119 |
| 19 | `R_19` | 126443605257, 122111579294, 42960731622, 5521030944, 331078344, 4100656 |
| 19 | `P_19` | 130965270477, 123659531220, 41184418943, 4903920444, 251464164, 2221081 |
| 21 | `R_21` | 12030499746751, 12330182836208, 4796354751404, 738531326288, 58868297768, 1454221328, 12571712 |
| 21 | `Z_21{1,2,3,7,10,12,13,15,16,17}` | 12302801976553, 12402025082712, 4631261764086, 675784771935, 48461006643, 1012015200, 5954940 |

**P8 — split table.** Every entry is reproduced exactly (`a13c_alpha_claims.out`):

| `n` | `2Δα_1 + 4Δα_2` | `Σ_(k≥3) 2^kΔα_k` | `ΔH` |
|---|---|---|---|
| 9 | −162 | +72 | −90 |
| 11 | −4,620 | +2,552 | −2,068 |
| 13 | −90,636 | +98,800 | +8,164 |
| 17 | −189,570,298 | +318,784,408 | +129,214,110 |
| 19 | −15,235,138,144 | +26,752,215,992 | +11,517,077,848 |
| 21 | −831,973,445,620 | +2,686,870,223,200 | +1,854,896,777,580 |

At `n = 21` the "flat" class is `Z_21{1,2,3,7,10,12,13,15,16,17}`. The other energy minimum, `Z_21{1,2,3,4,6,9,10,13,14,16}`, gives a different `ΔH`: `125547534942879 − 123713161714233 = 1,834,373,228,646`. THM-4751 §7's header "`H(energy minimum)`" therefore needs the class named at `n = 21` (P13).

**P9.** At `n = 13`: `8Δα_3 = 81,120` and `16Δα_4 = 17,680`. Without `α_4`, `ΔH = −90,636 + 81,120 = −9,516 < 0`.

**P10 — champions** (method Z, all classes):

| `n` | `α_1` | `α_2` | `α_3` | `α_4` | `α_5` | `α_6` |
|---|---|---|---|---|---|---|
| 9 | Emin | R | tie of 2 classes | — | — | — |
| 11 | Emin (`P_11`) | `Z_11{1,2,3,4,6}` | R | — | — | — |
| 13 | Emin | `Z_13{1,2,3,4,6,8}` | R | R | — | — |
| 15 | `Z_15{1,2,3,4,6,7,10}` (an Emin) | the same | R | R | `Z_15{1,3,4,6,7,10,13}` | — |
| 17 | Emin | Emin | R | R | R | — |
| 19 | Emin (`P_19`) | Emin | R | R | R | R |

All are unique except the `α_3` tie at `n = 9`. Every stated champion claim holds, including "`α_5` maximised elsewhere" at `n = 15`.

## 3. THM-137 and THM-161

**P11.**

| `n` | flat class | `α_2(R)` | `α_2(flat)` | `α_2(R)` vs `α_2(flat)` |
|---|---|---|---|---|
| 7 | `P_7` | 14 | 7 | greater |
| 9 | `Z_9{1,2,3,5}` | 354 | 336 | greater |
| 11 | `P_11` | 11,110 | 10,879 | greater |
| 13 | `Z_13{1,2,3,5,6,9}` | 436,748 | 439,452 | smaller |
| 15 | each of the 4 energy minima | 22,518,662 | 22,762,572 … 22,812,707 | smaller |
| 17 | `Z_17{1,2,3,5,7,8,11,13}` | 1,482,234,998 | 1,499,127,864 | smaller |
| 19 | `P_19` | 122,111,579,294 | 123,659,531,220 | smaller |
| 21 | `Z_21{1,2,3,7,10,12,13,15,16,17}` | 12,330,182,836,208 | 12,402,025,082,712 | smaller |

At `n = 21` the other energy minimum has `α_2 = 12,404,080,236,246`, also larger. For `k ≥ 3`, `Δα_k > 0` in every row from `n = 9` on. So "the interval has higher `α_2`" is true exactly for `n ≤ 11` in these tables, and the interval's advantage from `n = 13` lies in `α_k`, `k ≥ 3`, as claimed.

**P12 — the texts, quoted.**
* `THM-137-paley-crossover-mechanism.md`, Main Result 2: "Interval structural advantage: The interval's cycles are more INDEPENDENT (higher alpha_2), contributing more to the higher-order OCF terms."
* Its Layer 3 table has only `p = 7` (`α_2(P) = 7`, `α_2(I) = 14`): "Paley has MORE total cycles (alpha_1) but FEWER independent pairs (alpha_2)."
* THM-161, "Why Paley Maximizes H": "Paley maximizes EVERY c_k simultaneously (Savchenko's theorem for DRTs)". Its verification is `N_maximization_paley.out`, `p = 7, 11` only.
* THM-162 has the parallel line "Savchenko: c_k(Paley) = max c_k over all DRTs (2016–2024)".
* MISTAKE-660 quotes all three correctly.
* THM-161's "Key structural fact 5" ("H-class ordering follows c_11 (Hamiltonian cycles)", `p = 11`) is nearly true at larger `n` (`a12c_order.out`). Spearman`(H, c_n)` is 0.9991 at `n = 19` (2 discordant pairs of 435) and 0.9997 at `n = 23` (15 of 4371); the `H`-maximum and the `c_n`-maximum coincide (`R_n`). This is why HYP-474 fails at exactly the `n` where Paley stops being `H`-maximal.

**P13 — wording (low).**
* THM-4751 §7 / note §13 say "`R_n` has more packings of every size `k ≥ 3`, and the gap grows with `k`". This holds for the ratio `α_k(R)/α_k(flat)` (at `n = 19`: 1.043, 1.126, 1.317, 1.846 for `k = 3..6`). It fails for the `H`-contributions `2^kΔα_k`, which shrink after `k = 3` (`1.42·10^10`, `9.87·10^9`, …).
* "the switch … comes exactly when they fit" and "a tournament mechanism: … `deg I(Ω(P_q), x) = ⌊q/3⌋ ≤ 3`" are observations on `n ≤ 25`. The second is the threshold `q ≤ 11` restated. It should stay typed NUMEROLOGY-with-a-mechanism-sketch, not "mechanism" (pass-1 #21).

## 4. Rotational local maxima

**P14.** `a3_hcount ie` reproduces exactly (`a14_rotflip.out`):

| `n` | flip | `H` | ratio to `H(R_n)` |
|---|---|---|---|
| 27 | `{13,14}` | 425,337,994,159,252,763,931 | 0.996130 (the best) |
| 27 | `{1,26}` | 422,880,000,888,645,377,865 | 0.990374 |
| 29 | `{14,15}` | 86,820,266,011,846,253,815,571 | 0.996408 (the best) |
| 29 | `{1,28}` | 86,367,908,592,651,965,660,383 | 0.991216 |

`H(R_27)` and `H(R_29)` are from pass 1 (`a10_rot.out`). With the session's other neighbours, `R_27` and `R_29` are strict local maxima. Only the best and first flips are independently confirmed; the rest are single-source but not needed for strictness, since the best flip is below `H(R_n)`.

**P15.** `R_31`, flip `{15,16}`: `H = 20,284,103,532,402,676,770,852,133` (`a14_rot31flip.out`, 69,273,668 necklaces, 3 min). This equals the session's `rot31_H.txt`. The ratio to `H(R_31) = 20,351,875,726,012,821,643,802,839` (pass 1, three sources) is 0.996670, the best neighbour. So `R_31` is a strict local maximum, given the session's other 14 neighbours, all of which are lower.

## Required edits (old → new)

* **E1 (MISTAKE-660, first bullet).** After "At `p = 19` the odd lengths `5 ≤ k ≤ 17` still favour Paley; …", add: `At p = 23 the same holds for 5 ≤ k ≤ 21; at both primes P_q is the unique maximiser over all circulant classes for every odd 5 ≤ k ≤ p − 2 (audit pass 2, a12b_analysis.out).`
* **E2 (MISTAKE-660, THM-161 bullet).** `**THM-161**'s "Paley maximizes EVERY c_k simultaneously (Savchenko's theorem for DRTs)" overstates Savchenko, whose results compare doubly regular tournaments with `RLT_n`.` → `**THM-161**'s "Paley maximizes EVERY c_k simultaneously (Savchenko's theorem for DRTs)", THM-162's "Savchenko: c_k(Paley) = max c_k over all DRTs", and TANGENTS T134 / INV-053's "Savchenko (2024) proves c_m(DR_n) > c_m(RLT_n) for ALL m ≡ 1,2,3 mod 4" misstate Savchenko: arXiv:2403.07629 treats m = 8 only. Since P_19, P_23 are doubly regular and R_n = RLT_n, c_n(RLT_n) > c_n(DR_n) at n = 19, 23 refutes the T134 paraphrase at m = n.`
* **E3 (THM-4751 §7, split table).** `**The split of `ΔH = H(R_n) − H(energy minimum)`**` → `**The split of `ΔH = H(R_n) − H(flat)`**, flat = the energy minimum (at n = 21: Z_21{1,2,3,7,10,12,13,15,16,17}, one of the two energy minima; n = 15, with four energy minima, is omitted)`.
* **E4 (THM-4751 §7 and note §13).** `and the gap grows with `k`` → `and the ratio α_k(R_n)/α_k(flat) grows with k (at n = 19: 1.043, 1.126, 1.317, 1.846 for k = 3..6)`.
* **E5 (note §13, "Consequences").** `This gives characterization (1) of `{3, 7, 11}` in §4 a tournament mechanism:` → `This gives characterization (1) of `{3, 7, 11}` in §4 a mechanism sketch (an observation for n ≤ 25, and the threshold q ≤ 11 restated):`
* **E6 (verification pointers).** In THM-4751 §7, the note §13 and MISTAKE-660, add: `Independently verified: c_k and C by Held–Karp and by inclusion–exclusion (all classes, n = 19, 23), α_k by explicit enumeration (n ≤ 13) and by ranked zeta transforms (n ≤ 21) (lattice_bridges_audit_pass2_20261022.md).`

## Not done

* Explicit enumeration (method E) stops at `n = 13`; `n ≥ 15` uses method Z only.
* At `n = 21` only `R_21` and the two energy minima were computed. The champion claims are stated only for `n ≤ 19`, so this suffices.
