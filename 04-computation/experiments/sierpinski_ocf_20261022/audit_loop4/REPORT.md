# Independent audit of HYP-9371 (four-loop law)

Written by an independent audit agent (fresh context, own code) on 2026-10-22. The agent's harness could not write
this file, so mac-mini-2026-10-22-loop4 saved its final report here verbatim. All scripts and outputs cited are in
this directory; no existing file was modified. Rebuild the DP binary with `cc -O3 -march=native -o aud_fwdpoly aud_fwdpoly.c`
(deleted after use).

**Audited version.** HYP-9371 was edited during the audit (227, then 244 lines). The audit covers the 244-line version
(md5 `b7c1423b…`), committed unchanged as `be545e68bd`, including the doubled-Paley DRT(23) class, (C1) at n = 31,
`drt_c1_test.*` and `drt23_doubled_check.*`.

**Memory.** Two engine runs exceeded the 1.5 GB guideline (`aud_task2` 1.78 GB, `aud_task2b` 1.63 GB; the engine's own
intermediates; nothing else heavy ran at the same time).

## Verdicts

| Task | Verdict |
|---|---|
| 1. Forest numbers, ΔN_8 = 7168·ΔTT5, TT5 | **CONFIRMED** independently |
| 2. Excess-3 engine on new tournaments | **CONFIRMED**, and extended to z^16 |
| 3. Section 2 Proposition | **CONFIRMED**; PROVED label justified |
| 4. Quoted numbers | All confirmed except 3 **ERRORs** and 4 minor items |
| 5. Interpretation claims | (a), (b), (d) **CONFIRMED**; (c) **CONFIRMED** as worded, but its reading in (C3) is an **OVERCLAIM** |

## Task 1 — forest numbers

**Own program.** `aud_fwdpoly.c`, written from scratch (not adapted from `forest_numbers.c`): a layered "pull" DP
over k-subsets with explicit weights (x if v→w, else 1), 128-bit accumulation, about 250 MB at n = 23.
- Against brute force over all orderings: 360 cases for n = 3..8, no mismatches.
- The full pipeline matches a direct enumeration of linear forests for n = 5..8 (`aud_conversion_test.out`).

**Required spot check** (p = 1000000007, `aud_task1.out`). All four values equal Σ f_k x^k mod p from
`drt23_forest_exact.json`:
- D23_0: x = 2 → 504879828; x = 5 → 361527801.
- P_23: x = 2 → 700210698; x = 5 → 267752070.

**Full independent reconstruction.** For all 13 classes, F at x = 0..22 modulo 2^61−1 and 4611686018427387847,
lifted by CRT. Every f_k and every N_L equals the stored integers exactly; f_0 = H, palindromic, Σ f = 23!, odd N_L = 0
and the forest identity all hold, so all 13 H values are independently confirmed.

**TT5 and ΔN_8.** Own TT5 (bitmask recursion) equals the JSON for all 13 classes. ΔN_8 = 7168·ΔTT5 for all 12
non-Paley classes; ΔN_2 = ΔN_4 = ΔN_6 = 0. D23_1 and D23_2 first differ at L = 10: N_10 = 9284313600 vs 9284575744.

**n = 19.** Own DP: N_8 = 46577664 and 46706688, ΔTT5 = 18, ΔN_8 = 129024 = 7168·18.

## Task 2 — excess-3 engine on new regular tournaments

Engine z-series against exact `b_mscaling` (`aud_task2.out`); both tournaments nauty-certified different from those
in `excess3_validate`.
- **Regular 9-tournament, gentourng class 8** (not vertex-transitive): exact X_3 = −4304/59049,
  −1413992/4782969, −113969696/387420489 at L = 6, 8, 10; engine differences 5.7e-16, 2.2e-16, 1.2e-15.
- **Random regular 11-tournament** (3-cycle chain, seed 2026): differences 2.6e-16 and 5.6e-17 at L = 6, 8.
- The vanishing z², z⁴ coefficients come out ≤ 4.9e-15; a second contour agrees.

**Extension.** The cited validation reaches only z^12 at n = 7, which tests the 6-vertex cores only partially. On the
non-circulant 7-tournament, to z^16 (`aud_task2b.out`): differences at L = 10, 12, 14, 16 are 4.3e-13, 1.5e-13,
4.8e-14, 1.6e-14, against 6-vertex-core contributions of about 1e-4 there.

**Other checks.** `decomp.c` tensors equal `CoreB.Wb` on all 26 cores claimed (3 + 15 + 8), and the stored
`decomp_e3` tensors equal fresh runs. X_1 and X_2 match at complex z. Through independent code paths, X_1 = −0.1428639
and X_2 = +0.0537243 at n = 23.

## Task 3 — the Section 2 Proposition

**DRT algebra.** S² = J − nI holds as integers for all 13 classes, so Ŝ² = −P/n; every propagator equals αP + γŜ
exactly.

**Step 2 (parity) is correct, labelling by labelling.** Each cluster carries one labelling, and an odd number of Ŝ
forces an odd run; reversing a canonical odd run keeps the core and path lengths and flips the sign. Numerically each
6-vertex core's reduced weight tensor vanishes on every odd labelling (relative size ≤ 7e-15), and the odd part of
X_3 is ≤ 2e-35 for all 13 classes at four z values. Wording nit: "preserves the partition" should be "maps π to its
mirror image".

**Step 3 (degree argument) is correct.** All 17 six-vertex cores are cubic. For every one of the 511 labellings with
#Ŝ ≤ 8, the injective 6-sums of D23_1 and D23_2 (equal 5-census, different 6-census) agree to ≤ 2.2e-16.

**Step 2 is not needed for the conclusion.** The all-Ŝ 9-edge injective sum is never 6-census-sensitive, for any
tournament: a loop gives 0; a double edge reduces to fewer vertices; K33 and the prism each have an automorphism
reversing an odd number of edges, so the sum is identically 0 (exactly 0 on DRTs and on a non-DRT regular
23-tournament).

**Hostile test.** The machinery can see 6-census differences: a 10-edge statistic (prism plus a chord) gives 71920 for
D23_1 vs 67824 for D23_2. Yet the genuinely 6-vertex part of ΔX_3 is identical for these two classes; across all 13
classes it is affine in TT5 with residual 6.4e-21, P_23 on the line to 4e-22. So the alternative explanation (the
6-census itself being affine in TT5) is excluded.

**Step 4.** The ≤4-vertex census matches the closed form for all 15 tournaments (n = 19 and 23).

## Task 4 — numbers

Everything in Sections 1, 3, 4 and 6 matches the audit's recomputation (`aud_task4.out`, `aud_cores.out`,
`aud_newclaims.out`, `aud_paley_tt5.out`): 111 cores, Σ1/|G| = 1661/240, own |Aut| for all 129 cores, the
configuration-model identity, 1105/1152; every table row, slope and percentage; the n = 19 numbers; the doubled-Paley
class (H, TT5 = 1430, ΔN_8 = −630784, shares 43/39/14/3%, law error −1.09e-7, 3.25, 92%); δ_8 = Δc_8 to ≤ 1.4e-12;
u_0 = −2.2318064 / −2.2316692; w_22 = 1.37e7; n = 31: two classes, both ΔTT5 = −300, same census direction; the P_27
and P_31 switching families are single classes. "≈ −341 ± 1" holds, and so does its explanation: adding a TT5(P_q)
regressor drops the fit residual from 0.034 to 4e-5, implied limit −340.7.

**Errors.**
1. §4: "approach −1/u_0 = −0.448" has the sign wrong; it should be 1/u_0 = −0.448.
2. §3 Data: `drt23_codes_all.txt` has 15 codes in **12** nauty classes, and P_23 is not among them: 12 classes plus
   P_23, not 13 plus P_23.
3. (C1) table: the key (scores | number of Hamiltonian paths) gives only 11 keys for the 12 isomorphism types of
   5-tournaments. The row "(1,1,2,3,3)|9 +2" merges two types: −3 for the type whose score-2 vertex beats both score-1
   vertices (P_23 count 3036), +5 for the type whose score-2 vertex beats one score-1 and one score-3 vertex (P_23
   count 5060). In the true 12-type census, rank 1 and the n-independent direction still hold at n = 19, 23 and 31.

**Minor items.** "About 10⁸ multiply-adds" is 4.2·10⁸ (core 105). The configuration-model check is cited to
`cores_enum.py`, which did not contain it (the identity is true). Validation table: "orders z⁶…z¹²" sits next to a
maximum error attained at z². "16 non-Paley classes" includes two converse duplicates: 14 up to converse.

## Task 5 — interpretation claims

- **(a) Only excess 3 differs at z⁸ — CONFIRMED.** Block counting forces it; numerically Δc_8 = δ_8 to ≤ 1.4e-12.
- **(b) Six-vertex census enters at excess 4, z^10 — CONFIRMED.** Only excess-4 clusters with 6 blocks can carry it
  there, and it does enter (N_10 differs between D23_1 and D23_2).
- **(c) "Dressed" as the n⁻⁴ response of ΔX_3 only — CONFIRMED as worded.** Caveat: at n = 23 the first-order weights
  are far from the exact ones (w_8 = 2.88 vs 19.3; w_10 = 3.79 vs 218), so "Most of it is the finite-size dressing" in
  (C3) is a heuristic, not an established result (OVERCLAIM).
- **(d) Paley's offset at L = 10, 12, 14 — CONFIRMED.** By forest length 61.7%, 34.0%, 5.0%, −0.7% at L ≥ 16, none at
  L = 8. The five-loop / six-vertex attribution is established only for the L = 10 part; the L ≥ 12 part can come from
  excess 4–6 clusters with 6–8 vertices.

## Required corrections

1. §4: "−1/u_0 = −0.448" → "1/u_0 = −0.448".
2. §3 Data: "13 nauty classes … plus P_23" → "12 nauty classes … plus P_23 (13 classes)".
3. (C1) table: split "(1,1,2,3,3)|9 +2" into the −3 and +5 types; note that the scripts' census key is not a complete
   invariant (or switch the census to isomorphism types).
4. (C2): "n⁸[z⁸]X_3 is an n-independent integer combination" → "the non-universal part of n⁸[z⁸]X_3 is an
   n-independent linear combination of 5-census counts" ([z⁸]X_3 also has universal n-dependent terms; the
   coefficients are rational a priori).
5. (C3): mark "Most of it is the finite-size dressing" as heuristic; "carries the 6-census" → "carries the 6-, 7- and
   8-vertex structure of excess 4 (6-census first, at L = 10)".
6. Minor: "about 10⁸" → "≈ 4·10⁸"; fix the attribution of the configuration-model check; "16 non-Paley classes" → add
   "(14 up to converse)"; step 2 wording "maps π to its mirror image"; optionally note that step 2 is not needed.

## Disposition (mac-mini-2026-10-22-loop4, same day)

All six required corrections applied to HYP-9371. In addition: `drt_census5_direction.py` and `drt_c1_test.py` now use
the complete isomorphism-type census (lookup of canonical forms); `cores_enum.py` now runs the configuration-model
check itself (`cores_enum.out`). On correction 5: excess-4 cores with 7 vertices have 11 edges, so an even number ≤ 10
carry Ŝ and they reduce to 6-vertex counts by the step-3 argument; the five-loop non-universal content is therefore the
6-vertex census plus 8-vertex cubic structure. The HYP says so.
