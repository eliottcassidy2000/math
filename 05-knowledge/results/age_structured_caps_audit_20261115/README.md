# Age-structured caps, the audit of the sharp growth lemma, and six typed threads

Session klein-2026-11-15-age (2026-11-15/16), on the owner's directive "think age-structured potential for the sharp
lemma Z ≤ 38H; prioritise creative reasoning and past-session exploration; MRB and the logarithmic repair; {7, 21, 63},
{27, 63}, 27 → 53 → 105; (8, 15, 17) scaled by 3; Mersenne, Brier, Riesel, 78557, 509203, {3, 5, 7, 13, 17, 241}; 24,
Leech, Ramanujan". Canon entry: THM-5253. Updated: HYP-9803 (b) → PROVED, THM-5243 status.

## Status

- **AUDIT, SOUND** ([AUDIT.md](AUDIT.md)). The lemma `H(T) ≤ H(T+x) ≤ 39·H(T)` for span-3 ribbons, proved on
  2026-11-15/16 by another session (`chanpak_target_ribbon_audit_20261116/SHARP_GROWTH.md`, tail deletion) with the
  per-letter sharp constants 1, 9, 17, 19, 33, 33, 39, 31 (`chanpak_eight_braid_20261116/SHARP_EIGHT.md`), is correct;
  both checkers replay; the exhaustive maxima over all `8^n` words of length `n ≤ 9` reproduce the eight constants letter
  by letter, with `max Z_ℓ/H = C_ℓ − 1`. HYP-9803 (b) is therefore PROVED, and the "age-structured potential" proposed as
  its proof device is withdrawn. The sharpness fixture of the proof is the block `B` of THM-5243 (2), whose count
  `39 = 5·2³ − 1` THM-5243 (2′) proves in closed form.
- **PROVED, new** ([CAPS.md](CAPS.md)). Age-structured caps: every one of the 119 profile coordinates satisfies
  `x_s(T) ≤ κ(s)·H(T)` with `κ(s)` an explicit rational depending on the ages in the profile (`data/kappa_table.txt`);
  `max κ = 30`, two-fragment states `≤ 68/3`. Every raw transfer coordinate is at most `30·H` (THM-5252 (B) had 26,425).
  The proof is the tail deletion applied fragment by fragment over the last six old sites.
- **FINITE-EXACT** ([CAPS.md](CAPS.md) §1, §3). Observed maxima of `x_s/H` over words of length ≤ 9: 14, 14, 7 for
  two-, three-, four-fragment states; 16 states have `κ` = observed; three caps still grow between lengths 7 and 9
  (13/9 → 25/17, 8/5 → 44/27, 7/5 → 25/17), so some suprema are not finite-length maxima.
- **REFUTED, pre-registered** ([CAPS.md](CAPS.md) §4, [PREREGISTERED.md](PREREGISTERED.md)). No history-refined
  box cone (history length 0–3, caps from the data) is letter-invariant; value iteration stalls on dead directions
  (`K = 0`) or diverges (`K = 1`). Conceptual reason: the deletion proof uses the product inequality
  `H(T−S)·H(T[S]) ≤ H(T)`, information about sub-tournaments that no linear functional of the state carries.
- **Typed threads** ([THREADS.md](THREADS.md)). T1 MRB/logarithmic repair: ANALOGY with one exact new ingredient (the
  per-vertex log-gain ceiling `log 39`, saturated only by the sealed gadget) and one exact disanalogy. T2 {7, 21, 63},
  {27, 63}, 27 → 53 → 105: STRUCTURAL arithmetic (27 and 63 are the two ternary dressings of the forbidden 21; the dyadic
  ray `13·2^e + 1` and the ternary ray `7·3^j` are disjoint), tournament content = the cited inverse ports, OPEN beyond.
  T3 (8, 15, 17)·3, five pieces on 1…2024, 1024 = 10³ + 3·2³: decomposition and repair PROVED (cited); the chain
  24 → 1024 → 2024 is COINCIDENCE. T4 Riesel/Sierpiński/Brier wheels: PROVED (cited, replayed); new classical remark —
  the Riesel covering set is the set of primitive prime divisors of `2^d − 1` over `d | 24`, six-pronged because
  `2^6 − 1 = 63` is the Zsigmondy exception. T5 24/Leech/Ramanujan: NUMEROLOGY (MISTAKE-559 policy). T6 the potential:
  withdrawn as proof device, PROVED as a cap theorem, REFUTED as a cone certificate.

## What is new

The audit (an independent check plus the exhaustive corroboration), the age-structured cap theorem with constant 30,
the growing-caps phenomenon, the refutation of history-refined cones, and the Zsigmondy reading of the six-prong wheel.
Known and used: THM-002 (odd-cycle packing identity), THM-4880 (the 119-state transfer), THM-5243, THM-5252,
SHARP_GROWTH.md, SHARP_EIGHT.md, and the packets cited in THREADS.md.

## Files

- [`scripts/hratio5.c`](scripts/hratio5.c): exhaustive maxima of 138 features (coordinates, fragment-class sums,
  `Z_ℓ`, `H(Wℓ)`) relative to `H`, keyed by the last `K` letters, over all words of length ≤ L (sparse DFS, exact
  fractions; eight restricted runs by first letter). Receipts `hratio5_9_b*.out.gz` (L = 9, K = 3), `hratio5_10_b*.out.gz`
  (L = 10, K = 0).
- [`scripts/profile_caps.py`](scripts/profile_caps.py): the deletion bound `κ(s)` for every profile, exact, maximised
  over the 4,096 orientations of the last six sites. Receipts `profile_caps_frag1to3.out`, `profile_caps_frag4.out`.
- [`scripts/cone_test.py`](scripts/cone_test.py): letter-invariance of history-refined box cones by linear programming
  (scipy HiGHS), with value iteration. Receipts `cone_test_K0_K1.out`, `cone_test_K2.out`, `cone_test_K3_slice*.out`.
- [`scripts/clocks_check.py`](scripts/clocks_check.py): exact arithmetic for THREADS.md (sympy). Receipts
  `clocks_check.normal.out` = `clocks_check.optimized.out`.
- `data/states3.txt` (the 119 profiles in matrix order), `data/kappa_table.txt`, `data/caps_compare.json`.
- `receipts/replay_check_sharp_growth.out`, `replay_check_sharp_eight.out`: the other sessions' checkers, replayed.
- The transfer matrices are `chanpak_span3_growth_two_apex_seam_20261115/data/mats3.txt.gz` and `frag.txt`.
