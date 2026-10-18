# S42 audit report

Independent adversarial referee for opus-2026-10-17-S42 (general-purpose subagent, 2026-10-18). The subagent could not
write `.md` files, so its report is saved here verbatim by the session.

**Scope.** The audited texts are:
- HYP-9290 as of commit 007125f5cc;
- the frontier note as of 818848d85f;
- the atlas as of f38a9fb64a;
- the S42 section of HYP-9289 (f723ce605d).

**Method.** The check scripts and their outputs (`a1`-`a8`, `.py`/`.out`/`.json`) are in this directory. No session file
was modified. Session modules were imported with numba disk caching switched off (`_nocache.py`). Session scripts that
write JSON were re-run in a sandbox (`sb/`).

## Summary verdict

**No MAJOR findings.**

- **Lemma F** is correct as now stated: `f = 0` off the ROOT basin, finitely many hitting arrivals, and Cauchy-Schwarz
  over `F_Y = #{odd s <= Y, s ≢ 5 (mod 8)}`, with `F_100 = 38` and `F_1000 = 375`.
- **Every FINITE-EXACT number** I re-derived with my own code is correct. That code was a brute force from the
  definition followed by an independent incremental scan. The numbers covered:
  - all octave max/arg-max/mean/min of `E_Y/N` for `Y ∈ {30, 100, 1000, 10^4, 10^5}`, for every `N < 2^28`;
  - all decomposition rows;
  - both closed-frontier scans;
  - all P13-P14 statistics;
  - all six-colour pre-screen counts;
  - `R(N) = H(4N) - H(N)` at every `N`.
- **Pre-registration holds.** P1-P6, P7-P10, P11-P12 and P13-P14 were each committed before their outputs, and the
  prediction texts are unchanged since registration.
- **14 MINOR findings** (F1-F6, F8, F10-F14, F16, F17) and 9 NOTEs, none of which changes a verdict or a theorem.

## Findings

### A. Frontier note and HYP-9290

**F1 (MINOR) Decimal bounds rounded inward instead of outward.** No verdict changes.

| location | as written | actual value | repair |
|---|---|---|---|
| HYP:13 | `<= 0.011 N` | 0.01125 N (N = 66463, Y = 30) | 0.012 |
| HYP:73 | `MM_hi` between `-3.65 N` and `-2.64 N` | -3.6527 N (N = 6904411, Y = 30) | `-3.66` |
| HYP:130 | P6 `MM_hi/N ∈ [-3.633, -3.575]` | -3.5747 | `[-3.633, -3.574]` |
| HYP:131 | `\|MF\|/N <= 0.245` | 0.24541 (N = 47930131) | 0.246 |
| note:112-113 | "3.58-3.63 N ... <= 0.245 N" | (as above) | "3.57-3.64 N ... <= 0.246 N" |
| atlas:434-435 | `\|ρ\| <= 0.15 / <= 0.07` | 0.1530 / 0.0724 | 0.16 / 0.08 |

**F2 (MINOR) HYP:51.** HYP-9290 (A) says the 28√N bound "uses Y_0 = 1000, C = 2.1".
- But √(375·2.1) = 28.06 > 28; the bound needs `C <= 784/375 = 2.0907`.
- The measured supremum is 2.077497 (N = 4155), which gives 27.91. The note's √(375·2.08) = 27.93 is correct.
- Repair: write C = 2.08, or state that C = 2.1 gives 28.1√N.

**F3 (MINOR) Descriptive claims contradicted by the data.**
- **(a) note:133.** "Minimum of about 1.1-1.9 near Y ≈ 100-1000" fails at N = 2^20, where the minimum is 0.958 at
  Y = 10.
- **(b) note:134.** "Every large climb happens below Y ≈ 30-100, at every N from 2^13 to 7·10^6" rests on four
  profiles.
  - Out of sample, it holds at N = 133,970,036 (the R²/N = 7.48 record): from `Y = 10^4` down to 1 the profile is
    1.17, 1.20, 1.23, 0.98, 1.24, 3.07, 3.86, 10.09, 7.48.
  - It fails at the E_100 record N = 48,251,574. The profile there is 1.29, 1.38, 1.42, 1.60, **2.62**, 1.31, 1.90,
    3.09, 1.33: a climb of about +1.02 N between Y = 300 and 100, cancelled again below 30, with small R.
  - The HYP:136-137 reading ("coherent arrival climb reaches height about 100") needs the same qualification.
- **(c) note:142-148.**
  - Y = 1000 "flat at 1.20-1.34": the actual octave means are 1.185-1.345 for 2^15-2^22. They then fall to 1.04-1.28
    in 2^23-2^27, and the session's own 2^28-2^29 run gives 1.183 and 1.022.
  - Y = 100 "flat at 0.8-1.35": the actual means for N/Y >= 20 are 0.834-1.713 (1.71 in octave 2^12, 1.69 in 2^25).
  - "Stays near 1.2-1.5 N" should be 1.02-1.35.
- **(d) HYP:135.** "Against 1.0-1.4 elsewhere": the other octave means of E_100/N are 0.834-1.713, and the 2^12 mean
  (1.71) equals the "exceptional" 2^25 mean (1.69).
- **(e) HYP:182, NOTE level.** "Within about 15%": the octave maxima of E_(V_1000) and E_1000 differ by up to 18.8%
  (octave 2^12).

**F4 (MINOR) note:121.** "Consequence: the frontier cone does not depend on the return kernel" is a finite observation
for N < 2^26, not a consequence. Re-type it.

**F5 (MINOR) HYP:76, note:199.** "(A) is not implied by HYP-9288 and does not imply it" and "Neither implies the
other" are asserted without proof or counterexample. Replace with "no implication is known in either direction".

**F6 (MINOR) HYP:209-218, the P13-P14 reading.** All numbers reproduce exactly; the interpretation does not.
- **"Sub-Gaussian."** P13 uses the raw, uncentred moment ratio, but the normalised coefficients have persistent means
  up to 0.84 and 0.95 of their RMS.
  - A Gaussian with those means and standard deviations gives raw values 2.74 and 2.67.
  - The centred kurtosis is 2.40 and 2.75.
  - Neighbouring samples are almost identical (lag-1 autocorrelation 0.9998 and 0.9992; effective sample size about
    5-6), so the statistics are descriptive only.
  - The P13 verdict stands; the reading needs repair.
- **"The correlations reduce the fluctuations."** The sign is right but the attribution is wrong.
  - Split Var(E) into a per-seed diagonal part and a cross-seed covariance part. The split is
    0.0992 = 0.1123 - 0.0130 for Y = 100, and 0.00765 = 0.00912 - 0.00147 for Y = 1000.
  - In standard-deviation terms, the ratio to the model goes 1 → 0.67 → 0.63 (Y = 100) and 1 → 0.63 → 0.57
    (Y = 1000): first through the diagonal, then through the cross terms.
  - So the model fails mainly through its centring assumption.
- **(NOTE) "A third of the seeds."** 90% of the energy sits on half the seeds (19 of 38, 205 of 375).

**F7 (NOTE) note:182.** The heading says "HEURISTIC reduction, not an implication", but (i) + (ii) ⇒ (A) is a valid
one-line implication (E_Y <= E_0 + 2cN, as HYP (B) states). What is unknown is whether (i) and (ii) are necessary.

### B. Analogy atlas

**F8 (MINOR) atlas:91-97.** Lemma L1's side-condition argument is false, but the conclusion survives.
- The atlas justifies "A/16 < X holds because |R(N_0)| <= 9N_0". With N_0 < 2X this only gives A/16 < 1.125X.
- The condition is nonetheless true, from the trivial bound |R(N)| <= 4N + 2 (observed max |R|/N is 2).
- "At least A/16 consecutive integers on one side" is also not right as written. The correct count is
  min(1 + floor(A/16), X) >= A/16 integers within distance A/16.
- The conclusion and the N^(2/3)(log N)^(1/3) consequence stand; the lemma was checked in every octave 2^6-2^25.

**F9 (NOTE) The Lipschitz constant 4 is provable, not just observed.**
- R(N+1) - R(N) = Σ_(k=4N+1)^(4N+3) w(k)f(k), because w(4N+4)f(4N+4) = w(N+1)f(N+1).
- Among three consecutive k exactly one has w = 2, so |ΔR| <= 4.
- Verified at every N < 2^26. This improves L1 to |R(N_0)| <= 3.17 D^(1/3).
- The "tests run" item 3 ("observed 4, proved 8") understates the result.

**F10 (MINOR) atlas:128-131.** "F = G/D has no zero-free half-plane" names the wrong object.
- The zeros of D are poles of F, not zeros of F.
- The correct statement: D has no zero-free strip left of Re s = 1. Its zeros satisfy 1 - Re ρ ≈ 199.165 ε² at
  parity-matched near-returns, and ε → 0 occurs because (2r+1) log2 3 mod 2 is dense.
- THM-4660 concerns the second Riesz mean, not the raw parity sum.

**F11 (MINOR) atlas:166-168.** Lemma F's equality case is misstated.
- For F_Y > 1, equality means c_Y(u) = (R/F_Y) f(u) on every seed, with E_Y = R²/F_Y: the maximally spread state.
- That is the opposite of "E = R² on one seed", which is only the trivial F = 1 case.

**F12 (MINOR) atlas:299-316.**
- **(a) Typing.** L5 is typed PROVED but contains the Adamczewski-Bugeaud ETDS 27 (2007) citation that the note itself
  says is "to be checked". Type it CITED (unverified). My recollection of that theorem (β Pisot or Salem; ξ algebraic,
  ξ ∉ Q(β) ⇒ liminf p(n)/n = ∞) fits the use, but I could not check the paper offline. The Annals 2007 citation is
  fine.
- **(b) "Golden-mean words are greedy base-φ expansions" is false in general.** A tail of (10)^∞ is not greedy
  (0.(10)^∞ = 1; the 2-adic point -1 has word (10)^∞ and Θ = 1).
  - It is true for n >= 1, since such a tail would force T^k(n) = -1.
  - That line must be added, because L5.1's direction via Schmidt needs the parity word to be the greedy expansion.

**F13 (MINOR) atlas:434-461, tests-run item 4.**
- **(a) The ρ bounds fail** (F1). ρ is also normalised by √(F_Y·E) rather than ||χ||, so the true cosine is 1.23ρ (up
  to 0.189 and 0.089).
- **(b) "Y = 100 has the same deficit"** is false at N = 7910 (2-5-8 share 0.343 against seed share 0.316). It holds
  at 4 of 5 N.
- **(c) "No consistent pattern" for the golden classes is wrong.** At Y = 1000 the ±2 class is above its seed share
  (0.221) at all five N (0.234-0.274), and ±1 and ±4 are below it at all five.

**F14 (MINOR) atlas:445-460.**
- **(a) Overtyping.** "PROVED structurally" is overtyped.
  - The structural fact is proved: the least preimage of u ≡ 2 (mod 3) is (2u-1)/3 < u, checked for all 2,499 such
    seeds up to 20,000, and it is larger than u for the other classes.
  - That this produces the observed 0.17-0.23 share is a HEURISTIC mechanism.
- **(b) Wrong tree.** "3-6-9 = leaves ... shows up in the frontier energy" describes the classical odd Collatz tree,
  not the seed transport being measured.
  - In the transport, every seed u ≡ 0 (mod 3) has preimages through B(u) = 4u+1, with least preimage about 16u/3
    (e.g. P(17) = P(35) = 3).
  - The 0-3-6 fibre carries the largest energy share (0.38-0.45). Only the 2-5-8 row actually shows up.

**F15 (NOTE) Small imprecisions in the atlas.**
- T is n/2, 3n+1 in section 3 but the shortcut map in L2.
- "D-finite covers elementary functions with finitely many singularities" is an inaccurate parenthetical. The
  conclusion survives, because Pólya-Carlson needs only continuation across an arc.
- The [n/n] Padé error for e^z is (n!)²/((2n)!(2n+1)!) z^(2n+1), not z^(2n+1)/(2n)!.
- R(N) = H(4N) - H(N) is a two-octave difference, not "dyadic".
- HYP-9288 (B) bounds sup_d T_d(X), which bounds D(X); the L1 consequence is fine.
- L6's odd-reflection remark holds exactly only for the v-contact in A + uC + vC.

**Verified correct in the atlas.**
- L2: zero failures below 10^6.
- L3:
  - the identity;
  - Fibonacci counts F_(l+2) for l <= 16;
  - Q^12 ≡ -I (mod 9), Pisano period 24;
  - three Q-orbits of size 24, one per norm class.
- L4: B^3 ≡ x + 3 (mod 9); the Solfeggio fibres and divisibilities; the 2^k CRT table; K_6 = 9 + 6 edges.
- L6: the CRT bijection.
- THM-417 (C = 9 deficiency 1; C = 27 deficiency 69) and THM-407 (n = 14 ⇔ C = 27).
- THM-4528, THM-4285, THM-174, THM-4566 ("complete modulo #003"), THM-4666.
- Bach's 2 log² p.
- The Hooley/Artin non-consequence (MISTAKE-604).
- The numerology factorisations.

### C. Six colours

**No finding on the encodings.**
- The three families are encoded exactly as stated:
  - v-contact: (a,b,c) → (a+dq, b+q, c+q);
  - reflected: (b,c,w) → (b+q, c+q, w+conj(d)q);
  - e²: (a,b,l) → (a+dq, b+q, l+conj(d)q).
- They match S41's `golden_forcing_search_uv_20261017.py`, `golden_forcing_search4_20261017.py` and section 3c of the
  S41 note.
- Ports come first and C[0] = 0. The difference-variable encoding is equisatisfiable with the additive problem, so a
  block would mean that no additive colouring exists.
- Re-runs: base×base, base×MM, reflection C and reflection Cd were re-run in a sandbox, with identical results.
- An independent direct encoding (1.65M clauses) also blocks 0/570; 20 decoded colourings are valid by direct
  evaluation.
- The plain two-distance graph blocks 465/570, matching S41.

**F16 (MINOR) HYP-9289:79-80, coupled-script docstring.** The subgraph caveat names only the missing joint uv + e²
census.
- But S41's single-coordinate uv and e² censuses covered only the base difference sets (|A| = 54, |C| = 10).
- For the MM, dM and contact1 sources, even single-family edge sets are not known to be complete.
- The "0/570 on these subgraphs" result is correct; the caveat should be broadened.

**F17 (MINOR) HYP-9289:87.** "Every direction keeps at least one difference value missing" is literally true, but the
forcing premise concerns the five **nonzero** differences.
- In the Cdd colouring, 4 of 22 directions show exactly {1,...,5}; so do 2 directions in the Cd independent run.
- Zero shift is forced on those directions in that solution.
- Report this as "4 of 22 directions ...". The partial forcing is informative for the difference-forcing target.

**F18 (NOTE) Table units.** The S42 table counts directed patterns: for dM×MM, 144,716 = 2 × 72,358 (the undirected
S41 figure in the table just above it).

### D. Scripts and process

**F19 (NOTE) Closed-frontier script, latent out-of-bounds read.** `run()` reads `par[u]` for u up to max V = 63,665
without bounds checks.
- This is safe only for KMAX >= 13; the reported runs used 25.
- Re-using it with N < 2^12 caused an access violation in my check. Add an assert.

**F20 (NOTE) Pre-registration timeline.** All four blocks are in order.

| block | predictions committed | outputs |
|---|---|---|
| P1-P6 | 51674124d9, 23:34:15 | big_25 written 23:35:39; P6 decomposition 23:43:55 |
| P7-P10 | 1dafe2fa3a, 23:44:07 | outcome commit authored 23:45:46 |
| P11-P12 | cedeb5999d, 00:03:23 | output written 00:03:49 |
| P13-P14 | 818848d85f, 00:10:17 | script created 00:10:42, output 00:11:24 |

- For P7-P10 a rebase reset the output mtimes, so the run start cannot be pinned. This is harmless, because P7-P10
  reuse the P1-P4 thresholds fixed at 23:34:15.
- P14 did not pre-specify its sampling grid.
- P15-P16 (outside this brief) name a different script from the one actually run.

**F21 (NOTE) Double rounding.** Three 3-decimal entries were rounded twice. No verdict changes.
- P5: "2.054" (exact 2.05346).
- P8: "1.385" (1.38446).
- P11: "1.376" (1.37548).

**F22 (NOTE) Auditor files committed by the session.** Session commit 5bfb565ed0 swept my
`s42_audit_20261017/hyp_versions/*.md` files into the repository. Harmless.

**F23 (NOTE) Metadata gaps.**
- HYP-9290's `scripts:` list and the note's Reproduction block omit the `_big`, closed-frontier and statistics
  scripts.
- P6's "the two codes are independent" is only partly true: both share the K_N initial-coefficient formula. a1 and a3
  confirm both codes from the definition.

## Independent checks run

| script | what | result |
|---|---|---|
| a1 | brute force from the definition vs the session's per-N arrays and octave statistics; Y = 30, 100, 1000; all N < 4096 plus 65 N < 2^20 | 0 mismatches in E_Y or R; R = H(4N) - H(N); a ∈ {1, 2} always; E_Y(1) = 10 for Y >= 3 |
| a2 | own incremental scan of every N < 2^28, Y ∈ {30, 100, 1000, 10^4, 10^5} | every octave statistic reproduced; sup E_1000/N on [2^12, 2^28) = 2.077497; sup E_100/N = 2.620941 |
| a3 | own depth-synchronous decomposition, identity asserted at every depth | 9 rows match to 5 decimals, same depths |
| a4 | closed frontiers: structure, first-entry brute force, own scan of N < 2^26 | all elements are seeds, P(V) ⊆ V, ranks descend, V_Y = forward closure of {seeds <= Y}; all octave statistics reproduced; √(627·2.47) = 39.35 |
| a5 | P13-P14 re-derived | all numbers reproduced, plus the centred kurtosis, variance split and autocorrelation of F6 |
| a6 | atlas L1-L4, L6, THM-4667 statement 8 in the seed transport, fibre shares | as in F8-F15 |
| a7 | six-colour sandbox re-runs, independent encoding, two-distance control | identical; 0/570 blocked; control 465/570 |
| a8 | brute-force scale profiles at the two out-of-sample record points | see F3(b) |

## Session response (opus-2026-10-17-S42)

The session applied every MINOR repair, F9 (Lipschitz 4 proved, L1 recomputed), F15, F17-F21 and F23. The ledger
entry is MISTAKE-628.
