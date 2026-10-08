# Run-compressed pair chains, the orphan law, and the reducible locus

**Session.** mac-mini-2026-10-07-twoanchor (continuation).

**Owner prompt.** "keep pursuing possible next steps and explore around and synthesize new ones that you pursue as well". The previous part is [two anchors](twoanchor_reset2_friezes_20261007.md) (THM-4600–4602).

The concurrent codex sessions build directly on THM-4601: [general head phases](collatz_general_head_phases_20261007b.md) and [two-twos grounding](collatz_two_twos_grounding_20261007c.md). They ground specific Mersenne children. This note does not duplicate that. It takes the theory, the positive-integer side, and the frieze question further.

COLLATZ IS STILL OPEN.

## 0. Results and types

| Item | Statement | Type |
|---|---|---|
| [THM-4603](../../01-canon/theorems/THM-4603-run-compressed-pair-chains-drift-lemma-ladder-completeness-and-long-run-transitions.md) (1) | **Drift lemma.** If `u − c = 3^k(v − c′)` with c, c′ cycle points of words w, w′, the orbits run their patterns in lockstep and the debt drifts by `|w|/Σw − |w′|/Σw′` per Terras step. THM-4600 is the zero-drift case. | PROVED; 720 exact tests |
| THM-4603 (2) | **Ladder completeness.** Every equal-time merge, from any state, is a ladder collision: the last distinct odd values are U-preimages of the first common odd value, `z′ = 4^i z + (4^i−1)/3`. So the ladder compiler of THM-4601 (iv) (all i, child and source ladders) is complete. | PROVED; all 1651 observed merges |
| THM-4603 (3) | **Long-run transitions.** A long run of any periodic pattern sends any state to absorption, re-anchoring, an out-of-phase state, or a linear debt drift. From the universal state: a child ones-run re-anchors at −1 with debt −5; a source ones-run drifts at +1/2. | PROVED classification; FINITE-EXACT table (5 states × 55 cycles); actual-integer checks |
| [THM-4604](../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md) | **Collatz is the reducible locus.** Carries form a twisted 1-cocycle, trivialised on ⟨w⟩ by the cycle point c_w. Traces depend only on (\|w\|, Σw) and equal 2cosh(δ/2), where δ is the tuning error. Every pair of words lies on the Cayley cubic x²+y²+z²−xyz = 4 (tr[A,B] = 2). Markov / Conway–Coxeter structures live at tr[A,B] = −2 instead. | PROVED (elementary, standard); the frieze reading is DICTIONARY |
| [HYP-9242](../hypotheses/HYP-9242-orphan-law-deletion-orphans-decay-like-inverse-square-root-of-log-n.md) | **Orphan law.** The fraction of residual sources with no deletion partner decays like `(log n)^(−1/2)`. On the Mersenne line it is about 1/√K (orbits to K = 12800). Orphans are 10× enriched among long orbits, and about 40% are rescued by 3-adic predecessors. | NUMERICAL |

## 1. Ladder completeness (THM-4603 (2))

**Why it holds.**
* Just before two orbits first coincide, their last odd values z_u ≠ z_v are both U-preimages of the first common odd value Z.
* The odd U-preimages of Z are `(2^a Z − 1)/3` for a in a progression of step 2, so they form the ladder `z ↦ 4z + 1`.
* Hence every merge is a child ladder `z_v = 4^i z_u + (4^i−1)/3` (final letters c, c+2i) or a source ladder (the mirror).
* With the debt and equal-time conditions, the heads satisfy `|v| = |u| + k_0` and `Σu − Σv = ±2i`, plus the anchored identity `F_v(c*) = 4^i F_u(c*) + (4^i−1)/3`.
* The merge time is `max(Σu, Σv) + 1`, independent of the final letter.

**Consequence.** The audit's "the i = 1 grammar is incomplete" (MISTAKE-586) is now sharpened: the full ladder grammar is complete, at every anchor and every debt.

**Frequencies.** Among the 1651 merges observed from the universal state of THM-4601 (3000 sources, horizon 3000):

| ladder | i = 1 | i = 2 | i = 3 | i = 4 |
|---|---|---|---|---|
| child | 870 | 50 | 18 | |
| source | 672 | 33 | 7 | 1 |

So source ladders are almost as common as child ladders.

## 2. Drift and long-run transitions (THM-4603 (1), (3))

**The drift lemma.**
* Debts change at the rate of the difference of the odd densities of the two cycles the orbits are near.
* This explains, with one formula, the collapse of the reset children in THM-4601 (v): the source sits on +1 (density 1/2) while the child's limit sits on 0 (density 0), giving drift +1/2. It also explains re-anchoring (zero drift).

**From the universal state (3, 1−27).** These are the exact limit transitions (`transitions.out` has all 55 cycles × 5 states):

| long run | outcome |
|---|---|
| child ones-run | re-anchor at −1, debt −5, after 12 steps (limit pair (−53, −1); 40 actual sources) |
| source ones-run | drift +1/2 after an 8-step transient (child limit 25/27 → trivial cycle; 40 actual sources) |
| child run of 3 or 4 | re-anchor, debt 4 or 9 |
| child run of 5, 7, 8; (1,6), (3,5) | drift +0.27 … +0.47 |
| source run of 5, 7; (1,5), (3,4), (3,5) | drift −0.18 … −0.43: the debt is paid down |
| source (2,4), (2,6), (1,2,5), (1,5,2), (2,2,4), (2,3,3); child (1,1,4), (1,1,4,2) | absorption at fixed times 17–37 (periodic heads) |

**What this gives.** It identifies the heavy tail of the post-run problem: positive-drift runs, above all source ones-runs, which are new Mersenne-like excursions. Section 3 shows that positive-integer orphans reflect aggregate excess odd density rather than single long runs.

## 3. The orphan law (HYP-9242)

**Positive integers.**
* A source n is certified by the deletion child h_D exactly when they merge at equal time before 1. This is equivalent to equal odd-step counts `o(n) = o(h_D)` together with `σ_T(n) = σ_T(h_D) + D` (THM-4556 (ii)).
* An **orphan** has no such D.
* The debt walk is a martingale with tail `T^(−1/2)` (THM-4581), and it has time ∝ log n before the child reaches 1. This predicts an orphan fraction of `c/√(log n)`.

**Mersenne line** (orbits of 2^K − 1 computed for K ≤ 12800; equal odd count plus Terras difference):

| K window | odd K | orphans | fraction | fraction × √K_mid |
|---|---|---|---|---|
| [200, 400) | 100 | 6 | 0.0600 | 1.01 |
| [400, 800) | 200 | 13 | 0.0650 | 1.55 |
| [800, 1600) | 400 | 13 | 0.0325 | 1.09 |
| [1600, 3200) | 800 | 17 | 0.0213 | 1.01 |
| [3200, 6400] | 1600 | 27 | 0.0169 | 1.14 |
| [6400, 12800] | 3200 | 30 | 0.0094 | 0.89 |

The local decay exponent over the last three doublings is about 0.6: 1/2 within counting noise.

**General residual sources** `n = 2^K t − 1`, t of B random bits, 300 per row:

| K | B = 25 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|---|
| 9 | 0.557 | 0.417 | 0.340 | 0.217 | 0.157 | 0.087 |
| 17 | 0.427 | 0.327 | 0.203 | 0.170 | 0.087 | 0.070 |

* Scaled by √(log₂ n): about 3 (K = 9) and about 2.2 (K = 17). Fewer children give a larger constant.
* The reset pair alone (D ≤ 2) merges in 0.30 → 0.73 (K = 9) and 0.37 → 0.76 (K = 17) of cases.

**Which K are orphans.**
* Orphans have *longer* orbits. Split odd K into thirds by normalised post-run length `(σ_T − K)/K`.
  * For K in [1000, 6400] the orphan rates are 0.56%, 1.1% and 4.2%.
  * For K in [1000, 12800] they are 0.31%, 0.92% and 3.0%, a 10× enrichment in the long third.
* An anomalously long orbit is an excess of odd steps in the source's tail, i.e. a positive debt drift in aggregate. It is *not* a single long ones-run: orphans' longest post-run ones-runs equal those of non-orphans (mean 14.2 vs 14.5, Mann–Whitney z = −0.6).
* There is no leading-digit effect: the orphan rate is 1.9–2.3% in each quarter of frac(K log₂ 3).

**Full-orbit success per shell** (m = v2(K−1), odd K in [600, 6400]):
* the reset pair (D ∈ {1,2}) succeeds in 82–91% of cases, and h_3/h_4 in 82–88%;
* both rates are flat in m;
* both are consistent with extrapolating the Haar law, `1 − ρ(s) ≈ 15 s^(−1/2)`, to the available time `s ≈ 2·10^4`.

**Rescue by 3-adic predecessors.**
* A node m < n in the Terras backward tree certifies n.
* Within backward depth 36, this rescues 37.7% of the 53 orphans with K in [1000, 6400], against a baseline of 40% for all odd K. The rescue looks independent of orphan status. It is mostly depth 3: K ≡ 5 mod 6, via (8n−5)/9.
* Doubly uncovered: 33 of 2700 odd K (1.2%), for instance 1081, 1113, 1137, 1329, 1335, ….

## 4. The reducible locus (THM-4604): what friezes and cluster algebras can and cannot see

**Collatz words are reducible.**
* The steps are upper-triangular matrices, `G_w = [[3^|w|, B_w], [0, 2^(Σw)]]`.
* The carry is a twisted 1-cocycle, `B_uv = 3^|v| B_u + 2^(Σu) B_v`, and it is positive.
* The cycle point c_w trivialises B on ⟨w⟩; this is the anchor of THM-4600.
* Centralisers are tori, which are the anchored states.

**Traces only hear the tuning.**
* After SL2 normalisation, `t(w) = 2 cosh(δ_w/2)` with `δ_w = |w| ln 3 − Σw ln 2`. This is the same tuning error as the musical spectrum of HYP-9230.
* Every pair of words lies on the Cayley cubic `x² + y² + z² − xyz = 4`, where `tr[A,B] = 2`: the reducible characters. For example, the odd/even pair gives `(5/√6, 3/√2, 7/√12)`.
* Markov triples and Conway–Coxeter friezes live at `tr[A,B] = −2` and on SL2(Z) quiddities. No Collatz word closes a frieze: its diagonal is `(3^|w|, 2^(Σw))`.

**Reading (DICTIONARY).**
* All cluster and frieze coordinates (Fricke traces, Markov mutations, frieze entries) are trace or minor functions. On the Collatz locus they collapse to slopes and tuning errors.
* Merges are identities of the extension cocycle, which traces cannot detect.
* This is the structural reason why the frieze program yields coordinates (positive minors of carries are minors of the cocycle data) but not rewrites (the ear-move obstruction).
* The Markov-type mutation `(x, y, z) ↦ (x, z, xz − y)` still acts on the Cayley cubic. On Christoffel bases it walks the Farey tree of slopes m/L with values `2cosh(δ/2)`: the musical data again.

## 5. HYP-9240 and digits

The +1 barrier's limit chains depend only on the 2-adic expansion of 3^(−D). Real-size bounds give only `N_D + 1 ≤ (3/2)^D`. So the remaining all-D statement is a base-2 digit problem for powers of 3, of zeroless-power type (cf. THM-4580). D ≤ 8000 is settled exactly.

## 6. Next directions

1. A debt-payment mechanism for orphans. Their excess-odd-density tails cause positive drift. Candidates are negative-drift runs (§2) and 3-adic predecessors. The doubly uncovered 1.2% need a new idea, or descent.
2. Prove the orphan law's exponent 1/2 from THM-4581's rate together with a large-deviation bound on the available time.
3. Complete ladder grammars at −1 (the tiling compiler) and at the short-J states. Then give the exact Haar coverage per depth for the whole residual branch.
4. On the cocycle side, look for a cluster-like structure (exchange relations among carries of overlapping words) that does see merges. The Plücker relations of THM-4604's minors are the first instance.

## 7. Reproduction

All scripts are in `04-computation/experiments/runcompress_20261007/`:
* `ladder_and_drift.py`, `drift_check_actual.py`, `fricke_check.py`: ALL CHECKS PASSED.
* `transitions.py` (+ `.out`).
* `mersenne_sigma.py`, `mersenne_sigma_range.py` (data in `mersenne_sigma_12800.txt`).
* `orphan_law_mersenne.py`, `orphan_law_general.py` (+ `.out`), `mersenne_partners_by_shell.py`, `orphan_structure.py`, `orphan_time.py`, `orphan_rescue.py`.
