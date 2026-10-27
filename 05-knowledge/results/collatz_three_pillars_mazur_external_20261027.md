# Collatz strategies through the three pillars (integral coordinates, the preserved operation, positivity): Mazur's positive-density mechanism is the one atlas row carrying all three; the 2026 trusted-advances ledger completed with the three ProofAtlas Collatz theorems (T2); where the theorem sits (the typical-word constant `3/ln(4/3)`, the no-parameter renewal scale law holding to `5%` at `2^28`, an elementary almost-all floor `3.307 ln n`); the mechanism is sign-blind

**Session:** opus-2026-10-27-S64 (worktree `codex/session-fable-2026-10-27-collatz`), 2026-10-27.
**Owner's directive:** "The modulo-two phenomenon is much broader: for any odd-order symmetry action on `n` gaps,
the cancellation polynomial becomes `(1+t)^n` modulo two. That does not mean `n` independent toggles suddenly exist.
I'm organizing the proposed theory around three separate things: integral coordinates, operations that must be
preserved, and positivity needed to recover genuine geometric or combinatorial models. Look for ways to improve
Collatz strategies. Consider recent results we have pulled in from major AI labs and other smaller sources that have
provided advances we can trust on the major mathematical conjectures, all of the Millennium Prize problems, plus any
more of similar ilk." The attached PDF is Mazur's version 2.1 of 2026-09-06, digested on 2026-09-28
([digest](mazur_positive_density_20260928.md), audited, MISTAKE-548); it is re-read here through the three pillars.

**Concurrent work read before writing (same day, same directive):**
[kind-pasteur-2026-10-27-S6](mod2_collapse_collatz_strategy_audit_20261027.md) (11:41: the mod-2 collapse lemma in
session 65's pair-orbit language, with `n = #pairs - #orbits`; the
[trust-tiered ledger](../reference/TRUSTED-ADVANCES-2026-10-27.md); six strategy shapes);
[mac-mini-2026-10-27-potentials](dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md) (THM-4813,
the numerator coordinate `X = 3s+1` with positive carries) and its
[reflection II](../../07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md)
(potentials, symmetry budgets, positive coordinates and coboundaries: "the boundary term is open");
[klein-2026-10-27](external_preprints_mub6_local_smoothing_eth_merge_20261027.md) (external preprints as typed hooks);
[session 63](session63_integral_observers_root_arrangements_20261025/SYNTHESIS.md) (the owner's three parts:
integral carrier, operation bank, realizability cone; Theorem 1(3) of its
[cancellation result](session63_integral_observers_root_arrangements_20261025/cancellation/RESULT.md): for every
odd-order action on `n` gaps, `P(t) = (1+t)^n (mod 2)`, `d_2 = n`, and the characteristic-two depth produces no legal toggle); [THM-4764](../../01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md)
(klein: the depths at every prime). This note adds what those do not contain and credits them where they overlap;
it re-proves nothing of theirs.

**Status.** SYNTHESIS (the pillar typing of the strategy atlas, section 3) + PROVED (Proposition F, the elementary
almost-all floor, section 5.4; the trunk inequality, section 5.7; the Terras null case, section 2, a restatement of
classical facts) + FINITE-EXACT (every `tau(n)`, `n < 2^28`, exact integers; P5) + VERIFIED (the
pre-registered P1 and P2 passed: the no-parameter renewal model agrees to within `5%` in the fluctuation scale;
the pre-registered HYP-9425 (4), run after P3/P4 failed, passed: the independent word model reproduces the 36
quarter-dyadic floor counts) + OBSERVED (the pre-registered P3 and P4 FAILED on a period-two oscillation in the dyadic
index, which the word model then explained; the block distribution of the normalised deviation is right-skewed, so the
Gaussian shape is not verified; the measured growth exponent `0.70 +- 0.06` of the floor set against the model's
`0.738` in the limit and `0.63-0.68` at finite scales) + HEURISTIC
(the independent-3-adic exponent `E(C)`, the Lagarias-Weiss renewal law) + CITED at tier T2 (Mazur's three
theorems and the companion: Lean-kernel-checked with the standard axioms, statement read here and matching, no
independent replay, platform review pending) + CONJECTURE (the strong form of Mazur's theorem; the sign-blindness
transfer to `3x-1`) + HYP-9425. Collatz OPEN. Independently audited (section 10b: SOUND WITH CORRECTIONS, eleven corrections applied;
MISTAKE-685).

Scripts: `04-computation/experiments/collatz_fast_convergence_density_20261027/fast_convergence_density_20261027.py`
(P1-P5, pre-registered in the docstring and pushed as `4232963739` before the run; output beside it),
`sawtooth_probe_20261027.py` (post hoc), `fast_set_exponent_model_20261027.py` (the heuristic exponent),
`word_model_fine_resolution_20261027.py` (HYP-9425 (4), pre-registered in `9c744f56e2`, PASSED).

---

## 0. The answer in one screen

1. **The three pillars are proof inputs; the atlas's walls are conclusion controls.** A proof *uses* an integral
   coordinate (an actual integer and its exact orbit, not a residue class or a Haar-typical point), the operation
   (the specific map, through a transfer identity, a drift or exact carries) and positivity (nonnegative weights,
   one-sided carries, `n > 0`); the barrier atlas types what a *conclusion* survives (INTEGRAL: a rational cycle
   point; SHEET: `n -> -n`; DRIFT/UNIFORM: another multiplier). A proof that uses a pillar symmetrically under the
   corresponding hostile yields a conclusion blind to that wall: Mazur's proof uses positivity (one-sided carries,
   nonnegative weights) in a sign-symmetric way, so its conclusion is SHEET-blind, which is the atlas's typing.
   The owner's "modulo-two phenomenon" has an exact Collatz instance that is not even
   modular: Terras's theorem makes every parity word of length `k` the word of exactly one residue class mod
   `2^k`, so at the residue level all `2^k` toggles exist and the enumerator is `(1+t)^k` on the nose. The three
   pillars are precisely the three data the residue forgets: the integer, its sign, and the multiplier (section 2).
   Session 63's `P(t) = (1+t)^n (mod 2)` for odd actions on `n` gaps, THM-4764's depths, and kind-pasteur's
   `n = #pairs - #orbits` are three exact readings of the same sentence; they do not conflict.
2. **Mazur's argument is the one row of the strategy atlas whose proof uses all three pillars at every scale.**
   Integral coordinate: only actual odd positive inverse orbits are admitted ("a congruence alone is not a
   substitute", his section 2) and the count is reached through the source charge `x omega(w) <= M`, which is the
   one-sided sign of the carries (the Steiner-Eliahou product identity `1/x = 3^d 2^(-A) prod(1 + 1/(3x_j))`, digest
   Theorem C (c); THM-4813 (A) is its cocycle reading in the numerator coordinate).
   Preserved operation: the weighted transfer's exact mean identity (his Lemma 2.1) and the finite 3-adic mixing
   (2.3). Positivity: nonnegative weights, indicators dropped only inside upper bounds (his Proposition 5.3), a
   seed chosen because a maximum dominates a mean (Proposition 4.2). Every earlier row leaves one pillar unused and stops
   at the wall that pillar guards (section 3). This is a reading, typed SYNTHESIS; it explains why the method
   crosses from `x^0.84` to `cX` and what it cannot do (it averages over seeds, so it is blind to any single orbit).
3. **The ledger.** The concurrent trust-tiered ledger did not list Mazur, and said Tao's theorem "remains the only
   general structural theorem". Added there today (section 4): the three ProofAtlas theorems and the companion at
   tier T2 with the platform's own evidence record (axioms `propext`, `Classical.choice`, `Quot.sound`; no sorry;
   no independent replay; review pending), the erratum Mazur's pipeline found in Tao's Corollary 6.3 and Tao's fix
   (T1-level, acknowledged), two downstream AI-generated Lean repositories (T3), the documented Lean-kernel exploit
   by an AI "disproof" of Collatz (a risk, not a finding against Mazur), and the credibility anchor, Tao's digestion
   of Mazur's Sendov proof. Also added: AlphaEvolve (PNAS 2026, T1), AlphaProof Nexus, the Erdős unit-distance
   disproof's companion, the sunflower preprint, the two 2026 Fields medals as acceptance signals.
4. **Where the theorem sits.** Mazur's `523/50 = 10.46` is `0.3%` above the typical-word constant
   `C* = 3/ln(4/3) = 10.4282`. For all `n < 2^28` the density of `{tau(n) <= 10.46 ln n}` is `0.54` at every
   dyadic scale from `2^16` on (`0.555 -> 0.542`), and the fluctuation law of the Lagarias-Weiss renewal model with
   **no free parameter** (`tau ~ C* ln n + 16.42 sqrt(ln n) Z`) holds to `5%`: the block standard deviation of
   the normalised deviation is `0.949 -> 0.963` (`k = 20..28`, rising) against `1`, the mean `-0.038 -> -0.034`.
   So the theorem's set is the lower half of the population; under the model its density tends to `1/2`, not to
   `1`, and the only constant that can improve is `c` (from the reciprocal of a four-fold tower to `1/2`), not
   `C`. Unconditionally (Proposition F, elementary): below `3.307 ln n` the fast set has counting function
   `O(X^(1-c))`; the independent-3-adic model puts its exponent at `E(3.307) = 0.738` and the measured two-block
   exponent is `0.70 +- 0.06`; `E(C) = 1` exactly at `C = C*` (section 5).
5. **Two pre-registered predictions failed** (monotone sharpening in the dyadic index; floor-count block ratios
   below `2`): the dyadic sequences oscillate with period two in `k`, and at quarter-dyadic resolution the floor
   set's density varies by factors up to `2.8` inside one block. The structure is the carry geometry of the few
   admissible word shapes `(k, d)`, not a sawtooth: the test pre-registered as HYP-9425 (4) and run afterwards
   PASSED (the independent word model reproduces the 36 quarter-dyadic floor counts with log-correlation `0.992`,
   sign agreement `24/27`, total ratio `1.098`), so no 3-adic correlation is visible at this resolution.
6. **The mechanism is sign-blind.** The `3x-1` seeds `(2^(2j+1)+1)/3` permute `Z/3^q` exactly as Mazur's
   `(4^j-1)/3` do (P5, FINITE-EXACT `q <= 6`, PROVED by the same lifting-the-exponent identity), the harmonic-mass
   identity is sign-symmetric (digest Theorem A), and the time-constant argument uses only `S(y) <= 2y`. Conjecture:
   Mazur's Theorem 1.1 and companion hold for `3x-1` with the size inequality reversed (the carries change sign;
   the loss is a factor `1 + O(d/16^(b_0))` on the selected histories), giving the three `3x-1` cycles
   positive-density log-time basins, which no human method has reached either (section 5.7, HYP-9425 (3)).
7. **The trunk inequality** `max_(j >= 2) H_n(R_j) >= 4 H_(n+1)(1) - 3 H_n(1)`, `R_j = (4^j-1)/3` (PROVED, one
   line from the digest's parents identity): Mazur's seed averaging, restricted to the trunk, is the seed-`1` test
   of the digest in another coordinate, and a moving seed is exactly what his source-charge lemma forbids (5.7).
8. **Improvement moves** M1-M6 (section 6), each typed with map, preserved predicate, loss, sidecar and cheapest
   test, cross-referenced to kind-pasteur's six shapes and mac-mini's reflection; and one hygiene finding: the same
   external claim is typed differently across sessions (`#158`: CONDITIONAL, CITED, accepted premise; `#003`:
   "PROVED modulo", CONDITIONAL, PREPRINT CLAIM), so the ledger now carries the rule "one tier per source" (6, M6).

---

## 1. Inheritance pass and concept board

**Closest proved mechanism.** The digest's Theorem A (`3^n mu_n(y mod 3^n) = sum_(S^n x = y) 3^n 2^(-A_n(x))`,
the harmonic mass of a Syracuse tree layer, PROVED, audited) together with THM-4813 (A) (the numerator coordinate
`X = 3x+1`, in which the Syracuse step is `X -> 3X/2^v + 1` with a positive carry, so that the forcing cocycle is
the drift cocycle plus the coboundary of `log X`). Mazur's source charge is the same identity read as an inequality.
**Canonical hostiles.** `3x-1` (three positive cycles, the sign), `5x+1` (positive drift, the multiplier), and `Z_2`
(every periodic parity word is the word of a rational cycle point, Lagarias 1990 / THM-4715: the integer). **Corrected near
misses.** MISTAKE-548 (the tower was read one level short; a pruning "loss bound" was a calibration; the minimum of
a sequence is read off the whole computed range), MISTAKE-677 and MISTAKE-680 (same-day packets closed a directive
while a session worked; key numbers typed new were canon). **Least-used sidecar.** The parents identity of the
digest (Corollary A2's form `H_(n+1)(1) = 3 sum_(j >= 1) 4^(-j) H_n(R_j)`), used there only to extrapolate.

| live concept | representation | invariant / extremal | operation applied today | lost coordinate | cheapest next test |
|---|---|---|---|---|---|
| the three pillars (I, O, P) | wall typing of the barrier atlas | which hostile breaks a strategy | re-type S1-S13, K-L, Tao, Mazur (section 3) | none (a typing) | a strategy with all three that is not a count |
| Mazur's weighted histories `3^d 2^(-A)` | harmonic mass `H_n(seed)` | `z(M) > 0` for one seed | trunk inequality; sign-blind seeds (5.7) | the single orbit | replay P1-P3 of the digest for `3x-1` |
| the fast set `F_C = {tau <= C ln n}` | dyadic densities, `z`-statistics | `C* = 3/ln(4/3)`; `sigma = 16.42` | the `2^28` run (5.2); the word-model test (5.3) | the fine 2-adic structure | `2^32` by sieve (HYP-9425 (1)-(2)); the skewness of `z` |
| the floor `tau >= 3.307 ln n` a.a. | Terras bijection + Chernoff | exponent `1 - c` (proved), `0.738` (model) | Proposition F (5.4); `E(C)` (5.5) | the 3-adic integrality | the two-block exponent at `2^32` |
| the toggle spaces | residues mod `2^k` (Terras); gaps under odd actions | `(1+t)^k` exact; `(1+t)^n mod 2` | the null-case reading (section 2) | integer, sign, multiplier | none (exact) |
| the pillars as proof inputs vs the walls as conclusion controls | K/D per row | sign-symmetric use => SHEET-blind conclusion | the re-typing (section 3) | — | a proof using all three on one orbit |
| the external ledger | tiers T1-T3 | one tier per source | Mazur rows, AI-lab rows (section 4) | — | re-check T2/T3 before any promotion |

Meta-pattern cards used: "Read cited papers fully before claiming new" (the digest and its audit were re-read
before a single sentence was written; the concurrent packets were read first); "Pre-register the cheap tests"
(P1-P5 pushed before the run; two failed and are reported as failed); "Retype on contact with concurrent work".

---

## 2. The three pillars, exactly, and the null case

**Definitions (for the Collatz lane; pillars are proof inputs).** A proof *keeps* a pillar when it uses it as an
input and *drops* it when it never does. *Integral coordinate* (I): the proof works with an integer `n`, its size
and its exact orbit, not only with a residue class, a 2-adic point, a rational, or a Haar-typical point. *Preserved
operation* (O): the proof uses the specific map (its multiplier and shift) through a transfer identity, a drift or
exact carries, rather than only the free toggle space of parity words, which every `3x+k` and `5x+1` share.
*Positivity* (P): the proof uses `n > 0`, nonnegative weights or one-sided carries, for steps of the kind "a maximum
dominates a mean" and "drop an indicator inside an upper bound". The barrier atlas's walls are *conclusion controls*:
INTEGRAL asks whether the conclusion distinguishes `Z^+` from `Z_(2)`, SHEET whether it changes under `n -> -n`
(which conjugates `3x+1` to `3x-1`, THM-4474 (E)), DRIFT/UNIFORM whether it changes with the multiplier. A proof that
uses a pillar symmetrically under the corresponding hostile (Mazur's one-sided carries work with either sign) gives a
conclusion blind to that wall; a proof that never uses it does so a fortiori. [Corrected after the audit: the first
draft defined the pillars by what the conclusion changes under, which contradicted the typing of rows M and K-L.]

**The null case (PROVED; a restatement of Terras 1976, CITED in the atlas, and of THM-4474 (E)).** Fix `k >= 1` and
the one-division map `T(n) = n/2` (`n` even), `(3n+1)/2` (`n` odd). The parity vector
`v_k(n) = (n mod 2, T(n) mod 2, ..., T^(k-1)(n) mod 2)` depends only on `n mod 2^k`, and `n mod 2^k -> v_k(n)` is a
bijection onto `{0,1}^k`. Hence the enumerator of residue classes mod `2^k` by the number of odd steps is exactly
`(1+t)^k`: every one of the `2^k` toggles is realised by a residue class, with no symmetry quotient and no reduction
modulo two. The same bijection holds, with the same words, for `3x-1` (negation of the parity graph) and for every
`3x+k`, `k` odd, and the words of `5x+1` form the same free space. So the residue level sees none of the three
pillars: not the integer (a class contains infinitely many), not the sign, not the multiplier. *This is the Collatz
instance of the owner's sentence:* "`n` independent toggles exist" at the level of the toggle space, and the three
pillars are exactly the data the toggle space forgets.

**Three exact readings of the modular sentence, no conflict.** (a) Session 63, Theorem 1(3) (PROVED, audited):
for every odd-order action on `n` gaps, each subset orbit has odd size, so `P(t) = (1+t)^n (mod 2)` and `d_2 = n`,
while no invariant vector is added; THM-4764 gives the depth at every prime and the fifty-vertex witness. (b)
kind-pasteur-S6, Lemma 1 (PROVED, 5 checks): in session 65's pair-orbit defect polynomial,
`F_sigma / 2^#orbits = (1+z)^(#pairs - #orbits) (mod 2)`, never binomial integrally. (c) Terras, above: the
unquotiented toggle space is binomial exactly. In (a) and (b) the modular collapse is what remains of the free
toggle space after an odd symmetry quotient; in (c) there is no quotient. In all three the lawful operations are a
strict subspace: invariant subsets (a), nothing inside an orbit (b), and in (c) the integer's own orbit, which is one
point of each class at each level.

**What the pillars say about the two walls.** The coalescence atlas (2026-10-04) found that every proof style sees
either the DRIFT wall (pointwise divergence) or the SHEET/INTEGRAL walls (bounded complexity), never both at
unbounded complexity. In pillar language: a statement about *one* integer's divergence needs I (that integer), O
(its orbit under the actual map) and P (its sign) simultaneously and unaveraged; every density method averages I
away, every residue method drops P, every relaxation with choice drops O. That is the reading, not a theorem.

---

## 3. The strategy atlas re-typed by pillars

Columns: which pillars the method's proof uses (K), never uses (D), or restores through a sidecar (S), and the wall
the atlas records for its conclusion. A proof may use a pillar symmetrically and still have a conclusion blind to the
wall (row M). Rows S1-S13 are the coalescence atlas's (IDs and paths there); the last four rows are added.

| row | method | I | O | P | the wall, in pillar terms |
|---|---|---|---|---|---|
| S1 | density and mixing (Terras, Korec, Tao, GGM) | D (classes; Haar) | K | D (sheet-blind: GGM covers `3N-1`) | the measure-zero set `E_inf`; log density (Tao) |
| S2 | thin divergence (THM-4476/4499) | K (an orbit's own integers) | K | D (holds on both sheets) | the floor `o(X^(h*))`; the bound holds on both sheets (THM-4476's status), so it cannot see the sign |
| S3 | certificates, strategy cube (THM-4474/4479) | D (residues) | K per level | D (signs free; class (i) only forces `sigma(+-1)`) | the window `[0,1]`: the residue sees `0`, `10`, `1` alike |
| S4 | Diophantine cycle exclusion (Steiner, SdW, Hercher; THM-4484/4515/4516/4633) | K (the cycle element) | K (rotation) | K (positive cycle) | bounded complexity: a finite check per parameter |
| S5 | Fourier / spectral (THM-4519/4520; the profile) | D | K | D (signed characters; `E_n` vs `O_n`) | H1: everything conjugation-invariant sees only the `-1` resonance |
| S6 | rigidity (THM-4523) | K | K | D (`3n-1`, `5n+1` rigid too) | component membership is archimedean |
| S7 | conjugacies and carriers | moved into the target | — | — | reformulation |
| S8 | Diophantine criticality | diagnostic | — | — | every slack is a quantity the conjecture controls |
| S9 | potentials, ranks (THM-4482/4483/4507) | K | K | needed (a rank must use the whole integer) | the only known rank is the stopping time |
| S10 | relaxations with choice | K | D (choice) | exposed | removing DIM exposes SHEET |
| S11 | undecidability | n/a | n/a | n/a | binds uniform proof systems |
| S12 | tournaments, clocks, motifs | D | K | D (residue invariants identical for `3x+-1`) | sheet-blind by THM-4474 (E) |
| S13 | refuted external proofs | fail a control | | | the `3n-1` control |
| K-L | Krasikov-Lagarias 2003 (difference inequalities) | K (tree integers) | K to depth `k = 11` only | K (counts) | `x^0.84`: an operation truncated at a fixed depth gives a sublinear exponent |
| M | Mazur 2026 (weighted inverse histories) | K (admitted integral orbits; source charge) | K at every scale (Lemma 2.1 + (2.3)) | K (nonnegative weights, one-sided carries, max over mean; used sign-symmetrically, so the conclusion is SHEET-blind) | `cX` with `c` a reciprocal tower; blind to the single orbit |
| N | the numerator coordinate THM-4813 (A) | K | K (conjugacy) | K (positive carries) | the boundary term (mac-mini II): an identity, not a bound |
| C | the cycle budget THM-4633 | K | K | K | bounded complexity (the cycle half is of lonely-runner type) |

**Reading.** The rows whose proofs use all three pillars (S4, C, M, N) are counting or identity statements: a finite check per
parameter (S4, C), an averaged count (M), an exact identity whose boundary term is open (N). No row has all three
pillars *and* a pointwise unbounded conclusion; that is the pillar form of the atlas's two-wall verdict. The honest
frontier named there (the orbit-coupled local-time statement, THM-4506) is the one candidate that would keep I, O
and P on a single orbit.

---

## 4. Mazur's three theorems in the trusted-advances ledger

Statements are the ones read on the platform by S19 (2026-09-28) and today (2026-10-27, through a bounded web
reading whose sources are listed in the ledger); tiers are the concurrent ledger's T1/T2/T3 scale.

| item | statement | certificate and review | tier |
|---|---|---|---|
| **Positive-density log-time convergence** (v2.1, 2026-09-06; the attached PDF) | explicit rational `c > 0`, integer `X_0`: for `X >= X_0`, `#{0 < n < X : exists k <= (523/50) ln n, step^k(n) = 1} >= cX`, `step` the raw map; also with `262/25` | Lean 4, `CollatzPositiveDensity.lean`, commit `830b9d3f`: `positive_density (X) (hX : cutoff <= X) : densityConstant * X <= countBelow (523/50) X`, `densityConstant_pos`, `cutoff_pos`; the platform's evidence record: axioms `propext, Classical.choice, Quot.sound`, no sorry, no unexpected axioms, no `native_decide`, build passed (1,523 s), toolchain `v4.30.0-rc2`, 386 files / 61,169 lines, collected 2026-09-07. The platform itself asserts no independent replay, no human review ("accountable review not recorded"), publication pending; the manuscript reviews are by AI systems. Statement correspondence with Theorem 1.1: exact (count strictly below `X`, natural log, raw steps). | **T2** |
| **Positive lower density of predecessors** (companion, "accepted formalization" on the platform) | every target `a`, `3 ∤ a`, has a predecessor set of positive lower natural density; corollaries: a single non-convergent integer forces a positive-density non-convergent set; upper-density-one convergence forces universal convergence | theorems `predecessors_positive_lower_density`, `nonconvergence_positive_lower_density_of_counterexample`, `universal_reaches_one_of_arbitrarily_dense_convergence` (388 files / 61,800 lines, same supporting revision); same review status | **T2** |
| **Natural-density almost-bounded orbits in logarithmic time** (v2, 2026-07-21) | for every `f -> infinity`, a natural-density-one set of positive integers has an iterate below `f(n)` within `C_ND ln n` raw steps, `C_ND = 1509503/(5000 ln 2) < 436` (`145 ln N` Syracuse steps for odd starts): Tao's theorem moved from logarithmic to natural density with an explicit clock | Lean (main file `Erdos1135/ND/LogTime/Paper.lean`, commit `ca3dd0d6`), "accepted in ProofAtlas"; axiom profile not shown on the page; four AI reviews. Tao (2026-08-23, blog comment): "there are now some AI-generated solutions to the natural density version of my result" — a generic acknowledgment, not a reading | **T2** |
| **Certified `x^0.90` lower bounds for predecessor sets** (v2, 2026-07-17) | for each `t > 0`, `3 ∤ t`: `pi_t(x) >= x^(9/10)` for large `x`, and `>= c_t x^(901/1000)`; level-18 LP certificate (129,140,163 rows) plus an adaptive-potential growth system | Lean, but two computation certificates are recorded as named `native_decide` axioms (not kernel-clean); four AI reviews; the paper does not claim the exponents are new (it cites an unnamed `0.88` bound). Superseded by the positive-density theorem if that stands | **T2** (with the `native_decide` caveat) |
| **Erratum in Tao's paper found by the same pipeline** (2026-07-05/06) | the printed proof of Tao's Corollary 6.3 has a range error in (6.8) (the displayed sum exceeds `3^n` by about `2^912`); Tao: replace the `1/2` in (6.7) onward by `0.99` and the bounds go through; not a counterexample to the theorem | acknowledged and repaired by the author on his blog | **T1** (as an erratum) |
| Downstream: AxiomMath `CollatzPositiveDensity` (GitHub) | a cutoff-free refinement: for every `N`, at least `floor(N / 2^(2^(2^163)))` of the `n < N` reach `1` within `floor(10.46 ln n)` steps — a triple-tower reciprocal density against Mazur's four-fold tower, with no `X_0` | claims a comparator check; axioms not stated; no independence claimed | **T3** |
| Downstream: N. El Jaouhari, `ConjectureA` (Zenodo, 2026-09-29) | Lean proofs of Applegate-Lagarias Conjecture A and of the growth-exponent conjecture, resting on Mazur's vendored `predecessors_positive_lower_density` | standard axioms reported; "AI-assisted and not independently refereed" | **T3** |
| Risk note | an AI-generated Lean "disproof of Collatz" exploited a Lean kernel bug (reported with a quotation of Leo de Moura on social media, 2026) | — | the repository rule stays: a Lean-checked claim without an independent replay is T2 at best; ProofAtlas says the same of itself |
| Credibility anchor | Mazur's AI-generated proof of Sendov's conjecture (2026-08-05) was digested and re-formalised by Tao (about 15k lines against 90k; blog, 2026-08-12) | expert-digested, not refereed | **T1** for Sendov; it transfers no trust to the Collatz claims, only to the pipeline |

**What this changes and does not change.** (i) In the barrier atlas (section 8 addendum, dated today) the Tao row's
"the natural-density version is called plausible and is open" and the Krasikov-Lagarias row's `x^0.84` are each
annotated with the T2 claim that supersedes them; no status in canon moves (CITED stays CITED; THM-4517 is
compatible: Mazur's `c` decreases with the seed `M`, `c = beta/(32M)` in the qualitative conclusion of his Section 7
and `c = 3/(256 M* 2^m)` in the explicit (8.6) as read by the referee, `M*` a bound on the seed, as THM-4517 requires
of any root-uniform statement). (ii) The concurrent ledger's sentence "Tao (2019, logarithmic density) remains the only general
structural theorem" is amended there with a pointer to the new rows. (iii) The digest's Theorem C (positive density
of `{tau <= C ln n}` for any `C` forces `limsup n^(1/6) H_n(1) > 0`, passed through `n = 48`) remains the cheapest
confrontation with the claim; nothing found today bears on `lim H_n(1)`.

**Other 2025-26 items added to the ledger today** (details and URLs there): AlphaEvolve (PNAS 2026, peer-reviewed;
67 problems; the finite-field Kakeya construction proved by Deep Think and formalised by AlphaProof) — T1;
AlphaProof Nexus (9 of 353 attempted open Erdős problems, Lean-verified per the authors; preprint 2026-05) — T2;
Erdős #90 unit-distance disproof (already T1 in the ledger; Sawin's refinement and the nine-author remarks) — T1;
the Erdős–Szemerédi sunflower preprint of Bucić–He–Lee (2026-10-19; `[n]`-subset variant; AI credited with the
central idea) — T2 (eight days old); BB(5) paper accepted at STOC 2026 — T1; Hong Wang (Kakeya in `R^3`) and Yu
Deng (Hilbert's sixth) awarded 2026 Fields Medals — acceptance signals for T1 entries already present; the
Anthropic FLT formalisation (secondary reports of a comparator-checked artifact; written up here by a Claude model,
which is noted) — T2/T3; OpenAI's 2026-10-06 release of 719 manuscripts and the AHM statement — T3 as a collection,
with the already-typed repository readings unchanged.

---

## 5. Where the theorem sits: the typical-word constant and the run

### 5.1 The constant and the model

A Syracuse step with valuation `a` costs `1 + a` ordinary steps and multiplies by `3/2^a`. Under the Terras model
(`a` geometric, `P(a) = 2^(-a)`): `E[1 + a] = 3`, `E[a ln 2 - ln 3] = 2 ln 2 - ln 3 = 0.28768`, so the typical
ordinary time to descend by `ln n` is `C* ln n`, `C* = 3/(2 ln 2 - ln 3) = 3/ln(4/3) = 10.42818`; Mazur's
`523/50` is `C*` plus his rounding losses (the digest, section 1; his (7.9) is `kappa(L+h) + 1 = 7.25012 <=
(523/50) L = 7.25032` with `kappa = 34881/10000`, audit item 82). The renewal-reward central limit theorem for the
same model gives

```
tau(n) ~ C* ln n + sigma sqrt(ln n) Z,   sigma^2 = Var(c - C* delta)/E[delta] = 2 (1 - C* ln 2)^2/(2 ln 2 - ln 3) = 269.7,
sigma = 16.422,   c = 1 + a,   delta = a ln 2 - ln 3,
```

since `c - C* delta = const + a(1 - C* ln 2)` and `Var(a) = 2`. There is no free parameter. This is the
Lagarias-Weiss stochastic model (CITED; the repository has its delay constant `41.677`, in T-steps, in THM-4504; the
ordinary analogue is about `67`), stated here for the ordinary convention; the prediction for a dyadic block around `N` is the sigmoid
`G(C) = Phi((C - C*) sqrt(ln N)/sigma)`, which sharpens only like `sqrt(ln N)`. So the model says: Mazur's set at
`C = 10.46` is the lower half of all integers, with density `1/2 + o(1)`; for `C > C*` the density tends to `1`
and for `C < C*` to `0`, at the rate of a Gaussian tail in `sqrt(ln N)`.

### 5.2 The run (`fast_convergence_density_20261027.py 28`; all `1 <= n < 2^28`)

Exact 64-bit integers; `tau` by first descent below `n` with memoisation (6.3 s); the largest orbit value met is
`6,404,797,161,121,264 = 2^52.51`; the largest `tau` is `956` at `n = 226,588,897`. Block densities
`G_k(C) = #{n in [2^(k-1), 2^k) : tau(n) <= C ln n}/2^(k-1)` with the model's `Phi((C - C*) sqrt(ln N_k)/sigma)`,
`N_k = 2^(k-1/2)`, and the statistics of `z(n) = (tau(n) - C* ln n)/(sigma sqrt(ln n))`:

| `k` | `G_k(8)` | model | `G_k(10)` | model | `G_k(10.46)` | model | `G_k(13)` | model | `G_k(16)` | model | mean `z` | sd `z` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.3597 | 0.2934 | 0.5054 | 0.4618 | 0.5423 | 0.5028 | 0.7425 | 0.7176 | 0.9020 | 0.8939 | -0.038 | 0.949 |
| 24 | 0.3249 | 0.2753 | 0.5034 | 0.4581 | 0.5408 | 0.5031 | 0.7651 | 0.7363 | 0.9228 | 0.9146 | -0.037 | 0.956 |
| 28 | 0.2988 | 0.2593 | 0.5101 | 0.4547 | 0.5412 | 0.5034 | 0.7796 | 0.7529 | 0.9300 | 0.9307 | -0.034 | 0.963 |

Cumulative densities `F_X(C)`, `X = 2^k`: at Mazur's `C = 523/50`, `0.5550, 0.5513, 0.5465, 0.5430, 0.5426,
0.5392, 0.5408, 0.5400, 0.5404, 0.5418, 0.5399, 0.5429, 0.5420` for `k = 16..28` (the digest's `54.22%` at
`10^6` is the `k = 20` value `0.5426`); the median of `tau(n)/ln n` over `n < 2^k` is `9.665, 9.695, 9.735, 9.765,
9.785, 9.825, 9.835, 9.855, 9.875` for `k = 20..28`, rising towards `C*` at about `0.025` per doubling. The
deficit of the median (`C* - 0.55` at `k = 28`) is not a mean offset: the mean of `z` is `-0.035` (about `2.5`
steps), which would move the median by `0.13` only; the block distribution of `z` is right-skewed (audit: skewness
`+0.60 -> +0.66`, block median of `z` `-0.168 -> -0.146`). So the renewal model is confirmed in scale (sd within `5%`)
but not in shape at these `n`: the Gaussian profile is off by up to `0.066` absolute (`C = 8`, `k = 20`), and the
shape is OBSERVED, not verified. The full tables, `C = 3.307` to `30`, are in the output file (densities are read
off a histogram of `tau/ln n` with bins of width `0.01`, exact for the integer, `10.46` and `10.48` columns; the
`3.307` column reads `r < 3.31` (`4718` against the exact `4699` at `k = 20`), the floor counts are exact, medians
have resolution `0.01`).

**Pre-registered outcomes (thresholds fixed in commit `4232963739` before the run).**

| prediction | outcome | numbers |
|---|---|---|
| P1: block sd of `z` in `[0.85, 1.15]` and mean within `0.35` for `k = 20..28` | **PASS** | sd `0.949, 0.951, 0.953, 0.954, 0.956, 0.958, 0.960, 0.961, 0.963`; mean `-0.038 .. -0.034` |
| P2: `F_X(523/50)` in `[0.45, 0.60]` for `k = 16..28`; median within `1.0` of `C*` for `k >= 20` | **PASS** | `0.555 .. 0.542`; medians `9.665 .. 9.875` (`C* - 0.55` at `k = 28`) |
| P3: `G_k(8)` non-increasing and `G_k(13)` non-decreasing in `k = 20..28`; `G_28(8) <= 0.33`, `G_28(13) >= 0.70` | **FAIL** (the monotonicity clauses) | `G(8)` is non-increasing (`0.3597 -> 0.2988`) and the two end values hold (`0.2988`, `0.7796`), but `G(13)` is `0.7425, 0.7481, 0.7481, 0.7720, 0.7651, 0.7683, 0.7817, 0.7860, 0.7796`: two decreases, both at even `k` |
| P4: floor count `#{n in [2^(k-1), 2^k) : tau < 3.307 ln n}` has block ratio `< 2` for `k = 20..28` | **FAIL** | counts `3489, 4699, 8153, 12233, 24130, 32643, 67829, 91956, 167377, 239065` (`k = 19..28`); ratios `1.347, 1.735, 1.500, 1.973, 1.353, 2.078, 1.356, 1.820, 1.428`: alternating, `2.078` at `k = 25` |
| P5: `3x-1` seeds `(2^(2j+1)+1)/3` permute `Z/3^q`, `q <= 6`; `v_3(2^(2j+1)+1) = 1 + v_3(2j+1)` | **PASS** | `3, 9, 27, 81, 243, 729` classes hit |

The two failures are reported as failures. Both come from a period-two oscillation in the dyadic index that the
predictions did not anticipate; the repaired statements below are post hoc and typed OBSERVED, to be pre-registered
at the next scale (HYP-9425).

### 5.3 The sawtooth that is a lattice (post hoc, `sawtooth_probe_20261027.py`)

At quarter-dyadic resolution the floor set's density (per `10^4`) in the four quarters of `[2^27, 2^28)` is
`24.1, 19.8, 8.5, 18.8`; in `[2^26, 2^27)` it is `30.4, 28.7, 13.2, 27.5`; in `[2^24, 2^25)` it is `32.8, 40.7,
31.1, 57.1`: factors up to `2.85` inside one block, with a pattern that drifts across blocks rather than repeating
with a fixed period. Blocks of ratio `4` have floor-count ratios `2.77, 2.49, 2.79, 2.81, 2.54` (exponents
`log_4`: `0.736, 0.658, 0.739, 0.746, 0.673`), blocks of ratio `3` have `1.95, 2.04, 2.35, 2.27, 1.94, 2.40`. The
two-block dyadic exponents `log_4(count_k/count_(k-2))` are `0.612, 0.690, 0.783, 0.708, 0.746, 0.747, 0.652,
0.689` for `k = 21..28`, mean `0.703`. **Reading, tested (HYP-9425 (4): test pre-registered in commit `9c744f56e2`, then run: PASSED).** An integer of
the floor set with T-word of length `k` and `d` odd steps is `n = (2^k - c_w)/3^d`, `c_w = sum_(odd j) 3^(d_j') 2^j`
(`d_j'` the number of odd steps after `j`), so its position inside a block is fixed by the carry of its word, and
the admissible shapes at this `C` are few (the model's optimum has `d = 0.72 ln X`, `k = 2.59 ln X`, section 5.5).
The independent-3-adic word model (`word_model_fine_resolution_20261027.py`: 169 shapes `(k, d)`, 20,000 random words
of each shape ending with an even step, weight `C(k-1, d) 3^(-d)`, the floor condition applied to the model's own `n`)
reproduces the 36 quarter-dyadic floor counts with Pearson correlation `0.992` between `log(model)` and
`log(observed)`, the sign of `24` of the `27` consecutive-quarter changes, and a total of `711,515` against `648,085`
observed (ratio `1.098`; the model ignores first arrival at `1`, which overcounts). So the oscillation is the carry
geometry of the admissible word shapes, and the integrality events behave as independent `3^(-d)` events to within
`10%` in total and a factor `1.56` per cell at this resolution: no 3-adic correlation is visible. (A first draft of
this test, a proximity test to the lattice `{k - d log_2 3}`, was withdrawn before any run because that lattice is
dense at this scale.) The same structure is the repository's "expense is Diophantine" ledger (deadlines at the
semiconvergents of `log_2 3`, S9 of the receipts thread), seen from the fast side instead of the slow side.

### 5.4 Proposition F: the elementary almost-all floor (PROVED)

**Proposition F.** Let `tau` be the ordinary total stopping time and `C_F = (1 + (1 + log_2 3)/2)/ln 2 =
2.29248/ln 2 = 3.3074`. For every `C < C_F` there is `c(C) > 0` with
`#{1 <= n < X : tau(n) <= C ln n} = O(X^(1 - c(C)))`. In particular `{tau(n) <= C ln n}` has natural density zero
for every `C < 3.307`, and `tau(n) >= (C_F - epsilon) ln n` for all but `O(X^(1-c))` of the `n < X`.

*Proof.* Work in the one-division convention `T(n) = (3n+1)/2` (`n` odd), `n/2` (`n` even); an odd `T`-step is two
ordinary steps (`3n+1`, then one halving), an even one is one. Fix `X` and count `n in [X/2, X)`. Put
`k_0 = floor(log_2(X/2))`, so `2^(k_0) <= X/2 < 2^(k_0 + 1)`; the interval contains at most two integers from each
class mod `2^(k_0)`. By Terras's bijection the parity vector `v_(k_0)(n)` runs over `{0,1}^(k_0)` with each word taken
by at most two `n` of the interval. Let `d_0` be the number of odd entries of `v_(k_0)(n)`. For `0 < epsilon < 1/2`,
the words with `d_0 <= (1/2 - epsilon) k_0` number at most `2^(h_2(1/2 - epsilon) k_0)` (`h_2` the binary entropy),
so at most `2 X^(h_2(1/2 - epsilon))` of the `n` in the interval have such a word. For every other `n`:
`m_0 := T^(k_0)(n) = (3^(d_0) n + c)/2^(k_0)` with `c >= 0`, hence `m_0 >= 3^(d_0) n/2^(k_0) >= 3^(d_0)` (as `n >= X/2
>= 2^(k_0)`), and the orbit reaches `m_0` after exactly `k_0 + d_0` ordinary steps, without having reached `1` earlier (`T^j(n) >=
n/2^j >= 2^(k_0 - j) > 1` for `j < k_0`, since every step at least halves). From `m_0` to `1` the orbit
needs at least `log_2 m_0` halvings, because along any orbit segment from `m` to `1` with `A` halvings and `d`
odd steps `2^A = 3^d m + c' >= m`. Therefore
`tau(n) >= k_0 + d_0 + log_2 m_0 >= k_0 (1 + (1/2 - epsilon)(1 + log_2 3)) >= (log_2 n - 2)(C_F ln 2 - 2.585 epsilon)`,
using `k_0 >= log_2 X - 2 >= log_2 n - 2`. Given `C < C_F`, choose `epsilon` with `C < C_F - 2.585 epsilon/ln 2`;
then `tau(n) <= C ln n` fails for all large `n` outside the exceptional words. Summing over the dyadic intervals
below `X` gives `O(X^(h_2(1/2 - epsilon)))`, and `c(C) = 1 - h_2(1/2 - epsilon) > 0`. ∎

**Scope and prior art.** The ingredients are Terras's bijection and the entropy bound; the floor `C_F` is one stage
of the argument (the first `log_2 n` steps only). Iterating would need the distribution of the later parity digits of
`n`, which is exactly the 3-adic mixing that Tao's and Mazur's analyses supply; so the gap between `3.307` and `C*`
is the gap between Terras and Tao. I did not locate this statement in the repository (the trivial bound
`tau(n) >= log_2 n`, with equality at the powers of two, is the only all-`n` floor there) nor in the sources the
atlas cites; it is elementary and may well be folklore, and is typed PROVED with that flag rather than filed as a
theorem.

### 5.5 The independent-3-adic exponent `E(C)` (HEURISTIC, `fast_set_exponent_model_20261027.py`)

If the integrality `3^d | 2^k - c_w` of a word with `d` odd steps is modelled as an independent event of probability
`3^(-d)` (the Lagarias-Weiss / Applegate-Lagarias tree model), the count of `n in [X/2, X)` with `tau(n) <= C ln n`
is `X^(E(C) + o(1))` with

```
E(C) = max { kappa h(delta/kappa) - delta ln 3 : delta = (kappa ln 2 - 1)/ln 3 >= 0, kappa + delta <= C },
```

`h` the entropy in nats (`kappa = k/ln X`, `delta = d/ln X`, and `k ln 2 - d ln 3 = ln X + O(1)` for `n` at scale
`X`). The unconstrained maximum is `1` exactly, at `delta/kappa = 1/2`, `kappa = 6.952`, `tau/ln X = 10.428 = C*`:
the model puts the whole population at the typical word, and its exponent is `E(C) < 1` for `C < C*`, `1` for
`C >= C*`. Values: `E(2) = 0.421`, `E(2.5) = 0.584`, `E(3) = 0.690`, `E(3.307) = 0.738`, `E(4) = 0.820`, `E(6) =
0.940`, `E(8) = 0.986`, `E(9) = 0.996`, `E(10) = 0.9996`, `E(C*) = 1`. The measured two-block exponent at `C =
3.307` is `0.70 +- 0.06` (section 5.3); the audit's finite-`X` evaluation of the same model sum gives two-block
exponents `0.63` at `2^16`, `0.68` at `2^40` and `0.716` at `2^140`, so the fair comparison is the finite-scale value
`0.68`, not the limit `0.738`; at `C = 8` the model
says `X^(-0.014)`, invisible at `2^28` next to the Gaussian finite-size factor. Proposition F is the unconditional
shadow of this curve (`1 - c` instead of `E(C)`), and `E(C*) = 1` is the model's own statement that no positive
density is possible below `C*`: Mazur's constant is sharp in the model, and `523/50` is `0.3%` above the edge.

### 5.6 What this says about improving Mazur's theorem

* **`C` cannot be lowered** below `C*` without contradicting the model (and below `3.307` unconditionally). The
  set is the lower half of the population; the theorem's natural strong form is "density `1/2 + o(1)` at `C*`,
  `1 - o(1)` for every `C > C*`" (CONJECTURE, HYP-9425 (1)); the second half is weaker than the conjecture (it
  leaves a density-zero exceptional set) and would need Tao's machinery in natural density with the fixed target
  `1`, i.e. the companion theorem's mechanism pushed to density one, which no one claims.
* **`c` is the only improvable constant**: from a four-fold-tower reciprocal (digest, audited) to the model's `1/2`.
  The downstream AxiomMath refinement claims a triple-tower reciprocal with no cutoff (T3). The gap is the explicit
  mixing coefficient of Mazur's Lemma 8.1 (a three-fold tower); any route to a usable `c` must replace it, and the
  cheapest handle the repository owns is the exact rational 3-adic law to level 7 and the FFT law to level 18
  (digest), on which the mixing distances are `0.85 -> 0.35` for `m = 1..9` at `N = 18`.
* **The cheapest confrontation with the claim stays the seed-`1` test** (digest Theorem C, `limsup n^(1/6)
  H_n(1) > 0` forced; `H_n(1)` at `0.51-0.54` for `n = 41..48`, lower bounds), now with the trunk inequality below.

### 5.7 The trunk inequality and sign-blindness

**Trunk inequality (PROVED).** With `H_n(y) = 3^n mu_n(y mod 3^n)` (digest Theorem A) and `R_j = (4^j - 1)/3`,
the depth-one predecessors of `1` are the `R_j`, `j >= 1`, with weights `3 · 2^(-2j) = 3 · 4^(-j)` summing to `1`,
and `R_1 = 1`. The digest's parents identity `H_(n+1)(1) = 3 sum_(j >= 1) 4^(-j) H_n(R_j)` is therefore a convex
combination, so

```
max_(j >= 2) H_n(R_j)  >=  4 H_(n+1)(1) - 3 H_n(1)    for every n,
```

(the weights `3 · 4^(-j)`, `j >= 2`, sum to `1/4`, and a maximum dominates a weighted mean; the first draft's
weaker second step was a sign slip, corrected by the audit). At level `18` the digest's values
(`H_18(5) = 0.42, H_18(85) = 0.34, H_18(341) = 1.74, H_18(5461) = 0.07, H_18(21845) = 1.50`, i.e. `3/2` times its
`rho_18` list) satisfy it with room (`4 · 0.374 - 3 · 0.410 = 0.27` with the digest's lower bound for `H_19(1)`,
against `1.74`). *Reading.* Mazur chooses his seed among the `R_j` by a residue average
("some unit residue has value at least the mean") and then keeps it fixed; the inequality says the trunk alone
always contains a seed at least as good as the seed-`1` increment. What forbids a seed that moves with the
generation is his Lemma 6.1 (source charge), which needs one fixed `M` so that "the canonical orbit of a source
cannot reach `M` twice". So a Mazur-type proof of `liminf H_n(1) > 0` would have to keep the seed `1` and lose the
residue average; this is the exact obstacle, recorded, not solved.

**Sign-blindness (P5 PASS; CONJECTURE for the transfer).** For the `3x-1` map on positive integers (equivalently
`3x+1` on the negatives) the seeds `r_j = (2^(2j+1) + 1)/3` are odd, `3 r_j - 1 = 2^(2j+1)`, and
`v_3(2^(2j+1) + 1) = 1 + v_3(2j+1)` (lifting the exponent), so `r_(j+t) - r_j = 2^(2j+1)(4^t - 1)/3` is divisible by
`3^q` iff `3^q | t`: the `r_j` permute `Z/3^q` exactly as Mazur's `R_j` do (his Lemma 4.1). The harmonic-mass
identity holds for seeds of either sign (digest Theorem A) and the transfer-operator layer is residue arithmetic;
what changes is the sign of the carries: for `3x-1` the inverse step is `x_i = (2^a x_(i+1) + 1)/3`, so
`3^d x = 2^A M + Sigma'` with `Sigma' > 0` and the size identity reads `1/x = 3^d 2^(-A) prod(1 - 1/(3x_j))`, which
reverses the direction of Mazur's source charge `x omega(w) <= M` into `x omega(w) <= M prod(1 - 1/(3x_j))^(-1)`; on
his selected histories every intermediate value is at least `16^(b_0)`, so the extra factor is `1 + O(d/16^(b_0))`
(expected; not checked against every use of the inequality in the paper). `S(y) <= 2y` holds for both maps (audit).
Prediction: Mazur's Theorem 1.1 and the companion hold for `3x-1` with these seeds and the reversed size inequality,
so each of the three `3x-1` cycles (`{1}`, `{5,7,10}`, the eleven-cycle of `17`)
has a basin of positive lower density in logarithmic time; THM-4517 (3) measured those basins at `0.327 / 0.325 /
0.348` to `2^29`. Cheapest test: replay the digest's P1-P3 for `3x-1` (seeds, inverse histories, Lemma 2.1), then
read the formal development for its sign-dependent lemmas. This is the atlas's SHEET-B typing made concrete: the
mechanism cannot tell the sheets apart, which is also why it cannot see a single orbit.

---

## 6. Improvement moves, typed

Each: source → target; map; preserved predicate; lost information; sidecar; cheapest decisive test. Nothing below is
a result.

* **M1. Positivity on a finite compression is already realised for Collatz, and its price is a tower.**
  kind-pasteur's shape 1 proposes certifying a second-eigenvalue gap of the finite-section transfer operator mod
  `3^J` or `2^K`. Mazur's (2.3) with his Lemma 8.1 *is* that statement for the 3-adic Syracuse operator, with an
  explicit coefficient of three-fold-tower size, and the repository already has an exact dimension-free gap `1/2`
  for the *loosened* walk mod `3^k` (five-papers synthesis, Proposition 3.1, "audit owed"). Map: the Syracuse
  transfer on `Z/3^k` → its loosened relative; preserved: the Gauss-twisted circulant structure (THM-4520); lost:
  the geometric valuation weights. Test (one afternoon): the second singular value of the exact level-`k` Syracuse
  transfer (the rational law to level 7) against the loosened walk's `1/2`; if the gap decays in `k` the shape is
  dead for the actual operator and the note should say so, if it is bounded the explicit mixing coefficient can be
  replaced at low levels.
* **M2. Keep the seed `1`, lose the residue average.** The trunk inequality (5.7) and the obstacle it names. Map:
  Mazur's fixed seed → the trunk; preserved: the harmonic-mass bookkeeping; lost: the source-charge uniformity;
  sidecar: a multiplicity bound for histories that revisit the trunk (they cannot: an orbit meets the trunk once,
  THM-4517, so the sidecar is a per-trunk-point source charge `x omega(w) <= R_j` with the weights `3 · 4^(-j)`).
  Test: whether `sum_(j >= 2) 3 · 4^(-j) H_n(R_j)` stays bounded below along `n = 19..48` using the digest's depth
  decomposition (its existing outputs contain `H_18` on the depth-`d` layers of `1`, which are the layers of the
  `R_j`).
* **M3. The sign-blind transfer** (5.7): a theorem about `3x-1` for free, if the formal development has no
  sign-dependent lemma. Test: the P1-P3 replay, then the Lean file list. Typed CONJECTURE.
* **M4. Lean-first for the elementary layer** (kind-pasteur's shape 2, agreed): the digest's Theorem A
  (one induction) and Theorem C (partial summation plus the size identity) are the Collatz statements closest to
  Mazur's development; the AxiomMath repository's `countBelow` is a reference statement to match. Preserved:
  everything; lost: nothing but time.
* **M5. The fast side of "expense is Diophantine".** The floor set's lattice structure (5.3). Map: the slow-side
  deadline ledger (records at the semiconvergents of `log_2 3`) → the fast set's size lattice `{k ln 2 - d ln 3}`;
  preserved: the `(k, d)` bookkeeping; lost: the carry `c_w`; test: HYP-9425 (4).
* **M6. One tier per external source.** The repository sweep found the same external claim typed three ways:
  openai/math `#158` (the plane is not 5-colourable, Lean) is CONDITIONAL in the second reading, CITED in session
  61's source status, and an accepted premise in the finite-algebraic-witness note; `#003` (zero-free region) is
  "PROVED modulo `#003`" (oaimath3), CONDITIONAL (S42) and PREPRINT CLAIM (S31); three sessions (S20, oaimath3,
  codex-tiling) took OpenAI manuscripts as accepted premises on the owner's instruction of 2026-10-07. None of
  these is a dependency of the main Collatz theorems (THM-4581, THM-4627). Rule added to the ledger's maintenance
  section: a source gets one tier, recorded once in the ledger, and every citation names the tier; an
  owner-directed acceptance is recorded as such (T2 with "owner-accepted premise"), not as CITED.

---

## 7. What is new, what is inherited, what is not claimed

| item | status | where |
|---|---|---|
| the pillar typing of the atlas (section 3) | SYNTHESIS | this note |
| Terras null case; three readings of the modular sentence | PROVED (classical restatements); CITED (session 63, THM-4764, kind-pasteur Lemma 1) | section 2 |
| Mazur rows, erratum, downstream, risk, anchor; AI-lab rows | CITED at the stated tiers | ledger, section 4 |
| `tau(n)` for all `n < 2^28`; densities; `z`-statistics | FINITE-EXACT; model comparison VERIFIED (P1, P2) | 5.2 |
| P3, P4 | FAILED as pre-registered; the oscillation OBSERVED | 5.2-5.3 |
| HYP-9425 (4): the word model reproduces the quarter-dyadic floor counts | VERIFIED (pre-registered after the failures; PASSED: correlation `0.992`, signs `24/27`, total ratio `1.098`) | 5.3 |
| the Gaussian shape of the fluctuation | OBSERVED off (right-skewed, skewness `+0.6`); scale verified to `5%` | 5.2 |
| Proposition F | PROVED (elementary; prior-art flag) | 5.4 |
| `E(C)`, `E(C*) = 1` | HEURISTIC (model) | 5.5 |
| trunk inequality | PROVED (one line on the digest's identity) | 5.7 |
| `3x-1` seeds | PROVED + FINITE-EXACT (P5) | 5.7 |
| the sign-blind transfer (with the reversed size inequality); the strong form of Mazur | CONJECTURE (HYP-9425) | 5.6-5.7 |
| M1-M6 | directions; M6 a maintenance rule | section 6 |
| Collatz; `lim H_n(1) > 0` | OPEN | — |

Not claimed: any change of status for Mazur's theorems (CITED, T2), any Collatz convergence or cycle statement,
any novelty for Proposition F beyond "not located", any 3-adic mechanism for the dyadic oscillation (the word
model accounts for it).

---

## 8. Reproduction

```bash
python 04-computation/experiments/collatz_fast_convergence_density_20261027/fast_convergence_density_20261027.py 28 > .../fast_convergence_density_20261027.out   # 0.6 GB, ~3 min
python 04-computation/experiments/collatz_fast_convergence_density_20261027/sawtooth_probe_20261027.py > .../sawtooth_probe_20261027.out                      # ~2 min
python 04-computation/experiments/collatz_fast_convergence_density_20261027/fast_set_exponent_model_20261027.py > .../fast_set_exponent_model_20261027.out    # seconds
python 04-computation/experiments/collatz_fast_convergence_density_20261027/word_model_fine_resolution_20261027.py > .../word_model_fine_resolution_20261027.out  # ~4 min, seeded
```

Numba is used for the stopping times (falls back to pure Python, slowly, without it). The Mazur text was read from a
`pdftotext` extraction of the attached PDF in a scratch directory (not committed), cross-checked against the
digest's section 2 and audit.

---

## 9. Obligations

1. Independent audit of this note (section 10), the three scripts and HYP-9425 before any index promotion.
2. The `2^32` run by a 2-bit sieve with the pre-registered thresholds of HYP-9425 (1)-(2), and the skewness of `z`
   as the shape test the Gaussian profile still owes.
3. M1's one-afternoon test; M2's depth-decomposition test from the digest's outputs; M3's P1-P3 replay.
4. The ledger rows are to be re-checked for external verification before any promotion (ledger section 4).

## 10. Audit records

### 10a. Adversarial referee reading of Mazur's paper (independent subagent, 2026-10-27; scratch scripts `audit_mazur/`, not committed; the report is summarised here)

The referee was given the paper's text only (a `pdftotext` extraction), not the digest, and asked to break the
chain. **Verdict: nothing checked is wrong.** Checked by hand or computation: Lemma 2.1; the identification of
(2.3) with Tao's Proposition 1.14 (the oscillation `Osc_(m,n)`, read from the arXiv HTML of 1909.03562v7: in Mazur's
recursion `Y_q mod 3^m` depends only on the last `m` valuations and has exactly the law `mu_m`, so the finite
distribution statement is Tao's); (3.6), (3.7) given (2.3); the residue multiplicity bound (3.9)-(3.11) (a history at
a fixed generation is determined by its endpoint, since the forward orbit is unique and meets `M` once; equal
`(D, A)` pairs have constant weight and endpoints spaced `>= 3^q` in an interval of length `2^(b_0+1)/W`; the role of
`q <= 2 b_n` is exactly to make `3^q W <= 16^(b_0) (9/16)^(b_n)`); Lemma 4.1 (seeds, `q <= 7`, three windows) and
Proposition 4.2; Lemma 5.1; the structure of Proposition 5.2 (the stated exponential is weaker than what
concentration gives, hence valid); Proposition 5.3 (the overshoot-reduction map `F_j(x) = 4^j x + (4^j-1)/3` is
`x + j mod 3`, so it does move units to nonunits and the full-group means are necessary); Lemma 6.1 (the identity
`3^d x = 2^A M - Sigma'`, `Sigma' > 0`, on 35,450 enumerated inverse orbits; injectivity history → source because the
forward orbit is unique and `M` is visited once; the range `[X, 32X)`); Lemma 6.2 and the 8.5 thresholds; the whole of
Section 7 at 50 digits (`296/1029` is the two-term series lower bound for `2 artanh(1/7) = log(4/3)`; (7.8) reproduces
`604680379/126443520000000` exactly after the two substitutions; (7.9) holds with slack `2.0e-4`, the final ratio being
`10.4597 <= 10.46`); the 8.3/8.4/8.5 arithmetic ((8.7), `(200/201)^6`, `2^(-1/100)`, `175/256`, the geometric moment
bound, the strip bound). The concentration inequality for geometric sums was verified exactly for `t <= 40`
(worst ratio `0.478`). Not verified: the explicit mixing coefficient `C` and Lemma 8.1 (Lean only), the exact rounding
constants in (6.2), the full Section 8 recipe beyond spot checks, and the Lean definitions behind `countBelow`
(the page shows the theorem but not the definitions; the source file was not inspected). The referee's reading of
why Tao's method stops at logarithmic density and Mazur's reaches natural density: Tao runs forward from a
scale-invariant start and needs the 3-adic scale to grow (so the target grows), while Mazur grows the inverse tree
from the fixed target with the weight `3^d 2^(-A)` as a Radon-Nikodym factor, applies mixing at scales `b_j -> infinity`
and extracts integer sources at the leaves without a `1/n` weighting; the price is `c` proportional to `1/M` with
`M >= 16^(2^80)`.

**Two expository findings that matter for the repository.** (i) *The extraction drops superscripts and the letter
`l`.* Read literally, "`b_0 >= 325`", "`b_0 = 5125`" and "`b = 280`" make several displayed Section 7 inequalities
false (the (7.2) coefficient bound fails at `b_0 = 325` and at `5125`; "`(2n+1)L <= B/10^6`" needs `b_0 >= 1.4 10^6`;
Lemma 7.1 needs `b >= 512^5`); under the consistent reading `b_0 >= 32^5`, `b_0 = 512^5`, `b = 2^80` everything
checks, and the owner should confirm against the PDF rendering. With the literal `5125` the achievable constant would
be `10.467`, so Corollary 1.2's `10.48` would still hold and `10.46` would not. The digest of 2026-09-28 carried one
such artifact: its "a fixed odd seed `M >= 16 b_0`" is `M >= 16^(b_0)` (corrected there today with a dated
bracket; MISTAKE-684 records the mechanism). (ii) *The constant is tight to its stated rational bounds*: (7.9) has
slack `2e-4`, so `523/50` cannot be lowered within the paper's own bounds without re-deriving `kappa = 34881/10000`.

**Surfaced, not read:** arXiv:2608.27617 (R. A. Tavares, 2026-08-27, math.GM), "Wirsching's Positive-Predecessor-
Density Program: Proofs of Conjectures 1 and 3". From its abstract (fetched today): it proves two of the three
conjectures in Wirsching's reduction of uniform positive predecessor density to a chain of conditions and reduces the
chain to the single condition `(*3)`, which it does not claim; no Lean, no AI assistance, no mention of Mazur. Entered
in the ledger at T2 (abstract only, conditional, not refereed). Note the direction: Wirsching's program is
root-uniform ("uniform positive predecessor density"), which THM-4517 shows cannot hold with a root-independent
constant; a per-root positive density is what Mazur's companion claims.

### 10b. Independent audit of this note (2026-10-27; report and scripts in `audit_S64/` beside the experiment)

**Verdict: SOUND WITH CORRECTIONS.** The auditor recomputed every number with code written before reading the
session's scripts: all `tau(n)` for `n < 2^24` exactly and the `2^28` tables, maxima, sawtooth figures and P1-P5
outcomes reproduced; the digest's `54.22%` at `10^6` is `542194/999999`; `sigma = 16.4220` from the renewal-reward
central limit theorem (Smith 1955; Asmussen ch. VI) with `t = ln n`; `523/50` is `0.305%` above `C*`; Proposition F
checked line by line (entropy bound to `k = 400`); `E(C)` rederived (`E(2.5) = 0.5841`, `E(3.307) = 0.7383`,
`E(8) = 0.9861`, maximum `1` at odd share `1/2`); the trunk weights, `R_1 = 1` and the `3/2` conversion of the
`rho_18` list verified; the `3x-1` seed identities for `j < 2000` and the permutation for `q <= 6`; THM-4476 on
both sheets and Tao's sheet-blindness via GGM confirmed; the citations of session 63 and kind-pasteur faithful; the
word-model script's weight `C(k-1, d) 3^(-d)` and carry formula correct (its run not re-executed).

**Eleven corrections, all applied above:** (1) the pillars were defined by what the *conclusion* changes under,
which contradicted the typing of rows M and K-L (Mazur's conclusion is sheet-blind, so by that definition it drops
P); re-defined as *proof inputs*, with the walls as conclusion controls (sections 0.1, 2, 3); (2) the withdrawn
lattice-proximity test survived in the title of HYP-9425, in sections 1, 5.3 and 9 and in the INDEX line; replaced
by the word-model test and its PASS; (3) the title's "fluctuation law verified" overstated: scale verified to `5%`,
shape right-skewed (skewness `+0.60 -> +0.66`), Gaussian profile not verified; the median deficit is skewness, not
a mean offset; (4) the trunk inequality's second step was a sign slip (`4x >= x` fails for negative `x`); the sharp
form `max_(j >= 2) H_n(R_j) >= 4 H_(n+1)(1) - 3 H_n(1)` is stated; (5) Proposition F needed the explicit remark that
the orbit has not reached `1` before step `k_0`; the finite-`X` model exponents are the fair comparison; (6) the
sign-blind transfer must say the size inequality reverses for `3x-1` (and the claim about the time-constant layer
was unverifiable from the repository); (7) the form of `c` is now cited with its equation; (8) two paraphrases were
in quotation marks; (9) `41.677` is in T-steps, "every periodic parity word", the S2 wall sentence, the attribution
of the product identity to Steiner-Eliahou / digest Theorem C, "within a factor 2" removed; (10) the histogram-bin
resolution of the density tables is now stated; (11) the ledger's T2 definition did not cover
"platform-kernel-checked, no dispute, no expert reading"; the clause was added. Not checked by the auditor: the
ProofAtlas pages and Tao's blog (not re-fetched), the paper's text (section 10a's referee record), the level-18
`rho_18(R_j)` values, the literature beyond the auditor's own knowledge. MISTAKE-685 records the mechanisms. The
auditor also noted that the digest's P7 row prints the size inequality in the wrong direction for `3x+1`; corrected
there with a dated bracket (MISTAKE-684, addendum).
