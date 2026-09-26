# Procedurally generated approaches to Collatz: the choice ladder, the missing no-divergence mechanism, and a relaxed problem that almost closes

**Status: SYNTHESIS of session `collatz-procgen-20260922` (mac-mini,
2026-09-22). PROVED statements are those proved in the lane notes (hand
proofs and one core-Lean theorem). FINITE-EXACT statements have stated
bounds and, where marked, two independent code paths. CITED statements
come from primary sources read in the barrier-atlas lane. Typology entries
are modelling judgments. Collatz, E-SCC (Q1 and Q2), HYP-9120--9122 and
Althöfer's game remain OPEN.**

**Wave 3 (2026-09-23), summarized in section 2c.**
* **Named target.** The divergence half is Lagarias's **Periodicity
  Conjecture** (PC) restricted to positive integers.
* **Theorem S** (PROVED and independently audited): no rational number has
  an eventually Sturmian parity vector under any `3x+r` map, for any
  slope.
* **Q1 mirror** (PROVED): every escape from `-1` costs more than 1.
* **Duality:** the forward–backward duality is REFUTED; the shared
  numerators are a PROVED straddle.
* **Endgames:** in both halves of E-SCC the endgame reads **low** `p`-adic
  digits. It is not Erdős's top-digit problem.
* **Foundry v4:** every all-orbits target in the family has the same
  single missing mechanism.
* **Wave 4 (HARD class):**
  * Theorem D (PROVED; audited): periodic approximants reach exactly the
    words with `Dio > eta`.
  * Theorem Y (PROVED; audited): a zero-entropy word with `Dio = 1`
    falls to a 2-adic Tschakaloff–Padé argument.
  * The smallest open instance is the explicit cube-swap word `Y3`,
    whose Bernstein number is a cubic 2-adic theta value (HYP-9127).

**Wave 11 (2026-09-25), summarized in section 2i.** **THM-4474**: in the strategy cube, bounded-lookahead provability holds iff every parity-graph cycle has odd density `< log_3 2`. Collatz's window is `[0,1]` at every level, and the OPEN fraction is about 0.435, which is the random-sign no-extra-cycle probability. **THM-4475** (HYP-9136 PROVED): explicit provable trees exist at flip density `<= 2^(1-0.05L)`; the sharp exponent is HYP-9137, and the cube analogue is HYP-9138. Cycle-gate equidistribution fails off the critical line, and a sheet-aware random model repairs the heuristic.

**Wave 12 (2026-09-25, opus session), summarized in section 2j.** The owner's seed "multiplication : squares :: addition : doubles" was run as a transport axis (128 cells, 47 cards). PROVED: the parity graph mod `2^k` is the `+-sqrt` graph on `F_p^*` for Fermat primes; the multiplicative Collatz `sqrt X / rad(X) X^3` is Collatz on the diagonal `m^e` and divergent off it; `chi_{-4}` is the unique level-2 provable strategy; the orbitwise two-place identity `m_L 2^(d_L)/3^L = n prod(1 + b/(3 m_l))` gives `R(d) <= n` on every `3n-1` orbit. **HYP-9160** (no slow divergence) was then PROVED by **THM-4476** (thin divergence: every non-periodic orbit has `O(X^(0.95+eps))` elements below `X`; Terras count plus a landing pigeonhole), which also recovers the in-house no-bounded-strip theorem and gives `R(d) < n` strictly; still no divergence exclusion.

**Waves 14–15 (2026-09-26), summarized in section 2l.** Five new audited theorems:
* **THM-4480.** The price of provability is the peak-discounted undecided density, with exponent `1 - H(log_q 2)` for every `q`. For `q = 3` it is `rho_L exp(-Theta(L^(1/3)))`, so P1 is negative for arbitrary edits.
* **THM-4481.** An entropy law for sign flips: no provable `qn±1` strategy exists for `q >= 23`, and 5n+1 needs constant flip mass.
* **THM-4482.** Ranks are tensions: Collatz's rank defect is `log(3/2)`, at `-1`.
* **THM-4483.** Forced charges on backward trees: no finite-mass 2-adic potential rank exists, and nonnegative ones exist iff Collatz.
* **THM-4484.** Free and sporadic cycles; the Belaga–Mignotte off-by-one is resolved.
* **THM-4485.** Periodic edit price = feedback sets of expanding cycles; Golomb–Mykkeltveit fails for thresholds.

New hypotheses: HYP-9140 (pairing peak price, rationale corrected), HYP-9141 (5n+1 concentration), HYP-9142 (Robin inequality).

**Wave 13 (2026-09-26), summarized in section 2k.** **THM-4479** (HYP-9138 PROVED): flip exactly the undecided residues; the cube distance is `2^(-(1-h)k+O(log k))`, so the exponent is sharp. With the crossroads THM-4478 (HYP-9137), the exponent `1-h(log_3 2)` is now sharp in three settings, by three mechanisms: integer capacity, necklace packing, and moments (the last stopping at `2(1-h)`). Gersonides: four of the five known integer cycles of 3x+1 on Z are forced by `2-1, 3-2, 4-3, 9-8 = 1`, and `-17` is the sporadic one (`139`). `3+1 = 2^2` is at once the pairing ladder, the moment criticality and the trivial cycle. Positive drift makes arbitrary edits exponentially cheap (the peak-discounted price). The owner's Kuratowski–Tutte triple is read as excluded structures, a counting reason, and a sporadic third. Wave 14 lanes: mykk, drift, tension, peak.

**Wave 10 (2026-09-25), summarized in section 2h.** Codex's incoming work was audited (20 checks, no errors), and the square-sum graph turns out to be planar iff `N ≤ 24`, a Kuratowski event at 25. Kohl's Collatz group is a Tait-coloured graph: its Kempe chains are the doubling orbits and the rising runs, and its only closed chain is `{−1, −2}`. In the strategy square, Collatz is the only open corner among four sign strategies. The exact triple shape is Tutte's "dual pair + self-dual" `{F_7, F_7*, U_{2,4}}`. Natural boundary is KNOWN (Bell–Lagarias 2015); the Mahler/harmonic bridges are blocked by the controls.

**Wave 9 (2026-09-24), summarized in section 2g.** arXiv 2502.20642 (a claimed fixed-point proof of Collatz) is invalid (**THM-4471**): its general theorem fails on `x -> x+1`, and its table also "proves" `3n−1`. The owner's triangle sandwich is its Lemma 2.1, and the sign law is exactly the sandwich's two equality cases with `0` in the middle. The correct fixed-point theorem is Banach in `Z_2`: one gate per word, and 18 integral points for `p ≤ 24`. **THM-4472**: the owner's converse 4-tournaments are forward map versus inverse tree (time reversal), not the sheets. **THM-4473**: Collatz's digit chains are exact Markov laws (mod 10 entries `2^j/15`, `3 -> 5` certain, `9 -> 9 = 8/15`), unlike the vanishing prime-digit bias. Repunit primes are the prime fixed points of digit rotation, and `Q` has the new odd 2-cycle `{−1/5, 5/7}`.

**Wave 8 (2026-09-24), summarized in section 2f.** The owner's tree/mod-192 programme is exact, and residue types are provably side-blind (`x -> -x` transport), so any proof needs the sign law. **THM-4469** (Mahler bridge): no-divergence on adjacent supercritical block pairs is equivalent to a generalized Mahler Z-number statement beyond the Flatto–Lagarias–Pollington length (HYP-9134). **THM-4470** (pairing ladder): 3x+1 is the unique AM-fair consecutive pairing, its graph is two perfect difference systems, and density-zero flips falsify it (HYP-9135 and HYP-9136). Brackets: `{2,3,11}` fully mapped. Verdicts on the owner's analogy: graceful is an ANALOGY; square-sum with the brackets is REAL; square-sum with Collatz is NUMEROLOGY.

## 0. What was asked and how it was answered

The request was to generate new approaches to open problems like Collatz
*procedurally*, test the hypotheses, and look through past work for the
fringe ideas the repository has brushed, "to see how the missing insights
have been lurking." Six earlier Collatz waves (`collatz_mod6_*`,
`arithmetic_braids*`, `collatz_guards_*`, `collatz_blueprint_*`) had
mostly audited pasted blueprints. The pasted snippet itself had already
been audited
([level11_short](level11_short_20260922.md): the `1/15` law is PROVED; the
gluing and tournament numerology is typed).

This session built two instruments and pointed them at many siblings:

1. **A foundry.** It generates approach cards
   `problem x sub-target x mechanism x lens` over seven Collatz-type
   problems. Each mechanism is typed by the controls it is blind to (the
   `3n-1` sheet, `5n+1` drift, planted defects, rational 2-adic cycles,
   undecidability) and by whether it needs a thin exceptional set
   ([foundry](collatz_procgen_20260922_foundry.md)).
2. **An exceptional-set profiler.** It computes, for a descent game, the
   residue classes that have no descending certificate, for Collatz, its
   sheets, `5n+1`, and any relaxation with choice
   ([choice ladder](collatz_procgen_20260922_choice_ladder.md)).

## 1. Results

**The ladder of freedoms** (FINITE-EXACT counts; one dimension PROVED):

| system | freedom | exceptional set |
|---|---|---|
| Collatz `T` (either sheet) | none | dimension `h(log_3 2)=0.9500` (PROVED; not found in print, per the atlas); `1,037,374` classes mod `2^26` |
| `E_S`, `S={6 mod 8}` | extra `3n+1` only at rising-run entries | `3,238` mod `2^22` (Collatz `93,222`). Dimension `>= 0.0536` PROVED; about `0.17` by counts (HYP-9121 REFUTED, wave 7) |
| greedy fingerprint | 8 of 32 even classes mod 64, led by `54 mod 64` | `893` mod `2^20` (full choice `664`) |
| graph `E` (Q1 forward) | extra `3n+1` at every even | `908` mod `2^36`, `2454` mod `2^64`; dimension OPEN, evidence favours 0 (power law about `m^1.6`); infinite (Theorem F) |
| backward `E` (Q2) | reverse moves with choice | `157` mod `3^32`; `338` mod `3^42` (dimension lane) |
| Applegate--Lagarias semigroup | arbitrary wild multipliers | one class, `-1 mod 2^j` (CITED) |
| additive choice `3n+b`, `b in B`, `|B|>=2` | choice of sheet | **empty** at `2^22` |
| `5n+1`, no choice | none | positive measure `mu_5=0.17603` (PROVED, sibling ladder) |
| `5n+1` with `E`-choice | branch choice | falls `0.26->0.066` (`2^26`, exact) and about `4e-5` (`2^64`, MC): collapse CONJECTURED; `Haar<=0.0664` PROVED. **Choice is blind to the Collatz drift** (corrected; MISTAKES 2026-09-22) |
| `qn+1`, `E`-choice, `q>=41` | branch choice | positive Haar measure (PROVED, sibling-ladder Thm 3); heuristic threshold `q` about 10 |

**The relaxed problem E-SCC** (HYP-9120). `E` turns out to be Le--Smith's
*Loosened Collatz Graph* (arXiv 2109.01180, CITED; they state neither Q1,
nor Q2, nor strong connectivity).

* Q1 (every `n` reaches `1`) is implied by Collatz, so it holds below
  `2^71` (CITED).
* Q2 (`1` reaches every `m` prime to `3`):
  * **reduced by a Lean-checked mod-27 lemma to the two classes
    `1, 14 mod 27`.** Every other `m` descends within two reverse moves by
    a factor of at most `8/9`;
  * **verified below `10^18`** (wave 3,
    [endgame lane](collatz_procgen_20260922_q2_endgame.md)): an exhaustive
    DFS over the `Psi`-alive classes mod `3^38` gives every `m <= 10^18` a
    multiplicative descent. Earlier bounds were `7.87*10^17` (dimension
    lane, exact threads to depth 41) and `2.02*10^13` (lane one's DFS);
  * **PROVED conditional form** (endgame lane, Theorem 6.1): Q2 follows
    from `X_min`, "no integer `m >= 2` lies in `Bad_inf`", with no
    threshold and no HYP-9122. Assuming HYP-9122, it also follows from
    `X_T`, a statement about the digit-exhausted points
    `2^(floor(k log_2 3)) w` (Theorem 6.2);
  * the neighbourhood of `1` is handled by a PROVED escape lemma together
    with loops through `1` (every depth up to 41);
  * every escape from the neighbourhood of `1/2` costs more than 1. The
    canonical price is `2c(k-1)/3`, in `(1,2)`, with infimum 1 at the
    records of `log_2 3` (PROVED). Its exits are exactly the Collatz trunk
    `(4^i-1)/3`. **CORRECTED 2026-09-23:** this line previously read "at
    least `32/27`", which is the floor of one route only (MISTAKES
    2026-09-22);
  * the canonical chains obey a sharp budget `exp(c* D)`, with
    `c* = ln(128/81)/4 = 0.1144` per digit consumed by expensive links,
    for every hostile thread (PROVED, conditional on HYP-9122 at the
    precisions used). The backward exceptional set is infinite
    (Theorem F_b, the mirror of Theorem F).
* Q1's thread of `-1` (wave 3, [Q1 mirror](collatz_procgen_20260922_q1_mirror.md)):
  * every exit from `n = 2^m u - 1` costs exactly `3^(eta(m+1))`, which
    lies in `(1,3)`. At the rigid precisions `m = 5..9` it costs `3.8`
    to `12.8` (PROVED);
  * so `-1` behaves like Q2's `1/2`, not like Q2's `1`, and a see-saw is
    needed on both sides;
  * every `n = 7 mod 8` below `2^32` descends within 43 halvings
    (FINITE-EXACT);
  * hostile chains follow `u' = (3^A u + 1)/2^(m'+1)`, a Collatz-type map
    with memory.
* Consequence: Le--Smith's Conjecture 1 (every `n` prime to `3` lies on an
  `E`-cycle) holds below `10^18`.
* Past-work link: Le--Smith's Conjecture 2 says every nontrivial `E`-cycle
  uses an `E`-only arrow (`3n+1` at an even `n`). It is equivalent to
  Collatz having no nontrivial positive cycle. The 2026-09-17 session's
  census, in which all `74` simple `E`-cycles of length at most `40` in
  `[1,2000]` use such an arrow
  ([extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md)),
  is a window check of it, made before the paper was known.
* Hostile points (PROVED examples, FINITE-EXACT families):
  * forward: `-1` and `-1-2^i c/3^j`, e.g. `-13/9`, negative, with
    power-of-3 denominators;
  * backward: `1`, `1/2` and `1/2+3^j c/2^e`, positive dyadics.

  Each family accumulates at its base point in its own adic metric. They
  are mirror images under `2<->3` and the sign.

**The typology** (foundry plus atlas):

* The Collatz cycle half has four live mechanism types (Baker-type,
  carry anti-concentration, functional equations, order patterns).
* The **no-divergence half has none**, and neither do the `3n-1` sheet or
  rational periodicity.
* No published mechanism overcomes SHEET and DRIFT together at unbounded
  complexity (atlas, 22 primary results).
* Removing DIMENSION by choice exposes SHEET: Applegate--Lagarias's
  hostile class and `E`'s hostile rationals sit at the minus sheet's
  cycles.

**Where the sign enters** (PROVED, elementary). Positive plus-sheet cycles
contract multiplicatively (`3^L<2^K`), so they are invisible to
class-level descent certificates and live only at the thresholds
`B/(2^K-3^L)`. Positive minus-sheet cycles expand (`9/8`, `2187/2048`), so
their minima are hostile points. Hence the class-level method sees only
the divergence half of Collatz, and the sign enters exactly as the sign of
`2^K-3^L`.

**Negative controls.**

* The undirected Collatz game is equivalent to Collatz and does not
  collapse: max depth `59` and `0.2%` unresolved to `10^6`.
* The F2[x] analogue has no drift and no sheets, which is why it is
  provable (atlas).

## 2. The remaining lanes (all integrated)

* **Sibling dimension ladder**
  ([note](collatz_procgen_20260922_sibling_dimension_ladder.md)):
  * PROVED: `dim=h(log_3 2)` for `3n+b`, for every odd `b`. Exact counts to
    `m=10^5` fit `C 2^(hm) m^(-3/2)` with fitted exponent `0.94995552`
    against the proved `0.9499555`, and power `-1.4997`.
  * PROVED: positive measure for `q>=5`, with `mu_5=0.17603`,
    `mu_7=0.30075` to `10^-22` (Spitzer series).
  * PROVED: the F2[x] exceptional set is `{1}`.
  * The choice-game threshold is heuristically near `q=10`; positivity is
    proved for `q>=41`. This **corrects** the drift claim above.
  * Mahler's safe set has dimension `log_2(3/2)` (THM-3848).
  * Lagarias 2009: `dim E^(1)=log_3 2`, `dim E^(2)<=1/2`; he conjectures
    `dim E(Z_3)=0` (CITED).
* **Sign-specific order laws**
  ([note](collatz_procgen_20260922_order_laws.md)):
  * FINITE-EXACT: 6,016 grammar-generated order statements on both sheets
    to `10^7`; no sign-specific law beyond the sign law.
  * PROVED (direction lemma): a plus window can contradict its parity
    word's prediction only by decay going up, a minus window only by
    growth going down. So windows of at most 12 odd steps are sheet-blind,
    and the realized ordinal patterns are exactly the word-realizable
    ones.
  * The inherited "min `3 mod 4`, max `1 mod 4`" law is a normalization
    artifact.
  * Terras's `sigma=tau` holds to `10^7` on both sheets and is PROVED for
    `tau<=38` (plus).
  * Best plus-only statement: "every growth window grows" (the sign law).
    Its smallest minus witness is the transient near-cycle `165->163`
    (clock `19/12`).
* **Loops through 1 and escapes**
  ([note](collatz_procgen_20260922_loops_and_escapes.md)):
  * PROVED: loop ratio `>3/2` (a Catalan-type fact), so the 1-escape needs
    exactly `K0` halvings and is the only exit. HYP-9122 reduces to the
    record denominators of `log_2 3`.
  * FINITE-EXACT: loops exist to `s=6000`, so the thread of `1` is settled
    to depth `6001`.
  * PROVED: the height lower bound (no finite family); the `1/2` escape
    price is exactly `2c(k-1)/3`; seven dyadic hostile points.
  * The see-saw does not close. **The endgame of Q2 is read from the base-3
    digits of `2^K w`.** Wave 3 (section 2c) refines this: the digits read
    are the **low** (3-adic) ones, which the kappa formula
    `v_3(2^K - r) = 1 + v_3(K - kappa(r))` controls. Erdős's problem is
    about the top digits, and every 3-adic-window version of it is FALSE.

## 2c. Wave 3 (2026-09-23): the named target, a Sturmian theorem, the Q1 mirror, and one missing mechanism

* **Transversality foundry**
  ([note](collatz_procgen_20260922_transversality_foundry.md)).
  * **The target is named** (CITED). The "p-adic irrationality of the
    Bernstein series" of section 5 is exactly Lagarias's **Periodicity
    Conjecture** (1985, section 2.8; Bernstein–Lagarias 1996;
    Monks–Yazinski 2004). Its positive-integer instance, "a Bernstein
    number is never a positive integer", is Bernstein's 1994 form of no
    divergence (T1).
  * **Catalogue** (35 typed items). Every counting or dimension theorem is
    DEFECT-blind. Every every-element 2-versus-3 theorem excludes only a
    **zero-entropy** class. The only proved every-orbit results on
    positive-entropy classes of parity words use growth or capacity:
    * subcritical words, by Monks–Yazinski;
    * bounded critical discrepancy, by the in-house
      [discrepancy theorem](collatz_guards_20260921_discrepancy.md),
      extended to rationals as Proposition B.
  * **Theorem S** (PROVED; independently audited on 9 rows by a separate
    code path). No rational with odd denominator has an eventually
    Sturmian parity vector under any `3x+r` map, for every slope and every
    intercept. It also holds for `5x+1` at slopes `alpha < 0.804` and for
    Mahler's map, so no Z-number has a Sturmian carry word.
    * The proof is a 2-adic Liouville argument. The parity-vector map is
      an isometry, and the periodic extensions `u v^inf` of a word's own
      repetitions give rational approximants. A rational cannot be
      approximated that well.
    * Sturmian words repeat a block of length `q_n` for a stretch of order
      `q_(n+1)`, which suffices because `mu(mu-1) < phi`.
  * **Generator.** It produces 49 statements, with a matrix of statements
    against targets. The rows that matter most:
    * Erdős's problem needs the top ternary digits; the 3-adic-window
      forms B1 and B2 are FALSE;
    * the E-SCC Q2 endgame needs the low digits along a dynamically chosen
      `K` (the kappa formula);
    * for Mahler, the drift control `5/2` is a proved theorem.
  * **The HARD class.** Words that are supercritical, positive-entropy and
    non-repetitive are untouched by every proved mechanism. This is the
    precise missing piece of T1 and PC. The smallest open instance is
    candidate **C2**: bounded discrepancy around a supercritical slope.
* **Q1 mirror** ([note](collatz_procgen_20260922_q1_mirror.md)).
  * **Loops through `-1`** (PROVED):
    * their equation is `2^K + B = 3^a`;
    * every loop has ratio above 2 except the basic loop `MH`;
    * the records are the **lower** best approximations of `log_2 3`,
      against the **upper** ones for Q2;
    * the first three record cycles are exactly the three negative
      Collatz cycles, which is how the Catalan unit-gap clocks enter.
  * **Mirror of HYP-9122.**
    * FINITE-EXACT: it holds for `10 <= K <= 4000`, with explicit
      verified loops.
    * PROVED: it fails for `K = 5..9`, the rigid `-1` neighbourhood.
  * **The exit from `-1` costs more than 1 at every precision** (PROVED).
    The chain recursion and the landing law are PROVED as well. The
    endgame reduces exactly (PROVED) to the low binary digits of `3^A u`,
    with landing depth `v_2(A - lambda(u)) + 1`. Dupuy–Weirich 2016 is
    the averaged analogue (CITED; not read). The pointwise statement is
    OPEN.
  * **Duality REFUTED.** The shared numerators `793585` and
    `419868489953` pair a hostile forward point with a backward point
    that descends, at depths 53 and 212. For `N = 3^j + 2^(e-1)` the
    straddle identity `(|x|-3/2)(y-3/2) = -(rho-1)^2/(2 rho)` (PROVED)
    places exactly one partner above `3/2`. No letter-count-linear word
    map preserves non-descent (PROVED via Gelfond–Schneider).
  * **Census.** 117 backward dyadic points below `3/2` (`e <= 45`) are
    certified hostile, and the 4 above `3/2` descend. Whether any hostile
    point lies above `3/2`, on either side, is OPEN.
* **Q2 endgame** ([note](collatz_procgen_20260922_q2_endgame.md)).
  * **Structure** (PROVED). Every 3-adic unit lies on the thread of `1`
    or of `1/2`. The canonical escape `Psi` has exactly one expanding
    branch, the `1/2`-transfer at precision `k >= 3`, with price
    `2c(k-1)/3` in `(1,2)`. Every hostile thread's link factors through
    it (lift lemma); all 98 census points are `Psi`-preimages of `1` or
    `1/2`.
  * **Budget theorem** (PROVED and sharp). Chains cost at most
    `exp(0.1144 D)`. Two claims are REFUTED: "digits run out" (a
    37.6-digit `m` runs 15 expensive links over 60 digits) and "at most
    `m^0.104`" (`6082250` reaches `6.57`).
  * **Conditional theorems** (PROVED implications). Q2 follows from
    `X_min`, "no integer `>= 2` is hostile", with no thresholds. Assuming
    HYP-9122, Q2 also follows from `X_T`, on the digit-exhausted points
    `2^(floor(k log_2 3)) w`.
  * **Theorem F_b** (PROVED). `1/2 + 3^i/2^(K0(i-1)+1)` and
    `1/2 + 3^i/2^(K0(i-1)+2)` are hostile for every `i >= 3`, so the
    backward exceptional set is infinite.
  * **Diophantine inputs** (CITED). LTE, Yu, Senge–Straus/Stewart and
    Lagarias 2009 each control a piece of the endgame; none controls all
    of it.
  * **Numerics.** Q2 holds below `10^18` (FINITE-EXACT).
  * **Census points above `3/2`.** The four undecided census points above
    `3/2` all descend, by the mirror lane's certificates, replayed
    independently by the orchestrator.
* **HARD class, wave 4** ([note](collatz_procgen_20260922_hard_class.md); audited).
  * **Theorem D.** `Phi_T(w)` is irrational whenever the Adamczewski–Bugeaud
    Diophantine exponent exceeds the height rate: `Dio(w) > eta(w)`.
    Bugeaud–Kim's `Dio >= 2.50994` for every Sturmian and quasi-Sturmian
    word then covers every map with `eta < 2.50994`, **including every
    slope of `5x+1`**. The threshold is sharp for the method: an extremal
    word under `19x+1` escapes it.
  * **Theorem Y.** The square-swap word `Y` is explicit, supercritical
    (`beta = 9/10`), of zero entropy and of width 1, and has `Dio(Y) = 1`,
    so Theorem D does not apply. Its Bernstein number is
    `-1 - 512/(3^9-2^10) - (256/3^9) sum_k rho^(k^2)` with
    `rho = 2^10/3^9`, a 2-adic theta value. A 2-adic transcription of
    Zudilin's Tschakaloff–Padé construction proves it irrational: no
    rational has an eventually square-swap parity vector under any
    `3x+r` map.
  * **The smallest open instance** is the cube-swap word `Y3`, whose
    Bernstein number is `sum rho^(k^3)`, a cubic theta value (HYP-9127).
    Both mechanisms fail on it: `Dio = 1`, and there is no first-order
    q-difference equation.
  * **C2** (HYP-9123) is equivalent, for integers, to a *coupled* Z-number
    statement: no `L != 0` and strip word `w` with
    `frac(L 3^(a_s)/2^s) = frac(E_s(w))` for all `s`. It is PROVED on strip
    words with `Dio > mu` and on square-swap words.
  * **Correction to foundry v4.** HARD is not "positive entropy". It is
    "`Dio <= eta` and no functional equation behind the word". The foundry
    requirement is renamed accordingly.
* **Foundry v4** ([note](collatz_procgen_20260922_foundry.md), section 3b).
  * It adds the structural requirements HARD (first called ENTROPY; refined in wave 4) and CHAIN. In the
    relaxation both base points cost more than 1 to escape, so the
    hostile chains are again a no-divergence problem, for a Collatz-type
    map with memory.
  * **Every all-orbits target now has the same single unblocked real
    mechanism type** (sound certificate search), plus the placeholder
    "transversality on HARD words". The targets are Collatz and `3n-1`
    divergence, PC, both halves of E-SCC, Mahler, Erdős, and the `5n+1`
    existence question.
  * **The relaxation renormalizes the difficulty rather than removing
    it.**

## 2d. Wave 5 (2026-09-23): the cube-swap number, the trunk and RH, outside results, and AMM 12592

* **Cube-swap number** (HYP-9127; [note](collatz_procgen_20260923_cube_theta.md); audited). HYP-9127 is still OPEN.
  * **Theorem Q** (PROVED). Every single quadratic swap family (triangular, pentagonal, any `alpha k^2 + beta k`, with a sign twist) is irrational when `mu_bar < phi`, so for every `3x+r` map.
  * **Why cubes resist** (PROVED). They have no linear q-difference equation of any order (via Garoufalidis). Mahler's method is inadmissible, since the matrix is unipotent. A Subspace argument needs a growing number of S-unit terms, and no shift nearly preserves the cubes; this is exactly what separates it from Erdős 1062(ii).
  * **Reductions.**
    * A combinatorial 2-adic zero estimate (HYP-9130) implies HYP-9127.
    * So does a uniform S-unit gap `U(1/3)`. The n-term abc conjecture does not suffice.
  * **Certificates.**
    * The number is not a rational of height `<= 2^4999999`.
    * Lattice hunts sit exactly at the Dirichlet baseline.
    * The Hankel census shows no Padé miracle for cubes. Squares show one: their determinants are 4.4 times below random.
  * **Cheaper open instances.** A square swap just past `phi` (HYP-9131), and the near-critical cube `sum (2^19/3^12)^(k^3)`.
* **The trunk and the Riemann hypothesis** ([note](collatz_procgen_20260923_trunk_rh.md); audited).
  * **The trunk as a 3-adic object** (PROVED). `T(i) = (4^i-1)/3` is a 3-adic isometry of `Z_3` with fixed points exactly `{0, 1, -1/2}`, and `T(1/2) = -1`. The plus and minus trunks interleave as the Jacobsthal numbers.
  * **The trunk as an Iwasawa coordinate** (REAL). The trunk is `1/3` of Iwasawa's coordinate `4^s - 1` for the Kubota–Leopoldt `zeta_3`, which has no zeros.
  * **The two "1/2"s differ** (PROVED). The E-game's hostile `1/2` sits at the irrational exponent `i* = log(5/2)/log 4`, while `s = 1/2` maps to `-1`.
  * **Blindness.** Every trunk identity is blind to DRIFT and to SHEET.
  * **RH controls.** An RH foundry types the bridges against Davenport–Heilbronn and Epstein, whose off-line zeros were located numerically.
  * **The only precise leftover** is a resonance gap for the E-game's own dynamical zeta (HYP-9133).
* **Outside results** ([note](procgen_sources_20260923_interplay.md); audited).
  * **Erdős 1062(ii)**, which is Lean-accepted, runs on the same three-place S-unit engine as Theorems S and D. Its only difference is the height exponent `eta = 1`. The same exponent bookkeeping gives the no-go for `Y3` along natural-height approximants.
  * **Proposition T′** (conditional on the Subspace Theorem): for supercritical bounded-discrepancy words with `Dio > 1`, the 2-adic and real Bernstein values are not both rational.
  * **Corollary M** (PROVED modulo the p-adic Mahler–Manin theorem, BDGP 1996). The square-swap number is a 2-adic theta value on the Tate curve `q = rho^2`, whose `j` carries the 196884. So the square-swap and pronic-swap numbers are not both algebraic.
  * **The owner's moonshine exercise** (verified): `6 = 5+1` goes through `S_5 = PGL(2,5)` inside `S_6`, then `M_12`, `M_24`, Golay and Leech, to `196884 = 300 + 24 + 196560 = 1 + 196883`. Also, the snippet's level-11 recurrence is the eta product `eta(tau)^2 eta(11 tau)^2` attached to `M_24`'s order-11 elements.
  * **Seven hexagons** fit at side `5/sqrt3`, exactly (Morandi 2015, a record; optimality open).
* **AMM 12592, the owner's extractor question** ([note](amm12592_procgen_20260923_uniform_frontier.md); audited and **promoted as [THM-4467](../../01-canon/theorems/THM-4467-uniform-polya-capacity-gap-amm12592.md)**).
  * A Pólya-capacity argument on the `p <-> 1-p` quotient proves the first uniform gap, **`C* >= 11/8`**. Pure interval arithmetic gives `27/20`, and the method reaches `1.3775`. So a constant `1+eps` is impossible for `eps <= 3/8`.
  * Theorem B makes the golden constant a natural boundary of the lacunary series.
  * Super-blocks beat every separately balanced block at `N = 8, 64, 128, 256` (FINITE-EXACT) and approach about `1.570` numerically (HYP-9128).
  * Whether `C* < 3/2` is OPEN (HYP-9129). The proved window is `[1.375, 1.598]`.

## 2e. Waves 6–7 (2026-09-23): the owner's "prove the rest" — status of every session hypothesis

| hypothesis | status after waves 6–7 | where |
|---|---|---|
| HYP-9120 (E-SCC) | OPEN. It contains Q1 and Q2; the dimension route via HYP-9121 is void | sweep note |
| HYP-9121 (rising-run excursion) | **REFUTED**. In `E_{6 mod 8}` the orbit of `-1` is a 16-point trap of rate `2/3`, so `dim >= 0.0536` (PROVED) | sweep note §2 |
| HYP-9122 (loops through 1) | OPEN; FINITE-EXACT to `s <= 190535` (was 6000) | sweep note §3 |
| HYP-9123 (C2, supercritical strips) | OPEN; placement `PC ⟹ C2 ⟹ HYP-9127` | sweep note §6 |
| HYP-9124 (`X_min`, implies Q2) | OPEN; `⟸ HYP-9126` | sweep note §5 |
| HYP-9125 (loops through -1) | OPEN; FINITE-EXACT to `K <= 176249` (was 4000) | sweep note §3 |
| HYP-9126 (3/2 wall) | OPEN; `⟹ HYP-9124 ⟹ Q2`; the easy clocks are PROVED given the loops; all 589 generation-1 points with `i <= 300` descend | sweep note §4 |
| HYP-9127 (cube-swap number) | OPEN. Proved no-go for every natural determinant family (parabola lemma, Theorem NG, exact `e3(n)`); loophole: quartic cyclotomic content | theta round 2 |
| **HYP-9128** (super-blocks) | **PROVED, THM-4468**: `C* <= 159/100 < 1+log_5(phi^2)` | hyp9128 proof note |
| HYP-9129 (`C* < 3/2`?) | OPEN. Proved window `[1.377, 197/125]` (THM-4494, opus 2026-09-26: exact ratio in THM-4468's bottom regime); realizable states stall at about `1.567` | same |
| HYP-9130 (cube zero estimate) | OPEN; restricted forms only | theta notes |
| **HYP-9131** (square swap beyond `phi`) | **PROVED**: 2-adic Hankel (Bézivin) for `mu_bar < 7/4`, and `28/11` and `2.878` via KRVZ; every square swap under `5x+1` and `7x+1` is settled | theta-beyond-phi |
| HYP-9132 (transcendence of the square-swap number) | OPEN; **degree `>= 3` PROVED** (2-adic KRVZ non-quadraticity, cited inputs) | theta round 2 |
| **HYP-9133** (E-game zeta gap) | **PROVED** (Pringsheim + aperiodicity + Hurwitz) | trunk/RH §9 |

**Structural findings of these waves.**
* **One product-formula argument, two places.** THM-4467's Pólya argument and the H1 Hankel proof are the same argument; AMM runs it at the real place, Collatz at the 2-adic place ([bridges](procgen_bridges_20260923_lrc_amm_collatz.md)).
* **`phi` marks the limit of one-scale methods.** In every case it is beaten by coupling scales:
  * Bugeaud–Kim `2.51` in Theorem D;
  * `7/4`, `28/11` and `2.878` for theta values;
  * super-blocks at `1.59`, below the golden `1.598`, with about `1.570` numerically.
* **The owner's AM–GM principle, made exact.**
  * The Collatz drift is exactly the AM–GM gap `log(2/sqrt3)`, and `q = 3` is the only odd `q` with arithmetic-mean step factor 1.
  * THM-4467's termwise majorant is sharp (a negative result, §9 of the hyp9128 note). The real loss is cross-level cancellation, and HYP-9128's super-blocks exploit exactly that.
  * For LRC, "excised Bonferroni" certificates are exact on the tested 13-speed rows. Example: `{1..12, 5460}` has lonely measure `301/10296`, re-computed independently.
* **The owner's inspiration texts, dispatched** ([dispatch](procgen_numerology_20260923_snippet_dispatch.md), [brackets](procgen_numerology_20260923_odd_square_brackets.md)).
  * PROVED: `{2,3,11}` are exactly the primes that reach their own double inside their odd-square bracket.
  * The "every 5th" pulse is a real local-density effect: `13σ` mod 5, with similar effects mod 3 and 7.
  * REFUTED: the pentagonal-pulse lemma (at `j = 15`) and the Tower Packing Limit (at `M_4 = 2^127-1`).
  * REAL anchors:
    * Catalan parity at Mersenne indices is THM-4467's Lemma P;
    * Euler's pentagonal function is irrational at `2^10/3^9` (Theorem E);
    * `11` is a base-3 Wieferich prime, but no mechanism was found for `33^2 = 1089`.

## 2f. Wave 8 (2026-09-24): what the features represent, the owner's mod-192 tree, the implication atlas, and the pairing ladder

The owner asked four things:
* what each fundamental feature of Collatz *represents*, and how it can be reduced or re-expressed;
* for links to famous problems "so that one can be shown to prove the other";
* to study "incongruency", singletons and small sets (`{2,3,11}`, odd-square brackets);
* to replace "prove every number reaches 1" by the tree programme: double every odd, double the evens, send each number `2 mod 3` to `(2n-1)/3`, and decide where it goes mod 192.

The owner also offered an analogy: graceful tree : 3N+1 :: square-sum arrangements : a microcosm–macrocosm transition. Three lanes ran, each audited by the orchestrator.

**1. The owner's tree programme is exact, and provably needs one outside input** ([inverse tree mod 192](collatz_procgen_20260924_inverse_tree_mod192.md); audited: byte-identical rerun, independent Perron root, ladder cycle, Z-cycles and lattice sum).
* **Equivalence (PROVED).** Collatz holds iff the `D/E` tree of `{1,2}` (`D(x) = 2x`, `E(x) = (2x-1)/3` at `x = 2 mod 3`) is all of `Z_(>0)`.
* **The "4x fractal recursion"** is the identity `E D^2 = S E` with `S(p) = 4p+1`. The owner's "`(N-1)/3` versus `(4N-1)/3`" are consecutive rungs of the sibling ladder. When `N` is odd, the lower rung is even, and the step from it is exactly the extra arrow of the E-graph relaxation.
* **The automaton mod `3*2^k`** has one core component (the non-multiples of 3) with growth `4/3`. Multiples of 3 are leaf cones.
* **What mod 192 decides.** A residue mod 192 decides six parity steps and `E(n) mod 128`, but never `E(n) mod 3`, which needs `n mod 9`.
* **The no-go (PROVED, Theorem 6 / Corollary 7).** `x -> -x` carries every residue, Haar, integrality, `l`-adic and **size** datum of `3n+1` onto `3n-1`. Residue types are therefore *side-blind*: they prove for the positive half-line exactly what they prove for the negative half-line, where `3n+1` has three more cycles. Any proof must use the **sign law** (`2^p T^p(x) - 3^a x = c_w >= 0`), and after that, integrality at the gates `2^p - 3^a` and pointwise avoidance of the null set `Bad`.
* **The two-sided target (Theorem 15).** Collatz and 3n-1 no-divergence together are equivalent to a sheet-symmetric statement about all integers.
* **Small sets.** The HYP-9121 16-point trap is exactly the mirror image of the three `3n-1` cycles. `Bad_+ cap Z` in `[-10^7, 10^7]` is `{-17, -5, -1}`.
* **The drift is a residue.** The averaged tree series `1/(1 - 2^-s - (3/2)^s/3)` has residue `1/log(2/sqrt3) = 6.95` at `s = 1`: the reciprocal of the drift.
* **"Incongruency".** The owner's mod-3 refinement removes about half of the undecided classes but keeps the exponent `h(log_3 2) = 0.95` (PROVED).

**2. The one arrow out of the family: THM-4469 (Mahler bridge)** ([atlas](procgen_atlas_20260924_collatz_implication_atlas.md); PROVED + INDEPENDENTLY AUDITED).
* **The equivalence.** No-divergence on an adjacent supercritical block pair is *equivalent* to a generalized Mahler Z-number statement. The smallest instance: no positive integer's parity vector is eventually made of `0111101110` and `1101100111` **iff** no `xi > 0` has all `xi (2187/1024)^j in Z + [4726, 4727]/1163`.
* **Beyond FLP.** Proved Mahler-type exclusions (Flatto–Lagarias–Pollington, Dubickas, Bugeaud) stop at length `1/p`, and this instance needs `1.88/p`. So Collatz, T1 or PC implies new Mahler-type theorems, and a Mahler-side proof would settle a HARD slice of no-divergence (HYP-9134).
* **The rest of the atlas** (58 nodes, 83 edges):
  * no famous conjecture is known to imply Collatz, T1, NC or PC;
  * Lang–Waldschmidt gives cycle length `>= N^(1/2-eps)` in the verification bound;
  * abc gives prefix bounds (extending Rozier), but at the cycle gate abc is dominated by Baker;
  * CST together with T1 implies Collatz;
  * with Hercher and Barina's `2^71`, a nontrivial cycle has at least `137,528,045,312` odd terms (CITED);
  * Kohl: Collatz is equivalent to transitivity of a three-generator class-transposition group, and wildness is DRIFT- and SHEET-blind;
  * Collatz's 1932 permutation shares the cycle gates, and its signed carries cancel `13` and `7153`;
  * the pair-sum property is the martingale property, and Doob's theorem gives Terras's.

**3. The pairing ladder: THM-4470** ([pairings and transitions](procgen_brackets_20260924_pairings_transitions.md); PROVED + INDEPENDENTLY AUDITED).
* **Arithmetic-mean fairness.** `T(2i-1) + T(2i) = (2i-1) + 2i`: each consecutive pair keeps its sum while its product shrinks by `3/4 + 1/(8i-4)`. This is the owner's AM–GM principle in exact pairwise form.
* **Uniqueness.** Only `3n+1` on `{2i-1, 2i}` and `3n-1` on `{2i, 2i+1}` preserve pair sums: the sheets *are* the two consecutive pairings.
* **Graceful form.** The halving edges and the up edges each realise every difference exactly once, yet no subtree with at least 2 edges is literally graceful.
* **The family.** In the family of all pairings every member shares these properties, and `3n-1` is Collatz's antipodal corner.
* **Fragility.** A single flip creates a cycle at exactly 24 fragile pairs, all `<= 2308` (to `10^7`; HYP-9135).
* **Density-zero falsification.** A density-zero flip set yields a divergent orbit. So no pairing statistic decides tree-ness (DEFECT, made concrete).
* **Provability.** No periodic pairing is provable by bounded lookahead (the pair-0 obstruction is Applegate–Lagarias's `-1`). A provable (landing ⇒ down) pairing must flip at least 29% of the pairs, and window designs suggest the price tends to 0 (HYP-9136).

**4. Brackets and the analogy (same note).**
* **The escape set.** `{2,3,11}` is the prime escape set for every ratio `r in (25/13, 25/11]`. The complete list of `(k,p)` with `kp` in `p`'s bracket is `(2,2), (3,2), (4,2), (2,3), (3,3), (2,11)`. So every other prime is the only multiple of itself in its own bracket.
* **The microcosm.** Every Collatz move from `n >= 54` changes bracket.
* **Nothing special about 11.** Being a base-3 Wieferich prime is a 1-in-11 event, and the trunk match has probability 1/6: NUMEROLOGY.
* **Verdicts on the analogy.**
  * "Graceful tree : 3N+1" is an ANALOGY: the structure is exact but shared by divergent members.
  * "Square-sum : brackets" is REAL: both are square-density thresholds ending near 24–25.
  * "Square-sum : Collatz microcosm" is NUMEROLOGY: Collatz's exceptions are Diophantine, pinned to convergents of `log_2 3`, and recur.
  * The square-sum problem has an exact `x49` microcosm-to-macrocosm engine (`(49a+c) + (49b-c) = (7s)^2`, partitioning `[25, 1249]`). Collatz has none.

**What the features represent (the owner's first question, answered in one table).**

| feature | represents | sees | blind to |
|---|---|---|---|
| parity bit | a 2-adic digit (Terras bijection; conjugacy to the shift) | all residue statistics | the sign, integrality |
| multiplier `3` | the unique AM-fair multiplier (pair sums preserved) | the drift, as the AM–GM gap `log(2/sqrt3)` = 1/(tree residue) | which orbits realize it |
| `+1` vs `-1` | which consecutive pairing, i.e. the side of `0` | only the sign law and the gates' signs | every residue, density, size datum |
| `2 mod 3` | branch points of the inverse tree = images of up-moves | branching `4/3` | SHEET |
| `4p+1` ladder | `E D^2 = S E`, the owner's "4x recursion"; limit `-1/3` | the 3-adic rotation | order |
| trunk `(4^i-1)/3` | the sibling ladder of the root | Q2's exits | SHEET, DRIFT |
| `2^K - 3^L` | the incongruence of 2 and 3: cycle gates | cycles and fragile pairs | divergence |
| exceptional set | classes no congruence decides (`dim 0.95`, with or without the mod-3 refinement) | where certificates fail | whether a positive integer lies in it |

**The wave's answer to "reduce or re-express, and connect".**
1. Every re-expression the owner proposed (the tree, the mod-192 automaton, the pairing, the graceful form) is exact and PROVED equivalent or faithful. Each one is also provably blind to the sign, to density-zero defects, or to both.
2. The inputs a proof must add are now named exactly: the sign law, integrality at `2^p - 3^a`, and pointwise avoidance of a null set.
3. The one new two-way bridge to a classical problem family is Mahler's (THM-4469). There Collatz implies exclusions just beyond the Flatto–Lagarias–Pollington length, and the smallest such exclusion is a HARD slice of no-divergence.

## 2g. Wave 9 (2026-09-24): Brouwer and fixed points, the owner's triangle sandwich, four-vertex tournaments, and digits/repunits

The owner asked us to:
* think Brouwer's fixed point theorem;
* connect Rédei, Hamiltonian paths and "fixed point chain growth" with arXiv 2502.20642;
* hone in on the triangle sandwich `|d(x,z) − d(z,y)| ≤ d(x,y) ≤ d(x,z) + d(z,y)` and its squares, reading the two sides as positive and negative with the centre as 0;
* consider the claim that "3" and "±1" are the two 4-vertex tournaments swapped by reversing all arcs;
* look at the non-uniform last-digit transitions of consecutive primes, circular primes such as `{337, 373, 733}`, and repunit primes and their lengths.

Three lanes ran, each audited by the orchestrator.

**1. The paper is a claimed proof of Collatz, and it is invalid: THM-4471** ([fixed-point audit](collatz_procgen_20260924_fixed_points_kawasaki_audit.md)).
* **What the paper is.** arXiv 2502.20642 is Kawasaki's "A proof of the Collatz conjecture". The owner's sandwich is its Lemma 2.1, and the lemma is correct: it is the triangle inequality followed by AM–QM, and it discards the cross term `±2ab`.
* **The general theorem is false.** Its weighted-pseudocontraction fixed-point theorem (Theorems 2.1(5)–2.3(5)) fails:
  * `x -> x+1` on `(N, |x−y|)` satisfies every hypothesis with the paper's own constants;
  * the least counterexample has 3 points.
* **The gap.** The proof swaps a quantifier: it needs the second alternative at the pair `(Tp, p)`, but the hypothesis only supplies it at `(p, Tp)`. Collatz falls into the gap at every odd step.
* **Sheet control.** The paper's coefficient table, copied verbatim, also "proves" that `3n−1` reaches 1.
* **No contraction argument in `|x−y|` can work.** Up-steps stretch consecutive distances by up to `3/2`. Caristi's principle and a discrete contraction metric are each *equivalent* to Collatz, and Bessaga's contraction metric to its cycle half.
* **The owner's sandwich, made exact.**
  * The sign law is precisely the two equality cases, with `z = 0` in the middle: positives realize the upper case and negatives the lower.
  * The `b = 0` "central" sheet is Mahler's `3x/2`, and THM-4469 says the `±1` orbits shadow the central Mahler orbit with a confined carry.
  * The paper breaks exactly at the `+2ab` up-steps.
* **The correct fixed-point theorem is Banach's in `Z_2`, applied to the inverse branches.** Every parity word has exactly one fixed point, its cycle gate `c_w/(2^p − 3^a)`. This fixed-point chain grows like `2^p`, and for `p ≤ 24` exactly 18 of its points are integers: the five known cycles.
* **My Sharkovskii suggestion was REFUTED.** Every continuous extension has all periods on both sides. The sign does decide *stability* in Chamberland's extension.

**2. The two diamonds are time reversal, not the two sheets: THM-4472** ([four-vertex note](procgen_tourn_20260924_four_vertex_sheets_redei.md)).
* **The quadruple.** Take the AM-fair pair and its two images, with the map arcs and every other pair oriented by numerical order. The result is the 3-cycle over a sink `(0,2,2,2)`, with `H = 3`, iff `b·s_o > 0`; otherwise it is transitive.
* **The inverse tree gives the converse.** The inverse tree yields the source over a 3-cycle `(1,1,1,3)`, and the reflection `x -> s_o + s_e − x` identifies the two.
* **The owner's "±1" is the direction of time.** `3n+b` inverts to `(y−b)/3`.
* **Negation fixes the class,** and reversing all six arcs equals negation composed with time reversal.
* **Rédei's parity is the opposite of Collatz's.** Rédei leaves exactly one unpaired configuration. Collatz's `2^p` periodic points are paired off completely by a free involution that commutes with `T`.

**3. Digits, rotation and repunits** ([repunit note](procgen_repunit_20260924_digits_rotation_repunits.md)).
* **Lemke Oliver–Soundararajan, reproduced and extended.** Their tables are reproduced exactly and extended to `10^11`. The bias is real but decays, and their second-order term explains 86% of it.
* **Collatz's analogue is exact and does not decay.** For consecutive odd Syracuse terms:
  * mod 10, `3 -> 5` has probability 1, and every other entry is `2^j/15`, `j ∈ {0,1,2,3}`; for example `9 -> 9` has probability `8/15`;
  * mod 3 the terms are i.i.d. with law `(0, 1/3, 2/3)`;
  * mod 9 the stationary law is `(8,16,11,4,2,22)/63`.

  The owner's first snippet's "1/15" reappears as the unit of this matrix.
* **Rotation.** Rotating digits is multiplication by `b` modulo `b^k − 1`. Repunit primes are exactly its prime fixed points, and `{337, 373, 733}` is a free orbit. So "Brouwer: fixed point versus orbit" is REAL. The claim that the digit bias *generates* circular primes is NOT SUPPORTED.
* **The repunits inside Collatz.**
  * `T^k(2^k − 1) = 3^k − 1`, i.e. base-2 repunits run up to twice base-3 repunits.
  * The 2-adic repunit limits `ξ_b = 1/(1−b)` are moved by the Möbius map `b -> −(b+2)/(b−4)`.
  * The parity-vector map `Q` of Bernstein–Lagarias has odd fixed points `−1` and `1/3`: the repunits of bases 2 and −2.
  * `Q` has a **second odd 2-cycle `{−1/5, 5/7}`** alongside `{1, −1/3}`. It is exact, and its literature status is UNVERIFIED.
* **Lengths of all-ones primes.** As hidden Collatz structure this is NUMEROLOGY: 24 tests, none significant.
* **Corrected lead.** `−1/3` is not a fixed point of `T`, since `T(−1/3) = 0`.

**What wave 9 adds to the picture.**
1. The owner's inequality instinct is right in a precise sense. The sign lives in the cross term of the triangle sandwich, i.e. at the centre `0`. Every even observable (`|x|`, `x²`, `|x−y|`) discards it. That is why a metric fixed-point argument in `|x−y|` cannot be side-aware.
2. The right fixed-point theorem, Banach in `Z_2`, supplies every periodic point for free. It leaves exactly the integrality question at the gates.
3. The owner's tournament picture is real, but it encodes *time reversal*. The sheet swap is negation, and the full converse is their composite.

## 2h. Wave 10 (2026-09-25): the incoming codex work, Kuratowski/Tutte, and discrete ↔ continuous

The owner asked for three things:
* a deep synthesis of all new incoming work;
* a long session of creative proof angles;
* a comparison of `{Petersen, K_{3,3}, K_5}` under Kuratowski and Tutte with the session's triples, including bridges between discrete and continuous mathematics.

Three lanes ran, each audited by the orchestrator.

**1. Incoming work** ([incoming synthesis](collatz_procgen_20260925_incoming_synthesis.md)). Codex's session `collatz-bugs-20260925` produced 11 commits and 40 notes.
* **What it built.** Careful re-encodings of Collatz: tournaments, marked Pythagorean triples, Fano/E8 lattices, and "decoders" with sidecar data.
* **Audit.** Twenty independent spot checks found no numerical or logical disagreement.
* **Genuine contributions.**
  * Three repairs. Codex fixed the square-sum degree-2 forcing and gave a new `Q_24` proof, fixed THM-060 Type A, and fixed THM-4473's `k ≥ 2` boundary.
  * Limitation theorems for finite certificate searches. One of them shows that no finite-state encoder computes the parity-vector map, which closes the automata route to PC.
* **Flags.** The "plus nonlinear / minus linear" carry result is a binary-digit artifact, and no DRIFT controls were run.
* **New here: the owner's square-sum transition is a Kuratowski event.** `Q_N` is planar iff `N ≤ 24`. Its first `K_{3,3}` (branch vertices `{3,4,5,11,12,13}`) runs through the ear `11–25–24`, and all 10 Hamiltonian paths of `Q_25` must use that ear. The orchestrator re-verified both facts.

**2. Kuratowski, Tutte and Kohl's Collatz group** ([Tait/Kempe note](procgen_kuratowski_20260925_tait_kempe_triples.md)).
* **Kohl's graph is a Tait colouring.** Kohl's theorem says Collatz is equivalent to transitivity of `G_C = <a, b, c>`, three class transpositions. The Schreier graph of `G_C` is exactly the undirected Collatz graph on `Z \ 0(6)`, **properly 3-edge-coloured**: `a` is the up-edge, and doubling edges are coloured `b` or `c` by the sign of `m mod 3`. So Collatz is the statement that one Tait-coloured subcubic graph is connected.
* **Kempe chains are the session's objects.**
  * `<b,c>` chains are the orbits of the Banach contraction `D(x) = 2x`.
  * `<a,c>` chains are the **rising runs**, segments of `E(x) = (2x−1)/3` of length `2 v_2(y+1)` (the 2-adic distance to `−1`; re-verified to `2·10^5`).
  * `<a,b>` chains have at most 3 edges.
  * The only closed Kempe chain is the digon `{−1, −2}`, the Banach fixed point of `E`.
* **Petersen obstructs nothing inside a single sheet.** Every component of a functional graph has at most one cycle, so each sheet's graph is planar. Kuratowski graphs appear only when the two sheets are superimposed, in Althöfer's `3n±1` union graph `U_N`:

  | first `N` | event |
  |---|---|
  | 52 | `K_{3,3}`, on branch vertices `{5, 7, 11, 14, 20, 26}` (re-verified), which lie on the `3n−1` cycle's neighbourhood |
  | 68 | `K_5` |
  | 76 | no planar double cover |
  | 92 | Petersen-family minors (Colin de Verdière `μ ≥ 5`) |
  | 104 | Petersen |

  Every kernel of `U_N` checked (all `N ≤ 1000`, and `N = 2000, 4000`) is 3-edge-colourable, so there are no snarks.
* **The strategy square.** Choose the sign of `3n±1` by `n mod 4`:
  * always `+` is Collatz (OPEN);
  * always `−` is `3n−1` (PROVED intransitive);
  * the sign that always forces two or more halvings is PROVED transitive;
  * the sign that always forces exactly one halving is PROVED divergent.

  **Collatz is the only open corner.** Negation swaps Collatz and `3n−1` and fixes the two mixed strategies.
* **The exact triple shape is Tutte's, not Kuratowski's.**
  * The recurring shape is **"dual pair + self-dual"**, as in Tutte's regular-matroid obstructions `{F_7, F_7*, U_{2,4}}`. It is REAL for:
    * the sheets `{+1, −1 | 0}`;
    * the means `{GM, QM | AM}`, where `QM² + GM² = 2AM²`;
    * the 4-tournaments `{diamonds | TT, strong}`;
    * Kohl's generators `{b, c | a}`, since the sheet switch exchanges `b ↔ c`;
    * the strategy square.
  * "Twins + container" (`K_5`, `K_{3,3} ⊂` Petersen) is REAL only for the SHEET/DRIFT twin atoms (no container) and literally inside `U`.
* **An excluded-minor theorem for methods (PROVED).**
  * A sound method proving Collatz must separate all three controls.
  * The three controls are independent (a diagonal table of invariants).
  * SHEET and DRIFT are incomparable excluded minors of the rcwa order, and DEFECT lies outside RCWA.
  * The rcwa order is not a well-quasi-order, so there is no finite Kuratowski list. Connectivity has infinitely many excluded minors.

**3. Discrete ↔ continuous** ([natural boundary / Mahler / harmonic note](procgen_continuous_20260925_natural_boundary_mahler_harmonic.md)).
* **The natural boundary is KNOWN.** Collatz holds iff the basin series `B(z)` is rational, iff it is D-finite, iff it continues analytically across some arc of `|z| = 1` (Pólya–Carlson). This is Bell–Lagarias, Acta Arith. 170 (2015), Thms 1.1–1.3. The coordinator's lead was a rediscovery; the lane extended it to every odd `q`.
* **Mahler/Cobham.**
  * Separate pure 2-Mahler and 3-Mahler equations for `B` would force rationality (Schäfke–Singer), hence Collatz.
  * The real equation mixes in roots of unity, and those provably destroy the rigidity.
  * SHEET, DRIFT and DEFECT each block the route; for example, the planted-defect map satisfies a pure 2-Mahler equation and still diverges.
* **Trees and Tutte embeddings.**
  * The branching numbers are `≥ 1.2334` for `3n±1`, `1.1227` for `5n+1` and `1.0870` for `7n+1`, so simple random walk is transient on all these trees.
  * The Tutte (barycentric) embedding is the harmonic-measure transform.
  * All of these invariants are sign-blind and drift-blind.
* **Hex/Brouwer does not transfer.** Althöfer's game has infinite plays from every odd `n ≥ 3`.

**What wave 10 adds.**
1. The owner's triple is best read through Tutte's matroid form, "dual pair + self-dual". That shape is realized exactly by the sheets, the means, the tournaments, Kohl's generators and the strategy square.
2. Collatz is literally a connectivity statement about a Tait-coloured cubic graph. Its Kempe chains are the doubling orbits and the rising runs, and the 2-adic point `−1` is the only closed chain.
3. Every discrete ↔ continuous bridge tested is either already known (natural boundaries) or blind to a control (harmonic, tree and Mahler data).
4. The creative frontier is where the three controls meet. The strategy square isolates Collatz as the one undecided corner among its nearest neighbours.

## 2i. Wave 11 (2026-09-25): the exact boundary of provability

Three lanes pursued proof angles suggested by wave 10. Each was audited by an independent re-implementation.

**1. The strategy cube (THM-4474)** ([cube note](procgen_cube_20260925_strategy_cube.md)). This extends wave 10's mod-4 strategy square. Each strategy assigns a sign to the odd residues mod `2^k` (Althöfer's 3n±1 with the choice frozen into residues). All 65,814 strategies at levels 1–5 are classified, and three theorems are PROVED:
* **(A)** A bounded-lookahead descent proof exists iff every cycle of the finite parity graph has odd density `< log_3 2`.
* **(B)** A residue divergence proof exists iff some closed class has *all* cycles expanding.
* **(C)** The drift, the quantity every heuristic uses, is sandwiched between the extreme cycle densities.

What this gives:
* Provability is a property of the extreme cycles, not of the average. Collatz has the maximal window `[0, 1]` at every level, from the loops at `0` and `−1`. So no residue argument of either kind exists, at any modulus.
* Collatz's distance to the provable class, measured in flips, is `1, 2, 2, 4, 5, 9, 14, 23` for `k = 2..9` (Haar 0.5 → 0.09). Whether it tends to 0 is HYP-9138.
* **Macrocosm and microcosm.** The OPEN fraction stabilizes at about 0.435. That is exactly the probability of *no extra cycle* in the `k = ∞` random-sign model (extra-cycle counts `{0: 43.5%, 1: 47.6%, 2: 8.5%}`).
  * In the macrocosm (drift, exceptional dimension, density window), random strategies converge to Collatz's values.
  * The microcosm, the signs at the smallest integers, decides. 97.8% of the cycle witnesses have minimum `≤ 100`.
  * So Collatz's cycle-free positive side is typical, not miraculous.

**2. The price of provability (THM-4475; HYP-9136 PROVED)** ([price note](procgen_price_20260925_provability_price.md)).
* **Construction.** For every `L ≥ 8`, an explicit greedy member of the pairing family descends within `L` steps and is a tree. Its flip density is at most `2ρ_L ≤ 2^(1−0.05L)`, where `ρ_L` is Collatz's undecided density. Each flip is paid for by one Collatz-undecided number; the only delicate rescue sits on the 2-adic neighbourhood of the cycle `−5 → −7 → −10`.
* **Bounds.** The lower bound `2^(−0.774L)` is proved; the sharp exponent is HYP-9137.
* **Collatz between two dense regions.** Collatz therefore sits at density zero from *divergent* members (THM-4470) and at density `2^(−0.05L)` from *provable* members.
* **Audit.** The orchestrator's re-implementation, written from the prose alone, reproduces the flip densities exactly.

**3. Cycle gates and equidistribution** ([gates note](procgen_gates_20260925_gate_equidistribution.md)).
* **Census.** Exact censuses for `p ≤ 40` find only the known cycles of `3x±1` and `5x+1`.
* **No equidistribution.** The residues `c_w mod (2^p − 3^a)` are not equidistributed off the critical line.
* **The naive model fails.** The naive random model predicts `ln P + O(1)` positive cycles. The sheet-aware repair, which counts only cycles with least point `≥ 1`, matches all controls.
* **Priority.** Polynomial per-clock bounds are classical (Belaga).
* **Unresolved.** Two Belaga–Mignotte 3x+d counts are off by one.
* **Ranking.** LOW as a proof angle, but a better heuristic.

**What wave 11 adds.** The boundary of provability is now exact on two ladders, and both are governed by the same `0.95` exponent:
* in the strategy cube, provability is decided by the extreme cycle densities (THM-4474);
* in the pairing family, the price of provability is the undecided density `2^(−0.05L)` (THM-4475).

Collatz sits exactly on that boundary. It has the widest possible window, it is approachable by provable trees at vanishing cost, and it is equally approachable by divergent ones.

## 2j. Wave 12 (2026-09-25, opus session `collatz-squares-doubles-20260925`): the squares/doubles transport axis

Full note: [squares/doubles foundry](collatz_sqdbl_20260925_squares_doubles_foundry.md).
The owner's seed "multiplication : squares :: addition : doubles" was used as
a transport axis (exponential, discrete log, monoid extension, polarization,
quadratic characters, Jacobi symbols, local and global squares, Hilbert
symbols, the real value of the digit series, cubic theta polarization,
GM-fairness, heights, sums of squares) against eight Collatz ingredients:
128 cells, 47 curated cards.

* **Exact model of THM-4474 (PROVED, FINITE-EXACT).** For a Fermat prime
  `p = 2^k + 1` the parity graph mod `2^k` is the `+-sqrt` graph on `F_p^*`
  (`y -> +-sqrt y` on squares, `y -> +-sqrt(g y^3)` on non-squares); checked
  for `p = 17, 257, 65537`. The two Collatz loops are `y = 1` and `y = g^-1`.
* **Multiplicative Collatz (PROVED, FINITE-EXACT to `10^5`).**
  `M(X) = sqrt X` if `X` is a square, else `rad(X) X^3`, returns `X` to its
  radical iff `X = m^e` with `m` squarefree and `e` Collatz-convergent, and
  diverges off that diagonal (parity desynchronisation plus Terras
  injectivity). The square gate is a synchronisation constraint, not a
  mechanism.
* **Residue laws (PROVED).** `v_2(3n+s) >= 2` iff `s = chi_{-4}(n)`, `>= 3`
  iff also `chi_8(n) = -1`; hence `chi_{-4}` and `-chi_{-4}` are THM-4474's
  unique level-2 class-(i) and class-(ii) strategies. `(3/m_i) = (-1)^(v_i + [v_(i+1) = 1])`
  on both sheets (the reciprocity sign couples adjacent valuations);
  `(5/m_i) = (-1)^(v_i)`. Odd squares `(2t+1)^2` map to `3t(t+1)+1`, again a
  square iff `2t+1 = y_j` (`j` odd) of `x^2 - 3y^2 = 1`; no length-3 chains.
* **Two places, orbitwise (PROVED; the one sheet- and drift-aware cell).**
  `m_L 2^(d_L)/3^L = n prod_(l<L)(1 + b/(3 m_l))`, so the real value of the
  Bernstein series is `R(d) = b n (prod (1 + b/(3 m_l)) - 1)`. On the minus
  sheet `R(d) <= n` for every positive orbit, with equality iff
  `sum 1/m_l = infinity`; on the plus sheet `R(d) < infinity` iff the orbit
  diverges with summable reciprocals. The sign enters the divergence half
  exactly once, as this inequality. This extends Proposition T (hard-class
  lane) and the Eliahou identity (S4) from cycles and bounded-discrepancy
  words to every orbit.
* **HYP-9160 (OPEN).** No slow divergence: every divergent orbit has
  `sum 1/m_l < infinity`; equivalently the real and 2-adic values of a
  non-periodic word never coincide at a positive integer. It isolates the
  only corner of the divergence half where the real place gives an equality.
* **THM-4476 (PROVED, same session, after the owner asked for HYP-9160 at
  discrepancy `O(log l)`).** Thin divergence: every non-eventually-periodic
  orbit of `x -> x/2, (3x+b)/2` (`b` odd) has at most `C X^(h*+eps)` elements
  below `X`, `h* = h(log_3 2) = 0.95`, uniformly. Proof: Terras's class
  count (no-dip words have at least `rho k` odd letters) plus a pigeonhole
  (a dipping point lands on the same orbit below `X^(1-theta)`), then a
  bootstrap. Consequences: HYP-9160 holds on both sheets and for rationals;
  every divergent orbit is full-rate; `R(d) < n` strictly; `m_j > j^a`
  infinitely often for `a < 1.05268`; the in-house no-bounded-strip theorem
  is recovered and extended to log-bands `C < 0.02634`; the union of all
  cycles is thin; `R_2(d)` is irrational for log-drift words with
  `a < 1.05268`. Not for non-constant sign strategies (`-chi_(-4)` is the
  witness: the parity map is not a bijection). Not a divergence exclusion. See
  [thin divergence](collatz_thin_20260925_thin_divergent_orbits.md).
* **THM-4487 (2026-09-26, opus).** The dip spectrum: `#{n <= X : T^i(n) >= n^gamma, i <= log_2 n} = X^(h(gamma/log_2 3)+o(1))` for `gamma in (log_4 3, 1]`, exponent `1` below Korec's `log_4 3`; THM-4476's no-dip count is sharp (`X^(h(rho)+o(1))`), so `0.95` is optimal for its dichotomy and the frontier item 1 of the thin-divergence note is settled for the method. The constants `0.95, 0.05, 0.7925, 0.2075, 0.488, 1.0527` are one curve; the repo-wide census of recurrent numbers is typed in [the atlas](constants_atlas_20260926_recurrent_numbers.md).
* **Verdict.** No winning reframe in the sense of a proof route. The
  transport axis produced one exact reformulation, one hypothesis of the
  right size, one exact model, one new object with a theorem, and residue
  laws that need not be regenerated (cosmetic, degenerate and blocked cells
  are recorded, including the cubic theta polarization already excluded by
  HYP-9127's unipotent Mahler obstruction).

## 2k. Wave 13 (2026-09-26): one exponent in three settings, the crossroads work, Gersonides's cycles, and the Kuratowski–Tutte reading

**1. HYP-9138 PROVED = THM-4479** ([cube-distance note](procgen_cubedist_20260925_distance_to_provability.md)).
* **Construction.** Flip exactly the undecided residues `Bad_k`. The construction works by heredity:
  * a decided residue meets no flip before its first descent;
  * at an undecided residue (always `3 (mod 4)`) the flip is itself a `3/4` descent.
  
  So every orbit of `sigma_k` is a chain of contracting blocks, and `sigma_k` is class (i). Its critical density `F_k`, the best lower approximation of `log_3 2` with denominator `<= k`, is attained by an upper Christoffel word.
* **Lower bound.** The disjoint expanding necklaces of `B(2,k)` force `delta_k >= N_k >= 2^(hk)/(3k^2)`. So the Haar distance is `2^(-(1-h)k + O(log k))`, and the exponent is sharp.
* **Exact data.** `delta_k = 1, 2, 2, 4, 5, 9, 14, 23` for `k = 2..9`, and `40 <= delta_10 <= 44`. The approximants are transitive for `k <= 300`.
* **DRIFT.** The 5n±1 provable class is empty at levels 2–6. 5n+1's distance is `29/64` at `k = 7` and about `3.2/k` up to `k = 14`; whether it tends to 0 is OPEN.
* **Audit.** An independent re-implementation, written from the note's statements, confirms every theorem-level claim, and the pipeline rerun is identical up to timing.

**2. Incoming work, synthesized.** Sources: crossroads-20260926, collatz-reframe-20260925, bernoulli-boundary-20260925, and opus S5.
* **THM-4478 (crossroads) PROVES HYP-9137.** Arbitrary fixed-horizon edits, and the pairing family, cost `2^(-(1-h)L+o(L))`.
  * *Mechanism.* On undecided words `T^k(n) = w_k(n + h_k)` with an affine offset `h_k in [0, k/3]`. So a fixed endpoint, time and odd count admit at most `floor(k/3)+1` integer ancestors. A growth band keeps `2^(hL-o(L))` words, and a first-hit cut finishes the proof.
  * *Relation to THM-4477.* This bypasses THM-4477's distribution-only barrier `2(1-h)` without contradicting it.
  * *Refuted in the same work.* (i) The conditioned energy, since `E Q^2 >= 1/8` via a Bernoulli(3/4) martingale. (ii) Fixed-slope injectivity: `N = 2^153 u - 1` and `N - 4` merge at time 233 with 153 odd steps, shadowing `-1` and the `-5` cycle. Injectivity does hold for every `k <= 31`.
  * HYP-9139 is untouched, and no exponent depends on it any more.
* **Reframe.** It audits the owner's Kuratowski–Tutte paste on Berggren addresses and proves three things:
  * the block repetition lemma `R_w(n) = floor((v_2(E_w n) - 1)/S)`;
  * the reset family `n_H = (2^(H+3) - 13)/9`;
  * no rank `a log n + sum c_i v_2(n - beta_i)` over finitely many fixed centers can decrease at every odd step.
  
  [derived] The reset family's 2-adic limit `-13/9` is HYP-9120's proved hostile point, and the repetition count for `w = (1)` is the `<a,c>` Kempe chain length of wave 10.
* **Bernoulli boundary.** The first Bernoulli function's jump is an exact divisor indicator: `J_m(z) = 1/m + psi((z-1)/m) - psi(z/m) = 1_(m|z)`. No `a log n + h(n)` with bounded `h` decreases at every odd step. [derived] THM-4474's Lemma P potentials are exactly bounded residue corrections, so this obstruction is "Collatz is not bounded-lookahead provable", caused by the loop at `-1`. Lane `tension` is proving the equivalence.
* **THM-4476 (S5) audited SOUND. Updates to §2j:**
  * the log-band statement is superseded by the one-sided bound: no injective orbit has `Delta_j >= -a log_2 j - O(1)` with `a < 1/h*`;
  * Cor 8 speaks only to the no-divergence half.
* **Cross-links** [derived by the digest, unaudited]:
  * THM-4469's HYP-9134 pair (carries `4726, 4727`) meets the crossroads Thue–Morse and Rudin–Shapiro irrationality conditions (`2187 < 4^10`, `2187^7 < 2^110`). So TM- and RS-selected tapes over this pair are not the parity vectors of any rational.
  * Kohl's colour `b` on `{y, 2y}` (`y = 1 mod 3`) is Berggren's letter `C`.
  * The source's B and C graphs are Kohl's quotients by the `<b,c>` Kempe chains, so `B ∪ C` is Althöfer's union graph with the doubling rays contracted.

**3. One exponent, three settings, three mechanisms.** The price of bounded-lookahead provability now has the sharp exponent `1 - h(log_3 2) = 0.0500445` in every setting studied.

| setting | upper bound | lower bound | mechanism of the lower bound |
|---|---|---|---|
| pairing family | THM-4475 (`2 rho_L`) | THM-4478 | integer capacity (affine offsets) |
| arbitrary fixed-horizon edits | THM-4478 (`rho_L`) | THM-4478 | integer capacity |
| strategy cube (periodic) | THM-4479 (`\|Bad_k\|`) | THM-4479 | necklace packing (periodicity) |
| (distribution-only arguments) | — | THM-4477: stops at `2(1-h)` | moments, which lose a square |

**4. Orchestrator findings** ([note](procgen_wave13_20260926_orchestrator_findings.md), with two check scripts).
* **Gersonides's four cycles.** An integer cycle of `3x+1` on `Z` is *free* when every parity word of its shape is integral. The free cycles are exactly `{0}, {-1}, {1,2}, {-5,-7,-10}`, i.e. `2-1`, `3-2`, `4-3`, `9-8 = 1`. The proof is a shift argument plus Levi ben Gershon's theorem of 1343.
  * The fifth known cycle `{-17, ...}` is sporadic: its shape `(11,7)` has `3^7 - 2^11 = 139`, and exactly one of its 30 necklaces is integral.
  * The five densities `0, 1/2 | 1, 2/3, 7/11` are the first best lower and upper approximations of `log_3 2`.
  * `7/11` is the mediant of `2/3` (the `-5` cycle) and `5/8` (the critical Christoffel density of `sigma_k`).
* **`3 + 1 = 2^2` in three roles.** It is simultaneously:
  * THM-4470's pairing-sum preservation;
  * THM-4477's moment criticality `g_q(2) = (1+q)/4 = 1`;
  * the trivial cycle `1 -> 2 -> 1`.
  
  Each of the three holds iff `q = 3`. The continuous moment curve `g_3(s)` meets 1 exactly at the integer points `s = 1, 2`, which are the free cycles `{0}` and `{1,2}`.
* **Positive drift makes arbitrary edits cheap.** For `5n+1` the undecided density stays near `0.2`. Yet catching each bad orbit high up makes `eps_L(5)` exponentially small.
  * The digest's sharper form is the *peak-discounted density*: `eps_L(q) = rho^peak_L(q) = 2^(-L) sum_(Bad_L) 1/max_j slope_j`, up to `poly(L)`. The upper bound sends each bad source to 1 at its peak; the lower bound stratifies THM-4478's capacity by peak height.
  * If this holds, the exponent is `1 - H(log_q 2)` for every odd `q`, whatever the drift sign (`0.01391` for `q = 5`). For `q = 3` the peak discount is `exp(-Theta(L^(1/3)))`, which would answer THM-4478's question P1 negatively for arbitrary edits.
  * Lane `peak` is auditing and extending this. The periodic cube cannot see height and is polynomially sharp.

**5. The owner's triple `{Petersen, K_{3,3}, K_5}` under Kuratowski and Tutte (typed; full table in the findings note §3).**
* **The pattern.**
  * **Kuratowski/Wagner** characterize a property (planarity) by excluded substructures.
  * **Euler's formula** is the counting reason: each of the three graphs violates `E <= g(V-2)/(g-2)` at its girth `3, 4, 5`.
  * **Tutte** adds a sporadic third obstruction for a neighbouring property, the Petersen graph for 4-flows, which contains both Kuratowski graphs as minors. His regular-matroid list `{F_7, F_7*}` + `U_{2,4}` is a dual pair plus a self-dual.
* **The same pattern in Collatz.**
  * Bounded-lookahead provability is characterized by excluded substructures, namely expanding cycles (THM-4474 A).
  * The counting reason is entropy: `h(log_3 2)` sets the price in all three settings above.
  * The two smallest expanding obstructions, `{-1}` and `{-5,-7,-10}`, are *forced by an identity*: `3 - 2 = 1` and `9 - 8 = 1` (Gersonides/Catalan).
  * The third, `{-17, ...}`, is *sporadic* (`139 | c_w`). Its word `11110111000` contains the words `1` and `110`, it sits at the next Stern–Brocot approximant, and it is barely expanding (`3^7/2^11 = 1.068`).
  * Every class-(i) strategy at every level must break all three.
  * The level-2 cube is Tutte's shape: the `nu`-dual pair Collatz/3n−1 (both class (iv)), plus the self-duals `chi_(-4)` (the unique provable) and `-chi_(-4)` (maximally expanding).
* **Typing.** The classifications are PROVED; the correspondences are ANALOGY. The flow numbers `2, 3, 5` against the multipliers `2, 3, 5` are NUMEROLOGY.
* **Where the analogy stops, and what it isolates.** Kuratowski's list is finite, whereas the expanding cycles are an infinite family (necklaces `~ 2^(hk)`). The finite, integral part of the obstruction list is exactly the negative-side cycle problem: are `{-1}`, `{-5,...}` and `{-17,...}` the only expanding integer cycles? That is the 3x−1 cycle conjecture.
* **Pythagorean root.** `(3,4,5)` reads as `3 + 1 = 2^2 = 5 - 1`: the trivial cycle of `3x+1` (contracting) and of `5x-1` (expanding), a ±1 sandwich around `2^2`. The arithmetic is PROVED; the reading is ANALOGY.

**6. Discrete ↔ continuous, as of wave 13.**
* **Integer spacing versus moments.**
  * Distribution-only (Haar/moment) arguments lose a square, giving `2(1-h)` (THM-4477).
  * Actual integer spacing recovers it (THM-4478), and so does periodicity through necklaces (THM-4479).
* **Christoffel words.** These discrete lines of slope `F_k -> log_3 2` are the critical cycles of the provable approximants. The five integer cycles sit at the first best approximants of `log_3 2`.
* **The moment curve** `g_3(s)` meets 1 at the integers `s = 1, 2`, which are free cycles.
* **Potentials versus cycles.** This is LP duality. Real potentials are exactly the bounded Bernoulli-type corrections, and they exist iff there is no expanding cycle.
* **Two places.** Periodic (2-adic) modifications cannot see height; archimedean ones can. Under positive drift this separates polynomial from exponential price. It is the same split as S5's `R(d)` against `R_2(d)`.

**7. Wave 14 lanes (launched).**
* **`mykk`.** Golomb–Mykkeltveit for expanding cycles: is the minimum feedback set of the expanding cycles of `B(2,k)` equal to the number of expanding necklaces? This gives the exact price of periodic deletions, and DRIFT for deletions via Mykkeltveit's sine-weight construction.
* **`drift`.** Is 5n+1 in the Haar closure of the provable sign strategies?
* **`tension`.** Three questions:
  * periodic rank functions `<=>` class (i) `<=>` Lemma P, which connects Bernoulli-boundary's obstruction;
  * Christoffel maximizers and "Christoffel rigidity" (`rho_max <= F_k` for every provable level-`k` strategy?);
  * the `nu`-duality census.
* **`peak`.** The peak-discounted price for every `q`, the `L^(1/3)` second-order term, and the pairing family.

## 2l. Waves 14–15 (2026-09-26): the price is peak-discounted, an entropy law for flips, ranks are tensions, forced charges, and free cycles

Six lanes, all audited by independent orchestrator code and pipeline reruns. Five new canon theorems.

**1. THM-4480: the price of provability is peak-discounted** ([peak note](procgen_peak_20260926_peak_discounted_price.md)).
* **The theorem.** For every odd `q`, `rho^peak_L/M*_L <= eps_L(q) <= rho^peak_L`. Here `rho^peak_L = 2^-L sum_(Bad_L) 1/max_j slope_j`, and `M*_L = O(L^3)`.
  * The upper bound sends each undecided orbit to 1 at its peak.
  * The lower bound is a single-scale integer-capacity count.
* **Consequences.**
  * **The exponent is `1 - H(log_q 2)` for every multiplier,** whatever the drift sign. For 5n+1 it is `0.013911`, although 5n+1's undecided density stays at `0.176`.
  * For `q = 3`, `rho^peak/rho_L = exp(-Theta(L^(1/3)))`, a Brownian strip cost with sharp constant `kappa_3 = 2.108` modulo Mogul'skii. So THM-4478's P1 is **negative** for arbitrary edits.
* **The pairing family (HYP-9140) stays OPEN.**
  * The pairpeak lane's coupling lemma shows that a single pairing flip re-merges with probability `1/2`. That corrects the hypothesis's rationale.
  * The private price (per-source flips) is proved peak-discounted. Conjecture R (HYP-9142, a Robin-type barrier inequality) would make it `O(L) rho^peak`.
  * Consistency between sources is the remaining obstacle.

**2. THM-4481: an entropy law for sign flips** ([drift note](procgen_drift_20260926_positive_drift_provability.md)).
* **The law.** For every stationary law of every `qn±1` sign strategy: `1 - h(pi(odd)) <= merge entropy <= pi(R ∪ R*) <= pi(odd)`. The chain makes one bit per step; the parity sequence carries only `h(pi(odd))`, and bits are destroyed only where a flip and its partner merge.
* **Consequences.**
  * `rho_max >= 0.2271` for every strategy, so **no `qn±1` with odd `q >= 23` has a bounded-lookahead-provable sign strategy at any level**.
  * Provable 5n±1 strategies need constant flip mass `0.01391`, the same constant as THM-4480's exponent, now as a stationary mass.
  * 5n+1 is in the Haar closure of the provable class only if invariant densities concentrate on the flips (HYP-9141, OPEN).
* **Three drift regimes.** Negative drift (`q = 3`): exponentially cheap (THM-4479). Middle positive band (`5..21`): constant stationary cost. `q >= 23`: impossible.

**3. THM-4482: ranks are tensions** ([tension note](procgen_tension_20260926_ranks_christoffel_duality.md)).
* **The theorem.** A strategy is provable iff a rank `a log n + h` exists, with `h` periodic or merely bounded. The least rank defect is the maximum cycle mean, by LP duality between tensions and circulations.
* **For Collatz.** The defect is `log(3/2)`, from the loop at `-1`, so **no periodic or bounded correction beats `log n`**. The Bernoulli-boundary obstruction is exactly that loop. Finite banks of 2-adic valuation counters also fail, which contains the reframe's §5.
* **Corrected premises.** The session's two guesses about `sigma_k` failed:
  * its maximizers form a positive-entropy family, of which the Christoffel word is only the balanced one;
  * "Christoffel rigidity" (`rho_max <= F_k` for provable strategies) is REFUTED. It is replaced by a numerator bound that is exact to `k = 6`.

**4. THM-4483: what a Collatz rank must look like** ([rank note](procgen_rank_20260926_two_place_lyapunov.md)).
* **Forced charges.** Any rank `a log n + bounded h + sum c v_2(n - beta)` with a height-summable bank must charge every point of the backward tree of every expanding cycle `x` with `>= a chi(x)`.
  * Inside a shadow the counter pays for the height growth.
  * At the seam the entering point must already carry the charge.
* **For Collatz.** The tree of `-1` alone has `2^(j-1)` points at depth `j`, so no such bank works. Nonnegative banks need infinite height moment in every open set of `Z_2`.
* **The converse.** Nonnegative rational-center banks with strict descent exist **iff** Collatz holds: the rank can store each stopping time in the height of a private center. So valuation-bank ranks are exactly as hard as Collatz.
* **Adaptive centers** fail at seams (resets), never inside shadows.
* **One fact.** Four rank obstructions — bounded corrections, finite banks, bounded lookahead, and height-summable banks — are one fact: the expanding cycles and their backward trees.

**5. THM-4484: free and sporadic cycles** ([sporadic note](procgen_sporadic_20260926_free_and_sporadic_cycles.md)).
* **The shift criterion.** A cycle shape of `(qy+d)/2` is free iff `(2^p - q^a) | d`.
* **Gersonides for every `q`.** `|2^p - q^a| = 1` only for `a = 1` or `3^2 - 2^3`. So the wave-13 classification, four free cycles of 3x+1 plus the sporadic `-17`, is now audited canon.
* **Belaga–Mignotte.** The gates lane's off-by-one is RESOLVED. Two long primitive cycles lay beyond its scan: `d = 14303`, least element 101, period 2155; and `d = 17021`, least element 5, period 2140. All 11 table entries now match.

**6. The picture after waves 13–15 (discrete ↔ continuous).**
* **Height is the continuous coordinate the 2-adic world cannot see.**
  * Arbitrary edits exploit it: peak discount, exponent `1 - H(log_q 2)` for every `q` (THM-4480).
  * Periodic edits cannot: polynomially sharp for `q = 3` (THM-4479), and a constant entropy cost for `5 <= q <= 21` (THM-4481).
  * Ranks need it: bounded 2-adic data fails (THM-4482), finite-mass 2-adic potentials fail (THM-4483), and the only working 2-adic potentials encode heights (stopping times).
* **Every obstruction is an expanding cycle and its backward tree.** For Collatz these are the Gersonides cycles `-1`, `-5` and the sporadic `-17`, plus infinitely many 2-adic rational ones.
* **Kuratowski–Tutte.** The owner's triple reads as a finite obstruction list forced by an identity, plus a sporadic member. Kuratowski's list is finite; here the obstruction list is infinite, and the finite integral part of it is the `3x-1` cycle conjecture.

**7. Open, with numbers.**
* HYP-9140: pairing peak price; consistency is the obstacle.
* HYP-9141: 5n+1 concentration.
* HYP-9142: Robin inequality.
* `M_7` (`17/27` or `29/46`).
* The residual 42757 vs 42765 in Belaga–Mignotte's totals.
* The mykk lane (periodic deletions and Golomb–Mykkeltveit for expanding cycles) is still running; see the addendum when it lands.

**Addendum (mykk lane landed): THM-4485, the exact price of periodic edits** ([mykk note](procgen_mykk_20260926_expanding_cycle_feedback.md)).
* **The price.** A periodic edit (`G = v_0` on residues `R`) is provable iff `R` meets every expanding cycle of `B(2,k)`. So the exact price is `FVS_c(k)/2^k`.
* **Mykkeltveit.** His solution of Golomb's conjecture (the least feedback set of `B(2,k)` is the necklace count, selected by a continuous sine weight) gives `FVS <= Z(k) - 1`.
* **Drift.** 5n+1's deletion price is `~1/k` (`k price -> 1`), 3n+1's is `2^(-(1-h)k)`: positive drift changes the rate, not the limit.
* **The Golomb analogue for thresholds is REFUTED.** Cycles of length `k+1` pack more densely than necklaces. For `log_3 2`, `FVS/N >= 1.35` on a density-0.369 set of `k`.
* **Chain and new bounds.** `delta >= FVS^odd >= FVS >= nu >= N`, with new flip bounds `delta_11 >= 58` and `delta_12 >= 95`.
* **Three kinds of modification, one table** (THM-4485 §1).
  * Arbitrary edits see height and are exponentially cheap for every `q`.
  * Periodic deletions are exponential for `q = 3` and `~1/k` for `q = 5`.
  * Periodic flips cost at least deletions, because they merge only two orbits (THM-4481).

## 2m. Wave 17 (2026-09-26, opus session `collatz-exponent-atlas-20260926`): one order in four settings (THM-4495), the general dip spectrum, and the audited AMM bound

**1. THM-4495 (PROVED, self-contained): the no-descent count has exact order `2^(hk) k^(-3/2)`.** Reading the day's incoming work, the same object appeared three times under three names: THM-4479's `Bad_k` (residues mod `2^k` with no `k`-step descent, sandwiching the strategy-cube distance `N_k <= delta_k <= |Bad_k|`), THM-4485's chain `delta_k >= FVS^odd >= FVS >= nu >= N_k`, and THM-4487's `gamma = 1` dip count. All three lanes had the exponent `1 - h` and a polynomial gap (`k^2` in THM-4479, `log^(5/2)` in THM-4487). The gap closes:
* **Identity.** `k W_k = sum_(n=1)^k B_n W_(k-n)` with `B_n = sum_(j : 3^j > 2^n) C(n, j)`: the no-descent counts are determined by binomial tails. This is Spitzer's 1956 lemma for two letters, re-proved in four elementary steps (minimum decomposition, reversal, ladder blocks, rotation averaging) using only that `log_2 3` is irrational (no ties). Checked against a ballot DP to `k = 300`, integral to `k = 3000`, and `W_1..W_11` is THM-4479's table.
* **Order.** `0.26 * 2^(hk) k^(-3/2) <= N_k <= delta_k <= |Bad_k| = W_k <= 545 * 2^(hk) k^(-3/2)`. Lower bound: the `n = k` term of the identity plus Stirling (THM-4479's `2^(hk)/(3k^2)` lost `k^(1/2)` only through `C(k,m) >= 2^(kh(m/k))/(k+1)`). Upper bound: `sum_k W_k 2^(-hk) u^k = exp(G(u))` with `G`'s coefficients `<= 2.05 n^(-3/2)`, and the split-at-`k/2` convolution lemma with the `1/m!` of the exponential.
* **Consequences.** `delta_k = Theta(2^(hk) k^(-3/2))` (THM-4479 sharpened from rate to order); the whole THM-4485 chain in one constant window (its refuted Golomb analogue, ratio `1.35`, is a point inside it); `D_(+-)(X, 1) = Theta(X^h (log X)^(-3/2))` on both sheets (THM-4487's bracket closes at its lower end); `sum_(t<24) W_t = 367698` equals THM-4487's brute-force count at `2^24` exactly, on both sheets. Reading: `1 - h` is the Chernoff rate (THM-4476 section 1.8) and `k^(-3/2)` the ballot correction; one exponent and one polynomial in four settings. Collatz is untouched.
* **Not claimed.** THM-4476's `eps` is now a small power of the logarithm (Lemma 1.4b and addendum 1.6b of the thin-divergence note: `N(X) <= K X^(h*) (log_2 X)^a` for every `a > lambda*/h* - 1/2 = 0.0138`, the recursion with `theta = (1+eta) log_2 log_2 X/(h* log_2 X)`, `h(rho) <= h* + lambda* theta`, and a geometric binomial tail in the counting lemma, which also sharpens THM-4487's upper bound to `X^(E(gamma)) (log X)^(-1/2)`); a uniform ballot factor at a moving barrier would give `o(X^(h*))`. The `gamma < 1` polynomial (expected `log^(-1/2)`) is open. The observed constants (`W_k k^(3/2) 2^(-hk)` about `10-11`) are far inside `[0.26, 545]`.

**1b. THM-4498 (PROVED): the polynomial orders of the dip spectrum.** `D_b(X, gamma) = Theta(X^(E(gamma)) (log X)^(-1/2))` for `log_4 3 < gamma < 1` (geometric binomial tail above; below, a prepended odd block and Hoeffding's inequality without replacement show that half of the arrangements with `ceil(rho t')` odd letters have a bounded final rise), `Theta(X)` for `gamma <= log_4 3` (a fifth of all words work), and `Theta(X^(h*) (log X)^(-3/2))` at `gamma = 1` (THM-4495). The polynomial along the entropy curve is `1`, then `log^(-1/2)`, then `log^(-3/2)` at Terras's endpoint. Exact DP counts to `t = 3000` show the `t^(-1/2)` window settled for `gamma >= 0.94` and converging slowly near Korec's endpoint, where the geometric ratio `(1-rho)/rho` tends to `1`.

**2. Theorem 4 of the dip-spectrum note (PROVED, post-audit): the general dip spectrum is a constrained maximum entropy.** For every Conway map `g(x) = (p_i x + q_i)/m` with `max p_i > m` and every `gamma in (0, 1]`, the exponent of `#{n <= X : g^j(n) >= n^gamma, j <= floor(log_m n)}` is `E_g(gamma) = max{H_m(pi) : sum pi_i log_m(p_i/m) >= gamma - 1}`, attained by the tilted law `pi_i ~ (p_i/m)^lambda`; `E_g = 1` iff the uniform law meets the constraint, which for `3n+-1` is exactly Korec's `log_4 3`, and `E_g(1) = 1 - I(g)` is the thin-divergence exponent. The carries are handled multiplicatively (`y_j = M_j n (1 + O(n^(-gamma) log n))` along a no-dip orbit), which removes Theorem 1's hypothesis `gamma > log_2(3/2)`. Control: the `m = 3` map `x/3, (2x+1)/3, (4x+1)/3` counted exactly below `3^15` agrees with the carry-free word model to a few parts in ten thousand, and the finite-size gap to `E_g` seen in THM-4476's control (C2) is the multinomial of a 12-letter word.

**3. THM-4494 (AMM 12592, `C* <= 197/125 < log_2 3`) independently audited: one error found and repaired.** The exact binomial factors `(A_0 + j)/(R_0 - j)` are not monotone in `N` (the residue `ceil(cN) - cN` varies along `N = 16 * 4^k`); the majorant factors `(aN + j)/(bN - j)` are, and the repaired bottom margin is `109.592` (was `110.967`). Conclusion intact; MISTAKES entry; the header's output hashes, which matched no committed file, recomputed.

**4. Cross-lane reading.**
* THM-4484's sporadic cycle `-17` with gap `139 = 3^7 - 2^11` is the atlas's `139` row (shape `(11, 7)`, expansion `2187/2048`); the three "small-number coincidences" `139`, `2187/2048`, `1093 = (3^7 - 1)/2` are one pair of powers.
* The crossroads correction of the tournament note is accepted: interior insertion gaps are `01` transitions, `1 + #10` counts all slots including the ends, and the metric reading of the seed is not equivalent to the tournament axiom.
* The procgen reflection ("height is the coordinate the 2-adic world cannot see") describes exactly what THM-4476/4487/4495 use: the real place (`n` versus `X`) against the 2-adic word count. Those theorems show how far that combination reaches with the free window `floor(log_2 n)`, and Theorem 4's remark (iii) says why the window cannot be lengthened without new input (the residue word of length `tau t`, `tau > 1`, is not equidistributed among `n < 2^(t+1)`): the same obstruction the reflection names, seen from the counting side.

## 2b. The approach deck: every approach generated or considered, with its disposition

| # | approach (lens + mechanism) | barrier verdict | probe run | outcome |
|---|---|---|---|---|
| 1 | residue certificates for Collatz (Terras) | blind to SHEET, DEFECT; needs THIN | exceptional counts | the dimension `0.95` set remains; not a route alone |
| 2 | relax by branch choice (graph `E`), certificates + escapes | live for E-SCC | profiler, DP/DFS, Lean | Q2 reduced mod 27 (Lean), verified to `10^18`; Q2 follows from `X_min` (PROVED implication); gap: `1/2` chains (endgame: ternary digits of `2^K`) |
| 3 | partial choice `E_S` | diagnostic | profiler over `S` | `6 mod 8` does most of the work at finite levels; fingerprint led by `54 mod 64`. But its exceptional set keeps positive dimension (`>= 0.0536`, PROVED): the orbit of `-1` is a 16-point trap of rate `2/3` (wave 7) |
| 4 | additive / sign choice | diagnostic | zoo | empty exceptional set; one-player trivial |
| 5 | two-player sign choice (Althöfer game = Conway's Beans-Don't-Talk, Guy Problem 42) | live | game lane | no draws below `2^32` (FINITE-EXACT); ray/exit-parity law PROVED; P-density about `0.48` |
| 6 | undirected Collatz (grand-orbit moves) | equivalent to Collatz | BFS to `10^6` | no collapse (negative control) |
| 7 | sideways moves `n~4n+1` | inside #6 | argued | no power beyond Terras |
| 8 | Applegate--Lagarias see-saw transplanted to `E` | needs cheap escapes | cost bounds, endgame lane | escape prices `2c(k-1)/3` in `(1,2)` (infimum 1), but the price per digit reaches `0.1144` (sharp budget), so the see-saw does not close as is. Corrected from "`>=32/27`" |
| 9 | 1-escape via loops through `1` | live (Q2 near `1`) | loops lane | PROVED exact and optimal (ratio `>3/2`); loops exist to `s=6000`; HYP-9122 reduced to record denominators; no finite family (PROVED) |
| 10 | Hecke / `X_0(11)` density lens | DEFECT, SHEET | exact density | echo density `2/15`; statistic only |
| 11 | Catalan unit-gap clocks | cycle half only | table | inherited; not a divergence tool |
| 12 | fixed sheet gluing `3n+sgn(n)` | no new freedom | zoo | gains nothing |
| 13 | transversality: p-adic irrationality of Bernstein series | = Lagarias's Periodicity Conjecture (wave 3) | transversality foundry | named; proved on SUB, BCD, STURM; OPEN on the HARD class (C2 smallest) |
| 14 | sound certificate searches (rewriting/automata) | live (divergence half) | literature | YAH prize conjectures fail on negatives (atlas) |
| 15 | sign-specific order laws | the only sign-aware class | order-laws lane | nothing beyond the sign law among 6,016 generated laws (FINITE-EXACT to `10^7`); direction lemma PROVED; windows of at most 12 odd steps sheet-blind |
| 16 | exceptional dimension of `E` | structural | dimension lane (exact threads to `2^64`, `3^42`) | OPEN; evidence favours 0; infinite (Theorem F PROVED); 382 hostile `-p/3^j` certified |
| 17 | sibling ladder (`qx+1`, F2[x], Mahler, Erdős) | calibration | ladder lane | dimensions and measures PROVED; F2[x] set `{1}`; choice games drift-blind (threshold near `q=10`) |
| 18 | Tao-type Fourier / renewal | blind to SHEET, DEFECT, INTEGRAL | literature | GGM 2025: PROVED for `3N-1` (sheet-blind as a theorem) |
| 19 | bounded-modulus Lyapunov potentials | none in the model | inherited | REFUTED (earlier sessions) |
| 20 | carry anti-concentration mod `2^K-3^L` | cycle half | heuristic count | live for cycles; not run |
| 21 | Baker / continued fractions | cycle half | literature | Hercher: no m-cycles, `m<=91` (CITED) |
| 22 | functional equations (Berg--Meinardus) | DEFECT | none | typed only |
| 23 | measure rigidity (`x2 x3`) | DEFECT, INTEGRAL | none | blocked |
| 24 | E-cycle covering (Le--Smith Conj. 1, 2) | relaxation / cycle half | our verifications | Conj. 1 holds below `10^18`; Conj. 2 is equivalent to no positive Collatz cycle |
| 25 | periodic-approximant Liouville (Theorem R) | HARD-blind (reaches `Dio > eta` only), DRIFT-blind | Sturmian test, audit | **Theorem S PROVED** (no rational has an eventually Sturmian parity vector, `3x+r`, all slopes) |
| 26 | capacity plus ordered carry (bounded critical discrepancy) | HARD-blind (critical slope only) | inherited, Prop B | PROVED for positive integers (in-house) and rationals (Prop B, via an in-house sketch) |
| 27 | Q1 mirror: loops through `-1`, exit prices | relaxation, Q1 | mirror lane | exit costs `3^eta > 1` (PROVED); mirror of HYP-9122 to `K = 4000`; endgame = low binary digits of `3^A u` |
| 28 | forward-backward word duality | — | mirror lane | REFUTED (straddle identity PROVED; Gelfond–Schneider lemma) |
| 29 | supercritical strips (C2) | a HARD instance (HYP-9123) | wave-4 lane | OPEN; equivalent to a coupled Z-number statement; PROVED on `Dio > mu` and on square-swap words |
| 30 | Diophantine-exponent criterion (Theorem D) | reaches `Dio > eta` only | HARD lane | PROVED: Sturmian and quasi-Sturmian words for every map with `eta < 2.50994` (all slopes of `5x+1`) |
| 31 | q-series Padé (2-adic Tschakaloff) | needs a q-difference equation | HARD lane | **Theorem Y PROVED** (square-swap words, `Dio = 1`, zero entropy) |
| 32 | cubic 2-adic theta values | the smallest open instance | HARD lane | the cube-swap word `Y3`: OPEN (HYP-9127) |
| 33 | single-orbit density counting (Terras count + landing pigeonhole; opus S5/S6) | divergence half: a constraint, not an exclusion; SHEET-blind count with sign-specific consequence | THM-4476, THM-4487 | every non-periodic orbit has `O(X^(0.95+eps))` points below `X`, reciprocal sums converge, `R(d) < n` strictly; the count is sharp for the lemma (`X^(h(rho)+o(1))`), so `0.95` is the method's floor |

## 3. The snippet, dispatched

| pasted claim | verdict | where |
|---|---|---|
| root cycle length 3 versus 4 "breaks sheet symmetry" | TRUE as the Catalan unit-gap table: plus `(2,1)`, minus `(1,1),(3,2)`. It is cycle-half content only; class-level data cannot see the sheet because all odd `b` are conjugate | lead reflection 1; choice ladder §5b |
| extra freedom "matches" loops at `-7/4, -29/16` | SCOPE: shared numerals only | zsigmondy lane |
| Hecke `b_r`, the 4-block law, density `1/15` | TRUE for the level-11 eigenform. Its exact Collatz echo is that odd `n` with `v_2(3n+1)=3 mod 4` have density `2/15`, a density statistic blocked by DEFECT and SHEET | level11_short; this session |
| gluings `3n+-sgn(n)` | TRUE: two reflected copies. A fixed gluing gains nothing; a *choice* of sign trivializes descent (additive zoo), and a two-player choice is Althöfer's open game | choice ladder |
| `1,5,17` system | three cycles, of lengths 1, 2, 7 | inherited |
| `36=18+18`, `C(6,2)+6=21`, Fano orientations | arithmetic true, no map | level11_short |
| `11_B x=Bx+x` | TRUE | level11_short |
| Lean `forced_zero_density_limit` | FALSE: `sum (1/16)^j=16/15`, and `1/16=0` in `N` | this session |
| Lean `shift_memory_rule` | FALSE at `B=3, x=1` | this session |
| Lean `descent_condition` | the correct descent certificate | inherited |

## 4. Where the missing insight lurks

The instrument locates it. Choice (graph `E`) shrinks the `0.95`-dimensional
exceptional set by three orders of magnitude at every tested level, into a
far thinner set whose rational points lie over the other prime. Its exact
dimension is open. The
collapse is driven by freedom at rising-run entries. Collatz has no such
freedom, and its equivalent choiceful reformulation (the undirected game)
does not supply it. A proof of the divergence half must therefore
substitute for choice. It must see the drift, handle each orbit exactly,
use a non-uniform arithmetic input, and cope with a `0.95`-dimensional
exceptional set, which no known escape construction covers; the
`E`-relaxation's much thinner set is already hard (see item 3 of section 5). **Corrected emphasis (after the sibling-ladder lane).** What choice buys is
*generic*. The `5n+1` relaxation collapses as well (its mass falls to about
`4e-5` at `2^64`), and choice games switch only near `q=10`, so choice
relaxations are blind to both the sign and the drift. A relaxation that
becomes provable by choice therefore discards exactly the two features
that make Collatz special. The missing ingredient must be specific to
`3n+1` in both respects: it must see `log 3<2 log 2` and the sign of
`2^K-3^L` together. Choice is a diagnostic of *where* the difficulty
sits (rising runs), not a template for the proof.

**Where the sign can act (order-laws lane, PROVED direction lemma).**
A window contradicts its parity word's order prediction only at a *gate
crossing*: a decay window going up on the plus sheet, a growth window
going down on the minus sheet. Windows of at most 12 odd steps are
provably sheet-blind. The first minus crossings sit at clock `19/12`
(the transient `165->163`), and every plus crossing found rides the
orbit of `27`. So a sign-aware argument must act on windows whose clock
`K/L` is near a convergent of `log_2 3`: the inherited Pillai clocks,
which the snippet calls "unit-gap clocks". That is exactly where a
Baker-type Diophantine input enters, matching the transversality
target in section 5.

**Both halves end in the same kind of statement.** The divergence half
of Collatz points to p-adic irrationality of the Bernstein series (the
2-adic digits of an integer versus a dense odd-step set). The relaxed
backward half, after every certificate and escape, points to the base-3
digits of `2^K` (loops lane). The missing insight is a *2-versus-3 digit
transversality* theorem. This is the fringe idea the session kept
returning to from different directions (Mahler `3/2`, Erdős `2^n`, the
S596 two-block question).

**Refined by wave 3 (section 2c).**
1. **The divergence half is named.** It is Lagarias's Periodicity
   Conjecture on the positive integers. It is proved on three word
   classes:
   * subcritical words;
   * bounded critical discrepancy;
   * Sturmian words (Theorem S, this session).

   Every proved every-orbit technique works either by *growth/capacity*
   (the first two classes) or by *repetition*, i.e. zero entropy (the
   third). What is missing is a mechanism for **supercritical,
   positive-entropy, non-repetitive** words: the HARD class, whose
   smallest open instance is C2.
   **Refined by wave 4.**
   * Entropy is not the dividing line. Periodic approximants reach
     exactly `Dio > eta` (Theorem D).
   * The zero-entropy square-swap word (`Dio = 1`) falls to a Padé
     argument from a q-difference equation (Theorem Y).
   * The first word we could not settle is the cube-swap word `Y3`,
     whose Bernstein number is a cubic 2-adic theta value (HYP-9127).
   * So the missing mechanism must handle words with `Dio <= eta` and no
     functional equation behind them. The simplest test case is a single
     explicit 2-adic number.
2. **Both halves of the relaxed problem end in low `p`-adic digits, and
   the relaxation renormalizes rather than removes the difficulty.**
   * Q2's hostile chains are governed by `v_3(2^K w - h)`, via the kappa
     formula.
   * Q1's chains are governed by `v_2(3^A u + 1)`, via
     `v_2(A - lambda(u))`.
   * Each chain is again a Collatz-type map with memory. Its termination
     is a no-divergence statement of the same HARD type.
   * The Erdős attribution above is corrected. Erdős's problem concerns
     the **top** ternary digits, and every 3-adic-window version of it is
     FALSE. Both E-SCC endgames are pointwise versions of Dupuy–Weirich's
     *averaged* low-digit equidistribution.
3. **One missing mechanism for the whole family.** In foundry v4, every
   all-orbits target has the same single real unblocked mechanism type:
   Collatz and `3n-1` divergence, PC, both halves of E-SCC, Mahler,
   Erdős, and the `5n+1` existence question. So the "missing insight" is
   one kind of statement: a non-integrality or irrationality statement,
   2-adic or 3-adic, that reaches positive-entropy digit words.

The pasted snippet's themes map to exact places: the sign enters
as the sign of `2^K-3^L`, the Catalan clocks belong to the cycle half, and
the trunk `(4^i-1)/3` is the exit set of the hardest relaxed obstruction.

## 5. Frontier and next probes

1. Q2 near `1/2`: an amortized see-saw. The chain recursion
   `w_{t+1}=(2^(K_t+3)w_t-1)/3^(j_{t+1})` governs repeated hostile
   landings. **Wave-3 update ([endgame lane](collatz_procgen_20260922_q2_endgame.md)).**
   * This recursion is correct for its route, but it lands twice as high
     whenever `K0(k-1) = K0(k-2) + 1`.
   * The optimal (canonical) recursion is
     `w_(t+1) = (2^(K0(k_t - 1)+1) w_t - 1)/3^(k_(t+1))`.
   * Every 3-adic unit lies on the thread of `1` or of `1/2`. Every
     dyadic thread's link factors through the single expensive link, the
     `1/2`-transfer with `k >= 3`.
   * The budget `exp(0.1144 D)` is sharp.
   * The two remaining statements, `X_min` and `X_T`, are integer-avoidance
     statements of Collatz/Erdős type.

   **Why this is the real core (analysis).** Applegate--Lagarias close
   because their escape cost `1+2^(-j)` tends to `1`, so a chain of
   escapes costs a bounded product. In `E`, every escape from the `1/2`
   neighbourhood costs more than 1: the canonical price `2c(k-1)/3` has
   infimum 1, but the price per digit reaches `0.1144` ("at least `32/27`"
   was corrected on 2026-09-23). The post-escape value
   `x=1+3*2^(K-1)w (mod 3^D)` depends only on the loop's total halving
   count `K`. Multi-move exits (factor `3^(-t)`) need `x` in specific
   classes mod `3^(t+1)`, but a loop ratio below `27/8` (or `81/8`)
   leaves only one to three admissible values of `K`, too few to steer
   `2^(K-1)w` against an adversarial `w` modulo `ord_(3^t)(2)`. So chains
   of hostile landings can accumulate cost up to `exp(0.1144 D)` over `D`
   digits. That bound is sharp, and chains can outrun the digits of `m`.
   Bounding their length for integers is a Collatz-type
   digit-propagation question. The relaxation is far thinner than
   Collatz but keeps a genuine Collatz core.
2. Loops through `1` of every length with bounded ratio (HYP-9122). This is
   a carry-covering statement at the convergent clocks. It is a concrete
   instance of the S596 "two-block" analogy
   ([reflection](../../07-reflections/lrc-collatz-the-same-two-block-question-s596.md)):
   which ratios `2^K/3^s` can legal words through `1` realize at every
   length? That analogy was left untested there; here it has a precise
   finite form (every length `<=40` realized with ratio in `[1.517, 2.96]`).
3. The exact dimension of `Bad_inf(E)`. **Update: the dimension lane's
   exact threads to `2^64` favour dimension 0** (a power law about
   `m^1.6`, a countable set). It PROVED the family
   `-1-2^i/3^alpha(i)` hostile (Theorem F, so the set is infinite) and
   showed the "credit" below is a fixed slack `1/2-(|x|-1)` that every
   block spends, which undercuts the heuristic. The original heuristic,
   kept for provenance: it favoured positive dimension. A stretch near `-1` banks a
   multiplicative credit of at least `3/2` (the PROVED bound on paths from
   `-1`). The next segment then only has to avoid descending by more
   than that credit, a weaker and larger condition. Concatenating such
   blocks should give uncountably many hostile points, consistent with
   the observed growth of about `0.1` bit per level. If this holds, an
   E-SCC proof needs parametrized escape families, as the atlas notes,
   not finitely many lemmas.
4. Althöfer's `3n+-1` game (prize open to 2037). This is Conway's
   *Beans-Don't-Talk* and Guy's 1996 Problem 42 (CITED via the
   [game lane](collatz_procgen_20260922_althofer_game.md)). A public
   repository under the claimant's name labels its own proof OPEN.
   * FINITE-EXACT: **no draws below `2^32`** (every odd start has finite
     remoteness; heights up to `1.31*10^13`), matching OEIS A005694--A005698.
   * PROVED: ascending moves split the odd numbers into increasing rays,
     and a position is N iff the first point on its ray whose descending
     move reaches a P-position has even index.
   * The P-density is about `0.48`. **CORRECTION:** my earlier "about 28%
     P" was an artifact of the value cap.
   * Observed: the phase `frac(log_2 n)` predicts the value 84--100% of
     the time. Negation equivariance explains the `r<->-r` symmetry.
5. The divergence half of Collatz: the four-property job description above.
   **Wave 3 names it:** Lagarias's Periodicity Conjecture on the positive
   integers (T1, Bernstein 1994). It is proved on SUB, BCD and STURM
   (section 2c), and by Theorem D and Theorem Y also on every
   `Dio > eta` word and on the square-swap words. The open targets are C2
   (HYP-9123, a coupled Z-number problem) and, smallest of all, the
   cube-swap word `Y3` (HYP-9127): is
   `sum_(k>=1) (2^10/3^9)^(k^3)` irrational in `Q_2`?
   The corrected foundry leaves exactly two unblocked mechanism types:
   sound certificate searches, and transversality with a Diophantine
   input. The concrete transversality target (a restatement, Bernstein's
   formula, CITED via the atlas) is the following. A positive integer
   `n=-sum_l 2^(d_l)/3^l` in `Z_2`, with `d_l` the times of odd steps, can
   diverge only if `d_l` grows no faster than about `l log_2 3`, in a
   non-periodic way. The needed input is therefore a *p-adic irrationality*
   statement: such a series is never a positive integer. This is the
   Mahler--Baker-type ingredient the job description calls for. Periodic
   `d_l` give exactly the rational cycle points `B/(2^K-3^L)`.
