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
