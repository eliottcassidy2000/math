# Procedurally generated approaches to Collatz: the choice ladder, the missing no-divergence mechanism, and a relaxed problem that almost closes

**Status: SYNTHESIS of session `collatz-procgen-20260922` (mac-mini,
2026-09-22). PROVED statements are those proved in the lane notes (hand
proofs and one core-Lean theorem). FINITE-EXACT statements have stated
bounds and, where marked, two independent code paths. CITED statements
come from primary sources read in the barrier-atlas lane. Typology entries
are modelling judgments. Collatz, E-SCC (Q1 and Q2), HYP-9120--9122 and
Althöfer's game remain OPEN.**

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
| `E_S`, `S={6 mod 8}` | extra `3n+1` only at rising-run entries | `3,238` mod `2^22` (Collatz `93,222`) |
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
  * **verified below `7.87*10^17`** by exact 3-adic certificate threads to depth 41 (dimension lane; `2.02*10^13` by lane one's DFS);
  * the neighbourhood of `1` is handled by a PROVED escape lemma together
    with loops through `1` (every depth up to 41);
  * the neighbourhood of `1/2` costs at least `32/27` to escape (PROVED),
    and its exits are exactly the Collatz trunk `(4^i-1)/3`.
* Consequence: Le--Smith's Conjecture 1 (every `n` prime to `3` lies on an
  `E`-cycle) holds below `7.87*10^17`.
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
  * The see-saw does not close. **The endgame of Q2 is the base-3 digits of
    `2^K` (Erdős ternary type).**

## 2b. The approach deck: every approach generated or considered, with its disposition

| # | approach (lens + mechanism) | barrier verdict | probe run | outcome |
|---|---|---|---|---|
| 1 | residue certificates for Collatz (Terras) | blind to SHEET, DEFECT; needs THIN | exceptional counts | the dimension `0.95` set remains; not a route alone |
| 2 | relax by branch choice (graph `E`), certificates + escapes | live for E-SCC | profiler, DP/DFS, Lean | Q2 reduced mod 27 (Lean), verified to `7.87e17`; gap: `1/2` chains (endgame: ternary digits of `2^K`) |
| 3 | partial choice `E_S` | diagnostic | profiler over `S` | `6 mod 8` does most of the work; fingerprint led by `54 mod 64` |
| 4 | additive / sign choice | diagnostic | zoo | empty exceptional set; one-player trivial |
| 5 | two-player sign choice (Althöfer game = Conway's Beans-Don't-Talk, Guy Problem 42) | live | game lane | no draws below `2^32` (FINITE-EXACT); ray/exit-parity law PROVED; P-density about `0.48` |
| 6 | undirected Collatz (grand-orbit moves) | equivalent to Collatz | BFS to `10^6` | no collapse (negative control) |
| 7 | sideways moves `n~4n+1` | inside #6 | argued | no power beyond Terras |
| 8 | Applegate--Lagarias see-saw transplanted to `E` | needs cheap escapes | cost bounds | escape cost `>=32/27`, not `1+2^(-j)`: does not close as is |
| 9 | 1-escape via loops through `1` | live (Q2 near `1`) | loops lane | PROVED exact and optimal (ratio `>3/2`); loops exist to `s=6000`; HYP-9122 reduced to record denominators; no finite family (PROVED) |
| 10 | Hecke / `X_0(11)` density lens | DEFECT, SHEET | exact density | echo density `2/15`; statistic only |
| 11 | Catalan unit-gap clocks | cycle half only | table | inherited; not a divergence tool |
| 12 | fixed sheet gluing `3n+sgn(n)` | no new freedom | zoo | gains nothing |
| 13 | transversality: p-adic irrationality of Bernstein series | live (divergence half) | formulation only | the named missing ingredient |
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
| 24 | E-cycle covering (Le--Smith Conj. 1, 2) | relaxation / cycle half | our verifications | Conj. 1 holds below `7.87e17`; Conj. 2 is equivalent to no positive Collatz cycle |

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
digits of `2^K` (loops lane): Erdős's ternary problem, which the foundry
listed as a sibling with dimension `log_3 2`. The missing insight is a
*2-versus-3 digit transversality* theorem. This is the fringe idea the
session kept returning to from different directions (Mahler `3/2`,
Erdős `2^n`, the S596 two-block question).

The pasted snippet's themes map to exact places: the sign enters
as the sign of `2^K-3^L`, the Catalan clocks belong to the cycle half, and
the trunk `(4^i-1)/3` is the exit set of the hardest relaxed obstruction.

## 5. Frontier and next probes

1. Q2 near `1/2`: an amortized see-saw. The chain recursion
   `w_{t+1}=(2^(K_t+3)w_t-1)/3^(j_{t+1})` governs repeated hostile
   landings.

   **Why this is the real core (analysis).** Applegate--Lagarias close
   because their escape cost `1+2^(-j)` tends to `1`, so a chain of
   escapes costs a bounded product. In `E`, every escape from the `1/2`
   neighbourhood costs at least `32/27`. The post-escape value
   `x=1+3*2^(K-1)w (mod 3^D)` depends only on the loop's total halving
   count `K`. Multi-move exits (factor `3^(-t)`) need `x` in specific
   classes mod `3^(t+1)`, but a loop ratio below `27/8` (or `81/8`)
   leaves only one to three admissible values of `K`, too few to steer
   `2^(K-1)w` against an adversarial `w` modulo `ord_(3^t)(2)`. So chains
   of hostile landings can accumulate cost at least `32/27` per link.
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
