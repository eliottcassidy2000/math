---
source: collatz-procgen-20260922 (mac-mini)
status: REFLECTION (provenance, not truth). The results live in 05-knowledge/results/collatz_procgen_20260922_*.md; statuses there govern.
tags: [collatz, procedural-generation, exceptional-set, choice, relaxation-ladder, barriers, e-scc, applegate-lagarias, sheets, trunk, althofer]
---

# The missing Collatz insight lives where choice would be

**Prompt (user):** "procedurally generate new possible approaches for open
problems like collatz ... test your hypotheses and explore past work
extensively ... surprisingly unexpected fringe ideas we have brushed" (with
a pasted snippet on unit-gap clocks, signed sheets and a level-11 `1/15`
density).

## What the procedure was

Past sessions had audited pasted blueprints one claim at a time. This
session instead built a *generator* and an *instrument*:

* **The generator** (the foundry) types every approach by the controls it
  is blind to (`3n-1` sheet, `5n+1` drift, planted defects, rational
  2-adic cycles, undecidability) and by whether it needs a thin exceptional
  set. It turns "is this idea promising?" into "which control must this
  idea fail on, and does it?"
* **The instrument** (the exceptional-set profiler) counts residue classes
  with no descending certificate. It applies to Collatz, to both sheets, to
  `5n+1`, and to any relaxation with choice.

Neither is deep. What they buy is comparability: a hundred-odd
generated variants can be put on one scale.

## What the scale showed

Collatz's non-descending set has dimension `h(log_3 2)=0.95`. Three kinds
of freedom move it:

* branch choice (graph `E`) makes it three orders of magnitude thinner;
* sign choice (`3n+-1`) makes it empty;
* wild multipliers (Applegate--Lagarias) leave a single class.

Partial choice pins the effect to the entries of rising runs at finite levels. *Wave 7 correction:* the partial choice alone keeps a positive-dimensional exceptional set (`>= 0.0536`; HYP-9121 REFUTED), because the orbit of `-1` stays in a 16-point trap of rate `2/3` (`6 mod 8`;
the greedy fingerprint leads with `54 mod 64`). The equivalent
choiceful reformulation of Collatz, the undirected game, does *not*
collapse.

So the "missing insight" has a job description:

* **Where.** It must do, for positive integers and without choice, what
  choice does in the relaxations: control rising runs.
* **Which half.** It concerns the divergence half. The cycle half already
  has live mechanisms (Baker-type), and the literature lane confirms that
  no published all-orbits mechanism overcomes both the sheet and the
  drift.
* **Where the sign lives.** In the sign of `2^K-3^L`. Positive plus-sheet
  cycles contract multiplicatively, so class-level certificates never see
  them. Positive minus-sheet cycles expand, and their minima are hostile.

**Correction within the session.** I first read one level (`2^20`) of
the `5n+1` choice game as a "drift barrier". The sibling-ladder lane
showed that the fraction keeps falling (about `4e-5` at `2^64`). Choice
games tip only near `q=10`; positivity is proved for `q>=41`. So choice
is blind to the drift as well as the sign, and the job description
above must be read as "do what choice does, *but specifically for
`3n+1`*". The ledger entry is MISTAKES 2026-09-22.

## How the fringe ideas connected

* **Signed sheets (the snippet):** all class-level data are sheet-blind,
  because odd `b` are 2-adically conjugate. The negative cycles are exactly
  the hostile points of the relaxed plus-sheet problem.
* **The trunk `(4^i-1)/3`** (earlier sessions' Cipolla towers): it is the
  exit set of the 3-adic hostile point `1/2`.
* **Applegate--Lagarias's lone class `-1 mod 2^j`:** the one-point end of
  the ladder. Their near-free escape `1+2^(-j)` is exactly what `E` lacks
  near `1/2`. Every escape there costs more than 1: `2c(k-1)/3`, whose
  infimum is 1 but whose price per digit reaches `0.1144`. The earlier
  "at least `32/27`" was the floor of one route and is corrected. That is
  why E-SCC keeps a Collatz-type core.
* **Althöfer's `3n+-1` game:** the two-player version of the sign
  choice. The one-player version is trivial; the two-player version is
  open, with a prize.
* **The `1/15` density:** its exact Collatz echo, density `2/15` among odd
  `n`, is a density statistic, blocked by DEFECT and SHEET.

## What surprised me

1. **E is already in the literature** (Le--Smith's Loosened Collatz Graph).
   The reachability questions, Q1 and Q2, were not. A cheap search before
   the method saved a false novelty claim.
2. **Q2's hard core collapsed to a mod-27 statement.** Every class but two
   descends within two moves. It is short enough to kernel-check in core
   Lean in a second.
3. **Chained hostile landings are typical** (31% of shortest descents from
   the `1/2` class). A thin exceptional set does not mean an easy
   all-integers statement.
4. **My claims about the dimension of `Bad_inf(E)` moved twice.** I first
   said "dimension 0". I withdrew that when local minima to `2^40` fit
   about `0.1` per level. I then offered a "credit construction" for
   positive dimension. The dimension lane's exact threads to `2^64`
   favour `0` again (a power law about `m^1.6`). It proved the set
   infinite (Theorem F) and showed the credit is a fixed slack that every
   block spends. The question stays OPEN; the lesson is not to fit a
   dimension from a window of 24 levels.
5. **Two finite-truncation artifacts were caught by the lanes:** the
   `5n+1` "drift barrier" (one level) and the "28% P-positions" (a
   capped statistic). Both are logged in MISTAKES 2026-09-22.
6. **Both hard cores end in 2-versus-3 digits.** The divergence half
   points to p-adic irrationality of the Bernstein series, and the
   relaxed Q2 points to the base-3 digits of `2^K` (Erdős). *Wave 3
   corrected both halves of that sentence.*
   * The first is Lagarias's Periodicity Conjecture.
   * The second reads the **low** 3-adic digits, not Erdős's top digits.
     Q1 has the mirror image, the low binary digits of `3^A u`.

## Candidate card (not promoted)

**"Measure the exceptional set along a relaxation ladder."**
* Trigger: an all-orbits statement whose certificates are residue-class
  based.
* Action: compute exceptional-set growth for the system, a partial-choice
  family and a full relaxation. Read the missing freedom off the drop, and
  treat the hostile points of the relaxation as its escape obligations.
* Counterindication: equivalent reformulations with choice can keep the
  full dimension (the undirected game). A thin exceptional set can still
  carry an unbounded chain of escapes (Q2 near `1/2`). The ladder
  cannot separate `3n+1` from `5n+1`: choice games tip only near `q=10`.
  Use the ladder to locate difficulty, never as evidence that the relaxed
  statement carries the original's special features.
* Evidence: this session; Applegate--Lagarias's single-class obstruction;
  Caraiani's `5x+1` semigroup (per the atlas).

## Wave 3 (2026-09-23): the relaxation renormalizes the difficulty

Three lanes ran after the first close-out. They changed the picture in
three ways.

1. **The missing insight now has a name and a smallest instance.**
   * The divergence half is Lagarias's Periodicity Conjecture on the
     positive integers.
   * Every proved every-orbit result works either by growth and capacity
     (subcritical words; bounded critical discrepancy, an in-house
     theorem from 2026-09-21) or by repetition (zero entropy).
   * The transversality lane added the repetition case for Sturmian
     words. **Theorem S**: no rational has an eventually Sturmian parity
     vector. It is a 2-adic Liouville argument whose approximants are the
     word's own periodic extensions. I audited it by a separate code
     path.
   * What no technique reaches is the HARD class: supercritical,
     positive-entropy, non-repetitive words. The smallest open instance
     is bounded discrepancy around a supercritical slope.
2. **The relaxation does not escape Collatz. It renormalizes it.**
   * In E-SCC both base points cost more than 1 to escape:
     * `1/2` for Q2, at `2^eps`;
     * `-1` for Q1, at `3^eta`, which the mirror lane PROVED.
   * So hostile landings chain.
   * Each chain is a new Collatz-type map with memory, read off the low
     `p`-adic digits of `2^K w` or `3^A u`.
   * Excluding infinite chains is again a no-divergence statement of the
     same HARD type.
   * The foundry (v4) now gives every all-orbits target in the family the
     same single missing mechanism: Collatz, `3n-1`, the Periodicity
     Conjecture, both E-SCC halves, Mahler and Erdős.
3. **The duality I hoped for was a straddle.** The shared numerators
   `793585` and `419868489953` looked like a forward–backward duality. In
   fact each pairs a hostile point with a descending one. An exact identity
   `(|x|-3/2)(y-3/2) = -(rho-1)^2/(2 rho)` puts exactly one partner of
   each pair above `3/2`. The earlier observation had compared a certified
   census with an uncertified one.

**What this says about where the insight has been lurking.** It has not
been hiding in a clever relaxation, a symmetry or a choice. Each of those
moves the same obstruction somewhere else. It lurks in one kind of
statement that the literature and this repository have circled from many
sides without proving: a 2-adic or 3-adic non-integrality statement strong
enough to reach positive-entropy digit words. The Sturmian theorem shows the
method that works on the zero-entropy edge. C2 is where it has to be
pushed.

## Wave 4 (2026-09-23): the smallest open instance is one explicit number

The HARD-class lane corrected my wave-3 description in a useful direction.
1. **Entropy was the wrong axis.**
   * Periodic approximants settle exactly the words whose Diophantine
     exponent beats the map's height rate (Theorem D). Via Bugeaud–Kim,
     that covers every Sturmian word under every slope of `5x+1` too.
   * The zero-entropy square-swap word has exponent 1, so no Liouville
     argument touches it. Its Bernstein number is nonetheless a 2-adic
     theta value, and a 2-adic transcription of Zudilin's
     Tschakaloff–Padé construction proves it irrational (Theorem Y). I
     re-derived the linear form and checked the identity to `2^6000`.
2. **The frontier collapses to a single explicit 2-adic number.** For the
   cube-swap word both mechanisms fail. There is no repetition, and no
   first-order q-difference equation. Everything reduces to one question:
   is `sum_k (2^10/3^9)^(k^3)` irrational in `Q_2` (HYP-9127)?

   That is where the procedure ends up. It started with "procedurally
   generate approaches to Collatz" and went through relaxations, ladders,
   typed barriers, mirrors and a transversality catalogue. The missing
   mechanism turned out to be a p-adic transcendence question about lacunary
   q-series. It is the p-adic cousin of problems that are open even over
   `R` (cubic theta values). The pasted snippet's "Catalan clocks" and
   "theta-like densities" turn out to be closer to the core than its
   numerology suggested: the core is the arithmetic of `2^a/3^b` inside
   lacunary series.

## Wave 5 (2026-09-23): what the outside results did

The owner brought a new batch of papers and puzzles. The procedure typed
each one against the frontier.
* **The best-matched outside result is a Lean proof.** The accepted proof
  of Erdős 1062(ii) runs on exactly our engine: small nonzero S-unit forms
  at `{inf, 2, 3}`. The only thing between its series and a Collatz parity
  word is the height paid per bit. That one exchange rate, `eta`, is also
  why the cube-swap number resists.
* **The moonshine puzzle was not a detour.** The square-swap Bernstein
  number is a theta value on a 2-adic Tate curve, and the p-adic
  Mahler–Manin theorem then gives a transcendence-type corollary. The
  original snippet's level-11 recurrence is Mathieu moonshine: the eta
  product of `M_24`'s order-11 elements.
* **The "hostile 1/2 is RH's 1/2" hope failed cleanly.** The trunk really
  is a zeta coordinate, but of the 3-adic zeta, which has no zeros. The
  two halves are provably different points.
* **The fair-coin question produced a theorem.** A classical capacity
  theorem (Pólya 1928) applied after folding `p <-> 1-p` gives the first
  uniform gap for AMM 12592, `C* >= 11/8`. The repository's 44-theorem
  corpus had not found one. Folding by a symmetry to turn two-point
  integrality into one-point integrality is the same move as the session's
  E-game quotient, and it is a reusable card.

## Wave 8 (2026-09-24): every faithful re-expression is blind somewhere

The owner asked what each feature of Collatz *represents*, and offered four
re-expressions: the backward tree built by doubling and by `(2n-1)/3`, its
mod-192 automaton, a balanced up/down pairing, and a graceful-tree analogy.
All four turned out exact. The lanes proved each one faithful:
* the tree is equivalent to Collatz;
* the automaton is sound;
* the pairing identity holds for every pair;
* the two edge classes are perfect difference systems.

The same lanes also proved exactly where each one goes blind.

* **The sign.** `x -> -x` carries every residue, density, integrality and
  size fact about `3n+1` onto `3n-1`. So the tree and its automaton cannot
  tell the positive half-line from the negative one, where `3n+1` has
  three more cycles. The sign enters a proof only through the sign law.
* **Defects.** The pairing's balance and gracefulness hold for a whole
  family of maps. That family contains trees, maps with extra cycles, and a
  map that differs from Collatz on a density-zero set of pairs and
  diverges.

So the owner's "orderly, perfectly arranged" picture is right about Collatz's
local structure, and that local structure is provably not where the
conjecture is decided. This is the same lesson as the choice ladder from the
other side. Choice changes the exceptional set without changing the problem.
The pairing ladder keeps every statistic and changes the truth value.

**What was new.** A classical problem family turned out to be *equivalent* to
a slice of Collatz: Mahler's fractional parts of `xi (p/q)^n` (THM-4469). The
proved Mahler-type exclusions stop at interval length `1/p`, and the Collatz
slices need `alpha/(alpha-1)` times that. So the Collatz divergence half
sits just past the edge of what Flatto–Lagarias–Pollington-type arguments
reach. That edge is exactly where two free digits per step first become
possible. The theme of this session recurs once more: proved methods stop at
the threshold where Collatz-type freedom begins (Sturmian vs positive
entropy; `phi` vs coupled scales; `1/p` vs `1.88/p`).

**Corrections I made to my own prompts** (logged in MISTAKES 2026-09-24):
* a residue `n mod 3*2^k` decides the predecessor mod `2^(k+1)`, not `2^k`,
  and never its class mod 3;
* the bracket threshold formula I gave is the real-interval one, and for
  integers it overshoots by one on explicit windows;
* "tree growth counted by size" cannot separate the sheets, because `x -> -x`
  preserves `|x|`.

## Wave 9 (2026-09-24): the sign lives in the cross term

The owner pointed at a claimed proof of Collatz, arXiv 2502.20642. It is a
fixed-point argument in the metric `|x−y|`, and its key lemma is the
triangle sandwich `|a−b| ≤ c ≤ a+b`. The owner asked us to read the two
sides as positive and negative and the middle as zero.

The lemma is correct. The general theorem is false: the successor map
satisfies all of its hypotheses. It is still worth recording *why* such
an argument had to fail, because the owner's reading of the sandwich says
exactly why.

* **The sign lives in the cross term `±2ab`.** For a parity word `w`,
  `2^p |T^p x| = |3^a x + c_w|`. It sits on the upper side of the sandwich
  for `x > 0` and on the lower side for `x < 0`, with `0` as the middle
  point. So the sign law is precisely the sandwich's two equality cases.
* **Even observables throw it away.** The quantities `|x|`, `x²` and
  `|x−y|`, and the AM–QM step that the paper's lemma uses, all discard
  the cross term. The negation symmetry `x -> −x` (the mod-192 note) is an
  isometry of `|x−y|`. A metric argument in `|x−y|` therefore sees the
  positive half-line exactly as it sees the negative one, where Collatz
  has three extra cycles.
* **The breakdown is visible.** The paper's coefficient table, copied
  verbatim, also "proves" that `3n−1` reaches 1. Its error falls at every
  up-step, where positive orbits sit on the `+2ab` side.

**The fixed-point theorem that does fit Collatz is Banach's in `Z_2`.** The
inverse branches contract, and every parity word has exactly one periodic
point, its cycle gate. Existence is free; the problem is integrality.
Brouwer's one-dimensional form (the intermediate value theorem) gives all
periods on both half-lines. My guess that the sign decides the Sharkovskii
type was wrong.

**The owner's tournament picture found a different symmetry.** The two
4-vertex tournaments that swap under reversing every arc are the forward
map and the inverse tree of one AM-fair pair, i.e. time reversal. They are
not the two sheets. `3n+b` inverts to `(y−b)/3`, so the owner's "±1"
really does flip, but it flips with time, not with sign.

**Digits.** Collatz has an exact, non-decaying version of the
consecutive-prime digit bias. Its transition probabilities are
`1/15, 2/15, 4/15, 8/15`, and `3 -> 5` is certain. The difference from
primes is structural: halving counts are identically distributed at every
size, while prime gaps grow.

Repunit primes are the prime fixed points of digit rotation. In base 2 the
infinite repunit `...1111 = −1` is Collatz's own hostile fixed point.

## Waves 13–16 (2026-09-26): the obstruction has an exact anatomy, and height is the coordinate the 2-adic world cannot see

After waves 13–16 the obstruction to proving Collatz by the methods this
session can formalize has an exact anatomy. Every piece is a theorem.

**1. The obstructions are expanding cycles and their backward trees.**
* **Bounded-lookahead provability** holds iff no expanding cycle exists in the parity graph (THM-4474).
* **A rank** `a log n + h` with `h` bounded or periodic exists iff the same holds. For Collatz its least defect is `log(3/2)`, from the loop at `-1` (THM-4482).
* **Valuation ranks.** A rank that adds 2-adic valuation counters must charge every point of the backward tree of every expanding cycle (THM-4483). The tree of `-1` alone has `2^(j-1)` points at depth `j`.
* **The same list recurs.** The sign strategies' min-max density game (THM-4486) finds the same objects:
  * for `q = 3`, the free cycle `{1,2}` pins the game at `1/2`;
  * for `q = 5`, the sporadic cycle `(1,3,8,4,2)` fixes its value at exactly `2/5` for all `k >= 15`.

**2. The integral obstructions are four free cycles and one sporadic.** For 3x+1 on `Z`, four cycles are forced by the identities `2-1`, `3-2`, `4-3` and `9-8 = 1` (Gersonides), and one, `-17`, is sporadic, with gap `139` (THM-4484). This is the owner's Kuratowski–Tutte pattern:
* a characterization by excluded substructures (Kuratowski ↔ Theorem A);
* a counting reason (Euler's formula ↔ the entropy `h(log_3 2)`);
* obstructions forced by an identity, plus a sporadic one.

The correspondence is ANALOGY, but the pattern is exact. Where it stops is exact too: Kuratowski's list is finite, while the expanding cycles form an infinite, necklace-counted family. The finite integral part of that list is the `3x-1` cycle conjecture.

**3. Height is the continuous coordinate that 2-adic methods cannot see.** The price of making every orbit descend within a fixed horizon was measured in three settings:
* **arbitrary edits** (THM-4478, THM-4480): `2^(-(1-H(log_q 2))L)` for every multiplier, exponentially cheap even for 5n+1. The trick is to catch each orbit at its peak, where edits are sparse by height;
* **periodic deletions** (THM-4485): the feedback number of the expanding cycles of the de Bruijn graph, which is `~1/k` for 5n+1;
* **periodic sign flips** (THM-4479, THM-4481): no cheaper than deletions. A flip merges only two orbits, so an entropy law forces constant flip mass for positive drift and makes provability impossible for `q >= 23`.

Periodic modifications are residue classes. They live in `Z_2` and cannot see height, which is why positive drift costs them exponentially more. Ranks show the same thing: the only nonnegative valuation ranks that work are those that store each integer's stopping time in the height of a private center (THM-4483 D). So "2-adic plus real" is necessary and, in that class, exactly as hard as Collatz.

**4. Continuous methods lose exactly what integers keep.**
* **Moments.** Distribution-only (moment) arguments lose a square (THM-4477). Integer spacing recovers it (THM-4478), and so does periodicity (THM-4479).
* **Stationary laws.** Uniform-chain stationary laws cannot certify a density floor above `1/3` (THM-4486 F), while the true floors are worst-cycle phenomena. Averages cannot see the cycle that matters.
* **Strips.** The second-order `L^(1/3)` term of the price is a Brownian-strip constant, obtained by two independent sessions from exact sine eigenfunctions of the letter walk.

**What would count as a new idea.** Every method formalized here reduces Collatz to excluding expanding structures that live on the negative side of `Z_2`: `-1`, `-5`, `-17` and their backward trees, dense in `Z_2`. A proof must therefore use the one thing that separates positive integers from those trees: the sign, i.e. height seen from the real place. It must use it in a form that is not a bounded, periodic or finite-mass 2-adic potential. The theorems above rule those forms out one by one. What is left is an adaptive, unbounded, height-coupled quantity. The only known instance of such a quantity is the stopping time itself.

## Wave 18 (2026-09-26): methods reach their exact ceilings, and small graphs keep their arithmetic

This wave repeatedly found the *exact ceiling* of a method, not just a bound.
* **Landing multiplicity (THM-4506).**
  * The one-window recursion has exactly one input it can still improve, the averaged multiplicity. Its worst case is exactly `ceil((k-D)/log_2 3)`, and the recursion is saturated at `a*(mu)`.
  * Averages over depths, residue classes and parity words all leave `mu = 1`.
  * What would help is a local time of one orbit: how often one orbit revisits one dyadic shell before crashing below it. The concurrent crossings session reached the same object from the other side, as the length of a stay below `X` measured against the records inside it.
* **Sign strategies (THM-4508).**
  * A finite rule is a single 2-adic map with rational periodic points, so optimality at every level is one finite object.
  * The flip calculus says that only one local move ever pays: the pattern `(s,2),(s,2)`, the classes of `±1/3` for 7n±1.
  * For 5n±1 that single move lands on the optimum, the sporadic cycle `1,3,8,4,2`. For 7n±1 it stops at `2/5`, and after that the rules must grow without a visible pattern. Whether they ever reach `log_7 2` is the open question.
* **Fences (wave 19 lane).** A corner count at the junctions (tight exactly at T, Y and L) is the Euler identity of the problem. Combined with isoperimetry it caps the fence density at `0.5225`. The cap is an LP optimum made of regular pentagons, which do not tile, so the truth lies in `[1/2, 0.5225]`.

The small graphs behaved differently: they kept their arithmetic exactly.
* The square-sum threshold at 15, with ends 8 and 9, is a degree law plus one Pell identity. By Anglin's theorem that identity has exactly one solution, so the square case of the zigzag happens once.
* The Catalan identity `9 - 8 = 1`, the same one that forces the free Collatz cycle `-5,-7,-10`, is a coincidence of the pair (8,9) for square sums. It becomes a mechanism only where 8 and 9 are both sum targets.
* The Collatz alphabet `C_n` (sums that are powers of 2 or of 3) never closes into a Hamiltonian cycle. Its admissible windows are cut by the gaps `|2^p - 3^a|` in the order of the Beatty word of `log_2 3`.

**Discrete and continuous.** In every exact result of this wave, a continuous inequality closes a discrete count exactly, and the extremal case of the continuous inequality names the obstruction:
* optional stopping for the rise law;
* Lemma S's single dyadic shell for multiplicity;
* isoperimetry for fences;
* the contraction `g_s(y) = (4y - s)/q` for the max-halving skeleton.
