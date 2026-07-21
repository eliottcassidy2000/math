# Graffiti's WOWII-103 counterexample, and how to leverage it: TournamentGraffiti + the marginal-threshold refutation

**death-star-2026-07-20-S77** (HYP-8618). Owner: consider the counterexample to Written-on-the-Wall II
conjecture 103 (google-deepmind/formal-conjectures PR #4482) and think how similar ideas leverage into
the repo, or new problems we can attack. Verified the counterexample exactly and prototyped the leverage
(`graffiti_leverage_tournament_deathstar_S77.py`).

## The conjecture and its counterexample
**WOWII-103** (a Fajtlowicz/DeLaViña *Graffiti* machine-generated conjecture), for connected graphs $G$:
$$\alpha(G)\ \le\ \big\lfloor\, b(G)-\ln(\overline{\mathrm{ecc}}(G))\,\big\rfloor,$$
$\alpha$ = independence number, $b$ = largest induced bipartite subgraph, $\overline{\mathrm{ecc}}$ = average eccentricity.
**Counterexample** (CinnamonRolls1, Lean-verified by exhaustive $2^{11}$ enumeration): a **triangle with 4
leaves on each of two of its vertices** (11 vertices). Verified here: $\alpha=9$, $b=10$, $\overline{\mathrm{ecc}}=30/11$.

## The engine of the counterexample — a *marginal transcendental threshold*
The refutation works for a precise and beautiful reason:
$$\tfrac{30}{11}=2.7273\ >\ e=2.71828\quad\Longrightarrow\quad \ln(\tfrac{30}{11})>1\quad\Longrightarrow\quad
\big\lfloor 10-\ln(\tfrac{30}{11})\big\rfloor=8<9=\alpha.$$
The graph is tuned so its **average eccentricity sits a hair above Euler's $e$**, so the $\ln$ term barely
exceeds $1$ and the *floor* drops one integer below $\alpha$. The whole counterexample is one integer
boundary crossing driven by a transcendental constant. That is a **transferable technique**, not a
one-off.

## Leverage 1 — TournamentGraffiti: an automated conjecture/refute engine on the repo's invariant zoo
*Graffiti* (and its successor **TxGraffiti**, which has produced genuinely publishable graph theorems) is
just: an **invariant database** + a **candidate-inequality generator** + an **instant tester**. **This repo
already has the database** — for tournaments $n\le 7\!-\!8$ we compute $H$ (Ham-path count), $c_3$ (3-cycles
$=\tfrac13\mathrm{tr}A^3$), score sequences, char-poly / eigenvalues, the skew-determinant $\mathrm{disc}$, the
tournament pencil $M(t,u,v)$, SC/NS type, even-graph invariants of $E_n$, OCR, $\dots$ A `tournament_graffiti.py`
that generates inequalities among these and tests them exhaustively would surface (a) survivors = candidate
theorems and (b) refutations = counterexamples, exactly as TxGraffiti does — and it is a clean **engineering
deliverable** (per the equal-priority mandate; sibling to `mod_rank_library`, `tournament_tda`). The
prototype here already: confirms **Erdős–Moser** $\tau\ge\lfloor\log_2 n\rfloor+1$ survives $n\le 6$
($\tau$ = largest transitive subtournament), and auto-refutes naive floor+log guesses.

## Leverage 2 — the marginal-threshold refutation, aimed at the repo's own transcendental bounds
The repo is unusually rich in **transcendental constants inside invariant bounds**, which is exactly where
103-style counterexamples live:
- the "everything is the triangle" constants $\sqrt2,\ \pi,\ e,\ \gamma$ (fiber fraction $\sim1/\sqrt{\pi n}$,
  pseudo-doubling $2-\tfrac1{n-2}$, Burnside/Stirling $e$, $\Gamma(1/b)^b\sim b^be^{-\gamma}$);
- LRC's $1/(n+1)$ and the sharp measure-horn $1/(7L)$ (THM-1123);
- any repo bound of the shape "*integer invariant* $\le \lfloor f(\text{invariants})\rfloor$".
**Technique:** when a repo conjecture floors a transcendental $f$, hunt for the configuration where $f$
crosses an integer/rational boundary — the counterexample (or the sharp/hard case of a *proof*) is the
marginal one. Two immediate targets: my **$\{7,21\}$ Ham-path forbidden-gap** (S70) — is $21$ truly forbidden
or a not-yet-seen margin? — and the repo's own **width formula** $W(G_n)=\binom{n-2}{\lfloor(n-2)/2\rfloor}$,
which *holds $n=3\!-\!6$ and breaks at $n=7$* (CLAUDE.md): a **native 103** already in the ledger, the exact
"holds small, breaks larger" shape the prototype reproduces.

## Leverage 3 — the invariant analogues, which land on repo-native objects
Mapping WOWII-103's three invariants into tournament theory is fertile, and each lands on something the repo
already privileges:
- **$\alpha$ (independence) $\to \tau$ = largest transitive subtournament.** The acyclic/ordered analogue;
  Erdős–Moser is its Caro–Wei. (Verified survivor above.)
- **$b$ (largest induced bipartite) $\to$ largest induced *even* subtournament** — directly the repo's
  **first-class even graphs $E_n$**. And the counterexample's structural core, "one odd cycle, so delete one
  vertex to reach bipartite ($b=n-1$)", has the *exact* tournament mirror: **one 3-cycle** is the minimal
  non-transitive obstruction, and "almost-transitive" tournaments (a single 3-cycle) are my S75
  intransitivity-graded-by-cycle-length. The **odd-cycle $\leftrightarrow$ 3-cycle** correspondence is precise.
- **$\overline{\mathrm{ecc}}$ (eccentricity) $\to$ king-eccentricity** (kings reach all others in $\le 2$; every
  tournament has a king). A natural, unstudied-here tournament invariant to add to the zoo.

## Leverage 4 — the formal-conjectures repo as a two-way bridge (outreach + engineering)
`google-deepmind/formal-conjectures` is a curated Lean collection of open conjectures. Two concrete leads:
(i) **contribute** the repo's formalized artifacts (the LRC14 moment-floor route `LRCWitnessMomentFloor.lean`,
GMC(2) strata, tournament theorems) as formalized statements; (ii) **mine** the collection for conjectures
attackable by our parity/moment/tournament machinery. Added to the backlog.

## The deeper resonance
The 103 workflow — a machine conjecture, a small explicit counterexample, a *decide*-based exhaustive Lean
check — **is the repo's own epistemology** (CLAUDE.md: "check small cases exhaustively before generalizing";
MISTAKES.md is a ledger of holds-small-breaks-larger). So this is less a foreign import than a mirror: it
validates the exhaustive-small-$n$ + kernel-Lean style, and it hands us a ready-made *generator* (Graffiti)
and a ready-made *refutation heuristic* (marginal thresholds) for the invariant inequalities the repo keeps
stumbling into by hand.

## Honest scope
The counterexample is verified exactly (incl. the $30/11>e$ crux). TournamentGraffiti is a working prototype,
not yet the full engine; the four leverages are directions (one — the width-formula-as-native-103 — is
already a repo fact). No new theorem. Cross-links: S70 ($\{7,21\}$ forbidden-gap), S75 (intransitivity by
cycle length), S76 (invariant atlas / forbidden-gap lens), even graphs $E_n$ as first-class, the width
formula (CLAUDE.md), engineering-synthesis roadmap. `graffiti_leverage_tournament_deathstar_S77.py` (+out).
HYP-8618.
