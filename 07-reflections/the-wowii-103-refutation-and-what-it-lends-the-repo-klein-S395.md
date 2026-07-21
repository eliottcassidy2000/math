# The WOWII-103 refutation, and what its ideas lend the repo

*klein-2026-07-21-S395. Owner: consider the counterexample to Written-on-the-Wall-II Conjecture 103
(google-deepmind/formal-conjectures PR #4482) and how similar ideas can be leveraged here, or on
new problems.*

## What the counterexample is (verified from the PR)

**WOWII Conjecture 103** (DeLaViña's *Written on the Wall II* / Graffiti.pc): for a connected
graph `G`,

```text
   α(G)  ≤  ⌊ b(G) − log(ecc_avg(G)) ⌋
```

`α` = independence number (**the "hard"/extremal invariant**), `b` = largest induced bipartite
subgraph size, `ecc_avg` = average eccentricity (**the "easy"/structural invariants**), `log`
natural. **Disproved** by an explicit 11-vertex graph — *a triangle with four leaves on each of
two of its vertices* — with `α = 9`, `b = 10`, `ecc_avg = 30/11 ≈ 2.727`, so
`RHS = ⌊10 − 1.003⌋ = ⌊8.997⌋ = 8` and the inequality reads `9 ≤ 8`. Found by search, checked by
**exhaustive enumeration of all 2,048 vertex subsets** and **formally verified in Lean**
(`decide +native`).

## The three transferable ideas

1. **Shape: an easy invariant bounding a hard one, refuted by an explicit tuned witness.** The
   conjecture squeezes the extremal `α` under a cheap structural expression; the refutation is a
   small graph where the structural side lags by one integer. **This is the exact shape of the
   repo's own THM-1460/1580** (Hamiltonian-path count `H`, `#P`-hard, vs the poly-time
   arborescence count `Σa`, with THM-1580 already finding *where the easy invariant is loose* —
   the `n=7` residue). The repo has been doing WOWII-style analysis by hand; WOWII names it.

2. **Structure: an odd-cycle core + amplifier pendants, tuned to cross a floor.** The triangle is
   the sole non-bipartite obstruction (it keeps `b < n`), and the eight leaves inflate `α` to `9`
   while the `−log(ecc)` correction shaves the floor from `9` to `8`. **The triangle is a
   3-cycle** — the repo's intransitivity atom (THM-1805, THM-1840's both-signs atom). The
   template "one odd/3-cycle obstruction + independent amplifiers, tuned to break a tight integer
   bound" is directly buildable in the tournament/metagraph setting.

3. **Discipline: search → exhaustive check → formal (Lean) proof.** Three of these four steps are
   already the repo's culture (MISTAKES.md is a museum of "pattern breaks at n=6/7"; the LRC and
   TournamentH7 Lean threads exist). **The one piece the repo does not systematize is the
   automated-conjecture-generation front end** — Graffiti/WOWII *generates* the inequalities to
   test. Bolting a conjecture generator onto the repo's invariant zoo is the upgrade.

## Concrete leverage (targets, one demonstrated)

- **Demonstrated.** The repo's metagraphs `G_n`, `E_n` are ordinary connected graphs, so the
  whole WOWII list applies to them and has never been run. Computed the WOWII-103 invariants on
  `G_n`: `α(G_n) = 2, 5, 18` and `b(G_5) = 9`, `ecc_avg = 3/2, 5/2, 25/7` for `n = 4,5,6`. `G_n`
  **satisfies** 103 (e.g. `G_5`: `5 ≤ ⌊9 − log 2.5⌋ = 8`), so it is not a counterexample — but
  the machinery runs, and `α(G_n) = 2,5,18` and the metagraph's eccentricity profile are new
  data. (`04-computation/wowii_on_metagraph_klein_S395.py`.)

- **Off-by-one fragility of repo bounds.** WOWII-103 fails by *one* at a floor. The repo has
  tight integer bounds ripe for the same check: is THM-1790's detection depth **exactly** `d+1`
  (my lower bound) or does the upper bound miss by one? Is `H ≤ Σa` (THM-1460) ever *equal*, and
  where is it loosest? A WOWII-style explicit-witness pass could pin each.

- **Tournament analog of the counterexample template.** Build "3-cycle core + dominated
  amplifiers" tournaments and test tight tournament-invariant inequalities (H vs OCF vs the
  arborescence ranking THM-1750 vs the Pfaffian THM-1475). The 3-cycle is the natural
  non-transitive obstruction, exactly as the triangle is the non-bipartite one.

- **The H-spectrum as a WOWII "achievable-values" conjecture.** "The `H`-spectrum is odd `∖
  {7,21}`" (death-star HYP-8540) is a Graffiti-flavored achievability claim; the missing piece is
  a *construction* realizing every odd value — a generation problem WOWII-style search fits.

- **New problems.** The WOWII list is ~150 inequalities; formulating **directed / tournament
  analogs** (directed independence, 3-cycle-density vs score-spread) is a fresh problem family the
  repo is well-equipped to attack, and formalizing repo refutations in Lean (the WOWII discipline)
  raises their standard.

## The honest size of this

Not a theorem — an idea transfer. The one load-bearing recognition: **the repo's flagship
arborescence-vs-H result (THM-1460/1580) is a WOWII-shaped "easy invariant bounds hard invariant,
find the loose witness" statement, and the repo should adopt WOWII's automated-generate +
explicit-refute + Lean-verify pipeline for its whole invariant zoo.** The metagraph check
demonstrates the machinery ports; the off-by-one and 3-cycle-core targets are the concrete next
attacks.

*Files: `04-computation/wowii_on_metagraph_klein_S395.py` (+ `.out`). Backlog items filed.*
