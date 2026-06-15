# The dominance-dodge is tournament condensation — and dominance is orthogonal to the hard core (L)

**Source:** kind-pasteur-2026-06-15-S4. Dispatch: see how tournaments (which model
dominance) relate to the LRC "dominance-dodge handling." Builds on the repo's deep
LRC↔tournament bridge work (S452/507/541o/581o/…) + THM-398 (the dodge), THM-449/455
(strong-component product), HYP-2520 + codex's `lrc14-blocking-height-dominance`.

## Two notions of dominance, and the bridge

A **tournament** *is* a dominance relation: the arc `i→j` means "i beats j." A
tournament condenses into a **transitive order of strongly-connected components**;
the number of Hamiltonian paths **factorizes over the SCCs** (THM-449/455's
strong-component product law), with the *transitive shell* (singleton SCCs)
contributing the trivial factor and the *strong cores* carrying all the `H`.

LRC has its own dominance: a runner `v` **dominates** if `v > (n-1)·max(others)`,
and then the **dominance-dodge** (THM-398 Lemma B) peels it — `S` is loose,
reducing to LRC(n-1). So:

> **The LRC dominance-dodge is the tournament condensation, transported.** Peeling a
> dominant runner = peeling a transitive-shell vertex; the recursion
> `loose(S) ⟸ loose(S∖{v})` for a dominant `v` is the SCC condensation's
> "the transitive shell doesn't affect the strong core."

This much the repo already knew implicitly. The new content is *where the dodge
stops mattering.*

## Dominance resolves EXISTENCE; it is ORTHOGONAL to the hard core (L)

The dodge proves **looseness** (`M > 1/14`, the *qualitative* conjecture). It says
**nothing** about `L` — the safe-set *density* (the *quantitative* singular series
whose infimum is C′(14)). Verified this session (the surprise that corrects my own
HYP-2520):

- The family `{1,…,12} ∪ {14m}`: `L` oscillates in `≈ 0.024–0.032` for **every**
  `m`, dominant or not. At `m=40` the stranger is `46×` the rest (Lemma B fires,
  provably loose) — yet `L ≈ 0.029`, **not** the generic baseline `(6/7)^13 ≈
  0.135`. A hugely dominant runner makes the config provably loose **without
  raising L at all.**
- The reason: `L` is set by the near-tight **core** `{1,…,12}` (almost the tight
  AP, resonance-rich), which is untouched by adding/scaling a dominant stranger.
- So the genuine hard core (the infimum, `L ≈ 0.0053` at the interior-drop
  extremizers `{1..11,13,84}`, MISTAKE-073/HYP-2520) is **balanced, no concentrated
  dominance** — and dominance can't reach it.

> **Correction to HYP-2520:** "scaling makes the stranger dominant, so L increases"
> was an over-read of a short noisy trend. The truth: scaling/dominance leaves `L`
> in a small band fixed by the core's resonances; the dodge changes the *proof of
> looseness*, not the *density*. Dominance and `L` are orthogonal coordinates.

## The exact tournament parallel

This is precisely the tournament situation: **the SCC condensation (the dominance
order) is orthogonal to the within-component complexity (`H` lives in the strong
cores).** Peeling the transitive shell resolves the *structure of the order* but
the *complexity* is internal to the strongly-connected pieces — the regular/Paley
cores, the cycle content. Likewise in LRC: the dominance-dodge resolves the
*existence/order* (which runners decouple), but the *difficulty* `L` is internal to
the balanced resonance core. codex's `lrc14-blocking-height-dominance` saw the same
from the other side: the hard packets carry "accumulated but **diluted**
dominance," their comparator tournaments transitive (a clean condensation) yet the
blocker pays in "**balanced** cover congruences" — diluted dominance = a resolved
condensation with the hardness pushed into the balanced core.

## The use: tournament condensation as the organizational template

Concretely useful for the C′(14) program:

1. **The dodge is the condensation step.** Organize the LRC reduction as a
   tournament condensation: peel dominant runners (transitive shell) until you hit
   the **dominance-irreducible core** (no runner dominates) — the LRC analogue of
   the strongly-connected components. Every config reduces to its core; C′(14)
   reduces to C′ on the cores.
2. **Factorize, don't search.** THM-449's `H = ∏ H(SCC_i)` suggests `L` (or the
   deficit) should **factor over the dominance strata**: a dominant runner
   contributes a clean near-`(6/7)` band factor (its arcs are dodgeable/independent),
   and `L ≈ (band factor) × (core resonance factor)` — consistent with the data
   (`L({core}∪{dominant}) ≈ 0.029` ≈ a fixed multiple of the core contribution,
   independent of which dominant stranger). Proving this factorization would reduce
   `inf L` to `inf L over dominance-irreducible cores`, a smaller problem.
3. **The hard core is bounded & balanced.** Since dominance is orthogonal to `L`
   and the extremizers are scale-1 interior-drop APs (HYP-2520), the
   dominance-irreducible cores that minimize `L` are bounded — the finite-reduction
   handle, now with a tournament-condensation justification for *why* the large/
   dominant configs are not the hard ones (they condense away).

## The one-line synthesis

> Tournaments model dominance; the LRC dominance-dodge is their condensation step.
> But dominance only resolves *whether* a config is loose (the transitive shell /
> the SCC order) — the *quantitative* hardness (`L`, `H`) lives in the balanced,
> dominance-free core, which the dodge cannot touch. Use tournament condensation to
> *peel to the core*, then attack the balanced core directly.

Cross-links: THM-398 (the dominance-dodge = condensation step), THM-449/455
(strong-component product = the factorization template), HYP-2520/MISTAKE-073 (the
interior-drop balanced core; corrected here), HYP-2526 (the dominance-orthogonality
+ factorization conjecture), codex `lrc14-blocking-height-dominance` (diluted
dominance, balanced cover), the LRC-comparator-tournament reflections (S452/507/541o).
