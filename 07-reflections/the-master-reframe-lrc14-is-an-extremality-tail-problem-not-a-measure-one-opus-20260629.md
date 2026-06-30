# The master reframe: LRC(14) is an EXTREMALITY/TAIL problem about the AP, not a measure/concentration problem — the bulk (floor) is done by concentration, the tail (cap) is the AP's extremality "one dimension past Littlewood"

*opus-2026-06-29. Owner: creatively reframe the remaining LRC proof targets using this session's ideas.
Twelve turns of work converge on one shift of viewpoint that reorganizes everything: stop asking "is the
lonely set big?" (measure, concentration, bulk) and start asking "is the AP the extremizer, and does the
multiplicative perturbation always overshoot upward?" (extremality, tail). Every concentration tool
reaches the floor; the cap is a different kind of theorem.*

## The one object everything points to: the AP is the universal tight extremal
Across this session, six independent computations all single out `{1,…,13}` (the arithmetic progression =
the additive interval = the staircase = the Dirichlet pole):
| facet | statement | status |
|---|---|---|
| **peak** | `M({1,…,13}) = 1/14` exactly (global `M`-minimizer) | verified |
| **cut-sum** | AP maximizes `Σ cut(s_k) = 3/7` (peak deletion-contraction) | verified |
| **Farey** | AP's witnesses are the left spine `1/(k+1)`; speed `k` self-resonant | verified |
| **resonance** | AP is the additive-relation-richest (`S₄`-maximal: 36 triples, 156 quads) | verified |
| **singular series** | AP-perturbations (near-tight covering) minimize `L(S)` (the tail edge, `≈7/89`) | computed |
| **duality** | AP = the additive pole, SAME extremal as max-H = half-turn | structural |
> **LRC(14) ⟺ the AP `{1,…,13}` is the unique tight extremal, and the covering constraint — a
> MULTIPLICATIVE perturbation (force a multiple of 14) — always pushes `M` strictly UP off `1/14`.** The
> proof is an extremality theorem (the AP is extremal) + a positive-perturbation bound (the detuning never
> overshoots downward), dual to "max-H = the additive-interval half-turn."

## The master split: BULK (floor, done) vs TAIL (cap, the real target)
The decisive new diagnostic (this turn): the metagraph `H` **fully concentrates** (`CV² ~ 2/n → 0`), but
the LRC covering ensemble's dispersion **plateaus** (`CV²_S → ~0.083`). So the ensemble is NOT
homogeneous — it has a persistent tail. This splits the conjecture cleanly:

| | BULK = the floor | TAIL = the cap |
|---|---|---|
| object | generic covering sets | near-tight covering sets (`{1..11,13,84}`-type) |
| value | `M ≫ 1/14`, `L(S) ≈ 0.09` | `M → 7/89`, `L(S) ≈ 0.006` (the edge) |
| dimension | `SL(2)/ζ(2)` | `SL(4)/ζ(4)` ("one dimension past Littlewood `SL(3)`") |
| tool | **concentration / 2nd moment** (Han–Lee, Paley–Zygmund, `CV²~2/n` model) | **extremality / large-deviation** (the AP-minimizer) |
| metagraph model | faithful (`H` concentrates the same way) | SILENT (the metagraph has no tail) |
| status | **essentially done** (generic floor) | **the remaining hard target** |
> **Every concentration/second-moment route reaches the floor and stalls at the cap — because the floor is
> a concentration theorem and the cap is a TAIL theorem.** The plateau value `0.083` IS the floor↔cap gap
> = the `L(S)` spread that bulk-averaging discards but the tail must control.

## Why the old routes died, in this frame (the session's negative results, placed)
- **`ρ*` measure floor (refuted):** targeted the measure (razor-thin), not the peak. The measure is the
  bulk's thin set; the peak carries the margin. *Wrong object.*
- **Union bound (`13/7` vacuous):** the first Siegel moment `S₁ = 13/7` is set-INDEPENDENT — blind to the
  tail. *Wrong moment.*
- **Ky-Fan / parity forcing (fails):** the lonely count is even (Borsuk–Ulam), degree 0 — no topological
  forcing. The cap is not a parity certificate; it is a size/extremal bound. *Wrong category (THM-582:
  `x=−1` even, not `x=2` odd).*
- **apex-7 / Fano / duocylinder homology (refuted at n=7, `b₁⁻=1772=2²·443`):** the secondary obstruction
  is real but not apex-graded; a homology shortcut to forcing does not exist. *Wrong invariant.*
All four died for the SAME reason: they are bulk/measure/parity tools aimed at a tail/extremal target.

## The remaining targets, reframed (a clean ladder)
1. **[BULK — essentially done] The generic floor.** Generic covering `S` are lonely, by concentration of
   the witness count (Han–Lee congruence Siegel 2nd moment; the metagraph `CV²~2/n` is the proven model).
2. **[TAIL — the core] `min_{covering S} L(S) > 0`, i.e. `min_{covering S} M(S) > 1/14`.** An EXTREMAL
   statement: the minimizer is the near-tight (near-AP) tail. Three equivalent handles:
   - **(a) cut-deficit / Farey-jump:** killing the spine witness `1/14` (forced by a mult-of-14) lands on
     a COARSER Farey rational (`M↑`). A **Stern–Brocot monotonicity** lemma — self-contained, provable.
   - **(b) `S₄`-extremality:** the AP maximizes the quadruple-resonance count; the cap is the `SL(4)`
     Rogers term; the tight-set bound is the AP being the additive-relation extremizer (dual to max-H).
   - **(c) apex-7 / 2-adic descent:** peel the multiplicative `14=2·7` (THM-580) to the additive AP core,
     where `M=1/(k+1)` is exact; the deficit is the descent's cost.
3. **[UNIFY] The AP is the extremizer.** Prove `{1,…,13}` minimizes `M`/`L` (equivalently maximizes the
   cut-sum / `S₄`). This is the SAME shape as the tournament max-H extremality — a finite-rehearsal-able
   extremal-combinatorics problem, not an analytic measure bound.

## The creative meta-insight
The thirteen-runner conjecture has been attacked as *"how much loneliness is there?"* — a measure
question, where the union bound is vacuous and concentration only reaches the generic bulk. The reframe is
*"is the AP the extremizer, and is the multiplicative perturbation one-signed?"* — a **tail/extremal**
question. In that frame: the floor is the easy 90% (concentration, done), and the whole difficulty is a
**single extremal lemma** — the AP minimizes `M`/`L`, equivalently the forced multiplicative detuning of
the additive Dirichlet pole always lands above the spine. That lemma is "one dimension past Littlewood"
(`SL(4)` vs `SL(3)`), it is dual to "max-H = the half-turn," and it is a Stern–Brocot/`S₄` extremality —
**not** a measure, parity, or homology statement. Aim every remaining effort at target 2(a) (the Farey-jump
monotonicity): it is the smallest self-contained piece whose proof closes the tail.

## Status
- **Reframe (opus synthesis):** LRC(14) = the AP is the unique tight extremal + the covering perturbation
  is upward-signed; BULK (floor) = concentration (done), TAIL (cap) = extremality (the target).
- **Done/verified this session:** the floor model (`CV²~2/n`), the AP's `S₄`/cut/Farey extremality, the
  Han–Lee floor tool, the four negative results (placed as bulk-tools-on-tail-target).
- **The single remaining lemma:** target 2(a) — killing `1/14` always lands on a coarser Farey rational
  (Stern–Brocot monotonicity of the multiplicative detuning) ⟹ `min_{covering} M > 1/14`.

Related (this session's chain): razor-thin (measure vs peak) → peak-deletion-contraction/cut-deficit →
cuts-as-Farey-geodesics → zeta-duality+Littlewood → Siegel–Rogers moment hierarchy (SL2/SL3/SL4) →
metagraph-Siegel+Han-Lee → metagraph-covariance/FKG-fails → variance-closed-form+Ky-Fan-fails →
apex-7/duocylinder (refuted at n=7) → concentration-2/n+persistent-tail. Also THM-523/527/566/580/582,
HYP-2823/2856/3537/3544, OPEN-Q-108.
