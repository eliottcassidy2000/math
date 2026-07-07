---
source: kind-pasteur-2026-07-07-S68
status: SYNTHESIS (the history + the mechanisms + reusable structural knowledge) + one
  concrete bleeding-edge result (the anchored-gap subset lemma finitizes the 2-anchor (A')
  tail for all k). Owner: keep working the bleeding edge; synthesize how the proof changed,
  which angles worked and which did not and WHY, then abstract the mechanisms as structural
  knowledge.
tags:
  - lonely-runner
  - LRC14
  - synthesis
  - proof-mechanisms
  - structural-knowledge
  - finitization
  - anchor-floor
---

# The mechanisms of the LRC(14) proof: what worked, what failed, and why

**kind-pasteur-2026-07-07-S68 (HYP-5057).** Two things: a concrete step on the current
bleeding edge (the 2-anchor tail), and the synthesis the owner asked for — the arc of the
proof, which angles worked and which died, the *mechanism* behind each, and those mechanisms
lifted to reusable structural knowledge.

## Part 0 — the bleeding-edge step (mechanism in action)

The live open lemma is the per-`k` tail floor **(A') `μ_{1/7}(E) ≥ T_k`**, reduced (boxeph-S1)
to the **2-anchor tail** `PA_2(E) = P_x(max(gap@0, gap@1/2) > 1/7) ≥ T_k`, with the residual
"spread AP minimizes `PA_2`." I applied my S59 winning move — *domination by the AP →
finitization* — to this anchored object:

> **Anchored-gap subset lemma (proved; 0 / 27000 checks).** For integer sets `E ⊆ F`,
> `{frac(ex)} ⊆ {frac(fx)}`, so the gap containing a fixed anchor `a` is a union of ≥ 1 gaps
> of the `F`-config: `gap@a(E,x) ≥ gap@a(F,x)` **pointwise**, hence `PA_2(E) ≥ PA_2(F)`. With
> `F = {0..D}` (consecutive AP, `D = diam E`): `PA_2(E) ≥ PA_2(AP_{D+1})`, exact (shifted
> three-gap).

So the 2-anchor tail **finitizes on bounded diameter for EVERY `k = 8..13`** — diameter bites
`8, 9, 11, 15, 27, 68` (largest `n` with `PA_2(AP_n) ≥ T_k`), not just the `k=13` leg my
original diameter floor reached. The remaining spread-AP core has its `inf` at **bounded**
`(a,d)` (resonant dip at `d=2`; large `d` decorrelates to `≈ 0.72/0.56/0.37 ≫ T_k`), so it too
is a finite `(a,d)` check plus a decorrelation-past-`d₀` tail — the same shape as my S61
Part-A `V₀`. The winning mechanism ports cleanly; that is itself a data point for the synthesis.

## Part 1 — the arc of the proof (how the state changed)

1. **Route 2 (retired).** LRC(14) → (J–K reduction) → rank-2 torus → the 1-D Farey gap `(C)`.
   Lean built out. **Killed (MISTAKE-117):** Jain–Kravitz control the *accumulation points*
   `acc(S(n)) = S(n−1)`, not the *supremum* the LRC bounds. A counterexample would be an
   *isolated* point of `S(13)` above `acc`, which accumulation machinery cannot see.
2. **Route 1, first target (refuted).** Direct `M ≥ 1/14` via the good-period density with a
   `2/7` uniform floor. **Killed:** admissible zeros — the via-max criterion is
   sufficient-not-necessary. Corrected to the sharp `1/7` (boxeph); the object became
   `μ_{1/7}(E)`.
3. **The reductions that stuck.** Coarse/scale reduction (multi-scale → settled LRC(≤13),
   GREEN); the **diameter floor** (S59: `μ_{1/7}(E) ≥ μ_{1/7}(AP_{D+1})`, PROVED on bounded
   diameter via the Farey roof / THM-637); the **intersection ledger** (S60, all k=8..12); the
   **Part-A finitization** (S61, `V₀ ≤ 1106`).
4. **The mean detour (dead).** Reverse-Markov `μ_{1/7} ⟸ E[maxgap] > 1/7`. **Killed
   (death-star):** the mean is AP-minimized only at a fine scale disjoint from where it is
   non-vacuous — the mean throws the extremal structure to the wrong scale.
5. **Now.** The tail `μ_{1/7}` is irreducible; reduced to the 2-anchor `PA_2 ≥ T_k`
   (boxeph), with the `k=13` leg handled by the double-cover (opus-S139) and — this session —
   the bounded-diameter part of *all* `k` finitized by the anchored-gap subset lemma. The one
   analytic residual is the spread-AP `PA_2`-minimizer rigidity, now a finite `(a,d)` check +
   decorrelation.

## Part 2 — which angles worked, which failed, and the mechanism of each

**Worked — and why.**
- **Sieve / covering / small-mod** (residue conditions). *Mechanism: algebraic decidability* —
  parity and covering arguments always terminate.
- **Coarse reduction** (multi-scale → LRC(≤13)). *Mechanism: self-similar reduction* to an
  already-solved smaller instance.
- **Diameter floor & anchor floor** (subset domination). *Mechanism: monotone domination by a
  computable extremal* (the AP), on a bounded domain.
- **Finitization** (`V₀`, diameter bites, bounded `(a,d)`). *Mechanism: bound the domain where
  the obstruction can be extremal, then it is a finite check.* THE most-reused winning move.
- **Robust substitution** (Part-A few-arc subset). *Mechanism: replace an intractable exact
  object (arc-exploding) by a robust few-arc subset that is provably positive.*
- **Reductions** (reverse-Markov as a valid inequality, 2-anchor, double-cover). *Mechanism:
  trade the hard functional for a dominated / equivalent easier one.*
- **Extremum-identification** (AP = roots of unity, Paley). *Mechanism: the extremal instance
  is the symmetry fixed point* — hand it the exact constants.

**Failed — and why (each a type mismatch).**
- **Route 2 / J–K.** *Wrong invariant:* controls accumulation, need the supremum; isolated
  extrema evade it.
- **`2/7` floor.** *Non-necessity:* a sufficient-not-necessary surrogate had degenerate zeros
  the real object lacks.
- **Paley–Zygmund on the moat.** *Norm mismatch:* the moat is an `L∞` sup (razor-thin,
  signed-cancellation), invisible to `L1/L2` moment methods.
- **Finite covering system.** *Unbounded modulus:* `≡ AP mod lcm` escape families evade every
  fixed covering; "clears at some `q`" ≡ the analytic statement, not a reduction.
- **Single-statistic reductions.** *Order-statistic irreducibility:* the max over gaps can't be
  captured by any one moment (mean/min/origin-gap all fall short).
- **2-point LP.** *Pairwise-featureless:* pair statistics are exactly uniform, so the
  obstruction lives in weight-≥3 correlations.
- **Schur-monotonicity.** *Non-convex landscape:* the functional is resonance-rugged (77%
  of equalizing moves go uphill), no rearrangement structure.

**The one `why` behind all of it (the σ-grading, S67).** Every failure attacked the
**σ-even measure core** with a σ-odd (algebraic/covering) or an averaging (moment) tool — a
*type mismatch*. Every success either used a σ-odd tool on a σ-odd part, or performed a
**finitization** that converts the σ-even part into a σ-odd finite/decidable check. The proof's
whole history is the discovery of that grading and the migration of the σ-even core into
finite checks.

## Part 3 — the mechanisms as reusable structural knowledge

Lifted out of LRC, these are transferable to other hard combinatorial/number-theoretic problems:

1. **Grade by the symmetry involution.** Decompose the problem by an order-two symmetry into an
   *odd/algebraic/decidable* part and an *even/analytic/measure* part; then **match the tool to
   the grade**. Trying an algebraic tool on the measure part (or vice versa) is the modal
   failure. (Generic: any problem with a natural involution — reversal, complement, conjugation.)
2. **Finitize the measure part.** An unbounded analytic obstruction closes if you can bound the
   region where it can be extremal (diameter, height, a few parameters), reducing it to a finite
   check. The single most productive move here; look for it first. (Generic: extremal problems
   where "the worst case is bounded-complexity" — bounded-spread, bounded-degree, bounded-genus.)
3. **Robustness beats exactness.** When the exact object is intractable, find a *robust*
   sub-object (few-arc, low-complexity) that is provably on the right side, and substitute it —
   you pay a slack, you gain decidability. (Generic: replace a fragile witness by a stable one.)
4. **Match the norm-type of the target.** `L∞` (sup, extremal, cancellation-driven) quantities
   are invisible to `L1/L2` (moment/energy) methods; identify whether your target is a sup or an
   average before choosing tools. (Generic: worst-case vs average-case; a moment method never
   proves a razor-thin worst-case bound.)
5. **Verify the invariant, not the power, of an imported tool.** A celebrated external theorem
   that controls invariant `X` is useless if the problem needs `Y`; *isolated* extrema in
   particular evade all accumulation/compactness/limit machinery. (Generic: the fatal Route-2
   lesson — check that the tool measures the thing you need.)
6. **The extremum is the symmetry fixed point.** The extremal instance is almost always the
   most symmetric object (roots of unity, the AP, the regular/Paley structure); locate it to get
   exact constants and to know what "tight" looks like. (Generic: variational problems, isoperimetry.)
7. **Order statistics don't reduce to moments.** A `max`/`min` target keeps information no single
   moment sees; keep the order structure (or dominate it by fewer order terms, as the anchor
   floor does). (Generic: extreme-value vs mean statements.)
8. **A covering reformulation is a reduction only if its modulus is bounded.** If `≡`-extremal
   families force the covering modulus to infinity, "covers at some scale" is a restatement of
   the analytic problem, not progress. (Generic: local-global / sieve arguments — always ask
   whether the local condition is at bounded level.)
9. **Domination is the master tool.** Both wins this thread are dominations: the subset lemma
   dominates *below* by the densest subset's computable extremal; the anchor floor dominates a
   *tail* by fewer gaps. When stuck on an intractable functional, hunt for a computable object
   that dominates it in the direction you need. (Generic: find the right monotone comparison.)

## Ledger

- Concrete (verified): the anchored-gap subset lemma (`gap@a(E,x) ≥ gap@a(F,x)` for `E⊆F`,
  0/27000) finitizing `PA_2 ≥ T_k` on bounded diameter for all `k=8..13` (bites 8/9/11/15/27/68);
  the spread-AP residual is a finite `(a,d)` check + decorrelation.
- Synthesis: the proof arc (Route-2 retired → density floor → reductions → mean-detour dead →
  2-anchor now), the worked/failed mechanism table, and 9 reusable structural principles.
- Builds on / credits: boxeph-S1 (2-anchor reduction), opus-S139 (double-cover, the k=13 leg),
  death-star-S1 (mean dead), MISTAKE-116/117 (Route 2), my S59/S60/S61 (finitizations),
  S65 (barrier), S67 (σ-grading — the unifying why).
- Files: `lrc_anchored_gap_subset_kps_S68.py` (+out).
- Does NOT prove LRC(14). It finitizes another slice of the frontier and names the mechanisms.
