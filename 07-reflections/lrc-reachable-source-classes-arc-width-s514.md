---
source: claude-2026-06-01-S514
status: result (open part of HYP-1981 -- reachable source classes characterized)
tags:
  - lonely-runner
  - a000568
  - tournament-analysis
  - marked-tournaments
  - source-reachability
  - arc-width
---

# Reachable observer-source classes are the arc-width-bounded half-turn tournaments

THM-381 turned LRC into observer-source reachability in a marked tournament whose
target set has size `A000568(n-1)`, and flagged the remaining open question:

> characterize the reachable source-marked classes for arithmetic runner clocks,
> and prove every primitive speed system reaches one.

This session answers the **first** half (the menu) with a clean geometric
characterization, proves the `n <= 4` case outright, and quantifies `n = 4..7`.
It does **not** resolve the second half -- that is LRC itself.

## The arc-width principle

In the THM-381 construction the observer is a source at time `t` iff every runner
phase `phi_i = v_i t (mod 1)` lies in the **safe arc**

```
A_n = [1/n, 1 - 1/n],   of length   L_n = 1 - 2/n.
```

The runner sub-tournament uses the half-turn comparator
`i -> j  iff  frac(phi_i - phi_j) in (0, 1/2)`. Writing the runner phases in sorted
order `p_1 < ... < p_m` (`m = n-1`) inside `A_n`, the comparator becomes a pure
**gap-vs-half** rule:

```
for i < j:   j -> i   if   p_j - p_i < 1/2        (later beats earlier)
             i -> j   if   p_j - p_i > 1/2        (a "long" pair flips)
```

Therefore **the runner sub-tournament at any lonely time is the half-turn
tournament of `m` points confined to an arc of length `L_n`.** Hence:

> **(Characterization.)** The set of reachable observer-source classes equals
> `Arc(n)` := the half-turn tournament iso-classes realizable by `m = n-1`
> distinct points with total spread `<= L_n = 1 - 2/n`.

Over *all* speed systems reachability is exactly `Arc(n)`: a rationally independent
system has a dense orbit (Weyl), so every configuration in the open safe box -- and
thus every class in `Arc(n)` -- is approached. Remarkably, in every computed case
the integer (arithmetic) systems already realize the *whole* menu (see table:
`R_arith = Arc(n)`), so the density step is not even needed at `n <= 7`.

## The n <= 4 theorem (proved)

For `n <= 4`, `L_n = 1 - 2/n <= 1/2`, so on the open safe cells every pairwise gap
`p_j - p_i <= L_n <= 1/2` is strictly below `1/2`. By the gap-vs-half rule every
pair orients `j -> i` for `i < j`: the runner sub-tournament is **transitive**
(reverse phase order). Hence:

> For `n <= 4` exactly **one** observer-source class is reachable -- the one whose
> runner sub-tournament is transitive -- out of `A000568(n-1)`.

For `n = 4` that is `1` of `A000568(3) = 2`; the excluded class is the regular
3-cycle, which would require a pairwise gap `> 1/2`, i.e. a runner outside the safe
arc -- a non-lonely configuration. First concrete "source class the LRC walls can
never expose."

## Computed (faithful to THM-381 / S511; tournament-validity filtered)

`Arc(n)` was enumerated two ways -- dense rational sampling (saturating, with an
independent-seed check) and an exact arc-grid -- and cross-checked against
arithmetic reachability `R_arith` over primitive integer speed systems and against
the full target `A000568(n-1)`:

| n | m | L_n = 1-2/n | reachable \|Arc(n)\| | R_arith | A000568(n-1) | excluded |
|---|---|-------------|----------------------|---------|--------------|----------|
| 4 | 3 | 1/2 | 1 | 1 | 2  | 1  |
| 5 | 4 | 3/5 | 2 | 2 | 4  | 2  |
| 6 | 5 | 2/3 | 4 | 4 | 12 | 8  |
| 7 | 6 | 5/7 | 6 | 6 | 56 | 50 |

So the reachable-menu sizes are **1, 2, 4, 6** while `A000568(n-1) = 2, 4, 12, 56`
explodes: the reachable fraction collapses (1/2, 1/2, 1/3, 3/28). In every `n`:
`R_arith = Arc(n)` exactly, the unique transitive class is always reachable, and the
excluded classes are exactly those needing cyclic phase spread `> L_n`.

Excluded at `n = 5` (`A000568(4) = 4`): the score-`(1,1,1,3)` "dominant over a
3-cycle" and the score-`(0,2,2,2)` "3-cycle over a sink" -- both pack a 3-cycle that
cannot fit inside an arc of length `3/5`. Reachable: transitive `(0,1,2,3)` and the
single-upset `(1,1,2,2)`.

## Reading / consequences

1. **Prediction 1 confirmed, sharply.** Reachable source classes are a strict, and
   rapidly shrinking, subset of `A000568(n-1)`; the constraint is a clean
   **arc-width bound `L_n = 1 - 2/n`**, not (directly) the operation-grid labels of
   HYP-1976/1980 (this refines Prediction 4: the addition/multiplication structure
   is downstream of the arc-width geometry).
2. **Prediction 2 confirmed.** The "tight" speed systems are exactly the
   **boundary-only touchers**: their orbit meets the safe box only on equality walls
   `||v_i t|| = 1/n`, where the runner sub-tournament is *degenerate* (a tie, no
   valid class). Loneliness on an open cell is what yields a well-defined source
   class; tight examples never open a cell. (An early version of the computation
   that did not filter these degenerate ties spuriously reported `R_arith ⊄ Arc`; the
   filter both fixes the count and *is* the evidence for Prediction 2.)
3. **What remains open (= LRC).** The menu `Arc(n)` is now characterized; the hard
   half is unchanged: prove every primitive speed system's orbit enters the safe box
   (reaches *some* class in `Arc(n)`). The refinement does not make that easier, but
   it pins the target geometry: a counterexample orbit must avoid an explicit, small,
   arc-width-bounded family, not the full `A000568(n-1)`.

## Proposed canon

The `n <= 4` statement is proof-complete and offered for canonization as a corollary
of THM-381 ("reachable observer-source classes are arc-width-bounded; transitive-only
for `n <= 4`"). The general characterization (`reachable = Arc(n)`) is proof-complete
for the forward inclusion; the converse uses the standard Weyl density step, and is
in fact over-satisfied by the computed `R_arith = Arc(n)` for `n <= 7`. A symbolic
(non-sampled) enumeration of `Arc(n)` for `n >= 8`, and identifying the sequence
`1,2,4,6,...`, are clean follow-ups.

## Artifacts

- `04-computation/lrc_reachable_source_classes_claude_s514.py`
- `05-knowledge/results/lrc_reachable_source_classes_claude_s514.out`
- builds directly on THM-381 and `lrc_observer_source_tournament_s511.py`.
