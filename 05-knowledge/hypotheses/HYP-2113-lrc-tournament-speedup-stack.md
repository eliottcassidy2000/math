---
id: HYP-2113
status: OPEN roadmap; S583 ranks existing engines and proposes implementation stack
source: codex-2026-06-03-S583
related: [THM-381, THM-398, HYP-2112, HYP-2109, HYP-2108, HYP-2107, HYP-2106, HYP-2105, HYP-2024, THM-395]
---

# HYP-2113: LRC should use a layered tournament speedup stack

## Claim

The next LRC algorithms should not choose between exact interval search,
residue CRT, and tournament quotienting.  They should be stacked:

1. exact observer-source reachability;
2. certificate-gate DAG for Cprime and known exits;
3. fast tournament-class lookup by `A^2` fingerprints where certified;
4. SCC/good-cut protection-core peeling;
5. three-state `L/M/R` middle automata;
6. threshold-colored section/fiber pre-sieves;
7. exact interval or CRT only on the surviving middle cells.

This is a route to make the expensive exact checks rare and proof-relevant,
rather than merely faster brute force.

## Evidence from the repo

- THM-381 is exact: observer is lonely iff the marked observer is a source.
  This turns LRC into source-target reachability in a marked tournament movie.
- S535 shows raw phase classes mix good and bad states, but source-deleted,
  observer-source, gap-threshold, blocker-shadow, and apex-boundary fibers are
  pure in bounded audits and certify every tested speed set.
- S539 shows non-runner vertices matter: fixed sections, half-sections, and
  boundaries become useful tournament vertices once danger occupancy is colored.
- The `A^2` row fingerprint is a fast isomorphism lookup through verified
  tournament sizes `n<=9` and a strong cache key beyond that with fallback.
- THM-354/THM-395 turn SCC defect, good cuts, and backward wedges into a
  protection-core language.  This is the tournament analogue of LRC endpoint
  protectors and private pivots.
- HYP-2112 adds the exact circuit-to-gap functional `G(v)=Phi(C)`: on
  multiples of `n`, Cprime is equivalent to `Phi>0`, while `ker Phi` is the
  worry set.  This gives the Cprime gate a scalar terminal call instead of a
  raw residual search.
- S577's additive-circuit A/B probe is a downstream residual gate: a 3-term
  relation `v_c=v_a+v_b` makes the fold a shield on runner `c`'s clock, while
  sampled circuit-free rows sit well above `delta=1/(k+1)`.  Summand-graph
  folds should therefore be labelled exits in the certificate DAG.
- HYP-2109 keeps the middle state visible before projecting to a binary
  tournament, which matches HYP-2108's all-middle endpoint-cover residual and
  HYP-2107's live residue states.

## Proposed architecture

The proof search should use three layers.

**Layer A: exact target layer.**  Build the observer-source marked tournament
movie, but classify cells by source-target reachability rather than by raw
runner half-turn classes.

**Layer B: certificate quotient layer.**  Run the Cprime/certificate gate DAG:
no-multiple `t=1/n`, dominance, B', Lemma C, bounded CRT, endpoint-cover
positivity, the HYP-2112 `Phi` gap functional, additive fold/shield exits,
cheap-pair/sheaf exits.

**Layer C: middle-residual layer.**  Only unresolved cells are allowed to remain
terminal `M`.  Those cells carry endpoint-owner, interval-length, midpoint
phase, residue, or section labels.  A no-closed-middle lemma should then kill
the residual using HYP-2108 or HYP-2107.

## Assumption challenge

Tournament vertices should be chosen by the predicate being preserved:

- observer-source reachability uses observer-marked classes;
- Cprime residuals use safe components, endpoint owners, or residue states;
- endpoint-pressure uses proof obligations or owner clauses;
- section functors use fixed circle sections and section boundaries;
- certificate sheaves use germs and gluing failures.

The stack deliberately destroys raw continuous time early, but only after
recording labels needed by the next proof gate.  Dropping those labels too soon
is the known failure mode: raw phase tournaments mix lonely and non-lonely
states.

## First implementation target

Build `lrc_marked_source_speedup_stack.py`:

1. enumerate event cells for a speed row;
2. emit observer-source class and threshold-colored section class;
3. compute cheap TDA features: score, `c3`, SCC defect, terminal `M` count;
4. use `A^2` fingerprints as cache keys, with canonical fallback;
5. route cells through the Cprime/certificate gates and additive fold/shield
   exits;
6. print the remaining middle cells with endpoint-owner/residue labels.

The desired output is not just "safe/unsafe"; it is a residual certificate:
either a source target, a known positive-measure gate, a private pivot, or a
small labelled middle circuit ready for HYP-2112/HYP-2108/HYP-2107.
