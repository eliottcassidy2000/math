---
id: HYP-3563
title: RELATIONS, NOT THINGS -- the relational reframe extending the coboundary lens. The project's objects ARE relations (a tournament = the dominance relation; the LRC = the danger relation D(v,t)=[||vt||<1/14]; the metagraph = the arc-flip relation), and by klein THM-588 (no first-order/thing invariant, only the quadratic 3-cycle) the ONLY invariant is the relation's SELF-COMPOSITION D.D^T (the pair-correlation = the 2nd moment, THM-589/579). Essentiality = rank>1 (not separable); a COBOUNDARY = a rank-1/separable relation f(v)g(t) trivializable by a 'thing'/potential. A disproof = the obstruction relation is a coboundary; a proof = it is essential. The essentialities that forbid trivialization: the BILINEAR product v*t (LRC) and the CYCLIC 3-cycle (metagraph). The reference-collapse (klein-S5) = the change-of-base relation: S_n-collapse (clean, set-independent, CV(H)^2~2/n) vs Z_14-collapse (dirty, CV(N_R)^2 unbounded); the proof = MANUFACTURING the clean base relation (Gamma_0(N)).
status: SYNTHESIS + verified grounding (LRC safe relation full-rank=essential; the Gram D.D^T = the pair 2nd moment; metagraph 3-cycle the unique relational invariant). A reframe/ontology, not a new proof.
source: mac-mini-2026-06-29-S24
related:
  - HYP-3562  # the measure of the obstruction / the coboundary lens (this is its relational extension)
  - THM-588   # klein: no first-order invariant, only the quadratic 3-cycle (= forced relational)
  - THM-589   # the 2nd moment = the relation's self-composition (W(n), CV~2/n)
  - THM-579   # CV(N_R)^2 = the LRC pair-correlation (the relation's Gram)
  - HYP-3553  # Gamma_0(N) = the change-of-base relation (the clean reference-collapse)
external: arXiv:2507.05905 (Han-Lee); category theory (spans/correspondences); klein-S5 reference-collapse
results:
  - 04-computation/relations_not_things_danger_relation_macmini_20260629.py
  - 05-knowledge/results/relations_not_things_danger_relation_macmini_20260629.out
reflections:
  - lrc14-first-obstruction-cocycle-generation-codex-s259.md
---

# HYP-3563 -- relations, not things

## The reframe
Stop treating the project's nouns as objects; they are RELATIONS, and the proof is a statement about the
relation, not its endpoints:
| "thing" view | RELATION view |
|---|---|
| a tournament (a set with arcs) | the **dominance relation** `a -> b` (only invariant: cyclicity, THM-588) |
| a lonely time `t` | a full-safe COLUMN of the **danger relation** `D(v,t)=[||vt||<1/14]` (derived) |
| the metagraph (vertices) | the **arc-flip relation** between iso classes |
| the floor (a number) | the SELF-COMPOSITION `D.D^T` of the relation (the pair-correlation) |
The "things" (lonely point, SC tournament, the transitive order) are DERIVED from the relations; the
relations are primary.

## Why the project is FORCED relational (klein THM-588)
THM-588: the metagraph has NO first-order invariant (`mult(1)=0`) and EXACTLY ONE quadratic (the 3-cycle
count, `mult(2)=1`). The 3-cycle is a RELATION among three vertices -- the failure of transitivity. So
there is no "thing" invariant of individual arcs; the only invariant is a RELATION (a pair/triple). The
transitive tournament (a linear ORDER, the maximal "thing") is the trivial point (0 cyclicity); all
content is relational. VERIFIED: 3-cycle counts over iso classes = `{0,1,2}` (n=4), `{0..5}` (n=5),
`{0..8}` (n=6) -- the single relational invariant carries the whole metagraph.

## The second moment IS the relation, composed
The only invariant being quadratic means: study the relation's SELF-COMPOSITION. `D.D^T` (the Gram) is
the pair-overlap -- its diagonal the first moment (`#safe times`, a "thing" count) and its OFF-DIAGONAL
the pair-correlation (the relation, the 2nd moment, THM-588's quadratic). VERIFIED (LRC `{2,3,4}`): the
`3x3` Gram has the pair-overlaps `14,16,16` off-diagonal. So `CV(N_R)^2` (THM-579) and `Var(H)=W(n)`
(THM-589) are BOTH the relation composed with itself -- the only thing there is to bound.

## Essential vs coboundary = rank > 1 vs separable
A COBOUNDARY is a SEPARABLE (rank-1) relation `D(v,t)=f(v)g(t)` -- trivializable by two "thing"-functions
(a potential). An ESSENTIAL relation has rank `> 1`. VERIFIED: the LRC safe relation is FULL RANK
(`3,3,5` for the tested cases) -- essential, not a coboundary. **A disproof would need the danger
relation to be a coboundary (separable => covers everything); the BILINEAR product `v*t` inside `||vt||`
forbids separation** (the multiplicative/Littlewood structure, HYP-3551, is the essentiality). On the
metagraph the essentiality is the CYCLIC 3-cycle (non-transitivity); on the LRC it is the bilinear `v*t`.

## The reference-collapse is the change-of-base relation (klein-S5)
klein-S5: `CV(H)^2` (the metagraph 2nd moment, `S_n`-collapse) is clean / set-independent / `~2/n`;
`CV(N_R)^2` (the LRC 2nd moment, `Z_14`-collapse) is dirty / unbounded. The difference is the COLLAPSE
RELATION (the change of base): `S_n` (the transitive reference, no vanishing fiber) cleans the relation;
`Z_14` does not. So the proof is a RELATIONAL move: **manufacture the change-of-base relation that cleans
the LRC's second moment** -- which is exactly the `Gamma_0(N)` congruence (HYP-3553), the relation that
makes the moment depend only on `N` (set-independent = the essentiality is topological, HYP-3562).

## New frames (all relational)
- **Span/correspondence:** `D` is a span `Speeds <- D -> Time`; compose, transpose, pull back -- the
  calculus of correspondences replaces the calculus of points.
- **Composition = 2nd moment:** `D . D^op` is the only invariant; the proof is a bound on a composition.
- **Rank = essentiality:** rank-1 = coboundary = trivial; rank `> 1` = essential = the obstruction.
- **Morphisms are primary:** `R` (complement), `T` (vertex-add/Hecke), arc-flip GENERATE everything; the
  objects are their fixed points / orbits. The 2-adic involution `R` on three faces (floor/metagraph/
  witness, klein-S5) is ONE relation wearing three representations.

## What it buys
A clean ontology in which the proof target is sharp: not "a lonely point exists" (a thing) but "the danger
relation is essential" (rank `> 1`, not a coboundary), reduced by klein THM-588 to "the relation's
self-composition (the 2nd moment) is bounded under the right change-of-base relation (`Gamma_0(N)`)".
The disproof is the one thing an essential relation cannot be -- separable. Relations, not things.
