# THM-009: Per-Path Identity Failure Characterization at n = 6

**Type:** Theorem (computationally proved, theoretically clear)
**Certainty:** 4 -- VERIFIED + clear theoretical argument
**Status:** CANON
**Added by:** opus-2026-03-05-S1
**Source:** FINAL_FINDINGS.md (inbox contribution)
**Tags:** #per-path-identity #n6 #failure-characterization #mu-weights

---

## Statement

For n = 6, the per-path identity holds for path P' of Ham(T-v) if and only if:

    For every Type-II position j in P', mu(v, P'[j], P'[j+1]) = 1.

Equivalently: the identity FAILS for P' iff some Type-II position (a, b) has mu > 1, iff V \ {v, a, b} contains a directed 3-cycle in T-v.

This is a PERFECT BINARY SEPARATION (verified on all cases at n = 6):
- Path has some mu > 1 TypeII position => identity ALWAYS fails (100%)
- Path has no mu > 1 TypeII position => identity ALWAYS holds (100%)

---

## Proof Sketch

At n = 6, V\{v,a,b} has exactly 3 vertices. The restricted conflict graph for mu(v,a,b) contains odd cycles of T-v on those 3 vertices. The only possible odd cycle on 3 vertices is a directed 3-cycle.

- If the 3 vertices form a 3-cycle in T-v: mu(v,a,b) = I({one cycle}, 2) = 1 + 2 = 3.
- If they do not: mu(v,a,b) = I(empty, 2) = 1.

The per-path identity asks: #TypeII = sum_{TypeII at j} mu(v, P'[j], P'[j+1]).

When all mu = 1, both sides are #TypeII. When some mu = 3, the RHS exceeds the LHS.

---

## Verification Scripts

- `04-computation/n6_mu_distribution.py` — perfect binary separation at n=6: 100% correlation
- `04-computation/per_cycle_identity.py` — per-path identity failure analysis
- `04-computation/q009_proof_n5.py` — contrast with n=5 where all μ=1

## Consequences

Answers OPEN-Q-003: the structural property distinguishing passing/failing triples at n = 6 is whether any Type-II position's complement vertices form a 3-cycle.

The ~30% failure rate at n = 6 corresponds to the probability that a random triple of vertices in a random tournament on 3 vertices forms a 3-cycle (which is 2/8 = 25% for a single triple, but multiple Type-II positions increase the chance of at least one failure).
