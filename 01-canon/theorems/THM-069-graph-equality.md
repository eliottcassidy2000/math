# THM-069: Graph Equality for Cycle Deletion

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED (algebraic, trivial proof)
**Status:** PROVED
**Added by:** opus-2026-03-07-S34
**Tags:** #conflict-graph #cycle-deletion #claim-a #structural

---

## Statement

For any tournament T on n vertices, any vertex v, and any directed odd cycle C through v:

**Omega(T) \ N[C] = Omega(T-v)|_{avoid C\{v}}**

as induced subgraphs, where:
- Omega(T) is the conflict graph (vertices = directed odd cycles, edges = pairs sharing a vertex)
- N[C] is the closed neighborhood of C in Omega(T) (C and all cycles sharing a vertex with C)
- Omega(T-v) is the conflict graph of the tournament T-v
- "avoid C\{v}" restricts to cycles vertex-disjoint from all vertices of C except v

---

## Proof

Both sides equal the set D = {directed odd cycles D in T : V(D) cap V(C) = empty}.

**Left to right:** D in Omega(T) \ N[C] means D is a cycle of T that shares no vertex with C. So V(D) cap V(C) = empty, hence D in D.

**Right to left:** D in Omega(T-v)|_{avoid C\{v}} means:
1. D is a directed odd cycle of T-v (equivalently, a cycle of T not using v)
2. D is vertex-disjoint from C\{v}

From (2) and (1): D avoids all vertices of C\{v} and also avoids v, so V(D) cap V(C) = empty. Hence D in Omega(T) \ N[C].

The key observation: **v in C**, so any cycle disjoint from C automatically avoids v.

The adjacency relation (sharing a vertex) is identical in both subgraphs since the vertex sets are the same.

QED.

---

## Corollary: f(C) = 2 * mu(C)

Define:
- f(C) = sum_{S: C in S, S indep in Omega(T)} 2^|S| = 2 * I(Omega(T) \ N[C], 2)
- mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)

Since the two subgraphs are identical: **f(C) = 2 * mu(C)** for all C, v, T.

---

## Consequence for Claim A

The identity H(T) - H(T-v) = 2 * sum_{C ni v} mu(C) (Claim A) can be reformulated as:

**sum_{S indep in Omega(T), S ni some C ni v} 2^|S| = sum_{C ni v} f(C) - (higher-order IE)**

Since Claim B is proved and f(C) = 2*mu(C), this becomes:

**Higher-order inclusion-exclusion terms vanish.**

In other words, Claim B already implies Claim A if and only if the inclusion-exclusion over the set of cycles through v sums to exactly sum_{C ni v} f(C). This is equivalent to saying that independently picking cycles through v gives the correct count -- there is no overcounting from overlapping contributions.

---

## Verification

| n | Method | Pairs tested | Status |
|---|--------|-------------|--------|
| 5 | Exhaustive | 5120 | PASS |
| 7 | 20 random tournaments | 2664 | PASS |
| 9 | 10 random tournaments | 30,511 | PASS |

---

## Scripts

- `04-computation/fC_muC_n7_fast.py` -- Verification at n=7
- `04-computation/fC_muC_graphs.py` -- Graph equality proof and verification
- `04-computation/alpha0_claim_a.py` -- Alpha_0 cancellation
- `04-computation/claim_a_inc_exc.py` -- Inclusion-exclusion structure
