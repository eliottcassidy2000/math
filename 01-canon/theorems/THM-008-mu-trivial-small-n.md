# THM-008: mu(C) = 1 for All 3-Cycles at n <= 5

**Type:** Theorem (proved)
**Certainty:** 5 -- PROVED
**Status:** CANON
**Added by:** opus-2026-03-05-S1
**Source:** FINAL_FINDINGS.md (inbox contribution)
**Tags:** #mu #per-path-identity #small-n #trivial

---

## Statement

For any tournament T on n <= 5 vertices, any vertex v, and any directed 3-cycle C = (v, a, b) through v:

    mu(C) = 1

Consequently, the per-path identity holds trivially for n <= 5: it reduces to #TypeII = #TypeII.

---

## Proof

mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2).

C\{v} = {a, b}. The restricted conflict graph contains only odd cycles of T-v that are vertex-disjoint from {a, b}. Such cycles must use only vertices in V \ {v, a, b}, which has n - 3 vertices.

- n = 3: V\{v,a,b} has 0 vertices. No cycles possible. Conflict graph is empty. I(empty, 2) = 1.
- n = 4: V\{v,a,b} has 1 vertex. No odd cycle on 1 vertex. I(empty, 2) = 1.
- n = 5: V\{v,a,b} has 2 vertices. Minimum odd cycle needs 3 vertices. I(empty, 2) = 1.

In all cases, mu(C) = 1.

---

## Significance

This means the per-path identity for n <= 5 has **zero mathematical content** beyond the algebraic identity THM-004. The weighted sum sum_{TypeII at j} mu(v, P'[j], P'[j+1]) reduces to the unweighted count #TypeII when all mu = 1.

This resolves OPEN-Q-001: the "n=5 mystery" is no mystery at all. 5-cycles exist at n=5, but they are invisible to the per-path identity because that identity only involves 3-cycles, and all 3-cycle mu values are trivially 1.

---

## Verification Scripts

- `04-computation/n5_mu_analysis.py` — proves μ(3-cycle)=1 at n≤5 via structural analysis
- `04-computation/n6_mu_distribution.py` — comparison with n=6 where μ can vary
- `04-computation/q009_proof_n5.py` — verification of per-path identity trivial form at n=5

## Relationship to Other Results

- Explains why per-path identity holds at n <= 5 (MISTAKE-003 context)
- At n = 6, V\{v,a,b} has 3 vertices, which CAN form a 3-cycle, so mu can be 3. See THM-009.
- At n = 7, for 5-cycles C, V\{v,a1,a2,a3,a4} has 2 vertices, so mu(5-cycle) = 1 always.
