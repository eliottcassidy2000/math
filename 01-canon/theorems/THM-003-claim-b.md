# THM-003: Claim B (Algebraic Companion)

**Type:** Theorem
**Certainty:** 5 — PROVED
**Status:** PROVED
**Last reviewed:** SYSTEM-2026-03-05-S1
**Disputes:** none
**Tags:** #claim-b #independence-polynomial #conflict-graph #a-clique

---

## Statement

For every tournament T and every vertex v:

```
I(Ω(T), 2) − I(Ω(T−v), 2) = 2 Σ_{C∋v} μ(C)
```

where the sum is over all directed odd cycles C through v, and μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2).

---

## Proof / Proof Sketch

Proved via an **A-clique argument** on Ω(T).

Key steps:
1. Partition the independent sets of Ω(T) by whether they include a cycle through v or not.
2. Those not containing any cycle through v contribute exactly I(Ω(T−v), 2).
3. Those containing at least one cycle C through v: group them by their "top cycle" C through v. For each such C, the remaining independent set must avoid all cycles sharing a vertex with C\{v}, which contributes I(Ω(T−v)|_{avoid C\{v}}, 2) = μ(C).
4. Each such group has size μ(C), and each independent set containing C is counted once. The factor of 2 arises from the substitution x=2 and the structure of the clique argument.

Full proof in LaTeX paper §claim_strategies (subsec:scaffold).

## Verification Record

| n | Pairs (T,v) | Failures |
|---|-------------|---------|
| 4 | 256 | 0 ✓ |
| 5 | 5,120 | 0 ✓ |
| 6 | 196,608 | 0 ✓ |

## Verification Scripts

- `04-computation/q009_claim_b.py` — verifies Claim B identity I(Ω(T),2) − I(Ω(T−v),2) = 2·Σ μ(C)
- `04-computation/q009_claim_b_proof.py` — inductive proof verification across n=4,5,6
- `04-computation/q009_alternating_sum.py` — alternating sum identity (key lemma in Claim B)
- `04-computation/verify_all_theorems.py` — comprehensive verification including Claim B

## Notes & History

Claim B being proved while Claim A remains open is the central asymmetry of the paper. The OCF formula requires both: Claim B gives I(Ω(T),2) satisfies the right recurrence; Claim A gives H(T) satisfies the same recurrence; together they imply H(T) = I(Ω(T),2) by induction.
