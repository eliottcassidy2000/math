# CONJ-001: Claim A — The Central Open Problem

**Type:** Conjecture
**Certainty:** 3 — VERIFIED (exhaustive n≤6, random n=7; proof open)
**Status:** OPEN
**Last reviewed:** SYSTEM-2026-03-05-S1
**Disputes:** see 02-court/active/DISC-001-mu-bug-vs-verification.md
**Tags:** #claim-a #central-open-problem #vertex-deletion #hamiltonian-paths

---

## Statement

For every tournament T on n vertices and every vertex v:

```
H(T) − H(T−v) = 2 Σ_{C∋v} μ(C)
```

where the sum is over all directed odd cycles C through v, and μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2).

---

## What Is Known

**Proved for n ≤ 5:** The per-path identity (THM-005 + THM-004) implies Claim A for n ≤ 5, because for n ≤ 5 the only odd cycles through v are 3-cycles (the argument is in the paper; see OPEN-Q-001 for why this breaks at n=6).

**Verified exhaustively:**

| n | Pairs (T,v) | Failures |
|---|-------------|---------|
| 4 | 256 | 0 ✓ |
| 5 | 5,120 | 0 ✓ |
| 6 | 196,608 | 0 ✓ |
| 7 (random) | 3,500 | 0 ✓ |

⚠️ **Note:** There is an active discrepancy (DISC-001) regarding a μ computation bug in scripts 6-9. The verification numbers above are from the paper's main verification runs; the MASTER_FINDINGS document reports a bug in ind_poly_at_2_restricted() that may affect some verification logic. The paper's claims are stated to be unaffected. This needs formal resolution. See 02-court/active/DISC-001.

---

## Known Non-Proofs / Failed Approaches

| Approach | Status |
|----------|--------|
| H(T) = B_v + S_v + R_v (exact) | REFUTED (96/256 failures at n=4) |
| S_v + R_v = 2Σμ(C) | REFUTED (144/256 failures at n=4) |
| Per-path identity (3-cycle only) | Fails for n≥6 (2,758/9,126 at n=6) |
| Natural generalization (all cycles) | Overcounts at n=6 |
| Maximal-embedding-only formula | Fails at n=6 |

---

## Proof Strategies (from paper §claim_strategies)

Five organized approaches are listed in the paper. Details need to be extracted into a separate document. **Any Claude working on Claim A should start here.**

---

## Relationship to Other Results

- Claim A + Claim B (THM-003) + induction → OCF formula (THM-002)
- Claim A for n≤5 ← per-path identity ← THM-004 + THM-005
- Claim A for n≥6: per-path identity is insufficient (see OPEN-Q-001, OPEN-Q-004)
