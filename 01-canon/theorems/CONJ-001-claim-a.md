# CONJ-001: Claim A — PROVED (via Grinberg-Stanley OCF)

**Type:** Theorem (was Conjecture; now proved)
**Certainty:** 5 — PROVED for all n
**Status:** PROVED
**Last reviewed:** kind-pasteur-2026-03-05-S12
**Disputes:** RESOLVED (DISC-001 moot; OCF proved externally)
**Tags:** #claim-a #proved #vertex-deletion #hamiltonian-paths #grinberg-stanley

---

## Statement

For every tournament T on n vertices and every vertex v:

```
H(T) − H(T−v) = 2 Σ_{C∋v} μ(C)
```

where the sum is over all directed odd cycles C through v, and μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2).

---

## PROOF (via Grinberg-Stanley + Claim B)

**Claim A is a COROLLARY of OCF (THM-002) + Claim B (THM-003).**

**Step 1:** OCF states H(T) = I(Omega(T), 2) for all tournaments T.
OCF is now PROVED for all n — it follows from Theorem 1.39 & Lemma 6.5 of arXiv:2307.05569
(Grinberg & Stanley, 2023), restated as Corollary 20 in arXiv:2412.10572 (Irving & Omar, 2024),
combined with the observation that for tournaments, the complement D̄ = D^op (converse),
and ham(D^op) = ham(D) by path reversal. See THM-002 for full details.

**Step 2:** Claim B (THM-003) states:
I(Omega(T), 2) - I(Omega(T-v), 2) = 2 * sum_{C through v} mu(C).
This is proved via the A-clique argument.

**Step 3:** Since OCF holds for both T and T-v:
H(T) - H(T-v) = I(Omega(T),2) - I(Omega(T-v),2) = 2 * sum_{C through v} mu(C). QED.

---

## Verification Record (pre-proof)

| n | Pairs (T,v) | Failures |
|---|-------------|---------|
| 4 | 256 | 0 |
| 5 | 5,120 | 0 |
| 6 | 196,608 | 0 |
| 7 (random) | 3,500 | 0 |
| 8 (THM-015) | 134,217,728 | 0 |

All computational verifications are consistent with the now-proved result.

---

## Historical Failed Approaches (no longer relevant)

| Approach | Status |
|----------|--------|
| H(T) = B_v + S_v + R_v (exact) | REFUTED (96/256 failures at n=4) |
| S_v + R_v = 2Σμ(C) | REFUTED (144/256 failures at n=4) |
| Per-path identity (3-cycle only) | Fails for n>=6 (2,758/9,126 at n=6) |
| Natural generalization (all cycles) | Overcounts at n=6 |
| Maximal-embedding-only formula | Fails at n=6 |

These approaches attempted to prove Claim A directly. The actual proof goes through OCF
(Grinberg-Stanley) + Claim B (A-clique), bypassing per-path methods entirely.

---

## Relationship to Other Results

- **OCF (THM-002, PROVED)** + **Claim B (THM-003, PROVED)** → **Claim A (THIS, PROVED)**
- OCF proved by Grinberg-Stanley (arXiv:2307.05569, arXiv:2412.10572)
- OCF independently verified at n <= 8 by polynomial identity (THM-015, THM-018)
- Claim A for n<=5 also follows from per-path identity (THM-004 + THM-005)

## Key References

- Grinberg & Stanley, "The Rédei-Berge symmetric function of a directed graph", arXiv:2307.05569 (2023) — original OCF proof
- Irving & Omar, "Revisiting The Rédei-Berge Symmetric Functions via Matrix Algebra", arXiv:2412.10572 (2024) — Corollary 20 restates OCF via matrix algebra
