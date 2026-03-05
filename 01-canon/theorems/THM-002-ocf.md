# THM-002: Odd-Cycle Collection Formula (OCF)

**Type:** Theorem (conditional)
**Certainty:** 3 — VERIFIED (proof conditional on Claim A)
**Status:** VERIFIED, OPEN (proof incomplete)
**Last reviewed:** SYSTEM-2026-03-05-S1
**Disputes:** none
**Tags:** #ocf #independence-polynomial #conflict-graph #claim-a

---

## Statement

For every tournament T on n vertices (assuming Claim A, CONJ-001):

```
H(T) = I(Ω(T), 2) = Σ_{k≥0} α_k(Ω(T)) · 2^k
```

where Ω(T) is the conflict graph on directed odd cycles of T (vertices = odd cycles, edges = shared vertex), I(G,x) is the independence polynomial of G, and α_k counts independent sets of size k.

---

## Proof / Proof Sketch

The proof is by induction on n using the vertex deletion identity:

**Inductive step:** Pick vertex v. Then:
- H(T) = H(T−v) + (H(T) − H(T−v))
- Claim A gives: H(T) − H(T−v) = 2 Σ_{C∋v} μ(C)
- Claim B gives: I(Ω(T),2) − I(Ω(T−v),2) = 2 Σ_{C∋v} μ(C)
- By induction H(T−v) = I(Ω(T−v),2), so H(T) = I(Ω(T),2). □

**Gap:** The proof of Claim A (see CONJ-001). Everything else is proved.

## Verification Record

| n | Tournaments | Pairs (T,v) | Failures |
|---|-------------|-------------|---------|
| 4 | 64 | 256 | 0 ✓ |
| 5 | 512 | 5,120 | 0 ✓ |
| 6 | 32,768 | 196,608 | 0 ✓ |
| 7 | random | 3,500 | 0 ✓ |

## Notes & History

The formula H(T) = I(Ω(T), 2) is the central new result of the paper. It connects tournament counting to independence polynomials of conflict graphs, which in turn connects to the hard-core lattice gas model at fugacity 2 (see TANGENT T006).

The formula implies Rédei's theorem: I(G, 2) ≥ 1 (the empty set is always independent), and the parity follows from the induction.

**IMPORTANT — Closed-Form Clarification (kind-pasteur-2026-03-05-S5, DISC-002 resolved):**

The formula H(T) = I(Ω(T), 2) IS a valid closed-form identity where Ω(T) uses ALL directed odd cycles of T and I is the PLAIN independence polynomial (alpha_k = number of independent sets of size k, NO mu weights). MISTAKE-004, which claimed the formula was wrong, has been **RETRACTED** — it arose from a computational error in file.txt that used mu weights in place of the independence polynomial coefficients. See DISC-002 (now resolved) and the amended MISTAKE-004 entry for full details.

Exhaustive computational confirmation: H(T) = I(Omega(T), 2) verified for ALL 33,864 tournaments n<=6 (opus-2026-03-05-S2), 0 failures. Further confirmed by T_11: H(T_11) = 95095 matches the OCF expansion exactly (kind-pasteur-2026-03-05-S5).
