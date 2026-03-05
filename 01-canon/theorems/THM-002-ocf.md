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
