# THM-020: Real Roots of I(Omega(T), x)

**Type:** Theorem (n<=8), Conjecture (general n)
**Certainty:** 5 — PROVED for n<=8 via Chudnovsky-Seymour; computationally verified
**Status:** PROVED for n<=8. Open for n>=9.
**Added by:** kind-pasteur-2026-03-05-S13
**Tags:** #omega #independence-polynomial #real-roots #claw-free #chudnovsky-seymour

---

## Statement

**Theorem (n<=8):** For every tournament T on n<=8 vertices, all roots of the independence polynomial I(Omega(T), x) are real and negative.

**Conjecture (general n):** For every tournament T on n vertices, all roots of I(Omega(T), x) are real and negative.

---

## Proof for n<=8

**Step 1:** Omega(T) is claw-free for n<=8.

A claw (K_{1,3}) in Omega(T) requires a central odd cycle C0 adjacent to three pairwise non-adjacent cycles C1, C2, C3. "Non-adjacent" means vertex-disjoint; "adjacent to C0" means sharing a vertex. So C1, C2, C3 are pairwise vertex-disjoint odd cycles, each sharing at least one vertex with C0. Since each odd cycle has >= 3 vertices and C1, C2, C3 are pairwise disjoint, they use >= 9 vertices total. Therefore no claw exists for n<=8.

**Step 2:** By the Chudnovsky-Seymour theorem (J. Combin. Theory Ser. B, 2007), the independence polynomial of any claw-free graph has only real roots.

**Step 3:** Since all coefficients of I(Omega(T), x) are non-negative (they count independent sets) and I(Omega(T), 0) = 1 > 0, all real roots must be negative.

**Corollary:** I(Omega(T), x) > 0 for all x >= 0 and all tournaments T on n<=8 vertices. In particular, H(T) = I(Omega(T), 2) >= 1, giving an alternative proof of Redei's theorem for n<=8.

---

## Status at n>=9

At n=9, Omega(T) is NOT always claw-free (90% of random tournaments have a claw in Omega). The Chudnovsky-Seymour argument fails. However, real-rootedness may still hold for a different reason.

**Open question:** Does I(Omega(T), x) have all real roots for n>=9?

If true, this would imply:
- H(T) > 0 for all tournaments (alternative proof of Redei for all n)
- Log-concavity of independence polynomial coefficients for all n
- Unimodality of the coefficient sequence

---

## Verification Record

| n | Method | Real-root failures |
|---|--------|-------------------|
| 5 | exhaustive (1024) | 0 |
| 6 | 500 random | 0 |
| 7 | 1000 random (THM-019) | 0 |
| 8 | 5 random (full Omega) | 0 |

---

## Graph Properties of Omega(T)

| Property | n<=5 | n=6 | n=7 | n=8 | n=9 |
|----------|------|-----|-----|-----|-----|
| Chordal | YES | NO (72/2000 fail) | - | - | - |
| Perfect | YES | YES | YES | NO (53.8% fail) | - |
| Comparability | YES | YES (200/200) | NO (1/100 fail) | NO (58/100 fail) | - |
| Claw-free | YES | YES | YES | YES (trivial) | NO (90% fail) |
| Real roots | YES | YES | YES | YES | OPEN |

The hierarchy is: Chordal => Perfect => Comparability (no) ... Claw-free is independent.
Omega(T) is always claw-free for n<=8 (vertex counting), which suffices for real roots.

---

## Related Results

- **Chudnovsky-Seymour (2007):** Claw-free => all real roots for I(G,x).
- **Jerrum-Patel (2026):** Zero-free regions for subdivided claw-free graphs. May extend real-root results beyond strict claw-freeness.
- **THM-019:** Omega perfectness fails at n=8 but real roots persist.
- **Scott-Sokal theory:** Studies zeros of I(G,x) for hard-core model. Evaluation at x=2 (fugacity=2) is always in the zero-free region when all roots are real and negative.

---

## References

1. Chudnovsky, M., Seymour, P. "The roots of the independence polynomial of a clawfree graph." J. Combin. Theory Ser. B 97 (2007), 350-357.
2. Jerrum, M., Patel, V. "Zero-free regions for the independence polynomial on restricted graph classes." J. London Math. Soc. (2026).
3. Grinberg, D., Stanley, R.P. arXiv:2412.10572, Corollary 20 (proves H(T) = I(Omega(T), 2)).
