# THM-020: Real Roots of I(Omega(T), x)

**Type:** PROVED for n<=8, DISPROVED at n>=9 (THM-025 counterexample).
**Certainty:** 5 — PROVED for n<=8 via Chudnovsky-Seymour; DISPROVED at n>=9.
**Status:** PROVED for n<=8. DISPROVED at n>=9 (THM-025).
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

## Status at n>=9 — DISPROVED

**DISPROVED (opus-2026-03-06-S18, THM-025):** I(Omega(T), x) does NOT always have real roots at n=9. Explicit counterexample: tournament with score sequence [1,1,3,4,4,4,6,6,7], giving I(Omega, x) = 1 + 94x + 10x^2 + x^3. Newton's inequality fails at k=2 (100 < 141).

The n<=8 result (this theorem) remains valid. The claw-free bound n<=8 is sharp.

---

## Verification Record

| n | Method | Poly degree | Real-root failures |
|---|--------|-------------|-------------------|
| 5 | exhaustive (1024) | 1 | 0 |
| 6 | 500 random | 1-2 | 0 |
| 7 | 200 random (full Omega) | 1-2 | 0 |
| 8 | 50 random (3-cycle subgraph) | 2 | 0 |
| 9 | 50 random (3-cycle subgraph) | 2-3 | 0 |
| 10 | 50 random (3-cycle subgraph) | 3 | 0 |
| 12 | 30 random (3-cycle subgraph, ~55 cycles) | 3-4 | 0 |
| 15 | 10 random (3-cycle subgraph, ~110 cycles) | 5 | 0 |
| 20 | 1 random (3-cycle subgraph, ~298 cycles) | 6 | 0 |

**Ultra-log-concavity also verified** (0 failures) at all tested n. alpha_k/C(m,k) is log-concave where m = |V(Omega_3)|. (kind-pasteur-S16)

**Verification scripts:**
- `04-computation/real_roots_deep.py` (real-rootedness checks at various n)
- `04-computation/omega_claw_fast.py` (claw-freeness checks)
- `04-computation/verify_all_theorems.py` (THM-021 section, discriminant checks)

---

## Graph Properties of Omega(T)

| Property | n<=5 | n=6 | n=7 | n=8 | n=9 | n=10 |
|----------|------|-----|-----|-----|-----|------|
| Line graph | YES | NO (53% fail) | - | - | - | - |
| Chordal | YES | NO (72/2000 fail) | - | - | - | - |
| Perfect | YES | YES | YES | NO (53.8% fail) | - | - |
| Comparability | YES | YES (200/200) | NO (1/100 fail) | NO (58/100 fail) | - | - |
| Quasi-line | YES | YES | YES | NO (49% fail) | NO (90% fail) | NO (100%) |
| Claw-free | YES | YES | YES | YES (trivial) | NO (86% fail) | NO (98%) |
| S_{2,1,1}-free | YES | YES | YES | YES | YES (0/100 fail) | NO (92% fail) |
| Real roots | YES | YES | YES | YES | **NO** (THM-025) | - |

The hierarchy: Line graph => Claw-free => S_{2,1,1}-free. Real roots fail exactly when claw-freeness fails (n=9). The Chudnovsky-Seymour bound is SHARP for tournaments.

---

## Related Results

- **Chudnovsky-Seymour (2007):** Claw-free => all real roots for I(G,x).
- **Jerrum-Patel (2026):** Zero-free regions for subdivided claw-free graphs. May extend real-root results beyond strict claw-freeness.
- **THM-019:** Omega perfectness fails at n=8 but real roots persist.
- **Scott-Sokal theory:** Studies zeros of I(G,x) for hard-core model. Evaluation at x=2 (fugacity=2) is always in the zero-free region when all roots are real and negative.

---

## Potential Proof Approaches for n>=9

1. **~~Line graph hypothesis~~:** REFUTED at n=6. K5-e (Beineke forbidden subgraph) appears in 53% of n=6 tournaments. Heilmann-Lieb does not apply. (kind-pasteur-S14b)

2. **Subdivided-claw-freeness (PARTIALLY VIABLE):** Omega(T) is S_{2,1,1}-free at n=9 (0/100 failures) despite having claws (86%). But S_{2,1,1} appears at n>=10 (92%). The hierarchy is: claw-free (n<=8) -> S_{2,1,1}-free (n<=9) -> ??? (n>=10). Jerrum-Patel (2026, J. London Math. Soc.) prove all real roots for bounded-degree S_{a,b,c}-free graphs, but this requires a FIXED subdivided claw avoided for ALL n. The increasing chain of avoided subgraphs does not directly fit their framework.

3. **~~Quasi-line graphs~~:** FAILS at n=8 (49%). Not viable. (kind-pasteur-S14b)

4. **Direct structural argument:** Prove real roots from the specific combinatorial structure of Omega(T) (e.g., using the Grinberg-Stanley symmetric function framework). The independence number alpha(Omega) <= floor(n/3) is small, which limits the degree of I(Omega,x).

5. **Degree-dependent argument:** At each n, Omega(T) avoids some subdivided claw S_n and has max degree Delta_n. Both grow with n, but perhaps the Jerrum-Patel zero-free region grows fast enough to always include x=2.

---

## References

1. Chudnovsky, M., Seymour, P. "The roots of the independence polynomial of a clawfree graph." J. Combin. Theory Ser. B 97 (2007), 350-357.
2. Jerrum, M., Patel, V. "Zero-free regions for the independence polynomial on restricted graph classes." J. London Math. Soc. (2026).
3. Grinberg, D., Stanley, R.P. arXiv:2412.10572, Corollary 20 (proves H(T) = I(Omega(T), 2)).
4. Bezakova, I., Galanis, A., Goldberg, L.A., Stefankovic, D. arXiv:2404.07615 (2024). Hard-core model dichotomy.
5. Heilmann, O.J., Lieb, E.H. "Theory of Monomer-Dimer Systems." Commun. Math. Phys. 25 (1972), 190-232.
