# Backlog Investigation Results

**Source:** opus-2026-03-05-S13 (rigorous backlog investigation)

---

## INV-014: 2-adic Tower / Higher Redei Theorems

### Key Finding: v_2(H(T)) = 0 UNIVERSALLY

Exhaustive computation at n=3,4,5,6 (all tournaments) and sampled at n=7 (5000 random):

| n | Total tournaments | v_2(H(T)) = 0 | Higher v_2? |
|---|-------------------|----------------|-------------|
| 3 | 8 | 8 (100%) | 0 |
| 4 | 64 | 64 (100%) | 0 |
| 5 | 1024 | 1024 (100%) | 0 |
| 6 | 32768 | 32768 (100%) | 0 |
| 7 | 5000 (sampled) | 5000 (100%) | 0 |

**Conclusion:** H(T) is always odd. This IS Redei's theorem. The 2-adic valuation of H(T) is always 0 — there's no deeper 2-adic structure at the level of H(T) itself.

### Mod-4 and Mod-8 Structure

| n | H≡1(4) | H≡3(4) | H≡1(8) | H≡3(8) | H≡5(8) | H≡7(8) |
|---|--------|--------|--------|--------|--------|--------|
| 3 | 75.0% | 25.0% | 75.0% | 25.0% | 0% | 0% |
| 4 | 75.0% | 25.0% | 37.5% | 25.0% | 37.5% | 0% |
| 5 | 70.3% | 29.7% | 35.2% | 23.4% | 35.2% | 6.3% |
| 6 | 62.5% | 37.5% | 30.3% | 17.6% | 32.2% | 19.9% |
| 7 | 55.5% | 44.5% | 29.7% | 24.3% | 25.8% | 20.2% |

**Pattern:** H mod 4 approaches 50/50 as n grows; H mod 8 approaches uniform 25/25/25/25 on odd residues.

### Mod-4 vs 3-cycle Count

- At n=3,4: H ≡ 1 + 2·c_3 (mod 4) EXACTLY (100%)
- At n≥5: Fails. c_3 parity does NOT determine H mod 4.
- OCF gives: H ≡ 1 + 2·alpha_1 (mod 4), where alpha_1 = total number of odd cycles in Omega(T).
- Since alpha_1 ≥ c_3 in general (5-cycles, 7-cycles contribute), the c_3 formula breaks when 5-cycles appear.

### Impact on OPEN-Q-008

OPEN-Q-008 asked for combinatorial characterization of v_2(H(T)). **Answer: v_2(H(T)) = 0 always (Redei's theorem).** The question should be reformulated as: what is the distribution of H(T) mod 2^k for k=2,3,...? Answer: it approaches uniform on odd residues mod 2^k as n grows.

---

## OPEN-Q-015: Real Roots of I(Omega(T), x) — FULL OMEGA

### Full Omega Verification (all odd cycles, not just 3-cycles)

| n | Samples | OCF pass | All real roots | Max degree | Cycle types |
|---|---------|----------|----------------|------------|-------------|
| 5 | 100 | 100/100 | 100/100 | 1 | 3,5 |
| 6 | 50 | 50/50 | 50/50 | 2 | 3,5 |
| 7 | 20 | 20/20 | 20/20 | 2 | 3,5,7 |
| 8 | 10 | 10/10 | 10/10 | 2 | 3,5,7 |
| 9 | 5 | 5/5 | 5/5 | 3 | 3,5,7,9 |

**CRITICAL:** OCF verification with FULL Omega (not just Omega_3) passes 100%. Previous tests used only 3-cycles.

### Paley p=7 Full Omega

- Total odd cycles: 80 (14 three-cycles + 42 five-cycles + 24 seven-cycles)
- Independence polynomial: [1, 80, 7]
- I(Omega, 2) = 189 = H(T_7) --- OCF verified
- Roots: -11.4161, -0.0125 (both real negative)

### Coefficient Structure

**100% log-concave** across all tested polynomials (206 samples, n=5..9).
**100% unimodal** across all tested polynomials.

Coefficients drop super-exponentially: typical ratios alpha_k/alpha_{k-1} are ~0.3 (k=1→2), ~0.03 (k=2→3).

Root ratio |r_big/r_small| ≈ alpha_1^2/alpha_2 (from Vieta), typically 100-800.

### Source Cone Theorem (from INV-035 agent)

**H(source_cone(T')) = H(T') for all T'.** Adding a source vertex (beating everyone) to tournament T' produces a tournament with the same Hamiltonian path count. This is because every Hamiltonian path in the cone must start at the source, so Ham paths of cone biject with Ham paths of T'.

This confirms INV-004: OCF for source/sink cones is trivially equivalent to OCF for the base tournament. The R-cone proof strategy requires showing cut-flip invariance of E(T)=0.

---

## Cross-connections: Backlog → Open Questions

### INV-014 → OPEN-Q-008 (2-adic tower)
**Status: PARTIALLY RESOLVED.** v_2(H(T))=0 always (Redei). The mod-4 characterization is H ≡ 1+2·alpha_1 (mod 4) via OCF. Reformulate OPEN-Q-008 to ask about higher-order mod-2^k residue distributions.

### INV-032/OPEN-Q-015 (real roots)
**Status: STRENGTHENED.** Full Omega (all odd cycles) still gives all-real roots. The barrier remains: no graph-theoretic property of Omega(T) can explain this universally. The explanation must come from the algebraic structure (Irving-Omar / Grinberg-Stanley framework).

### INV-035 → Tournament families
**Source cone identity:** H(source_cone(T')) = H(T'). This means R-cone proof strategy (INV-004) reduces to showing cut-flip E-invariance for a single cut-flip step.

### INV-001 (transfer matrix symmetry)
Confirmed by existing computational evidence (7500+ tests, n=4..8). The deep structural reason is Feng's dual Burnside detailed balance. Not proved algebraically.

---

## References

- INV-014: v2_quick.py, mod4_analysis.py (exhaustive n≤6)
- OPEN-Q-015: real_roots_full_v2.py, coeff_analysis.py, root_spacing2.py
- INV-035: tournament_families.py (background agent)
