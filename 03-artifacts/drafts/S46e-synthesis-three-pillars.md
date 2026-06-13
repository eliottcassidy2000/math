# Three Pillars of Tournament Structure: OCF, Cumulants, and Path Homology

**Author:** opus-2026-03-07-S46e
**Date:** 2026-03-07

## Overview

This document synthesizes three apparently independent structural frameworks for tournaments, revealing deep connections between them.

### Pillar 1: The Odd-Cycle Collection Formula (OCF)
H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

- Proved by Grinberg-Stanley (arXiv:2412.10572)
- Counts vertex-disjoint odd-cycle collections
- Factor of 2 per odd cycle from T/T^op duality (Grinberg-Stanley Theorem 1.38)
- Encodes ALL Hamiltonian path information through the independence polynomial

### Pillar 2: The Cumulant-Cycle Hierarchy
kappa_{2k}(T) = Bernoulli_constant + (2/C(n,2k)) * t_{2k+1} + nonlinear lower terms

**Universal Coefficient Theorem (THM-117, PROVED):**
coeff(t_{2k+1} in kappa_{2k}) = 2/C(n, 2k)

The even cumulants of the forward-edge distribution form a graded hierarchy:
- kappa_2 (variance) probes 3-cycles via t_3
- kappa_4 probes 5-cycles via t_5 and alpha_2
- kappa_6 probes 7-cycles via t_7
- kappa_{2k} probes (2k+1)-cycles via t_{2k+1}

Each level introduces exactly one new family of odd cycles.

**Connection to OCF:** The factor 2 in 2/C(n,2k) is the SAME factor 2 from OCF (each directed odd cycle contributes 2 to H via I(Omega,2)). The 1/C(n,2k) comes from counting (2k+1)-vertex clusters in the moment expansion.

**Connection to Grinberg-Stanley:** The structural chain:
  U_T in N[p_1, 2p_3, 2p_5, ...]  -->  factor 2 per odd cycle
  --> 2*t_{2k+1} in forward path count  -->  2/C(n,2k) in cumulant

### Pillar 3: GLMY Path Homology
For a tournament T on n vertices, the GLMY path Betti numbers satisfy:
- beta_0 = 1 (connected)
- beta_{2k} = 0 for all k >= 1 (CONJECTURE, verified n <= 7)
- beta_1, beta_3 are mutually exclusive
- Euler characteristic chi in {0, 1}

**Even Betti Vanishing Conjecture (OPEN-Q-024):**
Tournaments have "odd-dimensional directed topology only."
Only beta_{2k-1} can be nonzero, mirroring the odd-cycle-only structure.

## The Structural Parallel

All three pillars express the same fundamental fact: tournaments have ODD-ONLY structure.

| Aspect | OCF | Cumulants | Path Homology |
|--------|-----|-----------|---------------|
| Fundamental object | I(Omega, 2) | kappa_{2k} | H_p^{path}(T) |
| Odd structure | Only odd cycles contribute | Only even cumulants nonzero | Only odd Betti numbers nonzero |
| Grading | alpha_k (size-k independent sets) | kappa_{2k} (2k-th cumulant) | beta_{2k-1} (odd Betti) |
| Factor of 2 | 2^{psi(sigma)} in Grinberg-Stanley | 2/C(n,2k) coefficient | Unknown |
| Alternating sum | H = 1 + 2*alpha_1 + 4*alpha_2 + ... | (moment-cumulant formula) | chi = 1 - beta_1 + beta_3 - ... |

## Key Results of This Session

1. **THM-117 (PROVED):** Universal coefficient 2/C(n,2k) in the cumulant hierarchy
2. **OPEN-Q-024 (NEW):** Even Betti vanishing conjecture for tournaments
3. **Mutual exclusivity:** beta_1 and beta_3 never both positive
4. **Grinberg-Stanley connection:** Factor 2 traced through all three pillars
5. **Celano-Sieger-Spiro clarification:** A_T(t) != F(T,x) (different statistics)
6. **Mod 3 amplification:** THM-086 implies n!*E[fwd^r] = 0 mod 3 for ALL r

## Open Questions

1. Prove even Betti vanishing (OPEN-Q-024)
2. Find formula for beta_1 in terms of cycle invariants
3. Is there a direct path-homology interpretation of I(Omega, 2)?
4. Does the mutual exclusivity of beta_1/beta_3 extend to beta_3/beta_5?
5. Can the cumulant generating function be expressed via the Ihara zeta?
6. Does the "topological phase transition" (beta_1 -> 0 -> beta_3) correspond to a threshold in cycle density?
