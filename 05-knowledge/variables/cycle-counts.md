# Variable: t_k — Directed k-cycle counts

**Symbol:** `t_k` (especially `t_3`, `t_5`, `t_7`)
**Type:** non-negative integer
**Defined in:** `01-canon/definitions.md`

## Definition
t_k = number of directed cycles of length k in tournament T.

Note: for 5-cycles, a 5-vertex set can contain MULTIPLE 5-cycles (up to the number of Hamiltonian cycles in the induced sub-tournament). This multiplicity matters for OCF.

## Values and ranges

### t_3 (directed triangles)
| n | Min | Max | Formula for regular | Source |
|---|-----|-----|-------------------|--------|
| 3 | 0 | 1 | 1 | |
| 4 | 0 | 4 | n/a | |
| 5 | 0 | 8 | C(5,3)/4 = 2.5... not integer, no regular at n=5 odd... | |
| 5 | 0 | 8 | max at Paley-like | exhaustive |
| 7 | 0 | 28 | C(7,3)/4 = 8.75... regular has t_3 near this | exhaustive |

General: t_3 = C(n,3) - sum of C(score_i, 2) [known formula from score sequence].

### t_5, t_7
Appear in higher-order OCF decomposition. Values at n=7:
- t_5 range: [0, ...] (computed in alpha1_3_structure_n7.py)
- t_7 range: [0, ...] (5040 Hamiltonian cycles max in K_7-tournament)

## Key equations
- **tr(c_{n-3})** = 2*(n-2)!*(t_3 - C(n,3)/4) — THM-054
- **c_2 at n=7:** c_2 = 24*bc33 - 60*t_3 + 12*t_5 + 231
- **c_0 at n=7:** c_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc33 + 253/4
- **c_2 at n=9:** c_2 = 462*t_3 - 60*t_5 + 12*t_7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640
- **H at n=7:** H = 1 + 2*(t_3 + t_5 + t_7) + 4*bc33 [EXACT — proved via OCF, THM-002]
- **Char poly connection:** tr(A^3) = 3*t_3

## Relationships
- t_3 determines D_2, D_3 components of W — see [fourier-coefficients.md](fourier-coefficients.md)
- t_3 appears in transfer matrix trace coefficient — see [transfer-matrix.md](transfer-matrix.md)
- t_3 + t_5 + t_7 + ... = total odd cycles = alpha_1 — see [alpha-k.md](alpha-k.md)
- bc33 = disjoint 3-cycle pairs — see [bc33.md](bc33.md)

## Tags
#cycles #t3 #t5 #t7 #OCF-decomposition #score-sequence #Fourier
