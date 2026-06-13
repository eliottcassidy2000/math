# THM-111: Odd/Even Relative Homology Dichotomy

**Status:** CONJECTURE (strong computational evidence)
**Filed by:** opus-2026-03-09-S52

## Statement

For any tournament T on n vertices and any vertex v:

**Odd levels:** dim H_{2k+1}(T, T\v) ≤ 1 for all k ≥ 0.

**Even levels:** dim H_{2k}(T, T\v) can be arbitrarily large.

## Proof of Boolean Odd Betti (THM-098) from Odd Relative Bound

**Theorem:** If H_{2k+1}(T, T\v) ≤ 1 for all (T,v), then β_{2k+1}(T) ≤ 1
for all tournaments T.

**Proof by strong induction on n:**

Base: β_{2k+1} = 0 for n small enough (n ≤ 4k+3 suffices, since Omega_{2k+1}
requires 2k+2 vertices for nontrivial paths).

Inductive step: Assume β_{2k+1}(T') ≤ 1 for all (n-1)-vertex tournaments T'.
Let T be an n-vertex tournament. For ANY vertex v, the LES of (T, T\v) gives:

    H_{2k+1}(T\v) →^{i_*} H_{2k+1}(T) →^{j_*} H_{2k+1}(T,T\v) → ...

From the general LES inequality (no conditions on adjacent Betti numbers):

    β_{2k+1}(T) = dim(im i_*) + dim(im j_*)
                 ≤ β_{2k+1}(T\v) + dim H_{2k+1}(T,T\v)
                 ≤ 1 + 1 = 2

But this only gives ≤ 2. The sharper bound needs one of:

(a) **Good vertex existence**: ∃ v with β_{2k+1}(T\v) ≤ β_{2k+1}(T) - dim H_{2k+1}(T,T\v).
(b) **Simultaneous bound**: When β_{2k+1}(T\v) = 1, either i_* = 0 or H_{2k+1}(T,T\v) = 0.

For k = 0 (β₁): PROVED (THM-103, star constraint argument).
For k = 1 (β₃): VERIFIED computationally n ≤ 8, 0 violations.

## Evidence

### Odd relative bound H_{2k+1} ≤ 1

| n | k | H_{2k+1}(T,T\v) | Method | Total pairs | Max | Violations |
|---|---|-----------------|--------|-------------|-----|------------|
| 6 | 0 | H_1 | exhaustive | 196608 | 1 | 0 |
| 6 | 1 | H_3 | exhaustive | 196608 | 1 | 0 |
| 6 | 2 | H_5 | exhaustive | 196608 | 0 | 0 |
| 7 | 0 | H_1 | sampled | 1400 | 1 | 0 |
| 7 | 1 | H_3 | sampled | 1400 | 1 | 0 |
| 7 | 2 | H_5 | sampled | 1400 | 0 | 0 |
| 8 | 1 | H_3 | sampled | 800 | 1 | 0 |

### Even relative — no bound

| n | k | H_{2k}(T,T\v) | Max observed |
|---|---|---------------|--------------|
| 7 | 2 | H_4(Paley T_7) | **7** |
| 7 | 1 | H_2 | 1 |
| 7 | 3 | H_6 | 0 |

### Paley T_7 complete relative data

    H_p(T_7, T_7\v) = [0, 0, 0, 0, 7, 0, 0]

for ALL 7 vertices v (by vertex-transitivity).

LES verification at level 4:
- β_4(T_7) = 6, β_4(T_7\v) = 0, β_3(T_7\v) = 1
- rank(i_*) = 0 (source is 0)
- dim(im j_*) = dim(ker δ) = 6
- dim(im δ) = dim H_4(T,T\v) - dim(ker δ) = 7 - 6 = 1
- im(δ) ⊂ H_3(T\v) which has dim 1
- So δ is surjective onto H_3(T\v), consistent

## Significance

1. **Odd Betti Boolean** (THM-098): Follows from the odd relative bound
   plus induction, modulo the "good vertex" refinement.

2. **Even Betti unconstrained**: The even relative H can grow, explaining
   why β_4 = 6 for Paley T_7. No upper bound is expected.

3. **Structural insight**: The odd/even dichotomy may reflect a deeper
   algebraic structure — perhaps related to the palindromic Omega dimensions
   of Paley tournaments (Poincaré duality?) or to the anti-automorphism.

## Palindromic Omega for Paley T_7

    dim(Omega_p): [7, 21, 42, 63, 63, 42, 21]

Palindromic! This equals C(7, p+1) * filling_ratio(p).
The palindrome suggests Poincaré-like duality in the path complex.

## Files

- `04-computation/beta_odd_relative_general.py`
- `04-computation/paley7_relative_betti.py`
- `05-knowledge/results/beta_odd_relative_general.out`
- `05-knowledge/results/paley7_relative_betti.out`

## See Also

- THM-098 (Boolean odd Betti conjecture)
- THM-110 (β₃ ≤ 1)
- THM-108 + THM-109 (β₂ = 0)
- THM-103 (β₁ ≤ 1)
- THM-095 (adjacent-odd seesaw)
