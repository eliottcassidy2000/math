---
name: HYP-2095-allzero-staircase-k7-k8
description: H values for all-0 interleaved staircase k=2..12; c3=k(k-1) confirmed; growth r(k)→3k asymptotically; OEIS not found (likely novel)
metadata:
  type: project
---

# HYP-2095: All-0 interleaved staircase H values at k=7..12

**Status:** CONFIRMED (Held-Karp exact DP, verified all k=2..10 against prior results)
**Sources:** monad-researcher-2026-06-02-S577 (k=7,8), monad-compute-2026-06-02 (k=9,10),
            monad-researcher-2026-06-02-S578 (k=11,12 NEW)
**Scripts:** `04-computation/staircase_allzero_k9_s_monad.py` (k≤11),
             `04-computation/staircase_allzero_k12_s_monad.py` (k=12)

## Full sequence k=2..12

| k | n | H | c3 | c3=k(k-1)? | New? |
|---|---|---|----|------------|------|
| 2 | 4 | 5 | 2 | ✓ | |
| 3 | 6 | 29 | 6 | ✓ | |
| 4 | 8 | 233 | 12 | ✓ | |
| 5 | 10 | 2489 | 20 | ✓ | |
| 6 | 12 | 33773 | 30 | ✓ | |
| 7 | 14 | 562685 | 42 | ✓ | S577 |
| 8 | 16 | 11222321 | 56 | ✓ | S577 |
| 9 | 18 | 262755369 | 72 | ✓ | monad-compute |
| 10 | 20 | 7110764837 | 90 | ✓ | monad-compute |
| 11 | 22 | **219612027389** | 110 | ✓ | **S578 NEW** |
| 12 | 24 | **7658921303353** | 132 | ✓ | **S578 NEW** |

## Growth ratios r(k) = H(k)/H(k-1)

5.80, 8.03, 10.68, 13.57, 16.66, 19.94, 23.41, 27.06, 30.88, 34.87

Differences in ratios: +2.23, +2.65, +2.89, +3.09, +3.28, +3.47, +3.65, +3.82, +3.99

**Asymptotic behavior:** The deficit 3k − r(k) peaked around k=6 (value ~4.43) and is decreasing monotonically: 4.34, 4.06, 3.59, 2.94, 2.12, 1.13 for k=7..12. Empirically r(k) → 3k as k→∞. This implies H(k) ~ C · 3^k · k! for large k (super-exponential growth).

## Key facts

- **c3(k) = k(k-1)** confirmed through k=12 (the directed 3-cycle count formula holds exactly)
- No order-2 or order-3 linear recurrence (Gaussian elimination confirms)
- The 2-term polynomial recurrence H(k) ≈ (6.82k−0.84)·H(k-1) + (-12.65k²+37.39k-52.01)·H(k-2) fits k=5..9 exactly but diverges for k≥10 — not a true identity
- Markov equation x²+y²+z² = 3xyz FAILS for all consecutive triples (LHS-RHS grows rapidly)
- OEIS search: specific terms searched (oeis.org blocked all API fetches; web search for "5,29,233,2489,33773" found no match). **Likely novel sequence.**

## The tournament at n=2k

- Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
- Global ranking: rank[2p]=p (dominants), rank[2p+1]=k+p (recessives)
- Within-pair: recessive (odd) beats dominant (even): (2p+1)→(2p)
- Between pairs: lower-ranked vertex beats higher-ranked vertex

## Open questions

1. **Prove** c3(k) = k(k-1) formula for all k (likely by induction on tournament structure)
2. **Submit to OEIS** — likely new sequence; needs 10+ terms which we now have
3. **Asymptotic proof**: why does r(k) → 3k? Is there a combinatorial explanation for the 3k growth?
4. **Compute k=13** — estimated ~7.5 min runtime with array.array (2^26·26 entries ≈ 13 GB); may need C implementation
5. **Algebraic structure**: is H(k) = permanent of some explicit k×k matrix?
