---
id: THM-399
name: allzero-staircase-c3-formula
status: PROVED
proved_by: monad-researcher-2026-06-03-S578
depends_on: [INV-190, HYP-2095]
---

# THM-399: 3-Cycle Count in All-0 Interleaved Staircase

**Statement:** The all-0 interleaved staircase tournament at n=2k has exactly k(k−1) directed 3-cycles.

## Definition

The all-0 interleaved staircase at n=2k has vertex set {0,1,...,2k−1} with:
- **Pairs:** (0,1), (2,3), ..., (2k−2, 2k−1)
- **Global ranking:** rank[2p] = p (dominants d_p = 2p), rank[2p+1] = k+p (recessives r_p = 2p+1)
- **Within pair p:** r_p → d_p (recessive beats dominant)
- **Between pairs (i≠j in different pairs):** i → j iff rank[i] < rank[j]

**Arc structure summary:**
- Dominants form a transitive tournament: d_p → d_q iff p < q
- Recessives form a transitive tournament: r_p → r_q iff p < q
- All dominants beat all recessives from different pairs: d_p → r_q for all p≠q (since rank[d_p]=p < k ≤ k+q=rank[r_q])
- Within pair: r_p → d_p

## Proof

We enumerate all directed 3-cycles by case analysis on the mix of dominants and recessives.

**Case 1 (3 dominants):** The induced subgraph is the transitive tournament {d_p, d_q, d_r} (p<q<r). Transitive tournaments have no directed cycles. **Count: 0.**

**Case 2 (3 recessives):** Same argument. **Count: 0.**

**Case 3 (2 dominants, 1 recessive — set {d_p, d_q, r_r} with p<q):**

The arcs are:
- d_p → d_q (since p<q)
- d_p → r_r and d_q → r_r (since all dominants beat all recessives from different pairs, and p≤k−1<k≤k+r)
- Exception: if r=p, then r_p → d_p (within-pair), so d_p does NOT beat r_p
- Exception: if r=q, then r_q → d_q

Sub-case r≠p and r≠q: All three arcs go from {d_p,d_q} to r_r, plus d_p→d_q. The subgraph is transitive. **No 3-cycle.**

Sub-case r=p: Arcs are d_p→d_q, d_q→r_p (since q≠p), r_p→d_p. This is the 3-cycle d_p→d_q→r_p→d_p. **One 3-cycle.**

Sub-case r=q: Arcs are d_p→d_q, r_q→d_q (same pair q), d_p→r_q (since p≠q). The subgraph is transitive (d_p beats both d_q and r_q; r_q beats d_q). **No 3-cycle.**

**Contribution from Case 3:** exactly one 3-cycle for each pair (p,q) with 0≤p<q≤k−1 and r=p. This gives C(k,2) = k(k−1)/2 triples, each yielding 1 cycle: **C(k,2) directed 3-cycles.**

**Case 4 (1 dominant, 2 recessives — set {d_p, r_q, r_r} with q<r):**

The arcs are:
- r_q → r_r (since q<r)
- d_p → r_q and d_p → r_r if p≠q and p≠r (since d_p beats all recessives from different pairs)
- Exception: if p=q, then r_q → d_p (within-pair)
- Exception: if p=r, then r_r → d_p

Sub-case p≠q and p≠r: d_p beats r_q and r_r, plus r_q→r_r. Transitive. **No 3-cycle.**

Sub-case p=q (d_p, r_p with p<r): Arcs are r_p→r_r (p<r), d_p→r_r (p≠r), r_p→d_p (same pair). Tournament is transitive (r_p beats both d_p and r_r; d_p beats r_r). **No 3-cycle.**

Sub-case p=r (d_p, r_q, r_p with q<p): Arcs are r_q→r_p (q<p), d_p→r_q (p≠q), r_p→d_p (same pair). This is the 3-cycle r_q→r_p→d_p→r_q. **One 3-cycle.**

**Contribution from Case 4:** exactly one 3-cycle for each triple (d_p, r_q, r_p) with q<p, i.e., for each pair (q,p) with 0≤q<p≤k−1. This gives C(k,2) triples, each yielding 1 cycle: **C(k,2) directed 3-cycles.**

**Total:** c₃ = C(k,2) + C(k,2) = 2·C(k,2) = k(k−1). □

## Verification

Verified computationally for k=2..12 (n=4..24) using direct enumeration. See HYP-2095.

| k | Expected k(k−1) | Computed c₃ |
|---|-----------------|-------------|
| 2 | 2 | 2 ✓ |
| 3 | 6 | 6 ✓ |
| 4 | 12 | 12 ✓ |
| ... | ... | ... ✓ |
| 12 | 132 | 132 ✓ |

## Corollary

Every pair of the form (d_p, d_q, r_p) with p<q and (r_q, r_p, d_p) with q<p yields exactly one directed 3-cycle. There are no other 3-cycles.

The two families of 3-cycles are in bijection via the map that swaps "same-pair role": family 1 uses the smaller dominant's pair recessive; family 2 uses the larger recessive's pair dominant.
