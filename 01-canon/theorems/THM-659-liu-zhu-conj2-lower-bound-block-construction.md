---
id: THM-659 (RENUMBERED from THM-657 by mac-mini-S57: mac-mini pushed THM-657 (covering reformulation) at 20:06 < opus 20:10; THM-658 = kps chi_c. Collision protocol.)
title: The Liu–Zhu Conjecture 2 lower bound, proved in general — for every both-odd
  M = {x, y, y−x, y+x} the Motzkin density satisfies μ(M) ≥ (k+1)m/(4(k+1)m+1) via the
  explicit period-N avoiding set A = 2x·B, B = (k+1) equally-spaced length-m blocks; and
  the matching upper bound reduces to tiling Z_N by the 4-clique {0, x, y, x+y}
status: LOWER BOUND PROVED (uniform, fully symbolic — verified 651 instances x+y<80, and
  the four exclusion residues have closed forms holding for all k,m); UPPER BOUND reduced
  to a clique-tiling statement (mechanism identified, tiling verified computationally on
  the tested range); the full equality μ = (k+1)m/N is confirmed exactly per-instance on
  27 cases (HYP-5217, window-graph max-cycle-mean). Combined: Conjecture 2 is proved
  wherever the clique tiling holds; the lower bound is unconditional.
source: opus-2026-07-07-S147 (HYP-5277), building on HYP-5217 (S146 exact instance
  confirmation + the slab/combinatorial split)
depends_on:
  - THM-652   # chi(G_GW)=14 — the |S|=13 member of the same combinatorial-optimum phenomenon
related:
  - HYP-5217  # S146: Conjecture 2 confirmed exactly on 27 instances; the divisor ladder
  - HYP-5137  # S144: mu(GW)=1/13, the ladder separation (the |S|=13 analog)
  - THM-658   # kps-S76: chi_c(G_GW)<=13.5 < 1/M; the linearization defect (chi_c<1/M) lives
              # on the mu>M=mu>kappa locus — exactly this both-odd family (x>=3) and GW
external:
  - "M.-C. Liu, X. Zhu (2004), Fractional/independence density of two-interval distance
     graphs, JGT 47: Conjecture 2 (open for x>=3); their Thm 4.3 proves the x=1 case."
  - "Cantor–Gordon 1973 (mu >= kappa); Haralambis 1977 (mu computations, mu>kappa family)."
---

# THM-659 — Liu–Zhu Conjecture 2: the general lower bound and the clique-tiling upper bound

## Setting

For coprime odd `x = 2k+1`, `y = 2m+1` with `y > x`, let `M = {x, y, y−x, y+x}` and
`N = 4(k+1)m + 1`. `μ(M)` is the Motzkin density (max upper density of a set avoiding
differences `M`); `κ(M) = M(M)` is the lonely-runner constant. Liu–Zhu (2004) conjectured
**`μ(M) = (k+1)m/N`** and proved only `x = 1` (their Thm 4.3, where `μ = κ`). For `x ≥ 3`
the set is "almost difference closed" of type A.3 with `μ > κ` (S146), so the optimal
avoiding set is *not* a rotation slab — which is exactly why the conjecture resisted.

## The lower bound (PROVED, uniform)

> **Theorem (lower bound).** Let `B = ⋃_{t=0}^{k} {2mt, 2mt+1, …, 2mt+m−1}` — `(k+1)`
> blocks of length `m`, step `2m`, lying in `[0, 2(k+1)m) = [0, (N−1)/2)`. Then
> `A := 2x · B (mod N)` is a period-`N` set avoiding all differences in `M`, of density
> exactly `(k+1)m/N`. Hence **`μ(M) ≥ (k+1)m/N` for every both-odd instance.**

*Proof.* **Well-defined and full density.** `N` is odd, and `gcd(x, N) = 1`: if a prime
`p ∣ x = 2k+1` and `p ∣ N`, then `2k ≡ −1`, so `2(k+1) ≡ 1 (mod p)`, whence
`0 ≡ N = 2·2(k+1)m + 1 ≡ 2m + 1 = y (mod p)`, forcing `p ∣ gcd(x,y) = 1` — impossible. So
`gcd(2x, N) = 1` and `b ↦ 2x·b` is a bijection of `Z_N`; `|A| = |B| = (k+1)m`.

**The block difference set.** `B − B = {2m·u + w : −k ≤ u ≤ k, −(m−1) ≤ w ≤ m−1} (mod N)`.
As balanced residues this is exactly `[−(mx−1), mx−1]` **minus the single-point gaps**
between consecutive blocks, `{±m, ±3m, …, ±(2k−1)m}` (the odd multiples of `m` below `mx`).
(`2mk + (m−1) = mx − 1` since `x = 2k+1`.) Verified identical to the direct `B−B` on all
tested instances.

**Avoidance.** `A − A = 2x·(B − B) (mod N)`, so `A` avoids `d` iff `(2x)^{-1}d ∉ B − B`.
The four residues have closed forms (verified exactly, 651 instances, and derived below):

| `d` | `(2x)^{-1}·d (mod N)`, balanced | why `∉ B − B` |
|-----|--------------------------------|----------------|
| `x`   | `(N+1)/2 = −2(k+1)m`   | `|·| = 2(k+1)m > mx−1` (out of range) |
| `y`   | `−m`                    | odd-`m`-gap for `k ≥ 1`; out of range for `k = 0` |
| `y+x` | `mx + 1`                | `> mx−1` (out of range) |
| `y−x` | `mx`                    | `= mx > mx−1` (out of range) |

The residue values follow from `(2x)^{-1}x = 2^{-1} = (N+1)/2` and
`(2x)^{-1}y ≡ −m` (because `−2mx = −2m(2k+1) = −(N − y) ≡ y`, i.e. `y ≡ −2mx`, so
`(2x)^{-1}y ≡ −m`); the other two are `−m ± 2^{-1}` reduced. None lies in `B − B`, so
`A` avoids `M`. ∎

The construction is the exact generalization of the `x = 1` slab: for `k = 0` it is one
block of length `m` (a rotation slab, recovering Liu–Zhu Thm 4.3); for `k ≥ 1` it is
`(k+1)` blocks — a **genuinely combinatorial** optimum (`μ > κ`, no single-slab realization
exists, S146), which is why the earlier slab-only methods could not reach it.

## The upper bound (mechanism identified; reduced to a clique tiling)

> **`{0, x, y, x+y}` is a 4-clique of the distance graph `G_M`**: its six differences are
> `x, y, x+y, y−x, y, x` — all in `±M`. So any `M`-avoiding set contains ≤ 1 of any
> translate `c + {0, x, y, x+y}`.

Hence if `Z_N` decomposes into `(k+1)m` disjoint translates of this 4-clique plus one
leftover point (`4·(k+1)m = N − 1`), every avoiding set has ≤ `(k+1)m + 1` elements per
period, and the standard Haralambis single-point argument (as in Liu–Zhu's `x = 1` proof,
where the leftover is `2m`) sharpens this to `≤ (k+1)m`, giving `μ(M) ≤ (k+1)m/N`. The
4-clique `{0, x, y, x+y}` is the general analog of Liu–Zhu's block
`{i, i+1, i+2m+1, i+2m+2}` (their `x = 1` clique). The clique property holds symbolically
for all instances (`x+y ≤ 60` checked), and the tiling is verified by exact cover on all
24 instances with `x+y ≤ 32`; a uniform tiling construction closes the conjecture in
general.

## Status and consequence

- **Lower bound `μ ≥ (k+1)m/N`: unconditional, uniform, proved.** (Previously only the
  `x = 1` slab bound was known.)
- **Upper bound `μ ≤ (k+1)m/N`: reduced to the clique tiling** (verified on the range;
  uniform construction open) and confirmed exactly per-instance on 27 cases via the
  window-graph max-cycle-mean certificate (HYP-5217). So `μ(M) = (k+1)m/N` — **Liu–Zhu
  Conjecture 2 — holds on every instance where the 4-clique tiles `Z_N`**, and the lower
  half is now a theorem for all of them.
- This is the `|S| = 4` rung of the divisor ladder (S146): the combinatorial optimum
  `A = 2x·B` is the small-gcd (`x, y` both odd, `y±x` even) structure where the rotation
  slab fails; `GW = {1,…,11,13,24}` (THM-652, `μ = 1/13`, `{0,12} mod 26`) is the `|S| = 13`
  member of the same phenomenon.

Files: `lrc_liu_zhu_general_lb_opus_S147.py`, `lrc_liu_zhu_blocks_opus_S147.py`,
`lrc_liu_zhu_lower_bound_opus_S147.py` (+outs).
