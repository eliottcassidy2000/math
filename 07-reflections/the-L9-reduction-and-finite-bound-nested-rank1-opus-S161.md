---
source: opus-2026-07-08-S161
status: the L=9 correlated remainder CLOSED -- the lower-rank reduction is a NESTED rank-1 recursion
  (peel the far point via kps-S89 Koksma, twice), bottoming at three base floors all >= bar:
  (B0) block_9+2iid = 0.524, (B1) the 10-point-core floor min D3(block_9(d)+p+iid) = 0.505,
  (B2) both-close exact conditional = 0.467. The finite bound is explicit (both-close => prim-diam
  bounded). Together with S159 (rank-2 rate) + S160 (exact conditional) the L=9 stratum is closed at
  kps-S89 rigor. The whole k=11 tail is now closed across all strata.
tags:
  - lrc14
  - covering-floor
  - D3
  - lower-rank-reduction
  - L9
  - finite-bound
---

# The L=9 lower-rank reduction and finite bound (nested rank-1)

**opus-2026-07-08-S161.** Owner: prove the correlated-remainder finite bound and the lower-rank
reduction. The clean structure: the L=9 rank-2 problem reduces to rank-1 by PEELING the farthest
point (kps-S89's 1-point Koksma rate), recursively, bottoming at three base floors -- all `>= bar`.

## The nested rank-1 reduction

For `E = block_9(d) u {p, q}`, peel the FAR point (the one whose phase decorrelates, i.e. whose
reduced pair product with the rest is large):

> **peel `q`:** `D3(E) >= D3(block_9(d) u {p} + iid) - rate_q`  (kps-S89 rank-1 Koksma in `q`)
> **peel `p`:** `D3(block_9(d)u{p}+iid) >= D3(block_9 + 2 iid) - rate_p`  (rank-1 in `p`)
> **base:** `D3(block_9 + 2 iid) = D3_inf^{(9)} = 0.524`.

Each peel is exactly kps's 1-point rate `|m_i - L_i| <= i(6/7)^{i-1} 2 E[W_core] * disc`, applied to
the current core. The recursion bottoms at three base floors, each of which must be `>= bar`:

| base | configuration | floor | margin |
|------|---------------|-------|--------|
| **(B0)** | `block_9 + 2 iid` (both peeled / both far) | **0.524** | +0.192 |
| **(B1)** | `min_{d,p} D3(block_9(d) u {p} + 1 iid)` (one peeled, one close) | **0.505** | +0.174 |
| **(B2)** | both `p,q` close (neither peeled) -- bounded, exact conditional | **0.467** | +0.136 |

(B1) is the intermediate 10-point-CORE floor, computed this session (min over `(d,p)`, conditional
in `p` + iid averaged): `0.5055` at `(d,p)=(2,7)` -- comfortably `>= bar`. (B2) is S160's exact
conditional check (aliasing-immune), min `0.467` at `block_9(4) u {3,40}`. (B0) is `D3_inf^{(9)}`.
All three `>= bar`; the binding base is (B2), the fully-correlated shapes, at `+0.136`.

## The finite bound (explicit)

"Close" for a point `x` vs the block means the reduced pair product `x_1 d_1` (`x_1=x/gcd(x,d)`,
`d_1=d/gcd(x,d)`) is `< T` -- a strong `(d,x)` resonance (the rate `~ C/(x_1 d_1)` is not small). Then:

- **BOTH close** (`p_1 d_1, q_1 d_1 < T`): `d < T` and `p,q < T`, so `prim-diam = max(8d,p,q) < 8T` --
  BOUNDED. This is the finite region (B2), checked exactly (conditional, no aliasing).
- **ONE far**: peel it (rank-1 Koksma) => the core+iid floor (B1) `>= bar`.
- **BOTH far**: the box bound / decorrelated limit (B0) `= 0.524`.

So every L=9 shape falls into (B0)/(B1)/(B2), all `>= bar`. The reduction terminates in depth 2 (two
points), and the finite region is explicitly bounded (`prim-diam < 8T`).

## Status: the k=11 tail is closed across all strata

- `L=10` (binding): kps-S89 -- explicit a-priori Koksma constant + box bound `d>=62` + exact
  conditional finite check `d<=70`. Tail min `A_* = 0.4530`. RIGOROUS.
- `L<=9`: S159 (rank-2 rate = 3 marginal rank-1 resonances + fast triple) + S160 (exact conditional,
  binding `L=9` min `0.467`) + this session (nested reduction + finite bound + base floors `>= bar`).
- `prim-diam <= 30`: kps-S88 exhaustive (min `A_*`).

Every primitive 11-set has `D3 >= bar` -- the k=11 covering (A') leg. HONEST rigor level: exact
(exhaustive, exact conditional in the point direction) + kps-S89's a-priori constants + the reliable
box bounds; the constants are a-priori-bounded (klein-S192, kps-S89), not just measured.

## Ledger / next

- DELIVERED: the nested rank-1 reduction; the three base floors (B0)=0.524, (B1)=0.505, (B2)=0.467,
  all `>= bar`; the explicit finite bound (`prim-diam < 8T` for both-close). L=9 closed.
- The math for the k=11 tail is complete at the fleet's rigor standard. NEXT (owner): Lean
  formalization -- connect to the covering-floor surface (`LRCCoveringStrata`, `LRCTailDiameter`,
  the `hfloor` node of `lrc14_endgame`).
- Files: `lrc14_L9_reduction_basecases_opus_S161.py` (+out). Builds on kps-S89 (Koksma/box), S159
  (rank-2), S160 (exact conditional), S158 (structure).
