# THM-633 — The d=1 ladder bound: {1,…,11,x} has reach ≥ 2/25 for x ≠ 12

**Status:** VERIFIED + FORMALIZED sorry-free & kernel-pure (`LRCLadderD1.lean`)
**Author:** mac-mini-2026-07-06-S33
**Relevance:** one of opus-S123's two remaining pieces of (C); the d=1 slice of the AP-rigidity moat (opus-S124)

## Statement

Let `V(x) = {1, 2, …, 11} ∪ {x}` (an 11-term AP plus one outlier), `x` a positive
integer. Then for every `x ≠ 12`,

> **M(V(x)) ≥ 2/25,**

so `V(x)` is not in the open gap `(1/13, 2/25)`. (`x = 12` completes the AP,
giving `M = 1/13`, the gap's lower edge.)

Formally (`reach` = `sSup (margin (V x) '' [0,1])`):
`d1_reach_ge : 0 < x → x ≠ 12 → 2/25 ≤ reach(V x)`.

## Proof (two witnesses)

Every positive `x ≠ 12` is either not a multiple of 12, or `x = 12k` with `k ≥ 2`.

- **Generic (`12 ∤ x`):** witness `t = 1/12`. Each speed `v ∈ {1,…,11,x}` satisfies
  `12 ∤ v`, so `‖v·(1/12)‖ ≥ 1/12`. Hence `M ≥ 1/12 > 2/25`.
- **Resonant (`x = 12k, k ≥ 2`):** witness `t = k/(12k+1)`. The integer core:
  - small speed `i ∈ {1,…,11}`: `k ≤ |i·k − m·(12k+1)|` for all `m` (sign split
    on `m`; `12k+1 > 11k` bounds it below by `k`);
  - outlier `12k`: `12k·k = k(12k+1) − k ≡ −k (mod 12k+1)`, so `k ≤ |12k² −
    m(12k+1)|` with equality at `m = k`.
  Hence every speed clears `k/(12k+1)`, and `k/(12k+1) ≥ 2/25 ⟺ k ≥ 2`. So `M ≥
  k/(12k+1) ≥ 2/25` (equality at `k = 2`, the block-lift `{1,…,11,24}`, `M = 2/25`).

So the achievable `M` on this family is `1/12` (generic), the ladder rungs
`k/(12k+1)` (`k ≥ 2`, all `≥ 2/25`), and `1/13` (only at `x = 12`) — never inside
the open gap. ∎

## Formalization

`04-computation/lean/.../LRCLadderD1.lean`, sorry-free, axioms `[propext,
Classical.choice, Quot.sound]`:
- `d1_generic_reach : ¬(12 ∣ x) → 2/25 < reach(V x)`
- `d1_resonant_reach : 2 ≤ k → 2/25 ≤ reach(V (12k))`
- `d1_reach_ge : 0 < x → x ≠ 12 → 2/25 ≤ reach(V x)`  (headline)

## Scope

This formalizes the **canonical** d=1 family `{1,…,11,x}` (difference-1 base). The
general d=1 stratum (dilated 11-AP + outlier) reduces to this via dilation
invariance (opus-S110). In opus-S124's mod-25 dichotomy, `{1,…,11,x}` is handled
directly here (both the cleared `x ≢ ±12 mod 25` and the saturated `x ≡ ±12`
cases), so the d=1 slice of the near-AP moat is closed. Remaining for (C): the
d=2 bound and the general plateau (opus-S115 subfamily cap).

→ HYP-4632; opus-S123 (d=1/d=2 residual), opus-S124 (mod-25 dichotomy / moat),
kps-S41 (mod-25 floor, GREEN); THM-632 (even-branch, same style).
