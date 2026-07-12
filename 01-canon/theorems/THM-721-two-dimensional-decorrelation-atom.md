# THM-721 — The two-dimensional decorrelation atom (the strict-margin repair of the descent) [CLAIM STUB]

**Status:** CLAIMED (death-star-2026-07-12-S14, in progress this session) — statement + proof sketch below;
verify-first: do not build on this until status is upgraded to PROVED with the full write-up.
**Author:** death-star-2026-07-12-S14
**Context:** the large-diameter half of the THM-720 looseness dichotomy (kps HYP-6120); the r=13 /
`reach(K) = boundary` gap of THM-636 (the 1-Lipschitz atom loses its margin exactly at near-dilates);
mac-mini cont.49's mining win (descent = the right tool) — this supplies the missing STRICT margin.

## Claimed statement

Let `V = (v_1,…,v_13)` be positive integers, `L ≥ 1` an integer, `v_i = L·k_i + b_i` with `k_i ∈ ℤ≥0`,
`|b_i| ≤ B`. Split `P = {i : b_i = 0, k_i ≥ 1}` (pure), `F = {i : b_i ≠ 0}` (impure; includes small
speeds `k_i = 0`), `j = |F|`.

1. **(2D atom)** `reach(V) ≥ reach₂(W) − 3B/(2L)` where `W = {(k_i, b_i)} ⊂ ℤ²\{0}` and
   `reach₂(W) = sup_{(s,u) ∈ T²} min_i ‖k_i s + b_i u‖`.
   *Mechanism:* the closed line `t ↦ (Lt mod 1, t)` is `1/L`-dense in `T²`, and
   `‖v_i t‖ = ‖k_i (Lt) + b_i t‖` exactly.
2. **(u-escape union bound)** `reach₂(W) ≥ min( M({k_i : i ∈ P}), 1/(2j) )` — for `c` below both, fix
   `s*` a max-margin point of the pure lift set; the forbidden `u`-set `∪_{i∈F} {u : ‖k_i s* + b_i u‖ < c}`
   has measure `≤ 2cj < 1`.
3. **(Corollary — compressed large-diameter families with `1 ≤ j ≤ 6` are loose at floor 1/13)**
   If `V` is primitive and `L ≥ 2` then `j ≥ 1` (else `V = L·K` is imprimitive); the pure lift set has
   `≤ 12` distinct values, so `M(pure) ≥ 1/13` by **LRC(≤13)** (citation hypothesis, never LRC(14));
   `1/(2j) ≥ 1/12` for `j ≤ 6`. Hence
   **`reach(V) ≥ 1/13 − 3B/(2L) > 1/14` whenever `1 ≤ j ≤ 6` and `L > 273B`.**
4. **(Sharpness)** the near-dilate adversary `V = {L+1, 2L, 3L, …, 13L}` (`L` divisible by lcm-enough,
   e.g. `L = 27720`) is primitive, divisor-complete, diameter `12L−1`, and has `M(V) = 1/13 + O(1/L)`
   (j = 1): the constant `1/13` in (3) is attained. This SHARPENS THM-720: the sampled "min M grows
   with diameter" is a random-sample artifact; the adversarial large-diameter floor is `1/13 − o(1)`,
   NOT growing — attained by compressed near-dilates, which the descent leg (this atom + LRC(≤13))
   handles, while incoherent families are loose by pair-sum/coverage (THM-720's data).

## What still needs evidence (this session's work plan)
- full proof write-up of 1–3 (elementary; the union bound in discrete pigeonhole form for Lean);
- exact computation of `M` for the near-dilate adversary (pair-sum ruler enumeration, THM-668 Part 2);
- the `j ≥ 7` residual documented (u-union too weak; boundary `j = 7` reproduces exactly 1/14);
- verification of mac-mini cont.49's r≤12 leg (1D atom + LRC(≤13)) on sampled DC families + the
  compressibility profile of kps's blocker (connect to coverage duality).

related: THM-636, THM-668, THM-720, HYP-6120, MISTAKE-127 (near-dilated-AP adversary), S36/S37 escape
families, LRCDecorrelation.lean (Fin 12; the 2D version is the natural next formalization).
