---
id: HYP-2220
title: The constructive route to LRC(14) — quantitative C′ on the 2n-1 shell (corrects S622)
status: OPEN (target); quantitative gap value VERIFIED n=4..8 (THM-415)
source: claudebox-2026-06-03-S623
depends_on:
  - THM-415   # min M = 2/(2n-1) for multiple-of-n configs
  - THM-398   # C′ ⟹ LRC(n)
related:
  - THM-412   # twisted-shell dodge (the S622 attempt, here corrected)
  - THM-401   # 2n-1 shell
  - THM-407   # opus-S599g shell-folding (13 shells → 3 strata mod 27)
---

# HYP-2220 — the corrected constructive route to LRC(14)

## What the small sizes taught us

Studying C′(n) exhaustively at small `n` (the lab being **n=5**, `2n-1 = 9 = 3²`, the baby of
n=14's `27 = 3³`) corrected the S622 picture and sharpened the target:

1. **The gap is `2/(2n-1)`, not just `> 1/n`** (THM-415): every multiple-of-n config has
   `M ≥ 2/(2n-1)`, the `s=2` Kravitz rung, with equality realised on the `2n-1` shell at
   band-distance 2.
2. **S622's `m ≤ 2n-1` dodge ceiling was a sampling artifact.** Exhaustively there is a real
   residual at every `n` (e.g. 225 configs at n=5) of loose configs whose optimum lives on shells
   `> 2n-1` — but all with `M > 2/(2n-1)` (deeper, narrower spikes). The `2n-1` shell is the
   **extremal** shell (where the minimum `2/(2n-1)` is attained), not a universal escape ceiling.
3. So the genuine difficulty of C′(n) is concentrated at the extremal: the multiple-of-n configs
   that are *tight on the `2n-1` shell* (`M = 2/(2n-1)`). Everything else is strictly looser.

## The route (target)

> **Prove `M(S) ≥ 2/(2n-1)` for every primitive multiple-of-n config `S`** ⟹ quantitative C′(n)
> ⟹ (via THM-398) **LRC(n)**. For `n = 14`: `multiple-of-14 ⟹ M ≥ 2/27`.

This is Kravitz's second-gap restricted to the multiple-of-n class. The constructive content is to
exhibit, for each such `S`, a witness `t = a/m` with `min_i ‖v_i t‖ ≥ 2/(2n-1)` — and the extremal
analysis says the binding witnesses live on `m = 2n-1`.

## The n=5 → n=14 lift (the 3-adic ladder)

- **n=5:** `2n-1 = 9 = 3²`. Shell tower `1 | 3 | 9`. Extremal `M = 2/9`. The twisted involution
  `±` on `(ℤ/9)*` (3 unit-pairs `{1,8},{2,7},{4,5}`); the inner shell `3·(ℤ/3)` is the "cross".
- **n=14:** `2n-1 = 27 = 3³`. Shell tower `1 | 3 | 9 | 27`. Extremal `M = 2/27`. Unit-pairs in
  `(ℤ/27)*` (9 of them); inner shells `3·(ℤ/9)`, `9·(ℤ/3)`, apex `0`.

Same prime-3 ramification, one level deeper — which is precisely why n=14 is the first even
frontier (cf. [[THM-407]]: doubling stops being shell-transitive at the odd-prime-power `2n-1`).
**Proving the n=5 extremal bound `M ≥ 2/9` constructively, on the `3²` tower, should template the
n=14 bound `M ≥ 2/27` on the `3³` tower.**

## To do
1. Prove `M ≥ 2/9` for every multiple-of-5 config (smallest prime-power shell; only 4 speeds —
   a finite case analysis on residues mod 9 should close it).
2. Lift the n=5 argument up the 3-adic tower to `M ≥ 2/27` at n=14.
3. Identify the extremal multiple-of-n configs (tight on the `2n-1` shell) for general n; is the
   minimiser set finite up to symmetry (cf. THM-411 tight-family finiteness)?
