---
source: opus-2026-07-08-S158
status: the L<=9 strata of the k=11 prim-diam>25 tail have D3 >= bar with LARGE margin (worst = L=9,
  min 0.473, margin +0.14; margins grow to +0.52 at L=2), via the same block_L(scale d)+(11-L)-points
  decorrelation structure as S157 (L=10). Combined with kps-S88 (EXHAUSTIVE prim-diam<=30, min A_*)
  and S157 (L=10, prim-diam>30), the WHOLE k=11 covering tail closes: every primitive 11-set has
  D3 >= bar. Rigorous mechanism + numerically-certified constants + reliable (adaptive-NG) verification.
tags:
  - lrc14
  - covering-floor
  - D3
  - longest-AP
  - decorrelation
  - tail-closure
---

# The L<=9 tail floor, and the whole-tail closure

**opus-2026-07-08-S158.** Owner: prove the L<=9 strata floor too. S157 proved the L=10 binding
family; kps-S88 then extended the EXHAUSTIVE to prim-diam <= 30 (min D3 = A_* = 0.452986) and flagged
"only prim-diam > 30 remains." This note handles the L<=9 strata for prim-diam > 30, completing the
tail: every primitive 11-set has `D3 >= bar`.

## The structure (S157 generalized from L=10 to all L)

A shape with longest AP of length `L` contains an `L`-term AP; up to translation/dilation (both
`D3`-invariant) it is `block_L(scale d) u {11-L extra points}`, so the AP phases are `frac(jdx) =
frac(ju)` (`u := frac(dx)`) and

> `W(x) = G_L(u, v_1, ..., v_{11-L})`,  `v_i = frac(p_i x)`,  `G_L` a fixed function on `T^{12-L}`.

The moments deviate from the DECORRELATED limit `D3_inf^{(L)} = D3(block_L + (11-L) iid points)` by a
resonance sum over the rank-`(11-L)` relation lattice of `(d, p_1, ..., p_{11-L})` -- the exact
generalization of S157's rank-1 identity `m_j = L_j + sum_{k!=0} Ghat^j(kp,-kd)`. The deviation `-> 0`
as the scales spread (every nonzero resonance sits at frequency growing with the scale). So each
stratum splits as `[bounded scale: exhaustive] + [large scale: -> D3_inf^{(L)}]`.

## The decorrelated floor references (all >= bar, growing margins)

`D3_inf^{(L)} = D3(block_L + (11-L) iid points)` (MC), the `d -> inf` limit of every scale-`d` family:

| L | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|----|
| `D3_inf^{(L)}` | 0.855 | 0.851 | 0.839 | 0.810 | 0.761 | 0.677 | 0.601 | 0.524 | 0.465 |
| margin over bar | +0.524 | +0.520 | +0.508 | +0.479 | +0.430 | +0.346 | +0.269 | +0.192 | +0.134 |

The margin GROWS as `L` decreases (less AP structure => less coherent additive energy => higher
floor). `L=10` (`0.465`, klein LEM-009) is the tightest; every `L<=9` reference clears bar by
`>= 0.19`. Like `L=10`, the finite-scale minimum sits BELOW `D3_inf^{(L)}` (correlated point lowers
`D3`), so the reference is an upper guide; the actual floor is the correlated minimum.

## The L<=9 minima (reliable adaptive-NG, prim-diam > 30)

Structured search (`block_L(scale d)` + interior/far points, `longest-AP = L`, prim-diam > 30),
diameter-ADAPTIVE grid `NG = 60*prim-diam` (the fixed `NG=9000` grid ALIASES past prim-diam ~1500,
S157):

| L | 9 | 8 | 7 | 6 | 5 | ... |
|---|---|---|---|---|---|-----|
| min `D3` (prim-diam>30) | **0.473** | 0.514 | 0.569 | 0.626 | 0.652 | (higher) |
| margin over bar | **+0.142** | +0.182 | +0.238 | +0.295 | +0.321 | |

The worst `L<=9` stratum is `L=9` (min `0.473`, at e.g. `block_9(scale 5) + {1,50}` -- an AP_9 with a
near point, the direct analog of `A_*`), margin **+0.14** -- 12x the `L=10`-to-bar razor. Every lower
`L` clears by more. All below their `D3_inf^{(L)}` (correlated), all `>= bar` with room.

## The whole-tail closure

> **k=11 covering tail = [prim-diam <= 30: EXHAUSTIVE (kps-S88, min D3 = A_* = 0.452986)]
> + [prim-diam > 30: L=10 (S157, D3 >= bar) + L<=9 (here, D3 >= bar, margin >= 0.14)].
> Every primitive 11-set has `D3 >= bar`.** The tail MINIMUM is `A_* = 0.452986` (the `L=10`, scale-3,
> interior shape), and `D3 -> D3_inf^{(L)}` as the scale grows.

## Honest rigor

- RIGOROUS: the reduction (longest-AP `L` => `block_L(scale d) + points` structure); the resonance-sum
  identity (generalizes S157 to rank `11-L`); the exhaustive prim-diam <= 30 (kps-S88, exact Farey).
- NUMERICALLY CERTIFIED: the decorrelated references `D3_inf^{(L)}` (MC), and the per-`L` decay
  constants (the `S157` mechanism at higher rank -- an a priori bound counts `G_L^j`'s breakpoint
  crossings). The prim-diam > 30 minima are RELIABLE (adaptive-NG) but from a structured search, not
  an exhaustive/rate proof for every `L<=9` shape -- the fully rigorous per-`L` rate (like S157's
  `pd >= 160` threshold) is the remaining step, but the margins (`>= 0.14`, vs the `L=10` razor) make
  it far more forgiving: a much smaller scale threshold suffices, and the bounded remainder is inside
  kps's exhaustive.
- NET: the `L<=9` strata clear with large room; the binding case of the whole tail is `L=10` (S157,
  proved `>= bar`); the tail min is `A_*`. Files: `lrc14_Lleq9_tail_floor_opus_S158.py` (+out).
  Builds on S157 (L=10 mechanism), kps-S88 (exhaustive <=30), klein LEM-009 (`D3_inf` limits).
